// Copyright (C) 2026 Jerome Milan.
//
// This file is part of TIFA.
//
// TIFA is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.

// End-to-end correctness tests for the generic tifa_factor() entry point.
//
// Each composite is built from known primes, then asked for a COMPLETE
// factorization; we verify the returned factors are all prime and multiply
// back (with multiplicities) to n. This exercises the size-based dispatch to
// factor_u64 (trial division / SQUFOF without it) / ECM and p-1 pre-passes /
// SIQS, including inputs with small factors, and with each hint of
// tifa_factor_with_hint().

#include "tifa_test.h"
#include "tifa_factor.h"
#include "mont_ecm.h"
#include "time_limit.h"

static void factor_and_check(const char* label, const mpz_t n) {
    mpz_array_t*    factors = alloc_mpz_array(8);
    uint32_array_t* multis  = alloc_uint32_array(8);

    ecode_t rc = tifa_factor(factors, multis, n, FIND_COMPLETE_FACTORIZATION);

    char msg[160];
    snprintf(msg, sizeof(msg), "%s: rc == COMPLETE_FACTORIZATION_FOUND", label);
    CHECK(rc == COMPLETE_FACTORIZATION_FOUND, msg);

    snprintf(msg, sizeof(msg), "%s: factors are prime and multiply back to n",
             label);
    CHECK(tifa_is_complete_factorization(n, factors, multis), msg);

    free_mpz_array(factors);
    free_uint32_array(multis);
}

// A hint changes the methods tried, never the factors found.
static void hinted_and_check(const char* label, const mpz_t n,
                             const tifa_factor_hint_t hint) {
    static const char* const names[] = { "random", "balanced", "unbalanced" };
    mpz_array_t*    factors = alloc_mpz_array(8);
    uint32_array_t* multis  = alloc_uint32_array(8);

    ecode_t rc = tifa_factor_with_hint(factors, multis, n,
                                       FIND_COMPLETE_FACTORIZATION, hint);
    char msg[160];
    snprintf(msg, sizeof(msg), "%s, %s: complete factorization", label,
             names[hint]);
    CHECK(   (rc == COMPLETE_FACTORIZATION_FOUND)
          && tifa_is_complete_factorization(n, factors, multis), msg);

    free_mpz_array(factors);
    free_uint32_array(multis);
}

// Modes other than FIND_COMPLETE_FACTORIZATION may pass a NULL multis array.
// An algorithm that splits n completely (SIQS) may say so in any mode.
static void some_factors_and_check(const char* label, const mpz_t n) {
    mpz_array_t* factors = alloc_mpz_array(8);

    ecode_t rc = tifa_factor(factors, NULL, n, FIND_SOME_FACTORS);

    char msg[160];
    snprintf(msg, sizeof(msg), "%s: rc == SOME_FACTORS_FOUND", label);
    CHECK(   (rc == SOME_FACTORS_FOUND)
          || (rc == COMPLETE_FACTORIZATION_FOUND), msg);

    int ok = (factors->length > 0);
    for (uint32_t i = 0; i < factors->length; i++) {
        if (   mpz_cmp_ui(factors->data[i], 1UL) <= 0
            || mpz_cmp(factors->data[i], n) >= 0
            || !mpz_divisible_p(n, factors->data[i])) {
            ok = 0;
        }
    }
    snprintf(msg, sizeof(msg), "%s: factors are proper divisors of n", label);
    CHECK(ok, msg);

    free_mpz_array(factors);
}

int main(void) {
    TIFA_TEST_BEGIN("tifa_factor (complete factorization)");

    mpz_t n;
    mpz_t p;
    mpz_init(n);
    mpz_init(p);

    // 15 = 3 * 5
    mpz_set_ui(n, 3UL);  mpz_mul_ui(n, n, 5UL);
    factor_and_check("15 = 3*5", n);

    // 8051 = 83 * 97
    mpz_set_ui(n, 83UL); mpz_mul_ui(n, n, 97UL);
    factor_and_check("8051 = 83*97", n);

    // 720 = 2^4 * 3^2 * 5 (repeated factors, tests multiplicities)
    mpz_ui_pow_ui(n, 2UL, 4UL);
    mpz_ui_pow_ui(p, 3UL, 2UL); mpz_mul(n, n, p);
    mpz_mul_ui(n, n, 5UL);
    factor_and_check("720 = 2^4*3^2*5", n);

    // 1234567 = 127 * 9721
    mpz_set_ui(n, 127UL); mpz_mul_ui(n, n, 9721UL);
    factor_and_check("1234567 = 127*9721", n);

    // ~2^40 = 1000003 * 1000033
    mpz_set_ui(n, 1000003UL); mpz_mul_ui(n, n, 1000033UL);
    factor_and_check("1000003*1000033 (~2^40)", n);

    // ~2^60 = 1000000007 * 1000000009 (factor_u64, or SQUFOF without it)
    mpz_set_ui(n, 1000000007UL); mpz_mul_ui(n, n, 1000000009UL);
    factor_and_check("1000000007*1000000009 (~2^60)", n);

    // 26 bits, both factors just below sqrt(2^26) = 8192
    mpz_set_ui(n, 8167UL); mpz_mul_ui(n, n, 8171UL);
    factor_and_check("8167*8171 (26 bits)", n);
    mpz_set_ui(n, 8179UL); mpz_mul_ui(n, n, 8191UL);
    factor_and_check("8179*8191 (26 bits)", n);

    // Small factors on numbers too large for the trial division branch
    mpz_set_str(n, "15318842556695263238", 10);
    factor_and_check("2*7659421278347631619 (64 bits)", n);
    mpz_set_str(n, "790617260706", 10);
    factor_and_check("2*3*7*13*19*76211419 (40 bits)", n);
    mpz_set_str(n, "15266532123449400042", 10);
    factor_and_check("2*3^2*619*1370178794062951 (64 bits)", n);
    mpz_set_str(n, "7407147867164031811", 10);
    factor_and_check("11*673377078833093801 (63 bits)", n);

    // Square cofactor once the small primes are removed: 4 * r^2
    mpz_set_str(n, "3308912282058341325890884", 10);
    factor_and_check("2^2*909520791689^2 (82 bits)", n);

    // Other powers: no congruence of squares splits a prime power
    mpz_set_str(n, "48703862963609784559", 10);
    factor_and_check("3651919^3 (66 bits)", n);
    mpz_set_str(n, "1329227995839317534787207543490874671", 10);
    factor_and_check("(2^40+15)^3 (121 bits)", n);
    mpz_set_str(n, "1329231983523305052739811905113505143624013", 10);
    factor_and_check("1000003*(2^40+15)^3 (141 bits)", n);
    mpz_set_str(n, "1766847064923008753571539193567776624376842367262499964598644648629358241", 10);
    factor_and_check("(2^40+15)^6 (241 bits)", n);

    // Primes, including those within the trial division bound
    static const unsigned long primes[] = { 3, 997, 8191, 10007, 1000003 };
    for (unsigned i = 0; i < sizeof(primes) / sizeof(primes[0]); i++) {
        mpz_set_ui(n, primes[i]);
        mpz_array_t* factors = alloc_mpz_array(4);
        CHECK(tifa_factor(factors, NULL, n, FIND_SOME_FACTORS)
                  == NO_FACTOR_FOUND, "a prime has no factor");
        free_mpz_array(factors);
    }

    // 128 bits, p - 1 is 2000-smooth: found by the p-1 pre-stage without
    // ECM, by SIQS with it
    mpz_set_str(n, "206010655745875572166197184221848916071", 10);
    factor_and_check("340730379077795843*604614875560711445197 (128 bits)", n);
    some_factors_and_check("FIND_SOME_FACTORS (128 bits)", n);

    mpz_set_ui(n, 8051UL);
    some_factors_and_check("FIND_SOME_FACTORS, 8051 = 83*97", n);
    mpz_set_str(n, "15318842556695263238", 10);
    some_factors_and_check("FIND_SOME_FACTORS, 2*7659421278347631619", n);

    // Random composites up to 64 bits (factor_u64; trial division and SQUFOF
    // without it)
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261014UL);
    uint32_t nincomplete = 0;
    for (uint32_t i = 0; i < 400; i++) {
        uint32_t bits = 8 + (i % 57);
        do {
            mpz_urandomb(n, state, bits);
            mpz_setbit(n, bits - 1);
        } while (mpz_probab_prime_p(n, 25));

        mpz_array_t*    factors = alloc_mpz_array(8);
        uint32_array_t* multis  = alloc_uint32_array(8);
        ecode_t rc = tifa_factor(factors, multis, n, FIND_COMPLETE_FACTORIZATION);
        if (   rc != COMPLETE_FACTORIZATION_FOUND
            || !tifa_is_complete_factorization(n, factors, multis)) {
            gmp_fprintf(stderr, "  incomplete: %Zd (rc = %d)\n", n, rc);
            nincomplete++;
        }
        free_mpz_array(factors);
        free_uint32_array(multis);
    }
    CHECK_EQ_UI(nincomplete, 0, "random 8-64-bit composites fully factored");

    // Primes of every size: above 64 bits, they used to go through p-1, p+1
    // and SIQS until these gave up.
    uint32_t nmisreported = 0;
    for (uint32_t bits = 20; bits <= 200; bits += 15) {
        do {
            mpz_urandomb(n, state, bits);
            mpz_setbit(n, bits - 1);
        } while (!mpz_probab_prime_p(n, 25));

        mpz_array_t*    factors = alloc_mpz_array(8);
        uint32_array_t* multis  = alloc_uint32_array(8);
        ecode_t rc = tifa_factor(factors, multis, n, FIND_COMPLETE_FACTORIZATION);
        if ((rc != NO_FACTOR_FOUND) || (factors->length != 0)) {
            gmp_fprintf(stderr, "  prime %Zd: rc = %d\n", n, rc);
            nmisreported++;
        }
        free_mpz_array(factors);
        free_uint32_array(multis);
    }
    CHECK_EQ_UI(nmisreported, 0, "primes of 20-200 bits: NO_FACTOR_FOUND");

    for (int h = TIFA_FACTOR_RANDOM; h <= TIFA_FACTOR_UNBALANCED; h++) {
        const tifa_factor_hint_t hint = (tifa_factor_hint_t) h;

        mpz_set_ui(p, 1UL);
        mpz_mul_2exp(p, p, 49);
        mpz_nextprime(p, p);
        mpz_set_ui(n, 1UL);
        mpz_mul_2exp(n, n, 50);
        mpz_nextprime(n, n);
        mpz_mul(n, n, p);
        hinted_and_check("balanced semiprime (100 bits)", n, hint);

        // What trial division leaves keeps the hint.
        mpz_mul_ui(n, n, 6UL);
        hinted_and_check("6 * balanced semiprime (103 bits)", n, hint);

        mpz_set_ui(p, 1UL);
        mpz_mul_2exp(p, p, 34);
        mpz_nextprime(p, p);
        mpz_set_ui(n, 1UL);
        mpz_mul_2exp(n, n, 85);
        mpz_nextprime(n, n);
        mpz_mul(n, n, p);
        hinted_and_check("35-bit factor (120 bits)", n, hint);

        mpz_set_str(n, "1329231983523305052739811905113505143624013", 10);
        hinted_and_check("1000003*(2^40+15)^3 (141 bits)", n, hint);

        for (uint32_t bits = 66; bits <= 156; bits += 10) {
            do {
                mpz_urandomb(n, state, bits);
                mpz_setbit(n, bits - 1);
            } while (mpz_probab_prime_p(n, 25));
            hinted_and_check("random composite (66-156 bits)", n, hint);
        }
        mpz_set_str(n, "206010655745875572166197184221848916071", 10);
        mpz_array_t* factors = alloc_mpz_array(4);
        ecode_t rc = tifa_factor_with_hint(factors, NULL, n,
                                           FIND_SOME_FACTORS, hint);
        CHECK(   ((rc == SOME_FACTORS_FOUND)
                  || (rc == COMPLETE_FACTORIZATION_FOUND))
              && (factors->length > 0)
              && mpz_divisible_p(n, factors->data[0]),
              "FIND_SOME_FACTORS with a hint (128 bits)");

        mpz_nextprime(n, n);
        reset_mpz_array(factors);
        CHECK(tifa_factor_with_hint(factors, NULL, n, FIND_SOME_FACTORS, hint)
                  == NO_FACTOR_FOUND, "a prime has no factor, with a hint");
        free_mpz_array(factors);
    }
    mpz_array_t* factors = alloc_mpz_array(4);
    CHECK(tifa_factor_with_hint(factors, NULL, n, FIND_SOME_FACTORS,
                                (tifa_factor_hint_t) 3) == FAILURE,
          "an unknown hint fails");
    CHECK_EQ_UI(factors->length, 0, "an unknown hint finds nothing");
    free_mpz_array(factors);

#if TIFA_HAVE_MONT_ECM
    // Above 256 bits, the ECM pre-pass on five to eight limbs: factors of 30
    // to 50 bits, the last cofactor prime (SIQS would take minutes or more).
    static const uint32_t prime_bits[][4] = {
        { 45, 235, 0, 0 }, { 40, 50, 300, 0 }, { 30, 40, 45, 385 }
    };
    static const char* const labels[] = {
        "45-bit factor (280 bits)", "40- and 50-bit factors (390 bits)",
        "30-, 40- and 45-bit factors (500 bits)"
    };
    for (uint32_t i = 0; i < 3; i++) {
        mpz_set_ui(n, 1UL);
        for (uint32_t j = 0; (j < 4) && (prime_bits[i][j] != 0); j++) {
            mpz_urandomb(p, state, prime_bits[i][j]);
            mpz_setbit(p, prime_bits[i][j] - 1);
            mpz_nextprime(p, p);
            mpz_mul(n, n, p);
        }
        factor_and_check(labels[i], n);
    }
#endif
    // Above 512 bits, ECM in GMP arithmetic: factors of 30 to 40 bits, the
    // last cofactor prime.
    static const uint32_t wide_bits[][3] = { { 35, 40, 560 }, { 30, 1000, 0 } };
    static const char* const wide_labels[] = {
        "35- and 40-bit factors (635 bits)", "30-bit factor (1030 bits)"
    };
    for (uint32_t i = 0; i < 2; i++) {
        mpz_set_ui(n, 1UL);
        for (uint32_t j = 0; (j < 3) && (wide_bits[i][j] != 0); j++) {
            mpz_urandomb(p, state, wide_bits[i][j]);
            mpz_setbit(p, wide_bits[i][j] - 1);
            mpz_nextprime(p, p);
            mpz_mul(n, n, p);
        }
        factor_and_check(wide_labels[i], n);
    }
#if TIFA_HAVE_MONT_ECM
    // The ECM pre-pass in 4 threads.
    static const uint32_t threaded_bits[][3] = { { 40, 160, 0 }, { 45, 50, 205 } };
    static const char* const threaded_labels[] = {
        "40-bit factor (200 bits), 4 threads",
        "45- and 50-bit factors (300 bits), 4 threads"
    };
    tifa_set_threads(4);
    for (uint32_t i = 0; i < 2; i++) {
        mpz_set_ui(n, 1UL);
        for (uint32_t j = 0; (j < 3) && (threaded_bits[i][j] != 0); j++) {
            mpz_urandomb(p, state, threaded_bits[i][j]);
            mpz_setbit(p, threaded_bits[i][j] - 1);
            mpz_nextprime(p, p);
            mpz_mul(n, n, p);
        }
        factor_and_check(threaded_labels[i], n);
    }
    tifa_set_threads(1);
#endif

    gmp_randclear(state);

    mpz_clear(n);
    mpz_clear(p);
    TIFA_TEST_END();
}
