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

// SIQS reliability: complete factorizations of semiprimes in the range where
// SIQS used to give up, plus inputs that exercise its retry and large-prime
// edge cases, and parameters that used to make it run forever. SIQS is
// deterministic for a given n, so failures reproduce.

#include "tifa_test.h"

#include "siqs.h"

static int siqs_factors(const mpz_t n) {
    mpz_array_t*    f = alloc_mpz_array(8);
    uint32_array_t* m = alloc_uint32_array(8);
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);

    ecode_t rc = siqs(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
    int ok = (rc == COMPLETE_FACTORIZATION_FOUND)
             && tifa_is_complete_factorization(n, f, m);

    free_mpz_array(f);
    free_uint32_array(m);
    return ok;
}

// SIQS with explicit factor base size, threshold and sieve half-width.
static ecode_t siqs_with(const mpz_t n, uint32_t fb, uint32_t th, uint32_t hw,
                         int* valid) {
    mpz_array_t*    f = alloc_mpz_array(8);
    uint32_array_t* m = alloc_uint32_array(8);
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);
    params.nprimes_in_base  = fb;
    params.nprimes_tdiv     = fb;
    params.threshold        = th;
    params.sieve_half_width = hw;

    ecode_t rc = siqs(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
    *valid = (rc != COMPLETE_FACTORIZATION_FOUND)
             || tifa_is_complete_factorization(n, f, m);

    free_mpz_array(f);
    free_uint32_array(m);
    return rc;
}

static void check_siqs(const char* label, const char* decimal) {
    mpz_t n;
    mpz_init_set_str(n, decimal, 10);
    CHECK(siqs_factors(n), label);
    mpz_clear(n);
}

static void random_prime(mpz_t p, gmp_randstate_t state, unsigned long bits) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);
        mpz_setbit(p, 0);
    } while (!mpz_probab_prime_p(p, 25));
}

int main(void) {
    TIFA_TEST_BEGIN("siqs reliability");

    // The multiplier (23) divides one of the primes of 'a' for many starts:
    // successive b's repeated and filled the matrix with duplicates.
    check_siqs("56 bits, multiplier prime in a", "45844815609288407");

    // Duplicate relations from different 'a' values: needs the retry.
    check_siqs("80 bits, duplicate relations", "1066778129149670260218683");

    // A large prime (96289) divides n.
    check_siqs("89 bits, large prime divides n", "325578952376929918515281297");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261014UL);

    mpz_t n, p, q;
    mpz_init(n);
    mpz_init(p);
    mpz_init(q);

    uint32_t nfailed = 0;
    for (uint32_t i = 0; i < 300; i++) {
        unsigned long bits = 48 + (i % 17);
        random_prime(p, state, bits / 2);
        random_prime(q, state, bits - bits / 2);
        mpz_mul(n, p, q);

        if (!siqs_factors(n)) {
            gmp_fprintf(stderr, "  siqs failed on %Zd\n", n);
            nfailed++;
        }
    }
    CHECK_EQ_UI(nfailed, 0, "48-64-bit semiprimes fully factored");

    mpz_clear(n);
    mpz_clear(p);
    mpz_clear(q);
    gmp_randclear(state);

    // Parameters that used to make SIQS run forever. A factor base of 46
    // primes gives too few distinct polynomials for this number: once its
    // duplicate relations are dropped, new ones are all duplicates too.
    // SIQS now notices and retries with a wider factor base and sieve.
    {
        mpz_t n;
        mpz_init_set_str(n, "387266287363586789317", 10);
        int valid;
        ecode_t rc = siqs_with(n, 46, 29, 9440, &valid);
        CHECK(rc == COMPLETE_FACTORIZATION_FOUND && valid,
              "stalling factor base: widened and factored");
        //
        // A threshold no sieve value reaches yields no candidate at all:
        // SIQS must give up rather than sieve forever.
        //
        rc = siqs_with(n, 78, 90, 11800, &valid);
        CHECK(valid && (rc == COMPLETE_FACTORIZATION_FOUND || rc == GIVING_UP),
              "threshold out of reach: returns");
        mpz_clear(n);
    }

    TIFA_TEST_END();
}
