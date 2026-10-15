// End-to-end correctness tests for the generic tifa_factor() entry point.
//
// Each composite is built from known primes, then asked for a COMPLETE
// factorization; we verify the returned factors are all prime and multiply
// back (with multiplicities) to n. This exercises the size-based dispatch to
// trial division / SQUFOF / p-1 / SIQS, including inputs with small factors.

#include "tifa_test.h"
#include "tifa_factor.h"

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

// Modes other than FIND_COMPLETE_FACTORIZATION may pass a NULL multis array.
static void some_factors_and_check(const char* label, const mpz_t n) {
    mpz_array_t* factors = alloc_mpz_array(8);

    ecode_t rc = tifa_factor(factors, NULL, n, FIND_SOME_FACTORS);

    char msg[160];
    snprintf(msg, sizeof(msg), "%s: rc == SOME_FACTORS_FOUND", label);
    CHECK(rc == SOME_FACTORS_FOUND, msg);

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

    // ~2^60 = 1000000007 * 1000000009 (dispatches to SQUFOF)
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

    // 128 bits, p - 1 is 2000-smooth: found by the p-1 pre-stage
    mpz_set_str(n, "206010655745875572166197184221848916071", 10);
    factor_and_check("340730379077795843*604614875560711445197 (128 bits)", n);
    some_factors_and_check("FIND_SOME_FACTORS, p-1 pre-stage (128 bits)", n);

    mpz_set_ui(n, 8051UL);
    some_factors_and_check("FIND_SOME_FACTORS, 8051 = 83*97", n);
    mpz_set_str(n, "15318842556695263238", 10);
    some_factors_and_check("FIND_SOME_FACTORS, 2*7659421278347631619", n);

    // Random composites over the trial division and SQUFOF ranges
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
    gmp_randclear(state);

    mpz_clear(n);
    mpz_clear(p);
    TIFA_TEST_END();
}
