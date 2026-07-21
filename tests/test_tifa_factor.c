// End-to-end correctness tests for the generic tifa_factor() entry point.
//
// Each composite is built from known primes, then asked for a COMPLETE
// factorization; we verify the returned factors are all prime and multiply
// back (with multiplicities) to n. This exercises the size-based dispatch to
// trial division / SQUFOF / SIQS (the ~2^60 case reaches SIQS).

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

    // ~2^60 = 1000000007 * 1000000009 (dispatches to SIQS)
    mpz_set_ui(n, 1000000007UL); mpz_mul_ui(n, n, 1000000009UL);
    factor_and_check("1000000007*1000000009 (~2^60)", n);

    mpz_clear(n);
    mpz_clear(p);
    TIFA_TEST_END();
}
