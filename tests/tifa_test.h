// Minimal, dependency-free test harness for the TIFA test suite.
//
// Each test is a standalone C program: it runs CHECK* macros and ends with
// TIFA_TEST_END(), which returns non-zero if any check failed. CTest uses the
// process exit code as the pass/fail signal.

#if !defined(_TIFA_TEST_H_)
#define _TIFA_TEST_H_

#include <stdio.h>
#include <gmp.h>

#include "array.h"

static int _tifa_test_failures = 0;
static int _tifa_test_checks   = 0;

#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        _tifa_test_checks++;                                                \
        if (!(cond)) {                                                      \
            _tifa_test_failures++;                                          \
            fprintf(stderr, "  FAIL [%s:%d]: %s\n", __FILE__, __LINE__,     \
                    (msg));                                                 \
        }                                                                   \
    } while (0)

#define CHECK_EQ_UI(actual, expected, msg)                                  \
    do {                                                                    \
        _tifa_test_checks++;                                                \
        unsigned long long _a = (unsigned long long)(actual);              \
        unsigned long long _e = (unsigned long long)(expected);            \
        if (_a != _e) {                                                     \
            _tifa_test_failures++;                                          \
            fprintf(stderr, "  FAIL [%s:%d]: %s (got %llu, expected %llu)\n",\
                    __FILE__, __LINE__, (msg), _a, _e);                     \
        }                                                                   \
    } while (0)

#define TIFA_TEST_BEGIN(name) printf("== %s ==\n", (name))

#define TIFA_TEST_END()                                                     \
    do {                                                                    \
        printf("  %d checks, %d failure(s)\n", _tifa_test_checks,           \
               _tifa_test_failures);                                        \
        return (_tifa_test_failures != 0) ? 1 : 0;                          \
    } while (0)

// Verify that (factors, multis) is a complete prime factorization of n:
// every factor is prime and the product of factors^multiplicities equals n.
// Every factor must have a multiplicity.
static inline int tifa_is_complete_factorization(const mpz_t n,
                                                 const mpz_array_t* factors,
                                                 const uint32_array_t* multis) {
    if (   factors == NULL || factors->length == 0
        || multis == NULL || multis->length != factors->length) {
        return 0;
    }
    mpz_t prod;
    mpz_t tmp;
    mpz_init_set_ui(prod, 1UL);
    mpz_init(tmp);

    int ok = 1;
    for (uint32_t i = 0; i < factors->length; i++) {
        if (mpz_probab_prime_p(factors->data[i], 25) == 0) {
            ok = 0;
        }
        unsigned long e = (unsigned long)multis->data[i];
        mpz_pow_ui(tmp, factors->data[i], e);
        mpz_mul(prod, prod, tmp);
    }
    if (mpz_cmp(prod, n) != 0) {
        ok = 0;
    }
    mpz_clear(prod);
    mpz_clear(tmp);
    return ok;
}

#endif
