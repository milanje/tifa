// Word-size primality and factorization (factor_u64.h), checked against GMP:
// known primes and strong pseudoprimes, edge values, and random integers of
// every size up to 64 bits.

#include <stdlib.h>

#include "tifa_test.h"
#include "factor_u64.h"

#if TIFA_HAVE_FACTOR_U64

static int gmp_is_prime(uint64_t n) {
    mpz_t z;
    mpz_init_set_ui(z, n);
    int r = mpz_probab_prime_p(z, 40) != 0;
    mpz_clear(z);
    return r;
}

// The factorization is complete, made of increasing primes, and multiplies
// back to n.
static int check_factorization(uint64_t n) {
    uint64_t primes[FACTOR_U64_MAX_PRIMES];
    uint32_t multis[FACTOR_U64_MAX_PRIMES];
    uint32_t nprimes = factor_u64(n, primes, multis);
    if (nprimes == 0) {
        return 0;
    }
    mpz_t prod;
    mpz_t pk;
    mpz_init_set_ui(prod, 1);
    mpz_init(pk);
    int ok = 1;
    for (uint32_t i = 0; i < nprimes; i++) {
        if (!gmp_is_prime(primes[i]) || (multis[i] == 0)) {
            ok = 0;
        }
        if ((i > 0) && (primes[i] <= primes[i - 1])) {
            ok = 0;
        }
        mpz_ui_pow_ui(pk, primes[i], multis[i]);
        mpz_mul(prod, prod, pk);
    }
    if (mpz_cmp_ui(prod, n) != 0) {
        ok = 0;
    }
    mpz_clear(prod);
    mpz_clear(pk);
    return ok;
}

#endif

int main(void) {
    TIFA_TEST_BEGIN("factor_u64");

#if TIFA_HAVE_FACTOR_U64
    char msg[160];

    // Primes, including the largest below 2^32 and 2^64.
    const uint64_t primes[] = {
        2, 3, 5, 997, 1009, 65537, 4294967291ULL, 2305843009213693951ULL,
        18446744073709551557ULL
    };
    for (uint32_t i = 0; i < sizeof(primes) / sizeof(primes[0]); i++) {
        snprintf(msg, sizeof(msg), "%llu is prime",
                 (unsigned long long) primes[i]);
        CHECK(is_prime_u64(primes[i]), msg);
    }
    // Composites, including strong pseudoprimes to several small bases
    // and Carmichael numbers.
    const uint64_t composites[] = {
        0, 1, 4, 9, 561, 1000001, 25326001, 3215031751ULL, 2152302898747ULL,
        3474749660383ULL, 341550071728321ULL, 3825123056546413051ULL,
        18446744073709551615ULL
    };
    for (uint32_t i = 0; i < sizeof(composites) / sizeof(composites[0]); i++) {
        snprintf(msg, sizeof(msg), "%llu is not prime",
                 (unsigned long long) composites[i]);
        CHECK(!is_prime_u64(composites[i]), msg);
    }

    // Edge values and structured inputs.
    const uint64_t specials[] = {
        2, 4, 1ULL << 63, 18446744073709551615ULL, 18446744073709551557ULL,
        1018081ULL,                              // 1009^2
        1018081ULL * 1009ULL,                    // 1009^3
        4294967291ULL * 4294967279ULL,           // two primes near 2^32
        4294967291ULL * 4294967291ULL,           // square of a prime > 2^31
        1000003ULL * 1000003ULL * 1009ULL,       // p^2 q
        3825123056546413051ULL,                  // spsp to bases 2..23
        4294967297ULL,                           // 2^32 + 1 = 641 * 6700417
        9007199254740881ULL                      // 2^53 - 111, prime
    };
    for (uint32_t i = 0; i < sizeof(specials) / sizeof(specials[0]); i++) {
        snprintf(msg, sizeof(msg), "factor_u64(%llu)",
                 (unsigned long long) specials[i]);
        CHECK(check_factorization(specials[i]), msg);
    }

    // Random integers and balanced semiprimes of every size from 2 to 64
    // bits; primality against GMP for each.
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261014UL);
    mpz_t z;
    mpz_t p;
    mpz_init(z);
    mpz_init(p);

    uint32_t nbad_prime = 0;
    uint32_t nbad_fact  = 0;
    for (uint32_t i = 0; i < 6000; i++) {
        const uint32_t bits = 2 + (i % 63);
        mpz_urandomb(z, state, bits);
        mpz_setbit(z, bits - 1);
        uint64_t n = (uint64_t) mpz_get_ui(z);

        if (is_prime_u64(n) != gmp_is_prime(n)) {
            fprintf(stderr, "  is_prime_u64(%llu) wrong\n",
                    (unsigned long long) n);
            nbad_prime++;
        }
        if (!check_factorization(n)) {
            fprintf(stderr, "  factor_u64(%llu) wrong\n",
                    (unsigned long long) n);
            nbad_fact++;
        }
        if (bits >= 8) {
            mpz_urandomb(z, state, bits / 2);
            mpz_setbit(z, bits / 2 - 1);
            mpz_nextprime(p, z);
            mpz_urandomb(z, state, bits - bits / 2);
            mpz_setbit(z, bits - bits / 2 - 1);
            mpz_nextprime(z, z);
            mpz_mul(z, z, p);
            if (mpz_sizeinbase(z, 2) <= 64) {
                n = (uint64_t) mpz_get_ui(z);
                if (!check_factorization(n)) {
                    fprintf(stderr, "  factor_u64(%llu) wrong\n",
                            (unsigned long long) n);
                    nbad_fact++;
                }
            }
        }
    }
    CHECK_EQ_UI(nbad_prime, 0, "is_prime_u64 agrees with GMP");
    CHECK_EQ_UI(nbad_fact, 0, "factor_u64 gives complete factorizations");

    mpz_clear(z);
    mpz_clear(p);
    gmp_randclear(state);
#endif

    TIFA_TEST_END();
}
