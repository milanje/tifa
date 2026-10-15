// Unit tests for the pure numeric primitives in lib/utils/src/funcs.c.
//
// The modinv_ui cases with ~2^31 and ~2^32 prime moduli are regression tests
// for the int-shift bug (funcs.c: `1 << n` with n == 31/32).

#include <limits.h>
#include <stdbool.h>

#include "tifa_test.h"
#include "funcs.h"

int main(void) {
    TIFA_TEST_BEGIN("funcs primitives");

    // most_significant_bit = floor(log2(n))
    CHECK_EQ_UI(most_significant_bit(1U),          0, "msb(1)");
    CHECK_EQ_UI(most_significant_bit(2U),          1, "msb(2)");
    CHECK_EQ_UI(most_significant_bit(255U),        7, "msb(255)");
    CHECK_EQ_UI(most_significant_bit(256U),        8, "msb(256)");
    CHECK_EQ_UI(most_significant_bit(0x80000000U), 31, "msb(2^31)");

    // powm: modular exponentiation
    CHECK_EQ_UI(powm(2U, 10U, 1000U), 24U, "2^10 mod 1000");
    CHECK_EQ_UI(powm(3U, 0U, 7U),      1U, "3^0 mod 7");
    CHECK_EQ_UI(powm(7U, 2U, 10U),     9U, "7^2 mod 10");

    // sqrtm: verify by squaring (2 and 4 are quadratic residues mod 7)
    {
        uint32_t r = sqrtm(2U, 7U);
        CHECK_EQ_UI((r * r) % 7U, 2U, "sqrtm(2,7)^2 == 2 mod 7");
        r = sqrtm(4U, 7U);
        CHECK_EQ_UI((r * r) % 7U, 4U, "sqrtm(4,7)^2 == 4 mod 7");
    }

    // is_square: returns the root for a perfect square, 0 otherwise
    CHECK_EQ_UI(is_square(1UL),     1UL,   "is_square(1)");
    CHECK_EQ_UI(is_square(144UL),   12UL,  "is_square(144)");
    CHECK_EQ_UI(is_square(10000UL), 100UL, "is_square(10000)");
    CHECK_EQ_UI(is_square(145UL),   0UL,   "is_square(145) == 0");

    // is_prime
    CHECK(is_prime(2U),           "is_prime(2)");
    CHECK(is_prime(97U),          "is_prime(97)");
    CHECK(is_prime(1000000007U),  "is_prime(1000000007)");
    CHECK(is_prime(4294967291U),  "is_prime(2^32 - 5)");
    CHECK(!is_prime(1U),          "!is_prime(1)");
    CHECK(!is_prime(4U),          "!is_prime(4)");
    CHECK(!is_prime(8051U),       "!is_prime(8051 = 83*97)");

    // gcd_ulint
    CHECK_EQ_UI(gcd_ulint(12UL, 18UL),    6UL,  "gcd(12,18)");
    CHECK_EQ_UI(gcd_ulint(1071UL, 462UL), 21UL, "gcd(1071,462)");
    CHECK_EQ_UI(gcd_ulint(17UL, 5UL),     1UL,  "gcd(17,5)");

    // modinv_ui: verify n * inv == 1 (mod p). The large moduli drive the
    // internal bit index n to 31 and 32 respectively -> exercises the shift fix.
    {
        struct { unsigned long n, p; } cases[] = {
            {3UL,         7UL},
            {2UL,         2147483647UL},  // p = 2^31 - 1 (Mersenne prime) -> n == 31
            {123456789UL, 4294967291UL},  // p = 2^32 - 5 (prime)          -> n == 32
            {987654321UL, 4294967291UL},
        };
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
            unsigned long p   = cases[i].p;
            unsigned long a   = cases[i].n;
            unsigned long inv = modinv_ui(a, p);
            // a, inv < p < 2^32, so the product fits in 64 bits.
            CHECK_EQ_UI((a * inv) % p, 1UL, "modinv_ui: n*inv == 1 mod p");
        }
    }

    // is_prime_mpz agrees with GMP on small and negative integers, around
    // 2^32 and 2^64, and on strong pseudoprimes to many bases.
    {
        mpz_t n;
        mpz_init(n);
        uint32_t ndisagree = 0;
        for (long i = -100; i < 70000; i++) {
            mpz_set_si(n, i);
            ndisagree += (is_prime_mpz(n) != (mpz_probab_prime_p(n, 50) != 0));
        }
        for (unsigned long e = 32; e <= 64; e += 32) {
            for (long d = -3000; d <= 3000; d++) {
                mpz_ui_pow_ui(n, 2, e);
                if (d < 0) {
                    mpz_sub_ui(n, n, (unsigned long) -d);
                } else {
                    mpz_add_ui(n, n, (unsigned long) d);
                }
                ndisagree += (is_prime_mpz(n) != (mpz_probab_prime_p(n, 50) != 0));
            }
        }
        const char* const spsp[] = {
            "3215031751", "2152302898747", "3474749660383", "341550071728321",
            "3825123056546413051", "318665857834031151167461",
            "3317044064679887385961981"
        };
        for (size_t i = 0; i < sizeof(spsp) / sizeof(spsp[0]); i++) {
            mpz_set_str(n, spsp[i], 10);
            ndisagree += is_prime_mpz(n);
        }
        CHECK_EQ_UI(ndisagree, 0, "is_prime_mpz agrees with GMP");
        mpz_clear(n);
    }

#if ULONG_MAX > 4294967295UL
    // modinv_ui with moduli above 2^32 (all wrong when the bit length of p
    // was taken on 32 bits), up to the largest prime below 2^63.
    {
        const char* const primes[] = {
            "4294967311", "1099511627689", "2305843009213693951",
            "9223372036854775783"
        };
        const unsigned long ns[] = { 2UL, 3UL, 1234567891011UL, ULONG_MAX };
        mpz_t p, n, inv;
        mpz_inits(p, n, inv, NULL);
        for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); i++) {
            mpz_set_str(p, primes[i], 10);
            for (size_t j = 0; j < sizeof(ns) / sizeof(ns[0]); j++) {
                mpz_set_ui(n, ns[j]);
                mpz_invert(inv, n, p);
                CHECK_EQ_UI(modinv_ui(ns[j], mpz_get_ui(p)), mpz_get_ui(inv),
                            "modinv_ui above 2^32 matches mpz_invert");
            }
        }
        mpz_clears(p, n, inv, NULL);
    }
#endif

    TIFA_TEST_END();
}
