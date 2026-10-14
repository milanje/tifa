// SIQS reliability: complete factorizations of semiprimes in the range where
// SIQS used to give up, plus inputs that exercise its retry and large-prime
// edge cases. SIQS is deterministic for a given n, so failures reproduce.

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

    TIFA_TEST_END();
}
