// tifa_cofactor (cofactor.h): smooth and non-smooth cofactors with two and
// three large primes, trial division when no factor base bound is known,
// prime powers, the large prime count, the give-up mode, and a randomized
// cross-check against complete factorizations by tifa_factor.

#include <stdlib.h>

#include "tifa_test.h"
#include "cofactor.h"
#include "funcs.h"
#include "tifa_factor.h"

// A random prime of exactly 'bits' bits.
static void random_prime(mpz_t p, gmp_randstate_t state, unsigned long bits) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);
    } while (!mpz_probab_prime_p(p, 30));
}

// The factors returned multiply back to n, are prime and below 2^lpb.
static int factors_ok(const mpz_t n, const mpz_array_t* f,
                      const uint32_array_t* m, uint32_t lpb) {
    if (f->length != m->length) {
        return 0;
    }
    mpz_t prod, t;
    mpz_init_set_ui(prod, 1);
    mpz_init(t);
    int ok = 1;
    for (uint32_t i = 0; i < f->length; i++) {
        if (!is_prime_mpz(f->data[i]) || (mpz_sizeinbase(f->data[i], 2) > lpb)) {
            ok = 0;
        }
        mpz_pow_ui(t, f->data[i], m->data[i]);
        mpz_mul(prod, prod, t);
    }
    if (mpz_cmp(prod, n) != 0) {
        ok = 0;
    }
    mpz_clear(prod);
    mpz_clear(t);
    return ok;
}

// The truth, from a complete factorization.
static int is_smooth(const mpz_t n, uint32_t lpb, const mpz_t fb,
                     uint32_t max_large) {
    mpz_array_t*    f = alloc_mpz_array(8);
    uint32_array_t* m = alloc_uint32_array(8);
    int smooth = 1;
    uint32_t nlarge = 0;
    if ((mpz_cmp_ui(n, 1) > 0) && is_prime_mpz(n)) {
        append_mpz_to_array(f, n);
        append_uint32_to_array(m, 1);
    } else if (mpz_cmp_ui(n, 1) > 0) {
        tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION);
    }
    for (uint32_t i = 0; i < f->length; i++) {
        if (mpz_sizeinbase(f->data[i], 2) > lpb) {
            smooth = 0;
        }
        if (mpz_cmp(f->data[i], fb) >= 0) {
            nlarge += m->data[i];
        }
    }
    if ((max_large != 0) && (nlarge > max_large)) {
        smooth = 0;
    }
    free_mpz_array(f);
    free_uint32_array(m);
    return smooth;
}

static cofactor_status_t run(const mpz_t n, const cofactor_strategy_t* s,
                             uint32_t lpb, int* ok) {
    mpz_array_t*    f = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    const cofactor_status_t status = tifa_cofactor(f, m, n, s);
    *ok = (status != COFACTOR_SMOOTH) || factors_ok(n, f, m, lpb);
    free_mpz_array(f);
    free_uint32_array(m);
    return status;
}

int main(void) {
    TIFA_TEST_BEGIN("cofactor");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261015UL);

    mpz_t n, p, q, r, fb;
    mpz_inits(n, p, q, r, fb, NULL);
    int ok;
    char msg[160];

    cofactor_params_t bad = { .lpb = 64 };
    CHECK(alloc_cofactor_strategy(&bad) == NULL, "lpb 64 refused");

    // ---- Two large primes: lpb 30, no factor below 2^20, at most 2.
    cofactor_params_t two = { .lpb = 30, .fb_bound = 1UL << 20,
                              .max_large_primes = 2 };
    cofactor_strategy_t* s2 = alloc_cofactor_strategy(&two);
    struct { unsigned pb, qb; cofactor_status_t want; const char* what; }
    cases2[] = {
        { 25, 29, COFACTOR_SMOOTH,     "p25 q29" },
        { 29, 30, COFACTOR_SMOOTH,     "p29 q30" },
        { 24, 36, COFACTOR_NOT_SMOOTH, "p24 q36 (q too large)" },
        { 31, 31, COFACTOR_NOT_SMOOTH, "p31 q31 (both too large)" },
    };
    for (uint32_t i = 0; i < sizeof(cases2) / sizeof(cases2[0]); i++) {
        random_prime(p, state, cases2[i].pb);
        random_prime(q, state, cases2[i].qb);
        mpz_mul(n, p, q);
        snprintf(msg, sizeof(msg), "2LP %s", cases2[i].what);
        CHECK(run(n, s2, 30, &ok) == cases2[i].want, msg);
        CHECK(ok, msg);
    }
    random_prime(n, state, 28);
    CHECK(run(n, s2, 30, &ok) == COFACTOR_SMOOTH && ok, "2LP prime below L");
    random_prime(n, state, 45);
    CHECK(run(n, s2, 30, &ok) == COFACTOR_NOT_SMOOTH, "2LP prime above L");
    random_prime(p, state, 22);
    random_prime(q, state, 23);
    random_prime(r, state, 24);
    mpz_mul(n, p, q);
    mpz_mul(n, n, r);
    CHECK(run(n, s2, 30, &ok) == COFACTOR_NOT_SMOOTH, "2LP three primes");
    mpz_set_ui(n, 1);
    CHECK(run(n, s2, 30, &ok) == COFACTOR_SMOOTH && ok, "n = 1");
    mpz_mul(n, p, p);
    CHECK(run(n, s2, 30, &ok) == COFACTOR_SMOOTH && ok, "2LP p^2");

    // ---- Three large primes, above 64 bits: lpb 31, no factor below 2^22.
    cofactor_params_t three = { .lpb = 31, .fb_bound = 1UL << 22,
                                .max_large_primes = 3 };
    cofactor_strategy_t* s3 = alloc_cofactor_strategy(&three);
    three.max_siqs_bits = 64;
    cofactor_strategy_t* s3g = alloc_cofactor_strategy(&three);
    uint32_t nwrong = 0, nbad = 0, nmissed = 0;
    for (uint32_t i = 0; i < 20; i++) {
        random_prime(p, state, 27 + i % 4);
        random_prime(q, state, 29 + i % 3);
        random_prime(r, state, 30 + i % 2);
        mpz_mul(n, p, q);
        mpz_mul(n, n, r);
        nwrong += (run(n, s3, 31, &ok) != COFACTOR_SMOOTH);
        nbad   += !ok;
        // Giving up is allowed, rarely; a wrong answer never is.
        const cofactor_status_t status = run(n, s3g, 31, &ok);
        nwrong  += (status == COFACTOR_NOT_SMOOTH);
        nmissed += (status == COFACTOR_GAVE_UP);
        nbad    += !ok;
    }
    CHECK_EQ_UI(nwrong, 0, "3LP smooth (exact and give-up modes)");
    CHECK(nmissed <= 1, "3LP smooth: at most one given up");
    CHECK_EQ_UI(nbad, 0, "3LP factors");
    random_prime(p, state, 25);
    mpz_mul(n, p, p);
    random_prime(q, state, 30);
    mpz_mul(n, n, q);
    mpz_mul(n, n, q);
    CHECK(run(n, s3, 31, &ok) == COFACTOR_NOT_SMOOTH, "3LP p^2 q^2 (four)");
    mpz_pow_ui(n, q, 3);
    CHECK(run(n, s3, 31, &ok) == COFACTOR_SMOOTH && ok, "3LP q^3");
    // A 65-bit composite with no factor below 2^22 has two prime factors,
    // so one is too large: decided without factoring.
    random_prime(p, state, 32);
    random_prime(q, state, 33);
    mpz_mul(n, p, q);
    CHECK(run(n, s3, 31, &ok) == COFACTOR_NOT_SMOOTH, "3LP 65-bit p q");
    // Two 40-bit primes: proven not smooth by SIQS, or given up.
    random_prime(p, state, 40);
    random_prime(q, state, 40);
    mpz_mul(n, p, q);
    CHECK(run(n, s3, 31, &ok) == COFACTOR_NOT_SMOOTH, "3LP p40 q40 exact");
    CHECK(run(n, s3g, 31, &ok) != COFACTOR_SMOOTH, "3LP p40 q40 give-up");

    // ---- No factor base bound: small primes are divided out first.
    cofactor_params_t none = { .lpb = 32 };
    cofactor_strategy_t* s0 = alloc_cofactor_strategy(&none);
    random_prime(p, state, 31);
    mpz_mul_ui(n, p, 32UL * 3 * 3 * 1009);
    CHECK(run(n, s0, 32, &ok) == COFACTOR_SMOOTH && ok, "small primes");
    mpz_mul_ui(n, n, 4294967311UL);
    CHECK(run(n, s0, 32, &ok) == COFACTOR_NOT_SMOOTH, "small primes, 33-bit");

    // ---- Random cross-check, exact mode: products of 1 to 4 primes of 21
    // to 40 bits, against complete factorizations.
    cofactor_params_t mix = { .lpb = 33, .fb_bound = 1UL << 20,
                              .max_large_primes = 3 };
    cofactor_strategy_t* sm = alloc_cofactor_strategy(&mix);
    mpz_set_ui(fb, 1UL << 20);
    uint32_t nmismatch = 0, nsmooth = 0;
    nbad = 0;
    for (uint32_t i = 0; i < 300; i++) {
        mpz_set_ui(n, 1);
        const uint32_t np = 1 + gmp_urandomm_ui(state, 4);
        for (uint32_t j = 0; j < np; j++) {
            random_prime(p, state, 21 + gmp_urandomm_ui(state, 20));
            mpz_mul(n, n, p);
        }
        const int want = is_smooth(n, 33, fb, 3);
        const cofactor_status_t got = run(n, sm, 33, &ok);
        if (((got == COFACTOR_SMOOTH) != want) || (got == COFACTOR_GAVE_UP)) {
            gmp_fprintf(stderr, "  %Zd: got %d, smooth %d\n", n, got, want);
            nmismatch++;
        }
        nbad    += !ok;
        nsmooth += want;
    }
    CHECK_EQ_UI(nmismatch, 0, "random products: exact answers");
    CHECK_EQ_UI(nbad, 0, "random products: factors");
    CHECK(nsmooth > 50 && nsmooth < 250, "random products: both outcomes");

    free_cofactor_strategy(s2);
    free_cofactor_strategy(s3);
    free_cofactor_strategy(s3g);
    free_cofactor_strategy(s0);
    free_cofactor_strategy(sm);
    mpz_clears(n, p, q, r, fb, NULL);
    gmp_randclear(state);

    TIFA_TEST_END();
}
