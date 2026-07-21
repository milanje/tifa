// Per-algorithm correctness tests: each algorithm is driven directly through
// its public interface on a composite suited to it, asking for a COMPLETE
// factorization and verifying the result. The ECM case sets b2 > 105 to run
// the second phase, exercising the precompute loop that held the abs() bug.

#include "tifa_test.h"

#include "cfrac.h"
#include "ecm.h"
#include "fermat.h"
#include "rho.h"
#include "p_minus_1.h"
#include "siqs.h"
#include "squfof.h"
#include "tdiv.h"

// Build n = a * b with a, b unsigned longs.
static void mul2(mpz_t n, unsigned long a, unsigned long b) {
    mpz_set_ui(n, a);
    mpz_mul_ui(n, n, b);
}

static void check(const char* label, ecode_t rc, const mpz_t n,
                  const mpz_array_t* factors, const uint32_array_t* multis) {
    char msg[160];
    snprintf(msg, sizeof(msg), "%s: rc == COMPLETE_FACTORIZATION_FOUND", label);
    CHECK(rc == COMPLETE_FACTORIZATION_FOUND, msg);
    snprintf(msg, sizeof(msg), "%s: valid complete factorization", label);
    CHECK(tifa_is_complete_factorization(n, factors, multis), msg);
}

int main(void) {
    TIFA_TEST_BEGIN("per-algorithm factorization");

    mpz_t n;
    mpz_init(n);

    // ---- SQUFOF ---- (~2^40 semiprime)
    {
        mul2(n, 1000003UL, 1000033UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        squfof_params_t params;
        set_squfof_params_to_default(&params);
        ecode_t rc = squfof(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("squfof", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- CFRAC ---- (~2^49 semiprime)
    {
        mul2(n, 15485863UL, 32452843UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        cfrac_params_t params;
        set_cfrac_params_to_default(n, &params);
        ecode_t rc = cfrac(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("cfrac", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- SIQS ---- (~2^60 semiprime, SIQS's intended range)
    {
        mul2(n, 1000000007UL, 1000000009UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        siqs_params_t params;
        set_siqs_params_to_default(n, &params);
        ecode_t rc = siqs(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("siqs", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Fermat (McKee) ---- close factors
    {
        mul2(n, 1000003UL, 1000033UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        fermat_params_t params;
        set_fermat_params_to_default(&params);
        ecode_t rc = fermat(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("fermat", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Pollard rho ---- (~2^37 semiprime)
    {
        mul2(n, 104729UL, 1299709UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        rho_params_t params;
        set_rho_params_to_default(n, &params);
        ecode_t rc = rho(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("rho", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Pollard p-1 ---- (10007 - 1 = 2*5003 is smooth)
    {
        mul2(n, 10007UL, 1000003UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_minus_1_params_t params;
        set_p_minus_1_params_to_default(n, &params);
        ecode_t rc = p_minus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_minus_1", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Trial division ---- 210 = 2*3*5*7, no mode parameter
    {
        mpz_set_ui(n, 210UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        (void)tdiv(f, m, n, 10U);   // try the first 10 primes
        CHECK(tifa_is_complete_factorization(n, f, m),
              "tdiv: 210 fully factors over the first 10 primes");
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- ECM ---- small factor; b2 > 105 forces the second phase
    {
        mul2(n, 101UL, 1000000007UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        ecm_params_t params;
        params.b1      = 2000U;
        params.b2      = 50000U;   // > 105 -> second phase runs (fixed abs path)
        params.ncurves = 200U;
        params.seed    = 1U;
        ecode_t rc = ecm(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("ecm", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    mpz_clear(n);
    TIFA_TEST_END();
}
