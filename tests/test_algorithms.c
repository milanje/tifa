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
#include "p_plus_1.h"
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

    // ---- CFRAC without the large-prime variation ---- regression test for a
    // NULL-hashtable dereference: the filter backends read filter->htable
    // unconditionally, but cfrac only allocates it when use_large_primes is
    // true, so use_large_primes=false used to segfault. A larger factor base
    // compensates for not recovering partial (one-large-prime) relations.
    {
        mul2(n, 15485863UL, 32452843UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        cfrac_params_t params;
        set_cfrac_params_to_default(n, &params);
        params.use_large_primes = false;
        params.nprimes_in_base  = 1024;
        params.nprimes_tdiv     = 1024;
        ecode_t rc = cfrac(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("cfrac no-large-primes", rc, n, f, m);
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

    // ---- Fermat with multipliers ---- q ~ 5p, 7p, 11p: the split is found on
    // k*n, and must be taken modulo n, not k*n.
    {
        const char* const inputs[] = {
            "22188730777", "12523140803", "1758490455136543", "2445747370601177"
        };
        for (uint32_t i = 0; i < 4; i++) {
            mpz_set_str(n, inputs[i], 10);
            mpz_array_t*    f = alloc_mpz_array(8);
            uint32_array_t* m = alloc_uint32_array(8);
            fermat_params_t params;
            set_fermat_params_to_default(&params);
            ecode_t rc = fermat(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
            check("fermat multiplier", rc, n, f, m);
            free_mpz_array(f);
            free_uint32_array(m);
        }
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

    // ---- Pollard p-1 backtracking ---- p - 1 and q - 1 are both 1000-smooth,
    // so stage 1 gives gcd = n; redoing it prime by prime splits n at 547.
    {
        mpz_set_str(n, "385402726458636771254457161", 10);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_minus_1_params_t params;
        set_p_minus_1_params_to_default(n, &params);
        params.b1   = 1000;
        params.b2   = 0;
        params.base = 2;
        ecode_t rc = p_minus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_minus_1 backtrack", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Pollard p-1 stage 2 ---- (10007 = 2*5003 + 1: b1 < 5003 < b2,
    // so stage 1 misses and stage 2 recovers the factor)
    {
        mul2(n, 10007UL, 1000003UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_minus_1_params_t params;
        set_p_minus_1_params_to_default(n, &params);
        params.b1 = 1000;    // below 5003: stage 1 misses
        params.b2 = 10000;   // above 5003: stage 2 catches it
        ecode_t rc = p_minus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_minus_1 stage2", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Pollard p-1 stage 2, batch rescan ---- two safe primes (10007,
    // 10079) whose (p-1)/2 primes 5003 and 5039 are adjacent and fall in the
    // same stage-2 batch, so the batch product is a multiple of n and the
    // per-prime rescan must isolate a factor.
    {
        mul2(n, 10007UL, 10079UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_minus_1_params_t params;
        set_p_minus_1_params_to_default(n, &params);
        params.b1 = 1000;
        params.b2 = 10000;
        ecode_t rc = p_minus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_minus_1 rescan", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Williams p+1 ---- (179 + 1 = 2^2*3^2*5 is smooth; (D/179) = -1 for
    // the default seed 6, so 179 is exposed through its p+1 group. b1 = 50 keeps
    // the 1000003 cofactor -- whose p+-1 are not 50-smooth -- out of stage 1.)
    {
        mul2(n, 179UL, 1000003UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_plus_1_params_t params;
        set_p_plus_1_params_to_default(n, &params);
        params.b1 = 50;
        params.b2 = 0;       // stage 1 only
        ecode_t rc = p_plus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_plus_1", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Williams p+1 stage 2 ---- (4051 + 1 = 2^2*1013: b1 < 1013 < b2, so
    // stage 1 misses and the large-prime continuation recovers the factor.)
    {
        mul2(n, 4051UL, 1000003UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_plus_1_params_t params;
        set_p_plus_1_params_to_default(n, &params);
        params.b1 = 1000;    // below 1013: stage 1 misses
        params.b2 = 10000;   // above 1013: stage 2 catches it
        ecode_t rc = p_plus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_plus_1 stage2", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Williams p+1 stage 2, batch rescan ---- two primes 4051, 4243 whose
    // p+1 large primes 1013 and 1061 fall in the same stage-2 batch, so the
    // batch product is a multiple of n and the per-prime rescan must isolate a
    // factor.
    {
        mul2(n, 4051UL, 4243UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_plus_1_params_t params;
        set_p_plus_1_params_to_default(n, &params);
        params.b1 = 1000;
        params.b2 = 10000;
        ecode_t rc = p_plus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_plus_1 rescan", rc, n, f, m);
        free_mpz_array(f);
        free_uint32_array(m);
    }

    // ---- Williams p+1 stage 2, many giant steps ---- 32323 + 1 = 2^2*8081,
    // and 8081 = 38*210 + 101, so the BSGS giant walk must advance to m = 38
    // (33 steps past m_min = 5) to catch it, exercising the giant-step
    // recurrence and a near-top baby residue (j = 101).
    {
        mul2(n, 32323UL, 1000003UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        p_plus_1_params_t params;
        set_p_plus_1_params_to_default(n, &params);
        params.b1 = 1000;
        params.b2 = 10000;
        ecode_t rc = p_plus_1(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
        check("p_plus_1 bsgs", rc, n, f, m);
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

    // ---- ECM, default bounds ---- a 40-bit factor of a 189-bit n: the
    // bounds start at the level for 20-bit factors and move up until found.
    {
        mpz_set_str(n, "665996203062451295188310324867037030192120874180607253493",
                    10);
        mpz_array_t* f = alloc_mpz_array(8);
        ecm_params_t params;
        set_ecm_params_to_default(n, &params);
        ecode_t rc = ecm(f, NULL, n, &params, FIND_SOME_FACTORS);
        CHECK(rc == SOME_FACTORS_FOUND, "ecm defaults: rc == SOME_FACTORS_FOUND");
        int found = 0;
        for (uint32_t i = 0; i < f->length; i++) {
            if (mpz_cmp_ui(f->data[i], 817225028183UL) == 0) {
                found = 1;
            }
        }
        CHECK(found, "ecm defaults: the 40-bit factor is found");
        free_mpz_array(f);
    }

    mpz_clear(n);
    TIFA_TEST_END();
}
