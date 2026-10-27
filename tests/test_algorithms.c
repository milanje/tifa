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

// Per-algorithm correctness tests: each algorithm is driven directly through
// its public interface on a composite suited to it, asking for a COMPLETE
// factorization and verifying the result. The ECM case runs stage 2 (b2 > b1);
// where the Montgomery engine handles n, it does not reach the GMP stage-2
// precompute that held the abs() bug.

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

    // ---- CFRAC with the trial division filters ---- they dropped every
    // negative smooth residue, and early abort split its base off by one
    // (reading past it with a three-prime base).
    {
        static const smooth_filter_method_t methods[] = {
            TDIV, TDIV_EARLY_ABORT
        };
        for (unsigned i = 0; i < 2; i++) {
            mul2(n, 15485863UL, 32452843UL);
            mpz_array_t*    f = alloc_mpz_array(8);
            uint32_array_t* m = alloc_uint32_array(8);
            cfrac_params_t params;
            set_cfrac_params_to_default(n, &params);
            params.filter_method = methods[i];
            ecode_t rc = cfrac(f, m, n, &params, FIND_COMPLETE_FACTORIZATION);
            check("cfrac trial division filters", rc, n, f, m);
            free_mpz_array(f);
            free_uint32_array(m);
        }
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

    // ---- Proper factors only ---- SQUFOF on 3*9626^2 (k*n a square for a
    // racer: it divided by 0 and hung), ECM on 110 (n divides 4u^3v: it
    // returned 110 * 1), SQUFOF on the prime 103 (the race returned 103).
    {
        static const struct { const char* n; int algo; } cases[] = {
            { "277979628", 0 }, { "110", 1 }, { "103", 0 }
        };
        for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
            mpz_set_str(n, cases[i].n, 10);
            mpz_array_t* f = alloc_mpz_array(8);
            if (cases[i].algo == 0) {
                squfof_params_t params;
                set_squfof_params_to_default(&params);
                squfof(f, NULL, n, &params, FIND_SOME_FACTORS);
            } else {
                ecm_params_t params;
                set_ecm_params_to_default(n, &params);
                ecm(f, NULL, n, &params, SINGLE_RUN);
            }
            bool ok = true;
            for (uint32_t j = 0; j < f->length; j++) {
                ok = ok && (mpz_cmp_ui(f->data[j], 1) > 0)
                        && (mpz_cmp(f->data[j], n) < 0)
                        && mpz_divisible_p(n, f->data[j]);
            }
            CHECK(ok, cases[i].n);
            free_mpz_array(f);
        }
    }

    // ---- Fermat, one run, 70 bits ---- the square roots mod p^2 were wrong
    // for p > 65535, so a single run found almost nothing above 64 bits.
    {
        mpz_set_str(n, "722091009885066927899", 10);
        mpz_array_t* f = alloc_mpz_array(8);
        fermat_params_t params;
        set_fermat_params_to_default(&params);
        CHECK(fermat(f, NULL, n, &params, SINGLE_RUN) == SOME_FACTORS_FOUND,
              "fermat finds 23107307611 * 31249465409 in one run");
        free_mpz_array(f);
    }

    // ---- ECM, p+1 and p-1 on small composites ---- their stages caught
    // every prime factor at once (gcd n) and gave up, for 18% of the odd
    // composites near 10^5 (ECM): they now redo the stage one prime at a
    // time. p-1 also goes through the primes in decreasing order (1247 =
    // 29.43, 4187 = 53.79: orders of 2 ending in the same prime), and tries
    // another base (2047 = 23.89: 2 has order 11 modulo both).
    {
        static const unsigned long ecm_ns[] = { 20003, 20077, 100000051 };
        static const unsigned long pp1_ns[] = { 323, 391 };
        static const unsigned long pm1_ns[] = { 1247, 2047, 4187 };
        mpz_array_t* f = alloc_mpz_array(8);
        for (unsigned i = 0; i < sizeof(ecm_ns) / sizeof(ecm_ns[0]); i++) {
            mpz_set_ui(n, ecm_ns[i]);
            ecm_params_t params;
            set_ecm_params_to_default(n, &params);
            f->length = 0;
            CHECK(ecm(f, NULL, n, &params, FIND_SOME_FACTORS)
                      == SOME_FACTORS_FOUND, "ecm splits a small composite");
        }
        for (unsigned i = 0; i < sizeof(pp1_ns) / sizeof(pp1_ns[0]); i++) {
            mpz_set_ui(n, pp1_ns[i]);
            p_plus_1_params_t params;
            set_p_plus_1_params_to_default(n, &params);
            f->length = 0;
            CHECK(p_plus_1(f, NULL, n, &params, FIND_SOME_FACTORS)
                      == SOME_FACTORS_FOUND, "p+1 splits a small composite");
        }
        for (unsigned i = 0; i < sizeof(pm1_ns) / sizeof(pm1_ns[0]); i++) {
            mpz_set_ui(n, pm1_ns[i]);
            p_minus_1_params_t params;
            set_p_minus_1_params_to_default(n, &params);
            f->length = 0;
            CHECK(p_minus_1(f, NULL, n, &params, FIND_SOME_FACTORS)
                      == SOME_FACTORS_FOUND, "p-1 splits a small composite");
        }
        free_mpz_array(f);
    }

    // ---- CFRAC on short periods ---- k*n = m^2 + 1 makes every residue 1,
    // and 1005679*1005681 has a short period without smooth residues: both
    // used to loop forever. Any outcome will do, as long as it ends.
    {
        static const char* const ns[] = { "1487813", "1011392262399" };
        for (unsigned i = 0; i < 2; i++) {
            mpz_set_str(n, ns[i], 10);
            mpz_array_t* f = alloc_mpz_array(8);
            cfrac_params_t params;
            set_cfrac_params_to_default(n, &params);
            cfrac(f, NULL, n, &params, FIND_SOME_FACTORS);
            bool ok = true;
            for (uint32_t j = 0; j < f->length; j++) {
                ok = ok && mpz_divisible_p(n, f->data[j]);
            }
            CHECK(ok, "cfrac on a short period: factors divide n");
            free_mpz_array(f);
        }
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

    // ---- ECM ---- small factor; b2 > b1 runs the second phase
    {
        mul2(n, 101UL, 1000000007UL);
        mpz_array_t*    f = alloc_mpz_array(8);
        uint32_array_t* m = alloc_uint32_array(8);
        ecm_params_t params;
        params.b1      = 2000U;
        params.b2      = 50000U;   // > b1: second phase runs
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
