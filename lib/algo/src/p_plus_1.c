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

/**
 * \file    p_plus_1.c
 * \author  Jerome Milan
 * \date    2026
 * \version 1.0
 */

#include <stdlib.h>
#include <inttypes.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "funcs.h"
#include "macros.h"
#include "factoring_machine.h"
#include "p_plus_1.h"
#include "tifa_factor.h"
#include "tuning.h"
#include "mont_ecm.h"

#define __PREFIX__      "p+1: "
#define __VERBOSE__     TIFA_VERBOSE_PP1
#define __TIMING__      TIFA_TIMING_PP1

#include "messages.h"

//
// Maximum number of times the smoothness bounds are grown (after the whole
// seed list has been tried at the current bounds) before giving up.
//
#define MAX_BOUND_LEVELS 4
//
// Factor by which the bounds grow on each bound-level increment.
//
#define B1_GROWTH_FACTOR 2UL
//
// Default parameters.
//
#define DEFAULT_SEED     6
#define DEFAULT_B1       50000
//
// Default stage-2 bound, as a multiple of B1 (clamped to LARGEST_PRIME).
//
#define DEFAULT_B2_RATIO 20UL
//
// Baby-step/giant-step stage-2 giant step D = 2*3*5*7 and its half D/2 (mirrors
// the ECM phase-2 constants). Every prime q > 7 is coprime to D, so it is
// representable as m*D +- j with a baby residue j in [1, D/2] coprime to D.
//
#define PP1_D            210
#define PP1_DO2          105

//
// Seeds tried in turn. A given seed P (with discriminant D = P^2 - 4) only
// exposes a prime factor p through its (p+1)-order subgroup when (D/p) = -1
// (and through its (p-1)-order subgroup when (D/p) = +1), i.e. roughly half of
// the primes per seed. Cycling several seeds raises the odds of catching a
// smooth factor. All values satisfy |P| >= 3, so D is non-zero and (for
// integer seeds) never a perfect square.
//
static const uint32_t PP1_SEEDS[] = { 6, 5, 8, 11, 7, 9 };
#define NSEEDS  (sizeof(PP1_SEEDS) / sizeof(PP1_SEEDS[0]))

//
// An ad-hoc structure holding the data variables used by Williams's p+1
// implementation.
//
struct struct_p_plus_1_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // p+1 parameters.
    //
    p_plus_1_params_t* params;
    //
    // Number of times the context was updated.
    //
    int nupdates;
    //
    // Index into PP1_SEEDS used once the initial (params->seed) run is done.
    //
    int seed_idx;
    //
    // Number of times the bounds have been grown.
    //
    int bound_level;
    //
    // Working variables. x carries V_E from stage 1 into stage 2.
    //
    mpz_t x;
    mpz_t g;
#if TIFA_HAVE_MONT_ECM
    //
    // Plan of the fixed-limb Montgomery implementation for the bounds
    // (mplan_b1, mplan_b2), used when n has at most MONT_ECM_MAX_LIMBS limbs.
    //
    mont_ecm_plan_t* mplan;
    uint32_t         mplan_b1;
    uint32_t         mplan_b2;
#endif
};
typedef struct struct_p_plus_1_context_t p_plus_1_context_t;

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_p_plus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_p_plus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_p_plus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_p_plus_1(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------
static void lucas_V_mod(
    mpz_ptr out,
    unsigned long m,
    const mpz_t seed,
    const mpz_t n
);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t p_plus_1(mpz_array_t* const factors, uint32_array_t* const multis,
                 const mpz_t n, const p_plus_1_params_t* const params,
                 const factoring_mode_t mode) {

    PRINT_INTRO_MSG("Williams's p+1");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_p_plus_1_context;
    machine.perform_algo_func   = perform_p_plus_1;
    machine.update_context_func = update_p_plus_1_context;
    machine.clear_context_func  = clear_p_plus_1_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------
void set_p_plus_1_params_to_default(const mpz_t n,
                                    p_plus_1_params_t* const params) {
    params->seed = DEFAULT_SEED;
    params->b1   = DEFAULT_B1;
    if (params->b1 > LARGEST_PRIME) {
        params->b1 = LARGEST_PRIME;
    }
    unsigned long b2 = DEFAULT_B2_RATIO * (unsigned long) params->b1;
    if (b2 > LARGEST_PRIME) {
        b2 = LARGEST_PRIME;
    }
    params->b2 = (uint32_t) b2;

    //
    // Override with machine-tuned values if a tuning file supplies a row for
    // this bit-size; otherwise the defaults above are kept as-is.
    //
    tuning_row_t row;
    if (tifa_tuning_get_row(
        "p_plus_1", TIFA_TUNING_NCOLS_P_PLUS_1,
        (uint32_t) mpz_sizeinbase(n, 2),
        &row)) {
        params->b1   = row.cols[0];
        params->b2   = row.cols[1];
        params->seed = row.cols[2];
        if (params->b1 > LARGEST_PRIME) {
            params->b1 = LARGEST_PRIME;
        }
        if (params->b2 > LARGEST_PRIME) {
            params->b2 = LARGEST_PRIME;
        }
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              Lucas sequence V_n(P, 1) modular arithmetic
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static void lucas_V_mod(mpz_ptr out, unsigned long m, const mpz_t seed,
                        const mpz_t n) {
    //
    // Computes out = V_m(P, 1) mod n, where P = seed and V is the Lucas
    // sequence V_0 = 2, V_1 = P, V_k = P*V_{k-1} - V_{k-2}. Assumes m >= 1.
    //
    // A Montgomery binary ladder is used: the pair (Vh, Vh1) = (V_h, V_{h+1})
    // is maintained and, for each bit of m from the second most significant
    // down to the least significant, h is mapped to 2h + bit via the
    // doubling/addition formulas (Q = 1):
    //     V_{2h}   = V_h^2 - 2
    //     V_{2h+1} = V_h*V_{h+1} - P
    //     V_{2h+2} = V_{h+1}^2 - 2
    // costing exactly two modular multiplications per bit.
    //
    mpz_t Vh, Vh1, t, P;
    mpz_init(Vh);
    mpz_init(Vh1);
    mpz_init(t);
    mpz_init(P);

    //
    // Private, fully reduced copy of the seed (used in the odd-step
    // subtraction), so that out may safely alias seed.
    //
    mpz_mod(P, seed, n);

    //
    // h = 1: Vh = V_1 = P, Vh1 = V_2 = P^2 - 2.
    //
    mpz_set(Vh, P);
    mpz_mul(Vh1, P, P);
    mpz_sub_ui(Vh1, Vh1, 2);
    mpz_mod(Vh1, Vh1, n);

    //
    // Index of the most significant bit of m (the MSB is consumed by the
    // h = 1 initialization above).
    //
    int top = -1;
    for (unsigned long mm = m; mm != 0; mm >>= 1) {
        top++;
    }

    for (int s = top - 1; s >= 0; s--) {
        if (((m >> s) & 1UL) == 0) {
            //
            // bit 0: (V_h, V_{h+1}) -> (V_{2h}, V_{2h+1}).
            //
            mpz_mul(t, Vh, Vh1);        // V_h * V_{h+1}
            mpz_sub(t, t, P);           //   - P            = V_{2h+1}
            mpz_mod(t, t, n);
            mpz_mul(Vh, Vh, Vh);        // V_h^2
            mpz_sub_ui(Vh, Vh, 2);      //   - 2            = V_{2h}
            mpz_mod(Vh, Vh, n);
            mpz_swap(Vh1, t);           // Vh1 <- V_{2h+1}
        } else {
            //
            // bit 1: (V_h, V_{h+1}) -> (V_{2h+1}, V_{2h+2}).
            //
            mpz_mul(t, Vh1, Vh1);       // V_{h+1}^2
            mpz_sub_ui(t, t, 2);        //   - 2            = V_{2h+2}
            mpz_mod(t, t, n);
            mpz_mul(Vh, Vh, Vh1);       // V_h * V_{h+1}
            mpz_sub(Vh, Vh, P);         //   - P            = V_{2h+1}
            mpz_mod(Vh, Vh, n);
            mpz_swap(Vh1, t);           // Vh1 <- V_{2h+2}
        }
    }

    mpz_set(out, Vh);

    mpz_clear(Vh);
    mpz_clear(Vh1);
    mpz_clear(t);
    mpz_clear(P);
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static ecode_t init_p_plus_1_context(factoring_machine_t* const machine) {
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context = (void*) malloc(sizeof(p_plus_1_context_t));
    p_plus_1_context_t* context = (p_plus_1_context_t*) machine->context;

    mpz_init_set(context->n, machine->n);
    mpz_init_set_ui(context->x, 1);
    mpz_init_set_ui(context->g, 1);

    context->params      = (p_plus_1_params_t*) machine->params;
    context->nupdates    = 0;
    context->seed_idx    = 0;
    context->bound_level = 0;
#if TIFA_HAVE_MONT_ECM
    context->mplan    = NULL;
    context->mplan_b1 = 0;
    context->mplan_b2 = 0;
#endif

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_p_plus_1_context(factoring_machine_t* const machine) {
    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    p_plus_1_context_t* context = (p_plus_1_context_t*) machine->context;

    if (context != NULL) {
#if TIFA_HAVE_MONT_ECM
        free_mont_ecm_plan(context->mplan);
#endif
        mpz_clear(context->n);
        mpz_clear(context->x);
        mpz_clear(context->g);
        free(context);
    }

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_p_plus_1_context(factoring_machine_t* const machine) {
    INIT_TIMER;
    START_TIMER;

    p_plus_1_context_t* context = (p_plus_1_context_t*) machine->context;
    p_plus_1_params_t*  params  = context->params;
    ecode_t ecode = SUCCESS;

    context->nupdates++;

    //
    // Cheap first: retry with the next seed at the same bounds (a full run at
    // the same cost, and the correct remedy for an overshoot). Only once the
    // whole seed list is exhausted do we (more expensively) grow the bounds.
    //
    if (context->seed_idx + 1 < (int) NSEEDS) {
        context->seed_idx++;
    } else {
        context->seed_idx = 0;
        context->bound_level++;

        if ((context->bound_level > MAX_BOUND_LEVELS)
            || (params->b1 >= LARGEST_PRIME)) {
            ecode = GIVING_UP;
            PRINT_UPDATE_GIVEUP_MSG;
        } else {
            unsigned long nb1 = (unsigned long) params->b1 * B1_GROWTH_FACTOR;
            if (nb1 > LARGEST_PRIME) {
                nb1 = LARGEST_PRIME;
            }
            params->b1 = (uint32_t) nb1;

            unsigned long nb2 = (unsigned long) params->b2 * B1_GROWTH_FACTOR;
            if (nb2 > LARGEST_PRIME) {
                nb2 = LARGEST_PRIME;
            }
            params->b2 = (uint32_t) nb2;
        }
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t perform_p_plus_1(factoring_machine_t* const machine) {
    p_plus_1_context_t* context = (p_plus_1_context_t*) machine->context;
    p_plus_1_params_t*  params  = (p_plus_1_params_t*)  machine->params;

    INIT_TIMER;
    START_TIMER;

    mpz_ptr n = context->n;
    mpz_ptr x = context->x;
    mpz_ptr g = context->g;

    ecode_t ecode = NO_FACTOR_FOUND;

    //
    // Effective seed for this run: the initial run uses params->seed, later
    // runs cycle through PP1_SEEDS (see update_p_plus_1_context).
    //
    unsigned long P = (context->nupdates == 0)
                        ? (unsigned long) params->seed
                        : (unsigned long) PP1_SEEDS[context->seed_idx];

    //
    // Cheap early check: a common factor between the seed and n is itself a
    // factor of n.
    //
    mpz_set_ui(g, P);
    mpz_gcd(g, g, n);
    if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(machine->factors, g);

        STOP_TIMER;
        PRINT_TIMING;
        return SOME_FACTORS_FOUND;
    }

    //
    // Cheap early check: a prime p dividing the discriminant D = P^2 - 4
    // degenerates the sequence mod p, but gcd(D, n) then splits n for free.
    //
    mpz_set_ui(g, P);
    mpz_mul(g, g, g);
    mpz_sub_ui(g, g, 4);
    mpz_gcd(g, g, n);
    if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(machine->factors, g);

        STOP_TIMER;
        PRINT_TIMING;
        return SOME_FACTORS_FOUND;
    }

#if TIFA_HAVE_MONT_ECM
    //
    // Fixed-limb Montgomery arithmetic when n has at most 4 limbs. When a
    // gcd is n, the code below redoes the run.
    //
    if (mont_ecm_handles(n) && (mpz_cmp_ui(n, P) > 0)) {
        if (   (context->mplan == NULL)
            || (context->mplan_b1 != params->b1)
            || (context->mplan_b2 != params->b2)) {
            free_mont_ecm_plan(context->mplan);
            context->mplan    = alloc_mont_ecm_plan(params->b1, params->b2);
            context->mplan_b1 = params->b1;
            context->mplan_b2 = params->b2;
        }
        const int found = mont_pp1(g, n, (uint32_t) P, context->mplan);
        if (found != -1) {
            if (found == 1) {
                append_mpz_to_array(machine->factors, g);
                mpz_divexact(g, n, g);
                append_mpz_to_array(machine->factors, g);
                ecode = SOME_FACTORS_FOUND;
            }
            STOP_TIMER;
            PRINT_TIMING;
            return ecode;
        }
    }
#endif

    //
    // Stage 1: compute x <- V_E(P, 1) mod n, where
    //     E = prod_{p prime <= b1} p^floor(log_p(b1)),
    // accumulated one prime power at a time. This relies on the composition
    // identity V_{m*k}(P, 1) = V_m(V_k(P, 1), 1) (valid because Q = 1): each
    // step reseeds the ladder with the current x and raises it to the next
    // prime power.
    //
    mpz_set_ui(x, P);

    uint32_t ip = 0;
    while ((ip < NFIRST_PRIMES) && (first_primes[ip] <= params->b1)) {
        unsigned long pr = first_primes[ip];
        unsigned long pe = pr;
        //
        // pe <- largest power of pr that is <= b1 (no overflow: pe*pr <= b1).
        //
        while (pe <= params->b1 / pr) {
            pe *= pr;
        }
        lucas_V_mod(x, pe, x, n);
        ip++;
    }

    //
    // Stage 1 result: if p | n has (D/p) = -1 and (p+1) | E, then V_E = 2
    // (mod p), so g = gcd(V_E - 2, n) reveals p. (When (D/p) = +1 the same
    // holds if (p-1) | E.)
    //
    mpz_sub_ui(g, x, 2);
    mpz_gcd(g, g, n);

    if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(machine->factors, g);

        STOP_TIMER;
        PRINT_TIMING;
        return SOME_FACTORS_FOUND;
    }

    //
    // Stage 2: if stage 1 did not overshoot (g != n) and a second bound is
    // set, look for a factor p such that p + 1 is a b1-smooth number times a
    // single prime q in (b1, b2]. Working in the reseeded sequence
    //     W_i := V_i(x, 1) = beta^i + beta^-i (mod p),   beta = alpha^E,
    // (x = V_E from stage 1), such a factor is exposed when beta^q = 1 (mod p),
    // i.e. W_q = 2 (mod p). We use a baby-step/giant-step continuation: every
    // prime q > 7 is written as q = m*D +- j with giant step D and baby residue
    // j in [1, D/2] coprime to D. Because
    //     R_m - SV[j] = (beta^mD + beta^-mD) - (beta^j + beta^-j) = 0 (mod p)
    // whenever beta^(m*D +- j) = 1, the single term R_m - SV[j] (with giant term
    // R_m = W_{m*D} and baby term SV[j] = W_j) catches both q = m*D + j and
    // q = m*D - j. ip already indexes the first prime > b1 (left by stage 1).
    //
    if ((mpz_cmp(g, n) != 0) && (params->b2 > params->b1)
        && (ip < NFIRST_PRIMES) && (first_primes[ip] <= params->b2)
        && (first_primes[ip] != 2)) {

        //
        // Baby residues: the j in [1, PP1_DO2] coprime to D (24 of them for
        // D = 210). Js lists them, njs counts them.
        //
        unsigned int Js[PP1_DO2 + 1];
        unsigned int njs = 0;
        for (unsigned long j = 1; j <= PP1_DO2; j++) {
            if (gcd_ulint(j, PP1_D) == 1) {
                Js[njs++] = (unsigned int) j;
            }
        }

        //
        // Baby table SV[j] = W_j = V_j(x, 1), stored only for coprime j. Built
        // from the Lucas recurrence W_0 = 2, W_1 = x, W_j = x*W_{j-1} - W_{j-2}
        // (mod n), carried across the non-coprime gaps by two rolling terms.
        //
        mpz_t* SV = (mpz_t*) malloc((PP1_DO2 + 1) * sizeof(mpz_t));
        for (unsigned int ij = 0; ij < njs; ij++) {
            mpz_init(SV[Js[ij]]);
        }

        mpz_t wprev, wcur, wtmp;
        mpz_init_set_ui(wprev, 2);       // W_0
        mpz_init_set(wcur, x);           // W_1
        mpz_init(wtmp);
        mpz_set(SV[1], wcur);            // 1 is coprime to D
        for (unsigned long j = 2; j <= PP1_DO2; j++) {
            mpz_mul(wtmp, x, wcur);      // x*W_{j-1}
            mpz_sub(wtmp, wtmp, wprev);  //  - W_{j-2}  = W_j
            mpz_mod(wtmp, wtmp, n);
            mpz_set(wprev, wcur);
            mpz_set(wcur, wtmp);
            if (gcd_ulint(j, PP1_D) == 1) {
                mpz_set(SV[j], wcur);
            }
        }
        mpz_clear(wprev);
        mpz_clear(wcur);
        mpz_clear(wtmp);

        //
        // Giant range: m in [m_min, m_max] brackets the primes in (b1, b2] via
        // m*D +- j (j <= PP1_DO2).
        //
        unsigned long m_min = ((unsigned long) params->b1 + PP1_DO2) / PP1_D;
        unsigned long m_max = ((unsigned long) params->b2 - PP1_DO2) / PP1_D + 1;
        if (m_min < 1) {
            m_min = 1;
        }

        if (m_max >= m_min) {
            unsigned long nrows = m_max - m_min + 1;

            //
            // Prime mask: mask[(m - m_min)*(PP1_DO2 + 1) + j] = 1 iff m*D +- j
            // is a prime in (b1, b2]. Built in a single pass over the tabulated
            // primes q in (b1, b2] (q = m*D + r, j = |r|); no primality tests.
            //
            unsigned char* mask = (unsigned char*)
                calloc(nrows * (PP1_DO2 + 1), sizeof(unsigned char));

            for (uint32_t i = ip;
                 (i < NFIRST_PRIMES) && (first_primes[i] <= params->b2); i++) {
                unsigned long q = first_primes[i];
                unsigned long m = (q + PP1_DO2) / PP1_D;
                long r = (long) q - (long) (m * PP1_D);
                unsigned long j = (r >= 0) ? (unsigned long) r
                                           : (unsigned long) (-r);
                if ((m >= m_min) && (m <= m_max) && (j >= 1) && (j <= PP1_DO2)) {
                    mask[(m - m_min) * (PP1_DO2 + 1) + j] = 1;
                }
            }

            //
            // Giant step G = W_D and the running giant terms R = W_{m_min*D},
            // Rprev = W_{(m_min-1)*D}. W_0 = 2 bootstraps the recurrence when
            // m_min == 1 (no point-at-infinity special case, unlike ECM).
            //
            mpz_t G, R, Rprev, Rnext, t, H;
            mpz_init(G);
            mpz_init(R);
            mpz_init(Rprev);
            mpz_init(Rnext);
            mpz_init(t);
            mpz_init_set_ui(H, 1);

            lucas_V_mod(G, PP1_D, x, n);
            lucas_V_mod(R, m_min * PP1_D, x, n);
            if (m_min == 1) {
                mpz_set_ui(Rprev, 2);    // W_0
            } else {
                lucas_V_mod(Rprev, (m_min - 1) * PP1_D, x, n);
            }

            //
            // Giant walk: accumulate H = prod (R_m - SV[j]) over the prime
            // (m, j) slots, advancing R by the three-term Lucas recurrence
            // W_{(m+1)D} = W_{mD}*W_D - W_{(m-1)D}.
            //
            for (unsigned long im = 0; im < nrows; im++) {
                unsigned char* row = &mask[im * (PP1_DO2 + 1)];
                for (unsigned int ij = 0; ij < njs; ij++) {
                    if (row[Js[ij]]) {
                        mpz_sub(t, R, SV[Js[ij]]);
                        mpz_mul(H, H, t);
                        mpz_mod(H, H, n);
                    }
                }
                if (im + 1 < nrows) {
                    mpz_mul(Rnext, R, G);
                    mpz_sub(Rnext, Rnext, Rprev);
                    mpz_mod(Rnext, Rnext, n);
                    mpz_set(Rprev, R);
                    mpz_set(R, Rnext);
                }
            }

            mpz_gcd(g, H, n);
            if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
                append_mpz_to_array(machine->factors, g);
                mpz_divexact(g, n, g);
                append_mpz_to_array(machine->factors, g);
                ecode = SOME_FACTORS_FOUND;
            } else if (mpz_cmp(g, n) == 0) {
                //
                // The whole product is a multiple of n: >= 2 factors were
                // caught at once. Re-walk the giant steps with a per-term gcd
                // to isolate a single proper factor.
                //
                lucas_V_mod(R, m_min * PP1_D, x, n);
                if (m_min == 1) {
                    mpz_set_ui(Rprev, 2);
                } else {
                    lucas_V_mod(Rprev, (m_min - 1) * PP1_D, x, n);
                }
                bool done = false;
                for (unsigned long im = 0; (im < nrows) && !done; im++) {
                    unsigned char* row = &mask[im * (PP1_DO2 + 1)];
                    for (unsigned int ij = 0; ij < njs; ij++) {
                        if (row[Js[ij]]) {
                            mpz_sub(t, R, SV[Js[ij]]);
                            mpz_gcd(g, t, n);
                            if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
                                append_mpz_to_array(machine->factors, g);
                                mpz_divexact(g, n, g);
                                append_mpz_to_array(machine->factors, g);
                                ecode = SOME_FACTORS_FOUND;
                                done = true;
                                break;
                            }
                        }
                    }
                    if (!done && (im + 1 < nrows)) {
                        mpz_mul(Rnext, R, G);
                        mpz_sub(Rnext, Rnext, Rprev);
                        mpz_mod(Rnext, Rnext, n);
                        mpz_set(Rprev, R);
                        mpz_set(R, Rnext);
                    }
                }
            }

            mpz_clear(G);
            mpz_clear(R);
            mpz_clear(Rprev);
            mpz_clear(Rnext);
            mpz_clear(t);
            mpz_clear(H);
            free(mask);
        }

        for (unsigned int ij = 0; ij < njs; ij++) {
            mpz_clear(SV[Js[ij]]);
        }
        free(SV);
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor(factors, multis, n, mode);

    //
    // The following code should be used to factor a number using _only_
    // Williams's p+1.
    //
    //p_plus_1_params_t params;
    //set_p_plus_1_params_to_default(n, &params);
    //return p_plus_1(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------
