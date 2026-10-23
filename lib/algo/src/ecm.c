// Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
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
 * \file    ecm.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "z2z8_curves.h"
#include "array.h"
#include "funcs.h"
#include "x_tree.h"
#include "macros.h"
#include "factoring_machine.h"
#include "ecm.h"
#include "mont_ecm.h"
#include "tifa_factor.h"
#include "tuning.h"

#define __PREFIX__   "ecm: "
#define __VERBOSE__  TIFA_VERBOSE_ECM
#define __TIMING__   TIFA_TIMING_ECM

#include "messages.h"
#include "tifa_alloc.h"

//------------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//------------------------------------------------------------------------------

//
// If PHASE_2_GLOBALS is set to non-zero, limits the number of mpz
// initializations and clearings by using declaring dedicated temporary
// variables in ecm_context_t. In pratice, this wins... only about 1%...
//
#define PHASE_2_GLOBALS 1

//
// An ad-hoc structure holding the data variables used by the ECM
// implementation.
//
struct struct_ecm_context_t {
    //
    // The number to factor
    //
    mpz_t n;
    //
    // k <- Prod_{pi<=B1} pi^{log(B1)/log(pi)}
    //
    mpz_t k;
    //
    // (x0 : z0) is the starting point on the curve using the Montgomery
    // representation.
    //
    mpz_t x0;
    mpz_t z0;
    //
    // (xk : zk) is the point obtained at the en d of the first phase.
    // (xk : zk) = [k] (x0 : z0) (uses the Montgomery representation).
    //
    mpz_t xk;
    mpz_t zk;
    //
    // Parameters of a Suyama curve: b.z.y^2 = x^3 + a.z.x^2 + x.z^2
    //
    mpz_t a;
    mpz_t b;
    //
    // Pointer to the found factors. Could be accessed differently, but
    // it's convenient to have a copy of the pointer here...
    //
    mpz_array_t* factors;
    //
    // Suyama's curve parameter
    //
    unsigned long int sigma;
    //
    // State of the generator choosing the curves (xorshift64*).
    //
    uint64_t rng;
    //
    // Curves run so far on fixed-limb arithmetic: the first ones come from
    // the Z/2 x Z/8 table (see mont_ecm_nth_curve).
    //
    uint32_t nmont_curves;
    //
    // Index of first "large prime" i.e. smallest prime greater than B1
    //
    uint32_t  first_lp_index;
    //
    // Temporary variable used in computations...
    //
    mpz_t tmp_add_1;
    mpz_t tmp_add_2;
    mpz_t tmp_add_3;
    mpz_t tmp_add_4;
    mpz_t tmp_add_5;
    mpz_t tmp_add_6;
    mpz_t tmp_mul_1;
    mpz_t tmp_mul_2;
    mpz_t tmp_mul_3;
    mpz_t tmp_mul_4;

#if PHASE_2_GLOBALS
    mpz_t tmp_phase2_1;
    mpz_t tmp_phase2_2;
    mpz_t tmp_phase2_3;
    mpz_t tmp_phase2_4;
    mpz_t tmp_phase2_5;
    mpz_t tmp_phase2_6;
    mpz_t tmp_phase2_7;
    mpz_t tmp_phase2_8;
    mpz_t tmp_phase2_9;
    mpz_t tmp_phase2_10;
    mpz_t tmp_phase2_11;
#endif

    //
    // A pointer to the ECM parameters. Could be accessed differently, but
    // it's convenient to have a copy of the pointer here...
    //
    ecm_params_t*  params;
    //
    // true if the second phase precomputations have been performed.
    //
    bool precomputations_done;
    //
    // Bounds used in phase 1 and 2.
    //
    unsigned long int B1;
    unsigned long int B2;
    //
    // Used in precomputations. See reference article for details.
    //
    unsigned long int D;
    unsigned long int DO2;
    unsigned long int m_min;
    unsigned long int m_max;
    unsigned int*  Js;
    unsigned int   njs;
    unsigned int*  gcd_table;
    unsigned int** prime_table;
    mpz_t* SX;
    mpz_t* SZ;
    //
    // Next level of ecm_levels to escalate to, last level worth trying (for
    // factors of half the size of n), and rounds run at that last level.
    //
    uint32_t next_level;
    uint32_t top_level;
    uint32_t rounds_at_top;
#if TIFA_HAVE_MONT_ECM
    //
    // Fixed-limb Montgomery curves for n of at most MONT_ECM_MAX_LIMBS limbs,
    // with their plan for the bounds (mplan_b1, mplan_b2).
    //
    bool             use_mont;
    mont_ecm_plan_t* mplan;
    uint32_t         mplan_b1;
    uint32_t         mplan_b2;
#endif
};
//------------------------------------------------------------------------------
typedef struct struct_ecm_context_t ecm_context_t;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                      CONSTANT DATA, SYMBOLS AND MACROS
//------------------------------------------------------------------------------
//
// Bounds by size of the smallest factor sought (fbits): the B1 minimizing
// the expected time to find such a factor (mean curves times time per
// curve, measured with mont_ecm_nth_curve on 120-bit numbers and B2 =
// ECM_B2_RATIO.B1) and the mean number of curves that takes. B1 cannot
// exceed the largest tabulated prime (stage 1 takes its primes from the
// table); above 60 bits, the bounds and curve counts are extrapolated.
//
typedef struct {
    uint32_t fbits;
    uint32_t b1;
    uint32_t ncurves;
} ecm_level_t;

static const ecm_level_t ecm_levels[] = {
    { 15,     50,  1 },
    { 20,     70,  2 },
    { 25,    100,  3 },
    { 30,    200,  5 },
    { 35,    560,  6 },
    { 40,    800, 13 },
    { 45,   1600, 16 },
    { 50,   4480, 15 },
    { 55,   8960, 20 },
    { 60,  17920, 18 },
    { 65,  36000, 28 },
    { 70,  72000, 32 },
    { 75, 144000, 40 },
    { 80, 288000, 50 },
    { 85, 576000, 64 }
};
#define NLEVELS (sizeof(ecm_levels) / sizeof(ecm_levels[0]))

#define ECM_B2_RATIO 50
//
// Rounds of curves run again at the last level (the one for factors of half
// the size of n) before giving up.
//
#define MAX_ROUNDS_AT_TOP 4
//
// Largest stage 2 bound, above the top level's: the GMP stage 2 table takes
// about two bytes per unit of b2.
//
#define MAX_B2 (1U << 26)

//------------------------------------------------------------------------------
//                    PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_ecm_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_ecm_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_ecm_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static void switch_to_level(factoring_machine_t* const machine,
                            const uint32_t level);
//------------------------------------------------------------------------------
static ecode_t perform_ecm(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t phase_1(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t phase_2(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t choose_curve(ecm_context_t* const context);
//------------------------------------------------------------------------------
static ecode_t perform_precomputations(ecm_context_t* const context);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------
static void multiply_point(
    mpz_t xk, mpz_t zk,  const mpz_t x0, const mpz_t z0,
    const mpz_t k, ecm_context_t* const context
);
//------------------------------------------------------------------------------
static void double_point(
    mpz_t x2a, mpz_t z2a, const mpz_t xa, const mpz_t za,
    ecm_context_t* const context
);
//------------------------------------------------------------------------------
static void add_points(
    mpz_t xc, mpz_t zc,
    const mpz_t xa, const mpz_t za,
    const mpz_t xb, const mpz_t zb,
    const mpz_t x0, const mpz_t z0,
    ecm_context_t* const context
);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t ecm(mpz_array_t* const factors, uint32_array_t* const multis,
            const mpz_t n, const ecm_params_t* const params,
            const factoring_mode_t mode) {

    PRINT_INTRO_MSG("elliptic curve method");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    //
    // A copy: the algorithm adjusts it, and the caller's may be shared.
    //
    ecm_params_t own_params = *params;
    own_params.b2 = MIN(own_params.b2, MAX_B2);

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) &own_params;
    machine.init_context_func   = init_ecm_context;
    machine.perform_algo_func   = perform_ecm;
    machine.update_context_func = update_ecm_context;
    machine.clear_context_func  = clear_ecm_context;
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
void set_ecm_params_to_default(const mpz_t n, ecm_params_t* const params) {
    //
    // Start with the bounds for the smallest factors: ecm() then moves up
    // the levels of ecm_levels as long as no factor is found.
    //
    unsigned int sizen = mpz_sizeinbase(n, 2);

    set_ecm_params_for_factor_size(0, params);

    //
    // Override the tabulated defaults with machine-tuned values if a tuning
    // file supplies a row for this bit-size. The RNG seed is not tuned (it is
    // set internally in init_ecm_context), so it is left untouched here.
    //
    tuning_row_t row;
    if (tifa_tuning_get_row(
        "ecm", TIFA_TUNING_NCOLS_ECM, (uint32_t) sizen,
        TIFA_TUNING_INTERP_ECM, &row)) {
        params->b1      = row.cols[0];
        params->b2      = row.cols[1];
        params->ncurves = row.cols[2];
    }
    return;
}
//------------------------------------------------------------------------------
static unsigned long int next_sigma(ecm_context_t* const context) {
    //
    // xorshift64*, top 24 bits (Suyama's sigma must be at least 6).
    //
    uint64_t x = context->rng;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    context->rng = x;
    return (unsigned long int) ((x * 2685821657736338717ULL) >> 40) + 6;
}
//------------------------------------------------------------------------------
void set_ecm_params_for_factor_size(const uint32_t fbits,
                                    ecm_params_t* const params) {
    uint32_t i = 0;
    while ((i + 1 < NLEVELS) && (ecm_levels[i].fbits < fbits)) {
        i++;
    }
    params->b1      = ecm_levels[i].b1;
    params->b2      = ECM_B2_RATIO * ecm_levels[i].b1;
    params->ncurves = ecm_levels[i].ncurves;
    params->seed    = 0;
}
//------------------------------------------------------------------------------
static void set_stage1_bound(ecm_context_t* const context,
                             ecm_params_t* const params) {
    if (params->b1 > LARGEST_PRIME) {
        //
        // Stage 1 takes its primes from the precomputed table; stage 2 then
        // starts right after it.
        //
        params->b1 = LARGEST_PRIME;
    }
#if TIFA_HAVE_MONT_ECM
    if (context->use_mont) {
        //
        // The Montgomery curves get their prime powers from their plan.
        //
        return;
    }
#endif
    //
    // Stage 2 takes giant steps of D = 210 from b1 + D/2: as on the
    // Montgomery curves, it needs b1 >= D/2, and b2 > b1.
    //
    if (params->b1 < 2) {
        params->b1 = 2;
    }
    if ((params->b2 > params->b1) && (params->b1 < 105)) {
        params->b1 = 105;
    }
    if (params->b2 <= params->b1) {
        params->b2 = 0;
    }
    //
    // Compute context->k = Prod_{pi<=B1} pi^{log(B1)/log(pi)} via a product
    // tree...
    //
    unsigned int ip    = 0;
    unsigned int power = 0;
    double       logb1 = log(params->b1);

    while ((ip < NFIRST_PRIMES) && (first_primes[ip] <= params->b1)) {
        ip++;
    }
    context->first_lp_index = ip;

    uint32_array_t* ppowers = alloc_uint32_array(ip);

    for (uint32_t i = 0; i < ip; i++) {
        power = (unsigned int) (logb1 / log(first_primes[i]));
        ppowers->data[i] = pow(first_primes[i], power);
    }
    ppowers->length = ip;

    mpz_tree_t* ptree = prod_tree_ui(ppowers);

    mpz_set(context->k, ptree->data[0]);

    free_mpz_tree(ptree);
    free_uint32_array(ppowers);
}
#if TIFA_HAVE_MONT_ECM
//------------------------------------------------------------------------------
static ecode_t mont_curve(factoring_machine_t* const machine) {
    ecm_context_t* const context = (ecm_context_t*) machine->context;
    ecm_params_t*  const params  = (ecm_params_t*)  machine->params;

    if (   (context->mplan == NULL)
        || (context->mplan_b1 != params->b1)
        || (context->mplan_b2 != params->b2)) {
        free_mont_ecm_plan(context->mplan);
        context->mplan    = alloc_mont_ecm_plan(params->b1, params->b2);
        context->mplan_b1 = params->b1;
        context->mplan_b2 = params->b2;
    }
    context->sigma  = next_sigma(context);

    ecode_t ecode = NO_FACTOR_FOUND;

    mpz_t g;
    mpz_init(g);
    const bool found
        = (context->nmont_curves < NZ2Z8_CURVES)
          ? mont_ecm_nth_curve(g, context->n, context->nmont_curves,
                               context->mplan)
          : mont_ecm_curve(g, context->n, (uint32_t) context->sigma,
                           context->mplan);
    context->nmont_curves++;
    if (found) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, context->n, g);
        append_mpz_to_array(machine->factors, g);
        ecode = SOME_FACTORS_FOUND;
    }
    mpz_clear(g);

    return ecode;
}
#endif
//------------------------------------------------------------------------------
static ecode_t init_ecm_context(factoring_machine_t* const machine) {

    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context             = (void*) tifa_malloc(sizeof(ecm_context_t));
    ecm_context_t* const context = (ecm_context_t* const) machine->context;
    ecm_params_t*  const params  = (ecm_params_t* const)  machine->params;
    
    context->factors = machine->factors;
    context->params  = params;

    mpz_init(context->a);
    mpz_init(context->b);
    mpz_init(context->x0);
    mpz_init(context->z0);
    mpz_init(context->xk);
    mpz_init(context->zk);
    mpz_init(context->k);
    mpz_init_set(context->n, machine->n);
    //
    // The curves follow from params->seed, or from n when it is 0: runs are
    // reproducible and share no state with other threads.
    //
    const uint64_t seed = (params->seed != 0)
                        ? params->seed
                        : mpz_fdiv_ui(context->n, 4294967291UL);
    context->rng = (seed << 1) | 1;
    context->nmont_curves = 0;

    mpz_init(context->tmp_add_1);
    mpz_init(context->tmp_add_2);
    mpz_init(context->tmp_add_3);
    mpz_init(context->tmp_add_4);
    mpz_init(context->tmp_add_5);
    mpz_init(context->tmp_add_6);

    mpz_init(context->tmp_mul_1);
    mpz_init(context->tmp_mul_2);
    mpz_init(context->tmp_mul_3);
    mpz_init(context->tmp_mul_4);

    context->sigma  = 5;

    context->precomputations_done   = false;
    context->rounds_at_top = 0;
    //
    // Escalate from the first level above the given bounds, up to the level
    // for factors of half the size of n.
    //
    const uint32_t fmax = (mpz_sizeinbase(context->n, 2) + 1) / 2;
    context->top_level = NLEVELS - 1;
    for (uint32_t i = 0; i < NLEVELS; i++) {
        if (ecm_levels[i].fbits >= fmax) {
            context->top_level = i;
            break;
        }
    }
    context->next_level = NLEVELS;
    for (uint32_t i = 0; i < NLEVELS; i++) {
        if (ecm_levels[i].b1 > params->b1) {
            context->next_level = i;
            break;
        }
    }

#if TIFA_HAVE_MONT_ECM
    context->use_mont = mont_ecm_handles(context->n);
    context->mplan    = NULL;
    context->mplan_b1 = 0;
    context->mplan_b2 = 0;
#endif

    set_stage1_bound(context, params);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_ecm_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    ecm_context_t* context = (ecm_context_t*) machine->context;

    mpz_clear(context->n);
    mpz_clear(context->a);
    mpz_clear(context->b);
    mpz_clear(context->x0);
    mpz_clear(context->z0);
    mpz_clear(context->xk);
    mpz_clear(context->zk);
    mpz_clear(context->k);

    mpz_clear(context->tmp_add_1);
    mpz_clear(context->tmp_add_2);
    mpz_clear(context->tmp_add_3);
    mpz_clear(context->tmp_add_4);
    mpz_clear(context->tmp_add_5);
    mpz_clear(context->tmp_add_6);

    mpz_clear(context->tmp_mul_1);
    mpz_clear(context->tmp_mul_2);
    mpz_clear(context->tmp_mul_3);
    mpz_clear(context->tmp_mul_4);

    if (context->precomputations_done) {
        for (unsigned int ij = 0; ij < context->njs; ij++) {
            unsigned int j = context->Js[ij];
            mpz_clear(context->SX[j]);
            mpz_clear(context->SZ[j]);
        }
        free(context->SX);
        free(context->SZ);

        unsigned int nrows = context->m_max - context->m_min + 1;
        for (unsigned int i = 0U; i < nrows; i++) {
            free(context->prime_table[i]);
        }
        free(context->prime_table);
        free(context->Js);
        free(context->gcd_table);

#if PHASE_2_GLOBALS
        mpz_clear(context->tmp_phase2_1);
        mpz_clear(context->tmp_phase2_2);
        mpz_clear(context->tmp_phase2_3);
        mpz_clear(context->tmp_phase2_4);
        mpz_clear(context->tmp_phase2_5);
        mpz_clear(context->tmp_phase2_6);
        mpz_clear(context->tmp_phase2_7);
        mpz_clear(context->tmp_phase2_8);
        mpz_clear(context->tmp_phase2_9);
        mpz_clear(context->tmp_phase2_10);
        mpz_clear(context->tmp_phase2_11);
#endif

    }
#if TIFA_HAVE_MONT_ECM
    free_mont_ecm_plan(context->mplan);
#endif
    free(context);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_ecm_context(factoring_machine_t* const machine) {
    //
    // Strategy: move up one level of ecm_levels (larger factors, larger
    // bounds) until the level for factors of half the size of n, then run
    // that level MAX_ROUNDS_AT_TOP more times before giving up.
    //
    INIT_TIMER;
    START_TIMER;

    ecm_context_t* context = (ecm_context_t*) machine->context;
    ecode_t ecode = SUCCESS;

    if (context->next_level <= context->top_level) {
        PRINT_UPDATE_LARGER_BOUNDS_MSG;
        switch_to_level(machine, context->next_level);
        context->next_level++;

    } else if (context->rounds_at_top < MAX_ROUNDS_AT_TOP) {
        context->rounds_at_top++;
        PRINT_UPDATE_MORE_CURVES_MSG;

    } else {
        PRINT_UPDATE_GIVEUP_MSG;
        ecode = GIVING_UP;
    }
    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static void switch_to_level(factoring_machine_t* const machine,
                            const uint32_t level) {
    //
    // Switch to the bounds of the given level and recompute what needs to be
    // recomputed...
    //
    ecm_context_t* context = (ecm_context_t*) machine->context;
    ecm_params_t*  params  = (ecm_params_t*)  machine->params;

    params->b1      = ecm_levels[level].b1;
    params->b2      = ECM_B2_RATIO * ecm_levels[level].b1;
    params->ncurves = ecm_levels[level].ncurves;

    if (context->precomputations_done) {
        //
        // _TO_DO_: Actually these clear_* are not needed here (except for
        //          prime_table). We can bypass this memory cleaning by
        //          adding a little bit of logic in perform_precomputations.
        //
        for (unsigned int ij = 0; ij < context->njs; ij++) {
            unsigned int j = context->Js[ij];
            mpz_clear(context->SX[j]);
            mpz_clear(context->SZ[j]);
        }
        free(context->SX);
        free(context->SZ);

        unsigned int nrows = context->m_max - context->m_min + 1;
        for (unsigned int i = 0U; i < nrows; i++) {
            free(context->prime_table[i]);
        }
        free(context->prime_table);
        free(context->Js);
        free(context->gcd_table);
    }
    context->precomputations_done = false;

    set_stage1_bound(context, params);
}
//------------------------------------------------------------------------------
static ecode_t perform_ecm(factoring_machine_t* const machine) {

    INIT_TIMER;
    START_TIMER;

    INIT_NAMED_TIMER(phase_1);
    INIT_NAMED_TIMER(phase_2);
    INIT_NAMED_TIMER(choose);

    ecm_context_t* context = (ecm_context_t*) machine->context;
    ecm_params_t*  params  = (ecm_params_t*)  machine->params;

    ecode_t exit_code = NO_FACTOR_FOUND;

    uint32_t ncurves = 0;

    while (ncurves < params->ncurves) {

        ncurves++;

#if TIFA_HAVE_MONT_ECM
        if (context->use_mont) {
            exit_code = mont_curve(machine);
            if (exit_code == SOME_FACTORS_FOUND) {
                break;
            }
            continue;
        }
#endif
        //
        // Choose a random curve using Suyama's parametrization...
        //
        START_NAMED_TIMER(choose);
        exit_code = choose_curve(context);
        STOP_NAMED_TIMER(choose);

        if (exit_code == SOME_FACTORS_FOUND) {
            //
            // Extremely unlikely!
            //
            break;
        }
        //
        // Perform first phase of ECM...
        //
        START_NAMED_TIMER(phase_1);
        exit_code = phase_1(machine);
        STOP_NAMED_TIMER(phase_1);

        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (params->b2 == 0) {
            //
            // No second phase is performed if b2 is set to 0.
            //
            continue;
        }
        if (!context->precomputations_done) {
            perform_precomputations(context);
        }
        //
        // Perform second phase of ECM...
        //
        START_NAMED_TIMER(phase_2);
        exit_code = phase_2(machine);
        STOP_NAMED_TIMER(phase_2);

        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
    }
    STOP_TIMER;

    PRINT_NCURVES_USED(ncurves);

    PRINT_CHOOSED_CURVES_IN;
    PRINT_NAMED_TIMING(choose);
    PRINT_PHASE_1_IN;
    PRINT_NAMED_TIMING(phase_1);
    PRINT_PHASE_2_IN;
    PRINT_NAMED_TIMING(phase_2);

    return exit_code;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor(factors, multis, n, mode);
}
//------------------------------------------------------------------------------
static ecode_t phase_1(factoring_machine_t* const machine) {

    ecm_context_t* const context = (ecm_context_t* const) machine->context;

    multiply_point(
        context->xk, context->zk, context->x0, context->z0, context->k, context
    );

    ecode_t exit_code = NO_FACTOR_FOUND;

    mpz_t g;
    mpz_init(g);

    mpz_gcd(g, context->zk, machine->n);

    if ((mpz_cmp_ui(g, 1) != 0) && (mpz_cmp(g, machine->n) != 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, machine->n, g);
        append_mpz_to_array(machine->factors, g);
        exit_code = SOME_FACTORS_FOUND;
    }
    mpz_clear(g);

    return exit_code;
}
//------------------------------------------------------------------------------
static ecode_t phase_2(factoring_machine_t* const machine) {
    //
    // This second phase implementation follows the standard continuation
    // and is implemented in a way reminiscent of the description
    // given in the article "Implementing the Elliptic Curve Method of Factoring
    // in Reconfigurable Hardware" by Kris Gaj et al.
    //
    ecm_context_t* const context = (ecm_context_t* const) machine->context;
    ecode_t exit_code = NO_FACTOR_FOUND;

    mpz_t g;
    mpz_t d;

    mpz_init_set_ui(d, 1);
    mpz_init(g);

    unsigned long int DO2   = context->DO2;
    unsigned long int D     = context->D;
    unsigned long int m_min = context->m_min;
    unsigned long int m_max = context->m_max;

    unsigned int* Js        = context->Js;
    unsigned int  njs       = context->njs;
    unsigned int** prime_table = context->prime_table;

    mpz_t* SX        = context->SX;
    mpz_t* SZ        = context->SZ;

#if PHASE_2_GLOBALS

    mpz_ptr xq0    = context->tmp_phase2_1;
    mpz_ptr zq0    = context->tmp_phase2_2;
    mpz_ptr xq     = context->tmp_phase2_3;
    mpz_ptr zq     = context->tmp_phase2_4;
    mpz_ptr xr     = context->tmp_phase2_5;
    mpz_ptr zr     = context->tmp_phase2_6;
    mpz_ptr xdiff  = context->tmp_phase2_7;
    mpz_ptr zdiff  = context->tmp_phase2_8;
    mpz_ptr xq_old = context->tmp_phase2_9;
    mpz_ptr zq_old = context->tmp_phase2_10;
    mpz_ptr tmp    = context->tmp_phase2_11;

    mpz_set(xq0, context->xk);
    mpz_set(zq0, context->zk);

    mpz_set(xq, xq0);
    mpz_set(zq, zq0);
#else
    mpz_t xq0;
    mpz_t zq0;

    mpz_t xq;
    mpz_t zq;

    mpz_t xr;
    mpz_t zr;

    mpz_t xdiff;
    mpz_t zdiff;

    mpz_t xq_old;
    mpz_t zq_old;

    mpz_init(xq0);
    mpz_init(zq0);

    mpz_init(xq);
    mpz_init(zq);

    mpz_init(xr);
    mpz_init(zr);

    mpz_init(xq_old);
    mpz_init(zq_old);

    mpz_init(xdiff);
    mpz_init(zdiff);

    mpz_t tmp;
    mpz_init(tmp);

    mpz_set(xq0, context->xk);
    mpz_set(zq0, context->zk);

    mpz_set(xq, xq0);
    mpz_set(zq, zq0);
#endif

    if (gcd_ulint(1, D) == 1) {
        mpz_set(context->SX[1], xq);
        mpz_set(context->SZ[1], zq);
    }
    double_point(xq, zq, xq, zq, context); // Q = [2] Q0

    mpz_set(xdiff, xq0);
    mpz_set(zdiff, zq0);

    for (unsigned long int j = 2; j <= DO2; j++) {
        if (gcd_ulint(j, D) == 1) {
            mpz_set(context->SX[j], xq);
            mpz_set(context->SZ[j], zq);
        }
        mpz_set(xq_old, xq);
        mpz_set(zq_old, zq);

        add_points(xq, zq, xq, zq, xq0, zq0, xdiff, zdiff, context);

        mpz_set(xdiff, xq_old);
        mpz_set(zdiff, zq_old);
    }
    mpz_set(xq, xq0);
    mpz_set(zq, zq0);

    mpz_set_ui(tmp, D);
    multiply_point(xq, zq, xq0, zq0, tmp, context);     // Q = [D] Q0

    mpz_set_ui(tmp, m_min-1);
    multiply_point(xdiff, zdiff, xq, zq, tmp, context); // DIFF = [m_min-1] Q

    mpz_set_ui(tmp, m_min);
    multiply_point(xr, zr, xq, zq, tmp, context);       // R = [m_min] Q

    for (unsigned int im = 0; im <= (m_max - m_min); im++) {

        unsigned int m = m_min + im;

        for (unsigned int ij = 0; ij < njs; ij++) {
            unsigned int j = Js[ij];
            if (prime_table[im][j] == 1) {
                mpz_mul(tmp, xr, SZ[j]);
                mpz_submul(tmp, zr, SX[j]);
                mpz_mul(d, d, tmp);
                mpz_mod(d, d, context->n);
            }
        }
        mpz_set(xq_old, xr);
        mpz_set(zq_old, zr);

        if (m == 1) {
            double_point(xr, zr, xr, zr, context);
        } else {
            add_points(xr, zr, xr, zr, xq, zq, xdiff, zdiff, context);
        }
        mpz_set(xdiff, xq_old);
        mpz_set(zdiff, zq_old);
    }
    mpz_gcd(g, d, context->n);

    if ((mpz_cmp_ui(g, 1) != 0) && (mpz_cmp(g, machine->n) != 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, machine->n, g);
        append_mpz_to_array(machine->factors, g);
        exit_code = SOME_FACTORS_FOUND;
    }

#if PHASE_2_GLOBALS
    mpz_clear(d);
    mpz_clear(g);
#else
    mpz_clear(d);
    mpz_clear(g);
    mpz_clear(xq0);
    mpz_clear(zq0);
    mpz_clear(xq);
    mpz_clear(zq);
    mpz_clear(xr);
    mpz_clear(zr);
    mpz_clear(tmp);
    mpz_clear(xdiff);
    mpz_clear(zdiff);
    mpz_clear(xq_old);
    mpz_clear(zq_old);
#endif

    return exit_code;
}
//------------------------------------------------------------------------------
static ecode_t choose_curve(ecm_context_t* const context) {
    //
    // Choose a curve and a starting point using Suyama's parametrization.
    //
    ecode_t exit_code = SUCCESS;

    context->sigma = next_sigma(context);

    mpz_ptr n = context->n;

    mpz_ptr v       = context->tmp_add_1;
    mpz_ptr u       = context->tmp_add_2;
    mpz_ptr u3v4    = context->tmp_add_3;
    mpz_ptr u3v4inv = context->tmp_add_4;
    mpz_ptr tmp     = context->tmp_add_5;
    mpz_ptr tmp2    = context->tmp_add_6;

    mpz_set_ui(u, context->sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5); // u = sigma^2 - 5

    mpz_set_ui(v, context->sigma);
    mpz_mul_ui(v, v, 4); // v = 4 * sigma

    mpz_mul(u3v4, u, u);
    mpz_mul(u3v4, u3v4, u);
    mpz_mul(u3v4, u3v4, v);
    mpz_mul_ui(u3v4, u3v4, 4); // u3v4 = 4 * u^3 * v

    mpz_gcd(tmp, u3v4, n);

    if (mpz_cmp_ui(tmp, 1) != 0) {
        //
        // We have a factor (very unlikely)!
        // Abort everything and save the found factor!
        //
        exit_code = SOME_FACTORS_FOUND;
        append_mpz_to_array(context->factors, tmp);
        mpz_divexact(tmp, n, tmp);
        append_mpz_to_array(context->factors, tmp);

        return exit_code;
    }

    mpz_invert(u3v4inv, u3v4, n); // inv(4 * u^3 * v) mod n

    mpz_sub(tmp, v, u);
    mpz_pow_ui(tmp, tmp, 3); // tmp = (v - u)^3

    mpz_mul_ui(tmp2, u, 3);
    mpz_add(tmp2, tmp2, v); // tmp2 = 3 * u + v

    mpz_mul(context->a, tmp, tmp2);
    mpz_mul(context->a, context->a, u3v4inv);
    mpz_sub_ui(context->a, context->a, 2);

    mpz_set_ui(context->z0, 1);

    mpz_pow_ui(context->x0, u, 4);
    mpz_mul_ui(context->x0, context->x0, 4);
    mpz_mul(context->x0, context->x0, u3v4inv);

    mpz_powm_ui(context->x0, context->x0, 3, n);

    mpz_pow_ui(v, v, 3); // v = v^3

    mpz_gcd(context->b, v, n);

    if (mpz_cmp_ui(context->b, 1) != 0) {
        //
        // We have a factor (very unlikely)!
        // Abort everything and save the found factor!
        //
        exit_code = SOME_FACTORS_FOUND;
        append_mpz_to_array(context->factors, context->b);
        mpz_divexact(context->b, n, tmp);
        append_mpz_to_array(context->factors, context->b);
        return exit_code;
    }
    mpz_invert(context->b, v, n);       // inv(v^3) mod n
    mpz_mul(context->b, context->b, u); // u * inv(v^3) mod n

    mpz_mod(context->a, context->a, n);
    mpz_mod(context->b, context->b, n);

    return exit_code;
}
//------------------------------------------------------------------------------
static void add_points(mpz_t xc, mpz_t zc,
                       const mpz_t xa, const mpz_t za,
                       const mpz_t xb, const mpz_t zb,
                       const mpz_t x0, const mpz_t z0,
                       ecm_context_t* const context) {
    //
    // Set (xc : zc) to (xa : za) + (xb : zb) knowing
    // (x0 : z0) = (xa : za) - (xb : zb).
    // Uses the Montgomery representation.
    //
    // _WARNING_: xc and x0 should point to different memory spaces!
    //
    mpz_ptr n  = context->n;

    mpz_ptr r = context->tmp_add_1;
    mpz_ptr s = context->tmp_add_2;
    mpz_ptr t = context->tmp_add_3;
    mpz_ptr u = context->tmp_add_4;
    mpz_ptr v = context->tmp_add_5;
    mpz_ptr w = context->tmp_add_6;

    mpz_sub(r, xa, za);
    mpz_sub(s, xb, zb);

    mpz_add(t, xa, za);
    mpz_add(u, xb, zb);

    mpz_mul(v, r, u);
    mpz_mul(w, s, t);

    mpz_add(xc, v, w);
    mpz_mul(xc, xc, xc);
    mpz_mul(xc, xc, z0);

    mpz_sub(zc, v, w);
    mpz_mul(zc, zc, zc);
    mpz_mul(zc, zc, x0);

    mpz_mod(xc, xc, n);
    mpz_mod(zc, zc, n);

    return;
}
//------------------------------------------------------------------------------
static void double_point(mpz_t x2a, mpz_t z2a,
                         const mpz_t xa, const mpz_t za,
                         ecm_context_t* const context) {
    //
    // Set (x2a : z2a) to [2] (xa : za).
    // Uses the Montgomery representation.
    //
    mpz_ptr n = context->n;
    mpz_ptr a = context->a;

    mpz_ptr plus  = context->tmp_add_1;
    mpz_ptr minus = context->tmp_add_2;
    mpz_ptr xz4   = context->tmp_add_3;
    mpz_ptr tmp   = context->tmp_add_4;

    mpz_add(plus, xa, za);
    mpz_sub(minus, xa, za);

    mpz_mul(plus, plus, plus);
    mpz_mul(minus, minus, minus);

    mpz_sub(xz4, plus, minus);
    mpz_mul(x2a, plus, minus);

    mpz_tdiv_q_2exp(plus, xz4, 2);
    mpz_add_ui(tmp, a, 2);
    mpz_mul(plus, plus, tmp);

    mpz_add(z2a, minus, plus);
    mpz_mul(z2a, z2a, xz4);

    mpz_mod(x2a, x2a, n);
    mpz_mod(z2a, z2a, n);

    return;
}
//------------------------------------------------------------------------------
static void multiply_point(mpz_t xk, mpz_t zk,
                           const mpz_t x0, const mpz_t z0,
                           const mpz_t k,
                           ecm_context_t* const context) {
    //
    // Set (xk : zk) to [k] (x0 : z0) using Montgomery's ladder algorithm.
    // Uses the Montgomery representation.
    //
    mpz_ptr qx = context->tmp_mul_1;
    mpz_ptr qz = context->tmp_mul_2;
    mpz_ptr px = context->tmp_mul_3;
    mpz_ptr pz = context->tmp_mul_4;

    mpz_set(qx, x0);
    mpz_set(qz, z0);

    double_point(px, pz, x0, z0, context);

    int s = mpz_sizeinbase(k, 2) - 2;

    while (s >= 0) {
        if (mpz_tstbit(k, s) == 1) {
            add_points(qx, qz, px, pz, qx, qz, x0, z0, context);
            double_point(px, pz, px, pz, context);
        } else {
            add_points(px, pz, px, pz, qx, qz, x0, z0, context);
            double_point(qx, qz, qx, qz, context);
        }
        s--;
    }
    mpz_set(xk, qx);
    mpz_set(zk, qz);

    return;
}
//------------------------------------------------------------------------------
static ecode_t perform_precomputations(ecm_context_t* const context) {
    //
    // Fills the tables needed for the second phase.
    //
    // These tables depends essentially on the bounds b1 and b2, so as soon
    // as the "optimal" parameters' values are determined, these tables can
    // be hardcoded.
    //
    // Ref: "Implementing the Elliptic Curve Method of Factoring in
    // Reconfigurable Hardware", K. Gaj et al., Cryptographic Hardware and
    // Embedded Systems - CHES 2006.
    //
    INIT_TIMER;
    PRINT_PERFORMING_P2_PRECOMP;
    START_TIMER;

#if PHASE_2_GLOBALS
    mpz_init(context->tmp_phase2_1);
    mpz_init(context->tmp_phase2_2);
    mpz_init(context->tmp_phase2_3);
    mpz_init(context->tmp_phase2_4);
    mpz_init(context->tmp_phase2_5);
    mpz_init(context->tmp_phase2_6);
    mpz_init(context->tmp_phase2_7);
    mpz_init(context->tmp_phase2_8);
    mpz_init(context->tmp_phase2_9);
    mpz_init(context->tmp_phase2_10);
    mpz_init(context->tmp_phase2_11);
#endif

    context->B2 = context->params->b2;
    context->B1 = context->params->b1;
    //
    // _WARNING_: This D value (representing a speed/memory trade-off) is from
    //            the aforementionned paper. Some tests are needed to assess
    //            the better value in our case.
    //
    context->D   = 210;
    context->DO2 = 105;

    context->m_min = (context->B1 + context->DO2) / context->D;
    context->m_max = (context->B2 - context->DO2) / context->D + 1;

    unsigned long int D   = context->D;
    unsigned long int DO2 = context->DO2;
    unsigned long int m_min = context->m_min;
    unsigned long int m_max = context->m_max;
    //
    // _NOTE_: More space than strictly required is allocated in order to keep
    //         a simple implementation without translating indices.
    //
    context->gcd_table = tifa_malloc((DO2 + 1) * sizeof(unsigned int));
    context->SX        = tifa_malloc((DO2 + 1) * sizeof(mpz_t));
    context->SZ        = tifa_malloc((DO2 + 1) * sizeof(mpz_t));
    context->Js        = tifa_malloc((DO2 + 1) * sizeof(unsigned int));
    context->njs       = 0;

    for (unsigned long int j = 1; j <= DO2; j++) {
        if (gcd_ulint(j, D) == 1) {
            context->gcd_table[j] = 1;
            mpz_init(context->SX[j]);
            mpz_init(context->SZ[j]);
            context->Js[context->njs] = j;
            context->njs++;
        } else {
            context->gcd_table[j] = 0;
        }
    }
    unsigned int nrows = m_max - m_min + 1;
    unsigned int ncols = DO2 + 1;

    context->prime_table = tifa_malloc(nrows * sizeof(unsigned int*));

    for (unsigned int i = 0U; i < nrows; i++) {
        context->prime_table[i] = tifa_malloc(ncols * sizeof(unsigned int));
    }
    unsigned int deltam = m_max - m_min;

    for (unsigned int im = 0; im <= deltam; im++) {
        unsigned int m = m_min + im;
        for (unsigned int j = 1; j <= DO2; j++) {
            if (is_prime(m * D + j)
                || is_prime((m * D >= j) ? (m * D - j) : (j - m * D))) {
                context->prime_table[im][j] = 1;
                continue;
            }
            context->prime_table[im][j] = 0;
        }
    }
    context->precomputations_done = true;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
