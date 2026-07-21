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
 * \file    p_minus_1.c
 * \author  Jerome Milan
 * \date    2026
 * \version 1.0
 */

#include <stdlib.h>
#include <inttypes.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "macros.h"
#include "factoring_machine.h"
#include "p_minus_1.h"
#include "tifa_factor.h"

#define __PREFIX__      "p-1: "
#define __VERBOSE__     TIFA_VERBOSE_PM1
#define __TIMING__      TIFA_TIMING_PM1

#include "messages.h"

//
// Maximum number of times the context is updated (i.e. B1 is grown and the
// stage is retried) before giving up.
//
#define MAX_NUPDATES     4
//
// Factor by which B1 grows on each update.
//
#define B1_GROWTH_FACTOR 2UL
//
// Default parameters.
//
#define DEFAULT_BASE     2
#define DEFAULT_B1       50000

//
// An ad-hoc structure holding the data variables used by the Pollard's p-1
// implementation.
//
struct struct_p_minus_1_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // p-1 parameters.
    //
    p_minus_1_params_t* params;
    //
    // Number of times the context was updated.
    //
    int nupdates;
    //
    // Working variables.
    //
    mpz_t a;
    mpz_t g;
};
typedef struct struct_p_minus_1_context_t p_minus_1_context_t;

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_p_minus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_p_minus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_p_minus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_p_minus_1(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t p_minus_1(mpz_array_t* const factors, uint32_array_t* const multis,
                  const mpz_t n, const p_minus_1_params_t* const params,
                  const factoring_mode_t mode) {

    PRINT_INTRO_MSG("Pollard's p-1");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_p_minus_1_context;
    machine.perform_algo_func   = perform_p_minus_1;
    machine.update_context_func = update_p_minus_1_context;
    machine.clear_context_func  = clear_p_minus_1_context;
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
void set_p_minus_1_params_to_default(const mpz_t n MAYBE_UNUSED,
                                     p_minus_1_params_t* const params) {
    params->base = DEFAULT_BASE;
    params->b1   = DEFAULT_B1;
    if (params->b1 > LARGEST_PRIME) {
        params->b1 = LARGEST_PRIME;
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static ecode_t init_p_minus_1_context(factoring_machine_t* const machine) {
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context = (void*) malloc(sizeof(p_minus_1_context_t));
    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;

    mpz_init_set(context->n, machine->n);
    mpz_init_set_ui(context->a, 1);
    mpz_init_set_ui(context->g, 1);

    context->params   = (p_minus_1_params_t*) machine->params;
    context->nupdates = 0;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_p_minus_1_context(factoring_machine_t* const machine) {
    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;

    if (context != NULL) {
        mpz_clear(context->n);
        mpz_clear(context->a);
        mpz_clear(context->g);
        free(context);
    }

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_p_minus_1_context(factoring_machine_t* const machine) {
    INIT_TIMER;
    START_TIMER;

    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;
    p_minus_1_params_t*  params  = context->params;
    ecode_t ecode = SUCCESS;

    context->nupdates++;

    if ((context->nupdates > MAX_NUPDATES) || (params->b1 >= LARGEST_PRIME)) {
        ecode = GIVING_UP;
        PRINT_UPDATE_GIVEUP_MSG;
    } else {
        //
        // Grow the smoothness bound and retry, clamping to the largest prime
        // available in the precomputed table.
        //
        unsigned long nb1 = (unsigned long) params->b1 * B1_GROWTH_FACTOR;
        if (nb1 > LARGEST_PRIME) {
            nb1 = LARGEST_PRIME;
        }
        params->b1 = (uint32_t) nb1;
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t perform_p_minus_1(factoring_machine_t* const machine) {
    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;
    p_minus_1_params_t*  params  = (p_minus_1_params_t*)  machine->params;

    INIT_TIMER;
    START_TIMER;

    mpz_ptr n = context->n;
    mpz_ptr a = context->a;
    mpz_ptr g = context->g;

    ecode_t ecode = NO_FACTOR_FOUND;

    //
    // Cheap early check: a common factor between the base and n is itself a
    // factor of n.
    //
    mpz_set_ui(g, params->base);
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
    // Stage 1: compute a <- base^E mod n, where
    //     E = prod_{p prime <= b1} p^floor(log_p(b1)),
    // accumulated one prime at a time (so no huge exponent is materialized).
    //
    mpz_set_ui(a, params->base);

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
        mpz_powm_ui(a, a, pe, n);
        ip++;
    }

    //
    // g = gcd(a - 1, n). A nontrivial gcd is a factor of n.
    //
    mpz_sub_ui(a, a, 1);
    mpz_gcd(g, a, n);

    if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(machine->factors, g);
        ecode = SOME_FACTORS_FOUND;
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
    // Pollard's p-1.
    //
    //p_minus_1_params_t params;
    //set_p_minus_1_params_to_default(n, &params);
    //return p_minus_1(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------
