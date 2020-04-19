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
 * \file    rho.c
 * \author  Jerome Milan
 * \date    Tue Apr 22 2014
 * \version 2014-04-22
 */

#include <stdlib.h>
#include <inttypes.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "macros.h"
#include "factoring_machine.h"
#include "rho.h"
#include "tifa_factor.h"

#define RHO_DEBUG 0
#if RHO_DEBUG
  #include <stdio.h>
  #define RHO_GMP_PRINT(...) gmp_printf(__VA_ARGS__)
  #define RHO_PRINT(...)     printf(__VA_ARGS__)
#else
  #define RHO_GMP_PRINT(...) /* intentionally left empty */
  #define RHO_PRINT(...)     /* intentionally left empty */
#endif

#define __PREFIX__      "rho: "
#define __VERBOSE__     TIFA_VERBOSE_RHO
#define __TIMING__      TIFA_TIMING_RHO

#include "messages.h"

#define MAX_NUPDATES              3
#define MAX_NITERS_GROWTH_FACTOR  4

#define DEFAULT_A          1
#define DEFAULT_B          1
#define DEFAULT_ACC_NITERS 100
#define DEFAULT_MAX_NITERS 100000
#define DEFAULT_RHO_METHOD 1

//
// An ad-hoc structure holding the data variables used by the Pollard's rho
// implementation.
//
struct struct_rho_context_t {
    //
    // Number to factor
    //
    mpz_t n;
    //
    // Rho parameters
    //
    rho_params_t* params;
    //
    // Number of times the ocntext was updated
    //
    int nupdates;
    //
    // Working variables
    //
    mpz_t tortoise;
    mpz_t hare;
    mpz_t tortoise_bck;
    mpz_t hare_bck;
    mpz_t g;
    mpz_t diff;
    mpz_t prod;


};
typedef struct struct_rho_context_t rho_context_t;

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_rho_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_rho_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_rho_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_rho(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_rho_floyd(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_rho_brent(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------
inline static void iteration_function(rho_context_t* context, mpz_t fx, mpz_t x);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t rho(mpz_array_t* const factors, uint32_array_t* const multis,
            const mpz_t n, const rho_params_t* const params,
            const factoring_mode_t mode) {

    PRINT_INTRO_MSG("Pollard's rho");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_rho_context;
    machine.perform_algo_func   = perform_rho;
    machine.update_context_func = update_rho_context;
    machine.clear_context_func  = clear_rho_context;
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
void set_rho_params_to_default(const mpz_t n MAYBE_UNUSED, rho_params_t* const params) {
    params->a          = DEFAULT_A;
    params->b          = DEFAULT_B;
    params->acc_niters = DEFAULT_ACC_NITERS;
    params->max_niters = DEFAULT_MAX_NITERS;
    params->method     = DEFAULT_RHO_METHOD;
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
static ecode_t init_rho_context(factoring_machine_t* const machine) {
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context       = (void*)malloc(sizeof(rho_context_t));
    rho_context_t* context = (rho_context_t*)machine->context;

    mpz_init_set(context->n, machine->n);

    mpz_init_set_ui(context->tortoise, 1);
    mpz_init_set_ui(context->hare, 1);
    mpz_init_set_ui(context->tortoise_bck, 1);
    mpz_init_set_ui(context->hare_bck, 1);
    mpz_init_set_ui(context->g, 1);
    mpz_init_set_ui(context->diff, 0);
    mpz_init_set_ui(context->prod, 1);

    context->params = (rho_params_t*)machine->params;

    context->nupdates = 0;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_rho_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    rho_context_t* context = (rho_context_t*)machine->context;

    if (context != NULL) {
        mpz_clear(context->n);
        mpz_clear(context->tortoise);
        mpz_clear(context->hare);
        mpz_clear(context->tortoise_bck);
        mpz_clear(context->hare_bck);
        mpz_clear(context->g);
        mpz_clear(context->diff);
        mpz_clear(context->prod);
        free(context);
    }

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_rho_context(factoring_machine_t* const machine MAYBE_UNUSED) {
    INIT_TIMER;
    START_TIMER;

    rho_context_t* context = (rho_context_t*)machine->context;
    ecode_t ecode = SUCCESS;

    if (context->nupdates == MAX_NUPDATES) {
        ecode = GIVING_UP;
        PRINT_UPDATE_GIVEUP_MSG;

    } else {
        PRINT_UPDATING_RHO_CONTEXT;

        rho_params_t* params = (rho_params_t*)context->params;

        params->a += 2;
        if (IS_EVEN(params->a)) {
            params->a -= 1;
        }
        params->b += 2;
        if (IS_EVEN(params->b)) {
            params->b -= 1;
        }
        params->max_niters *= MAX_NITERS_GROWTH_FACTOR;
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t perform_rho(factoring_machine_t* const machine) {
    rho_params_t* params  = (rho_params_t*)machine->params;
    switch (params->method) {
        case FLOYD_CYCLE_FINDING:
            return perform_rho_floyd(machine);

        case BRENT_CYCLE_FINDING:
        default:
            return perform_rho_brent(machine);
    }
}
//------------------------------------------------------------------------------
static ecode_t perform_rho_brent(factoring_machine_t* const machine) {
    //
    // This implementation blindly follows Brent's algorithm described in
    // section 7 of his paper "An improved Monte Carlo factorization algoritm",
    // BIT Numerical Mathematics 1980, Volume 20, Issue 2, pp 176-184.
    //
    // The variable names are exactly the ones from the article.
    //
    rho_context_t* context = (rho_context_t*) machine->context;
    rho_params_t*  params  = (rho_params_t*)  machine->params;

    PRINT_PERFORMING_BRENT_CYCLE;
    INIT_TIMER;
    START_TIMER;

    mpz_t q;
    uint32_t r = 1;

    uint32_t i = 1;
    uint32_t imax = params->max_niters;
    uint32_t m = params->acc_niters;
    uint32_t k = 0;

    mpz_ptr n  = context->n;
    mpz_ptr x  = context->tortoise;
    mpz_ptr y  = context->hare;
    mpz_ptr ys = context->hare_bck;
    mpz_ptr g  = context->g;
    mpz_ptr diff = context->diff;

    ecode_t ecode = NO_FACTOR_FOUND;

    mpz_set_ui(g, 1);
    mpz_init_set_ui(q, 1);

    while ((mpz_cmp_ui(g, 1) == 0) && (i < imax)) {
        mpz_set(x, y);
        for (uint32_t j = 0; j < r; j++) {
            iteration_function(context, y, y);
        }
        i += r;
        k = 0;
        do {
            mpz_set(ys, y);
            uint32_t mmin = MIN(m, r - k);
            for (uint32_t j = 0; j < mmin; j++) {
                iteration_function(context, y, y);
                mpz_sub(diff, x, y);
                mpz_abs(diff, diff);
                mpz_mul(q, q, diff);
                mpz_mod(q, q, n);
            };
            i += mmin;
            mpz_gcd(g, q, n);
            k = k + m;
            if (i >= imax) {
                break;
            }
        } while ((k < r) && (mpz_cmp_ui(g, 1) == 0));

        r = 2 * r;
    }

    if (mpz_cmp(g, n) == 0) {
        while (i < imax) {
            iteration_function(context, ys, ys);
            i++;
            mpz_sub(diff, x, ys);
            mpz_abs(diff, diff);
            mpz_gcd(g, diff, n);
            if ((mpz_cmp_ui(g, 1) != 0) && (mpz_cmp(g, n) != 0)) {
                append_mpz_to_array(machine->factors, g);
                mpz_divexact(g, machine->n, g);
                append_mpz_to_array(machine->factors, g);
                ecode = SOME_FACTORS_FOUND;
                break;
            }
        }

    } else if (mpz_cmp_ui(g, 1) != 0) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, machine->n, g);
        append_mpz_to_array(machine->factors, g);
        ecode = SOME_FACTORS_FOUND;
    }

    mpz_clear(q);

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t perform_rho_floyd(factoring_machine_t* const machine) {
    rho_context_t* context = (rho_context_t*) machine->context;
    rho_params_t*  params  = (rho_params_t*)  machine->params;

    PRINT_PERFORMING_FLOYD_CYCLE;
    INIT_TIMER;
    START_TIMER;

    uint32_t i = 1;
    uint32_t imax = params->max_niters;
    uint32_t jmax = params->acc_niters;

    mpz_ptr n    = context->n;
    mpz_ptr xi   = context->tortoise;
    mpz_ptr x2i  = context->hare;
    mpz_ptr xib  = context->tortoise_bck;
    mpz_ptr x2ib = context->hare_bck;
    mpz_ptr g    = context->g;
    mpz_ptr diff = context->diff;
    mpz_ptr prod = context->prod;

    ecode_t ecode = NO_FACTOR_FOUND;

    if (jmax > 1) {
        while (i <= imax) {
            mpz_set_ui(prod, 1);
            for (uint32_t j = 0; j < jmax; j++) {
                iteration_function(context, xi, xi);
                iteration_function(context, x2i, x2i);
                iteration_function(context, x2i, x2i);

                mpz_sub(diff, x2i, xi);
                //
                // Accumulate the product of the (x2i - xi) and only perform
                // a gcd computation every jmax iterations.
                //
                mpz_mul(prod, prod, diff);
                mpz_mod(prod, prod, n);
            }
            i = i + 3 * jmax;

            mpz_gcd(g, n, prod);
            if (mpz_cmp(g, n) == 0) {
                //
                // The product of some (x2i - xi) is equal to n. Backtrack
                // to the start of the last product accumulation and perform
                // gcd computations at every steps.
                //
                mpz_set(xi, xib);
                mpz_set(x2i, x2ib);
                i = i - 3 * jmax;
                break;

            } else if (mpz_cmp_ui(g, 1) != 0) {
                append_mpz_to_array(machine->factors, g);
                mpz_divexact(g, machine->n, g);
                append_mpz_to_array(machine->factors, g);
                ecode = SOME_FACTORS_FOUND;
                goto RETURN;
            }
            mpz_set(xib, xi);
            mpz_set(x2ib, x2i);
        }
    }

    while (i <= imax + 3 * jmax) {
        //
        // Perform a gcd computation at each iteration.
        //
        iteration_function(context, xi, xi);
        iteration_function(context, x2i, x2i);
        iteration_function(context, x2i, x2i);
        i += 3;

        mpz_sub(diff, x2i, xi);
        mpz_gcd(g, n, diff);
        if ((mpz_cmp_ui(g, 1) != 0) && (mpz_cmp(g, n) != 0)) {
            append_mpz_to_array(machine->factors, g);
            mpz_divexact(g, machine->n, g);
            append_mpz_to_array(machine->factors, g);
            ecode = SOME_FACTORS_FOUND;
            break;
        }
    }

  RETURN:

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
    // Pollard's Rho.
    //
    //rho_params_t params;
    //set_rho_params_to_default(n, &params);
    //return rho(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------
inline static void iteration_function(rho_context_t* context, mpz_t fx, mpz_t x) {
    mpz_mul(fx, x, x);
    mpz_mul_ui(fx, fx,context->params->a);
    mpz_add_ui(fx, fx, context->params->b);
    mpz_mod(fx, fx, context->n);
}
//------------------------------------------------------------------------------

