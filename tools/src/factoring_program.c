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
 * \file    factoring_program.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <string.h>

#include "tifa_config.h"
#include "funcs.h"
#include "macros.h"
#include "tdiv.h"
#include "factoring_program.h"

#include "tool_utils.h"
#include "common_funcs.h"

#define __TIMING__ 1
#include "../../lib/utils/include/timer.h"

//-----------------------------------------------------------------------------
//
// Default size of arrays.
//
#define PROG_DFLT_ARRAY_LENGTH 8
//
// Set to non-zero to perform benchmarks averaged over NITERATIONS iterations.
//
#define DO_AVERAGED_BENCHMARK  0
//
// Number of iterations to perform when performing averaged benchmarks.
//
#define NITERATIONS            100
//-----------------------------------------------------------------------------
void print_factors(mpz_array_t* const, uint32_array_t* const);
//-----------------------------------------------------------------------------
ecode_t run_program(factoring_program_t* const program) {

    ecode_t rcode = NO_FACTOR_FOUND;
    //
    // Reads arguments given on the command line (if any) and initializes
    // all the needed parameters.
    //
    program->process_args_func(program);

    //
    // Of course, check if program->n is prime before proceeding.
    //
    if (MPZ_IS_PRIME(program->n)) {
        gmp_printf("\nOOPS: %Zd is (probably) prime!\n", program->n);
        return 0;
    }
    printf("Integer to factor:\n");
    printf("------------------\n");
    gmp_printf("\t%Zd\n", program->n);

    //
    // Factors & multiplicities for program->n.
    //
    mpz_array_t*    factors = alloc_mpz_array(PROG_DFLT_ARRAY_LENGTH);
    uint32_array_t* multis  = alloc_uint32_array(PROG_DFLT_ARRAY_LENGTH);

    ecode_t ecode = NO_FACTOR_FOUND;

    if (program->nprimes_tdiv != 0) {
        printf("\n");
        //
        // Proceed with the trial divisions...
        //
        ecode = tdiv(factors, multis, program->n, program->nprimes_tdiv);

        switch (ecode) {
            case COMPLETE_FACTORIZATION_FOUND:
                printf("\n");
                rcode = COMPLETE_FACTORIZATION_FOUND;
                goto clear_tdiv_and_return;

            case SOME_FACTORS_FOUND:
                //
                // Adjust the factoring algorithm's default parameters to account
                // for the smaller unfactored part if default paramters were
                // asked for.
                //
                if (program->argc == 1) {
                    program->set_params_to_default_func(program);
                }
                break;

            default:
                break;
        }
    }
    //
    // Found factors (if any) cannot account for the complete factorization
    // of program->n.
    //
    mpz_t unfactored;

    printf("\nParameters used:\n");
    printf("----------------\n");
    program->print_params_func(program);

    if (ecode == NO_FACTOR_FOUND) {
        mpz_init_set(unfactored, program->n);
    } else {
        mpz_init_set(unfactored, factors->data[factors->length - 1]);
        factors->length--;
        multis->length--;
        rcode = PARTIAL_FACTORIZATION_FOUND;
    }
    if (program->verbose || program->timing
            || TIFA_VERBOSE_TDIV || TIFA_TIMING_TDIV) {
        printf("\n");
    }

    if (program->nprimes_tdiv != 0) {
        printf("Integer to factor after trial division:\n");
        printf("---------------------------------------\n");
        gmp_printf("\t%Zd\n\n", unfactored);
    }

    unsigned int shift = 0;
    //
    // First, check if the unfactored part is not a perfect square.
    //
    while (MPZ_IS_SQUARE(unfactored)) {
        mpz_sqrt(unfactored, unfactored);
        shift++;
        if (mpz_cmp_ui(unfactored, 1) == 0) {
            break;
        }
    }

    if (MPZ_IS_PRIME(unfactored)) {
        //
        // We have found the complete factorization!
        //
        append_mpz_to_array(factors, unfactored);
        append_uint32_to_array(multis, 1 << shift);

        rcode = COMPLETE_FACTORIZATION_FOUND;
        goto clear_sqt_and_return;
    }
    //
    // Factors & multiplicities for the unfactored part of program->n.
    //
    mpz_array_t*    progfa = alloc_mpz_array(PROG_DFLT_ARRAY_LENGTH);
    uint32_array_t* progmu = alloc_uint32_array(PROG_DFLT_ARRAY_LENGTH);

    if (program->verbose || program->timing) {
        printf("%s trace:\n", program->algo_name);
        int len = strlen(program->algo_name);
        for (int i = 0; i < len; i++) {
            printf("-");
        }
        printf("-------\n\n");
    }

#if DO_AVERAGED_BENCHMARK
    INIT_TIMER;
    uint32_t nfails = 0;

    for (unsigned int i = 1; i < NITERATIONS; i++) {
        START_TIMER;
        ecode = program->factoring_algo_func(
                    progfa,
                    progmu,
                    unfactored,
                    program->params,
                    program->mode
                );
        STOP_TIMER;
        if (ecode == NO_FACTOR_FOUND) {
            nfails++;
        }
        reset_mpz_array(progfa);
        reset_uint32_array(progmu);
    }
    START_TIMER;
    ecode = program->factoring_algo_func(
                progfa,
                progmu,
                unfactored,
                program->params,
                program->mode
            );
    STOP_TIMER;

    if (ecode == NO_FACTOR_FOUND) {
        nfails++;
    }

    if (program->verbose || program->timing) {
        printf("\n");
    }
    printf("Averaged time: %8.6f seconds\n\n", (GET_TIMING/(float)NITERATIONS));
    printf("Failure rate : %4.1f %% \n\n",
            100 * ((float)nfails/(float)NITERATIONS));

#else
    ecode = program->factoring_algo_func(
                progfa,
                progmu,
                unfactored,
                program->params,
                program->mode
            );
    if (program->verbose || program->timing) {
        printf("\n");
    }

#endif

    switch (ecode) {
        case COMPLETE_FACTORIZATION_FOUND:
            append_mpz_array(factors, progfa);
            append_uint32_array(multis, progmu);
            rcode = COMPLETE_FACTORIZATION_FOUND;
            break;
        case SOME_FACTORS_FOUND:
        case SOME_COPRIME_FACTORS_FOUND:
        case PARTIAL_FACTORIZATION_FOUND: {
            //
            // Simplify the list of found factors by computing and
            // keeping only coprime factors.
            //
            uint32_t     blen = progfa->length;
            mpz_array_t* base = alloc_mpz_array(blen * (blen + 1));

            find_coprime_base(base, unfactored, progfa);
            ins_sort_mpz_array(base);
            append_mpz_array(factors, base);

            free_mpz_array(base);

            rcode = PARTIAL_FACTORIZATION_FOUND;
            break;
        }
        default:
            if (factors->length != 0) {
                //
                // We found some factors via trial divisions. Add the
                // unfactored part in the list of factors so that it will be
                // printed.
                //
                append_mpz_to_array(factors, unfactored);
                append_uint32_to_array(multis, 1);
            }
            break;
    }

    free_mpz_array(progfa);
    free_uint32_array(progmu);

  clear_sqt_and_return:

    mpz_clear(unfactored);

  clear_tdiv_and_return:

    print_factors(factors, multis);

    switch (rcode) {
        case COMPLETE_FACTORIZATION_FOUND:
            printf("Factorization is complete.\n\n");
            break;
        case PARTIAL_FACTORIZATION_FOUND:
            printf("Factorization is (potentially) not complete.\n\n");
            break;
        case NO_FACTOR_FOUND:
        default:
            printf("Factorization failed.\n\n");
            break;
    }
    free_mpz_array(factors);
    free_uint32_array(multis);

    return rcode;
}
//-----------------------------------------------------------------------------
void print_factors(mpz_array_t* const factors, uint32_array_t* const multis) {
    printf("Found %2u factor(s):\n", (unsigned int)factors->length);
    printf("-------------------\n");

    uint32_t nmult = 0;
    if (multis == NULL) {
        nmult = 0;
    } else {
        nmult = multis->length;
    }

    for (uint32_t i = 0; i != nmult; i++) {
        //
        // _WARNING_: The multis array can be shorter than the factors array
        //            since the factors obtained by trial divisions always come
        //            with their corresponding multiplicities. This assumes
        //            that these factors are given at the beginning of the
        //            factors array.
        //
        gmp_printf("\t%Zd ^ %"PRIu32"\n", factors->data[i], multis->data[i]);
        fflush(stdout);
    }
    for (uint32_t i = nmult; i != factors->length; i++) {
        gmp_printf("\t%Zd\n", factors->data[i]);
        fflush(stdout);
    }
    printf("\n");
}
//-----------------------------------------------------------------------------
