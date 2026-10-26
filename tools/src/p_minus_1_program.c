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

#include <tifa_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>
#include <time.h>

#include "tool_utils.h"
#include "tdiv.h"
#include "p_minus_1.h"
#include "factoring_program.h"
#include "common_funcs.h"

//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s\n", program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\t<number_to_factor>\n", program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\t<nprimes_tdiv> <number_to_factor>\n",
            program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s <nprimes_tdiv> <b1> <b2> <base> <number_to_factor>\n\n",
            program->argv[0]);
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {

    p_minus_1_params_t* params = (p_minus_1_params_t*) program->params;

    int    argc = program->argc;
    char** argv = program->argv;

    uint32_t* nprimes_tdiv = &(program->nprimes_tdiv);
    uint32_t* nfactors     = &(program->nfactors);

    print_hello_msg(program->algo_name);

    switch (argc) {
    case 1: {
        //
        // No argument provided: proceed in interactive mode, i.e. get the
        // number to factor and then use default values for the p-1 parameters.
        //
        char str_buffer[MAX_NDIGITS];
        PRINT_ENTER_NUMBER_MSG();
        char* str_factor_me = read_number(str_buffer, MAX_NDIGITS);
        printf("\n");

        if (!is_a_number(str_factor_me, MAX_NDIGITS)) {
            PRINT_NAN_ERROR(str_factor_me);
            exit(-1);
        }
        mpz_init_set_str(program->n, str_factor_me, 10);
        set_p_minus_1_params_to_default(program->n, params);
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = 2;

        break;
    }
    case 2: {
        //
        // Only one argument provided: the number to factor. Use default values
        // for the p-1 parameters.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[1], 10);
        set_p_minus_1_params_to_default(program->n, params);
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = 2;

        break;
    }
    case 3: {
        //
        // Two arguments provided: the <nprimes_tdiv> parameter and the number
        // to factor. Use default values for the p-1 parameters.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        *nprimes_tdiv = strtoul(argv[1], NULL, 10);

        if (!is_a_number(argv[2], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[2]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[2], 10);
        set_p_minus_1_params_to_default(program->n, params);
        *nfactors = 2;
        break;
    }
    case 6: {
        //
        // Five arguments provided: the <nprimes_tdiv> parameter, the stage-1
        // bound <b1>, the stage-2 bound <b2>, the <base>, and finally the
        // number to factor.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        *nprimes_tdiv = strtoul(argv[1], NULL, 10);

        if (!is_a_number(argv[2], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[2]);
            exit(-1);
        }
        params->b1 = strtoul(argv[2], NULL, 10);

        if (!is_a_number(argv[3], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[3]);
            exit(-1);
        }
        params->b2 = strtoul(argv[3], NULL, 10);

        if (!is_a_number(argv[4], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[4]);
            exit(-1);
        }
        params->base = strtoul(argv[4], NULL, 10);

        if (!is_a_number(argv[5], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[5]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[5], 10);

        *nfactors = 2;
        break;
    }
    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(-1);
    }
    return;
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program,
                                  const mpz_t n) {
    set_p_minus_1_params_to_default(n,
                                    (p_minus_1_params_t*) program->params);
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    p_minus_1_params_t* params = (p_minus_1_params_t*) program->params;
    printf("\t  b1 : %"PRIu32"\n", params->b1);
    printf("\t  b2 : %"PRIu32"\n", params->b2);
    printf("\tbase : %"PRIu32"\n", params->base);
}
//------------------------------------------------------------------------------
static ecode_t p_minus_1_func(mpz_array_t* const factors,
                              uint32_array_t* const multis, const mpz_t n,
                              const void* const params, factoring_mode_t mode) {
    return p_minus_1(factors, multis, n,
                     (const p_minus_1_params_t*) params, mode);
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    p_minus_1_params_t params;
    factoring_program_t program;

    program.argc = argc;
    program.argv = argv;

    program.verbose = TIFA_VERBOSE_PM1;
    program.timing  = TIFA_TIMING_PM1;

    program.algo_name = "Pollard's p-1";
    program.params    = (void*) &params;
    program.mode      = FIND_SOME_FACTORS;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = p_minus_1_func;
    program.set_params_to_default_func = set_params_to_default;

    ecode_t ecode = run_program(&program);

    mpz_clear(program.n);

    return ecode;
}
//------------------------------------------------------------------------------
