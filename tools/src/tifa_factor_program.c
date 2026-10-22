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
 * \file    tifa_factor_program.c
 * \author  Jerome Milan
 *
 * Complete factorization with tifa_factor, which chooses the methods and their
 * parameters.
 */

#include <tifa_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"
#include "tifa_factor.h"
#include "factoring_program.h"
#include "common_funcs.h"

//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s <number_to_factor>\n", program->argv[0]);
    fprintf(stderr, "Or, to be asked for the number:\n");
    fprintf(stderr, "%15s\n\n", program->argv[0]);
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {

    int    argc = program->argc;
    char** argv = program->argv;

    print_hello_msg(program->algo_name);

    //
    // tifa_factor does its own trial divisions.
    //
    program->nprimes_tdiv = 0;

    switch (argc) {
    case 1: {
        char str_buffer[MAX_NDIGITS];
        PRINT_ENTER_NUMBER_MSG();
        char* str_factor_me = read_number(str_buffer, MAX_NDIGITS);
        printf("\n");

        if (!is_a_number(str_factor_me, MAX_NDIGITS)) {
            PRINT_NAN_ERROR(str_factor_me);
            exit(-1);
        }
        mpz_init_set_str(program->n, str_factor_me, 10);
        break;
    }
    case 2:
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            print_usage(program);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[1], 10);
        break;

    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(-1);
    }
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program) {
    (void) program;
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    (void) program;
    printf("\tnone: tifa_factor chooses the methods and their parameters\n");
}
//------------------------------------------------------------------------------
static ecode_t tifa_factor_func(mpz_array_t* const factors,
                                uint32_array_t* const multis, const mpz_t n,
                                const void* const params,
                                factoring_mode_t mode) {
    (void) params;
    return tifa_factor(factors, multis, n, mode);
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    factoring_program_t program;

    program.argc = argc;
    program.argv = argv;

    program.verbose = false;
    program.timing  = false;

    program.algo_name = "tifa_factor";
    program.params    = NULL;
    program.mode      = FIND_COMPLETE_FACTORIZATION;
    program.nfactors  = 2;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = tifa_factor_func;
    program.set_params_to_default_func = set_params_to_default;

    ecode_t ecode = run_program(&program);

    return ecode;
}
//------------------------------------------------------------------------------
