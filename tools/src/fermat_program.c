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
 * \file    fermat_program.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <tifa_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"
#include "tdiv.h"
#include "fermat.h"
#include "factoring_program.h"
#include "common_funcs.h"

//
// The number of arguments accepted is either 0, or FERMAT_F_MAX_ARGC-1
//
#define FERMAT_F_MAX_ARGC 3

//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s <nprimes_tdiv>", program->argv[0]);
    fprintf(stderr, " <number_to_factor>\n\n");
    fprintf(stderr, "Or:\n\n");
    fprintf(stderr, "%15s <number_to_factor>\n\n", program->argv[0]);
    PRINT_USAGE_WARNING_MSG();
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {

    int    argc = program->argc;
    char** argv = program->argv;

    fermat_params_t* params = (fermat_params_t*) program->params;
    uint32_t* nprimes_tdiv  = &(program->nprimes_tdiv);

    print_hello_msg(program->algo_name);

    switch (argc) {
    case 1: {
        //
        // No argument provided: proceed in interactive mode, i.e. get the
        // number to factor.
        //
        char str_buffer[MAX_NDIGITS];
        PRINT_ENTER_NUMBER_MSG();
        char* str_factor_me = fgets(str_buffer, MAX_NDIGITS, stdin);
        printf("\n");
        chomp(str_factor_me, MAX_NDIGITS);

        if (!is_a_number(str_factor_me, MAX_NDIGITS)) {
            PRINT_NAN_ERROR(str_factor_me);
            exit(-1);
        }
        mpz_init_set_str(program->n, str_factor_me, 10);

        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        set_fermat_params_to_default(params);

        break;
    }
    case FERMAT_F_MAX_ARGC: {
        //
        // Read parameters on the command line as they are provided by
        // the factorize.pl script... The script already checks for the validity
        // of the parameters, but let's check one more time while we're at it.
        //
        for (int i = 1; i < argc; i++) {
            if (!is_a_number(argv[i], MAX_NDIGITS)) {
                PRINT_NAN_ERROR(argv[i]);
                print_usage(program);
                exit(-1);
            }
        }
        char** endptr = NULL;
        *nprimes_tdiv = strtoul(argv[1], endptr, 10);
        mpz_init_set_str(program->n, argv[2], 10);

        break;
    }
    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(-1);
    }
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program) {
    set_fermat_params_to_default((fermat_params_t*) program->params);
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    printf("\tnprimes_tdiv : %u\n", program->nprimes_tdiv);
}
//------------------------------------------------------------------------------
static ecode_t fermat_func(mpz_array_t* const factors,
                           uint32_array_t* const multis, const mpz_t n,
                           const void* const params,
                           factoring_mode_t mode) {
    return fermat(factors, multis, n, params, mode);
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    fermat_params_t params;
    factoring_program_t program;

    program.argc = argc;
    program.argv = argv;

    program.verbose = TIFA_VERBOSE_FERMAT;
    program.timing  = TIFA_TIMING_FERMAT;

    program.algo_name = "Fermat";
    program.params    = (void*) &params;
    program.mode      = FIND_SOME_FACTORS;
    program.nfactors  = 2;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = fermat_func;
    program.set_params_to_default_func = set_params_to_default;

    ecode_t ecode = run_program(&program);

    mpz_clear(program.n);

    return ecode;
}
//------------------------------------------------------------------------------

