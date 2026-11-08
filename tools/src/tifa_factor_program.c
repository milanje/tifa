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
#include <string.h>

#include "tool_utils.h"
#include "tifa_factor.h"
#include "tifa_options.h"
#include "factoring_program.h"
#include "common_funcs.h"

//------------------------------------------------------------------------------
static const char* const hint_names[] = { "random", "balanced", "unbalanced" };
//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s [--prove] [--threads <n>] "
                    "[random|balanced|unbalanced] <number_to_factor>\n",
                    program->argv[0]);
    fprintf(stderr, "Or, to be asked for the number:\n");
    fprintf(stderr, "%15s [--prove] [--threads <n>]\n\n", program->argv[0]);
    fprintf(stderr, "balanced: two prime factors of about the same size;\n");
    fprintf(stderr, "unbalanced: a prime factor below a third of the size;\n");
    fprintf(stderr, "random (the default): nothing known.\n\n");
    fprintf(stderr, "--prove: also prove the prime factors, and print their "
                    "certificates in\nJSON, as tifa.check_certificate (Python) "
                    "takes them. A factor not proven\nprime makes the "
                    "factorization partial.\n\n");
    fprintf(stderr, "--threads <n>: run SIQS and the ECM search before it in "
                    "up to n threads\n(SIQS: one per 100 primes of its factor "
                    "base at most, so none below\nabout 100 bits).\n\n");
    print_exit_statuses();
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {

    int    argc = program->argc;
    char** argv = program->argv;

    tifa_options_t* options = (tifa_options_t*) program->params;

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
            exit(TIFA_EXIT_USAGE);
        }
        mpz_init_set_str(program->n, str_factor_me, 10);
        break;
    }
    case 2:
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            print_usage(program);
            exit(TIFA_EXIT_USAGE);
        }
        mpz_init_set_str(program->n, argv[1], 10);
        break;

    case 3: {
        bool known = false;
        for (int i = TIFA_FACTOR_RANDOM; i <= TIFA_FACTOR_UNBALANCED; i++) {
            if (0 == strcmp(argv[1], hint_names[i])) {
                options->hint = (tifa_factor_hint_t) i;
                known = true;
            }
        }
        if (!known) {
            fprintf(stderr, "Unknown hint: %s\n", argv[1]);
            print_usage(program);
            exit(TIFA_EXIT_USAGE);
        }
        if (!is_a_number(argv[2], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[2]);
            print_usage(program);
            exit(TIFA_EXIT_USAGE);
        }
        mpz_init_set_str(program->n, argv[2], 10);
        break;
    }

    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(TIFA_EXIT_USAGE);
    }
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program,
                                  const mpz_t n) {
    (void) program;
    (void) n;
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    const tifa_options_t* options = (const tifa_options_t*) program->params;
    printf("\thint: %s (tifa_factor chooses the methods and their parameters)\n",
           hint_names[options->hint]);
}
//------------------------------------------------------------------------------
static ecode_t tifa_factor_func(mpz_array_t* const factors,
                                uint32_array_t* const multis, const mpz_t n,
                                const void* const params,
                                factoring_mode_t mode) {
    return tifa_factor(factors, multis, n, mode,
                       (const tifa_options_t*) params);
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    factoring_program_t program;
    tifa_options_t      options = { .hint = TIFA_FACTOR_RANDOM };

    bool prove   = false;
    bool bad_arg = false;
    int  nargs   = 1;
    for (int i = 1; i < argc; i++) {
        if (0 == strcmp(argv[i], "--prove")) {
            prove = true;
        } else if (0 == strcmp(argv[i], "--threads")) {
            char* end = NULL;
            const unsigned long n = (i + 1 < argc)
                                  ? strtoul(argv[i + 1], &end, 10) : 0;
            if ((end == NULL) || (*end != '\0') || (n < 1) || (n > 1024)
                || (argv[i + 1][0] == '-')) {
                bad_arg = true;
                break;
            }
            options.threads = (uint32_t) n;
            i++;
        } else {
            argv[nargs++] = argv[i];
        }
    }
    argv[nargs] = NULL;

    program.argc = nargs;
    program.argv = argv;

    program.verbose = false;
    program.timing  = false;
    program.prove   = prove;
    prove_options   = &options;

    program.algo_name = "tifa_factor()";
    program.params    = &options;
    program.mode      = FIND_COMPLETE_FACTORIZATION;
    program.nfactors  = 2;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = tifa_factor_func;
    program.set_params_to_default_func = set_params_to_default;

    if (bad_arg) {
        fprintf(stderr, "--threads takes a number of threads, from 1 to "
                        "1024.\n\n");
        print_usage(&program);
        return TIFA_EXIT_USAGE;
    }
    const int status = run_program(&program);

    mpz_clear(program.n);

    return status;
}
//------------------------------------------------------------------------------
