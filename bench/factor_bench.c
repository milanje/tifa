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
 * \file    factor_bench.c
 * \brief   factor-bench: mean time of tifa_factor(), ecm(), p_minus_1() or
 *          p_plus_1() on numbers n = p.q with p of a given size.
 */

#include "bench_utils.h"
#include "array.h"
#include "ecm.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "tifa_factor.h"

static void usage(FILE* const out) {
    fprintf(out,
        "usage: factor-bench [--seed <s>] <factor|ecm|pm1|pp1> <bits> <pbits>\n"
        "                    [<count> [<seconds>]]\n"
        "\n"
        "Times tifa_factor() (a complete factorization), ecm() (some factors),\n"
        "p_minus_1() or p_plus_1() (a single run) on <count> numbers n = p.q\n"
        "of <bits> bits, p of <pbits> bits, run again and again for at least\n"
        "<seconds> of CPU time, and prints the mean time per run. Defaults: 30\n"
        "numbers, 3 seconds; --seed sets the seed of the numbers (default 1).\n"
    );
}

int main(int argc, char** argv) {
    unsigned long seed = 1;

    if ((argc > 1) && is_help(argv[1])) {
        usage(stdout);
        return EXIT_SUCCESS;
    }
    int i = 1;
    if ((argc > 2) && (0 == strcmp(argv[1], "--seed"))) {
        if (!parse_ulong(argv[2], &seed)) {
            usage(stderr);
            return 2;
        }
        i = 3;
    }
    if ((argc - i < 3) || (argc - i > 5)) {
        usage(stderr);
        return 2;
    }
    const char* const algo = argv[i];
    unsigned long bits  = 0;
    unsigned long pbits = 0;
    unsigned long count = 30;
    double secs = 3.0;

    bool ok =    (   (0 == strcmp(algo, "factor")) || (0 == strcmp(algo, "ecm"))
                  || (0 == strcmp(algo, "pm1")) || (0 == strcmp(algo, "pp1")))
              && parse_ulong(argv[i + 1], &bits)
              && parse_ulong(argv[i + 2], &pbits)
              && (pbits >= 2) && (bits >= pbits + 2) && (bits <= 1024);
    if (ok && (argc - i > 3)) {
        ok = parse_ulong(argv[i + 3], &count) && (count > 0);
    }
    if (ok && (argc - i > 4)) {
        ok = parse_double(argv[i + 4], &secs) && (secs >= 0);
    }
    if (!ok) {
        usage(stderr);
        return 2;
    }
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    mpz_t* const ns = malloc(count * sizeof(mpz_t));
    for (unsigned long k = 0; k < count; k++) {
        mpz_init(ns[k]);
        random_semiprime(ns[k], state, bits, pbits);
    }
    mpz_array_t*    const factors = alloc_mpz_array(8);
    uint32_array_t* const multis  = alloc_uint32_array(8);

    unsigned long runs = 0;
    const double t0 = cpu_seconds();
    do {
        for (unsigned long k = 0; k < count; k++) {
            factors->length = 0;
            multis->length  = 0;
            if (0 == strcmp(algo, "ecm")) {
                ecm_params_t params;
                set_ecm_params_to_default(ns[k], &params);
                ecm(factors, NULL, ns[k], &params, FIND_SOME_FACTORS);
            } else if (0 == strcmp(algo, "pm1")) {
                p_minus_1_params_t params;
                set_p_minus_1_params_to_default(ns[k], &params);
                p_minus_1(factors, NULL, ns[k], &params, SINGLE_RUN);
            } else if (0 == strcmp(algo, "pp1")) {
                p_plus_1_params_t params;
                set_p_plus_1_params_to_default(ns[k], &params);
                p_plus_1(factors, NULL, ns[k], &params, SINGLE_RUN);
            } else {
                tifa_factor(factors, multis, ns[k],
                            FIND_COMPLETE_FACTORIZATION);
            }
            runs++;
        }
    } while (cpu_seconds() - t0 < secs);

    printf("%lu bits: %.2f us per run\n", bits,
           1e6 * (cpu_seconds() - t0) / (double) runs);

    for (unsigned long k = 0; k < count; k++) {
        mpz_clear(ns[k]);
    }
    free(ns);
    free_mpz_array(factors);
    free_uint32_array(multis);
    gmp_randclear(state);

    return EXIT_SUCCESS;
}
