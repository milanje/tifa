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
 * \brief   factor-bench: mean time of tifa_factor(), ecm(), p_minus_1(),
 *          p_plus_1() or squfof() on numbers n = p.q with p of a given size.
 */

#include "bench_utils.h"
#include "array.h"
#include "ecm.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "squfof.h"
#include "tifa_factor.h"
#include "time_limit.h"

static void usage(FILE* const out) {
    fprintf(out,
        "usage: factor-bench [--seed <s>] [--threads <n>] <algo> <bits> <pbits>\n"
        "                    [<count> [<seconds>]]\n"
        "\n"
        "Times tifa_factor() (a complete factorization), ecm() or squfof()\n"
        "(some factors), p_minus_1() or p_plus_1() (a single run) on <count>\n"
        "numbers n = p.q of <bits> bits, p of <pbits> bits (or random\n"
        "composites of <bits> bits if <pbits> is 0), run again and again for at\n"
        "least <seconds> of CPU time, and prints the mean time per run and the\n"
        "number of runs that found no factor. Defaults: 30 numbers, 3 seconds;\n"
        "--seed sets the seed of the numbers (default 1). <algo> is factor,\n"
        "factor-balanced, factor-unbalanced (tifa_factor_with_hint()), ecm,\n"
        "pm1, pp1 or squfof. --threads runs tifa_factor() in up to <n> threads\n"
        "(tifa_set_threads()), timed by the wall clock.\n"
    );
}

int main(int argc, char** argv) {
    unsigned long seed = 1;

    if ((argc > 1) && is_help(argv[1])) {
        usage(stdout);
        return EXIT_SUCCESS;
    }
    unsigned long threads = 1;
    int i = 1;
    while ((i + 1 < argc) && (0 == strncmp(argv[i], "--", 2))) {
        const bool ok = (0 == strcmp(argv[i], "--seed"))
                      ? parse_ulong(argv[i + 1], &seed)
                      :    (0 == strcmp(argv[i], "--threads"))
                        && parse_ulong(argv[i + 1], &threads)
                        && (threads >= 1) && (threads <= 1024);
        if (!ok) {
            usage(stderr);
            return 2;
        }
        i += 2;
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

    tifa_factor_hint_t hint = TIFA_FACTOR_RANDOM;
    if (0 == strcmp(algo, "factor-balanced")) {
        hint = TIFA_FACTOR_BALANCED;
    } else if (0 == strcmp(algo, "factor-unbalanced")) {
        hint = TIFA_FACTOR_UNBALANCED;
    }
    bool ok =    (   (0 == strcmp(algo, "factor"))
                  || (hint != TIFA_FACTOR_RANDOM) || (0 == strcmp(algo, "ecm"))
                  || (0 == strcmp(algo, "pm1")) || (0 == strcmp(algo, "pp1"))
                  || (0 == strcmp(algo, "squfof")))
              && parse_ulong(argv[i + 1], &bits)
              && parse_ulong(argv[i + 2], &pbits)
              && ((pbits == 0) ? (bits >= 4)
                               : ((pbits >= 2) && (bits >= pbits + 2)))
              && (bits <= 1024);
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
        if (pbits == 0) {
            random_composite(ns[k], state, bits);
        } else {
            random_semiprime(ns[k], state, bits, pbits);
        }
    }
    mpz_array_t*    const factors = alloc_mpz_array(8);
    uint32_array_t* const multis  = alloc_uint32_array(8);

    tifa_set_threads((uint32_t) threads);
    double (*const now)(void) = (threads > 1) ? wall_seconds : cpu_seconds;

    unsigned long runs  = 0;
    unsigned long fails = 0;
    const double t0 = now();
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
            } else if (0 == strcmp(algo, "squfof")) {
                squfof_params_t params;
                set_squfof_params_to_default(&params);
                squfof(factors, NULL, ns[k], &params, FIND_SOME_FACTORS);
            } else {
                tifa_factor_with_hint(factors, multis, ns[k],
                                      FIND_COMPLETE_FACTORIZATION, hint);
            }
            if (factors->length == 0) {
                fails++;
            }
            runs++;
        }
    } while (now() - t0 < secs);

    printf("%lu bits: %.2f us per run (%lu of %lu runs without a factor)\n",
           bits, 1e6 * (now() - t0) / (double) runs, fails, runs);

    for (unsigned long k = 0; k < count; k++) {
        mpz_clear(ns[k]);
    }
    free(ns);
    free_mpz_array(factors);
    free_uint32_array(multis);
    gmp_randclear(state);

    return EXIT_SUCCESS;
}
