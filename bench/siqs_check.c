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
 * \file    siqs_check.c
 * \brief   siqs-check: factors balanced semiprimes of random sizes with siqs()
 *          and checks every result.
 */

#include "bench_utils.h"
#include "siqs.h"
#include "array.h"

static void usage(FILE* const out) {
    fprintf(out,
        "usage: siqs-check [--seed <s>] <lo> <hi> <count>\n"
        "\n"
        "Factors <count> balanced semiprimes of random sizes from <lo> to <hi>\n"
        "bits with siqs(), one after the other, and prints how many were not\n"
        "completely factored and how many were given a wrong factor. Exits with\n"
        "1 if any; --seed sets the seed of the numbers (default 1).\n"
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
    unsigned long lo    = 0;
    unsigned long hi    = 0;
    unsigned long count = 0;

    const bool ok =    (argc - i == 3)
                    && parse_ulong(argv[i], &lo)
                    && parse_ulong(argv[i + 1], &hi)
                    && parse_ulong(argv[i + 2], &count)
                    && (lo >= 20) && (lo <= hi) && (hi <= 300);
    if (!ok) {
        usage(stderr);
        return 2;
    }
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    mpz_t n, prod;
    mpz_inits(n, prod, NULL);

    mpz_array_t*    const factors = alloc_mpz_array(4);
    uint32_array_t* const multis  = alloc_uint32_array(4);

    unsigned long incomplete = 0;
    unsigned long wrong      = 0;

    for (unsigned long k = 0; k < count; k++) {
        const unsigned long bits = lo + gmp_urandomm_ui(state, hi - lo + 1);
        random_semiprime(n, state, bits, bits / 2);

        siqs_params_t params;
        set_siqs_params_to_default(n, &params);

        factors->length = 0;
        multis->length  = 0;
        const ecode_t ecode = siqs(factors, multis, n, &params, SINGLE_RUN);

        bool bad = false;
        for (uint32_t j = 0; j < factors->length; j++) {
            const mpz_srcptr f = factors->data[j];
            bad = bad || !mpz_divisible_p(n, f)
                      || (mpz_cmp_ui(f, 1) <= 0) || (mpz_cmp(f, n) >= 0);
        }
        if (ecode == COMPLETE_FACTORIZATION_FOUND) {
            mpz_set_ui(prod, 1);
            for (uint32_t j = 0; j < factors->length; j++) {
                const uint32_t e = (j < multis->length) ? multis->data[j] : 1;
                for (uint32_t r = 0; r < e; r++) {
                    mpz_mul(prod, prod, factors->data[j]);
                }
            }
            bad = bad || (mpz_cmp(prod, n) != 0);
        } else {
            incomplete++;
        }
        if (bad) {
            gmp_printf("wrong factors for %Zd\n", n);
            wrong++;
        }
    }
    printf("%lu-%lu bits: %lu numbers, not complete %lu, wrong %lu\n",
           lo, hi, count, incomplete, wrong);

    mpz_clears(n, prod, NULL);
    free_mpz_array(factors);
    free_uint32_array(multis);
    gmp_randclear(state);

    return ((incomplete == 0) && (wrong == 0)) ? EXIT_SUCCESS : EXIT_FAILURE;
}
