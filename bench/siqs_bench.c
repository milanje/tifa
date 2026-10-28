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
 * \file    siqs_bench.c
 * \brief   siqs-bench: mean time of siqs() on balanced semiprimes of a size or
 *          of a range of sizes, with its parameters scaled or shifted.
 */

#include "bench_utils.h"
#include "siqs.h"
#include "array.h"

static void usage(FILE* const out) {
    fprintf(out,
        "usage: siqs-bench [options] <bits>|<lo>-<hi> [<count> [<seconds>]]\n"
        "\n"
        "Times siqs() on <count> balanced semiprimes of <bits> bits, or of\n"
        "random sizes from <lo> to <hi> bits, factored again and again for at\n"
        "least <seconds> of CPU time, and prints the mean time per run.\n"
        "Defaults: 30 numbers, 3 seconds.\n"
        "\n"
        "options:\n"
        "  --seed <s>             seed of the numbers (default 1)\n"
        "  --fb-mul <x>           factor base (and trial division) times x\n"
        "  --hw-mul <x>           sieve half width times x\n"
        "  --threshold-delta <d>  sieve threshold plus d\n"
        "  --large-primes <0|1>   large prime variation off or on\n"
        "  --nrelations <n>       relations beyond the size of the factor base\n"
    );
}

int main(int argc, char** argv) {
    unsigned long seed = 1;
    double fbmul = 1.0;
    double hwmul = 1.0;
    long thdelta = 0;
    long large_primes = -1;
    long nrelations = -1;

    if ((argc > 1) && is_help(argv[1])) {
        usage(stdout);
        return EXIT_SUCCESS;
    }
    int i = 1;
    while ((i < argc) && (0 == strncmp(argv[i], "--", 2))) {
        if (is_help(argv[i])) {
            usage(stdout);
            return EXIT_SUCCESS;
        }
        unsigned long u = 0;
        const char* const opt = argv[i];
        const char* const val = (i + 1 < argc) ? argv[i + 1] : "";
        bool ok;
        if (0 == strcmp(opt, "--seed")) {
            ok = parse_ulong(val, &seed);
        } else if (0 == strcmp(opt, "--fb-mul")) {
            ok = parse_double(val, &fbmul) && (fbmul > 0);
        } else if (0 == strcmp(opt, "--hw-mul")) {
            ok = parse_double(val, &hwmul) && (hwmul > 0);
        } else if (0 == strcmp(opt, "--threshold-delta")) {
            ok = parse_long(val, &thdelta);
        } else if (0 == strcmp(opt, "--large-primes")) {
            ok = parse_ulong(val, &u) && (u <= 1);
            large_primes = (long) u;
        } else if (0 == strcmp(opt, "--nrelations")) {
            ok = parse_ulong(val, &u) && (u <= 1024);
            nrelations = (long) u;
        } else {
            ok = false;
        }
        if (!ok) {
            usage(stderr);
            return 2;
        }
        i += 2;
    }
    if ((argc - i < 1) || (argc - i > 3)) {
        usage(stderr);
        return 2;
    }
    unsigned long lo = 0;
    unsigned long hi = 0;
    unsigned long count = 30;
    double secs = 3.0;

    const char* const dash = strchr(argv[i], '-');
    bool ok;
    if (dash != NULL) {
        char lostr[32];
        const size_t len = (size_t) (dash - argv[i]);
        ok = (len < sizeof(lostr));
        if (ok) {
            memcpy(lostr, argv[i], len);
            lostr[len] = '\0';
            ok = parse_ulong(lostr, &lo) && parse_ulong(dash + 1, &hi);
        }
    } else {
        ok = parse_ulong(argv[i], &lo);
        hi = lo;
    }
    ok = ok && (lo >= 20) && (lo <= hi) && (hi <= 300);
    if (ok && (argc - i > 1)) {
        ok = parse_ulong(argv[i + 1], &count) && (count > 0);
    }
    if (ok && (argc - i > 2)) {
        ok = parse_double(argv[i + 2], &secs) && (secs >= 0);
    }
    if (!ok) {
        usage(stderr);
        return 2;
    }
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    mpz_t*         const ns  = malloc(count * sizeof(mpz_t));
    siqs_params_t* const par = malloc(count * sizeof(siqs_params_t));

    for (unsigned long k = 0; k < count; k++) {
        const unsigned long bits = lo + gmp_urandomm_ui(state, hi - lo + 1);
        mpz_init(ns[k]);
        random_semiprime(ns[k], state, bits, bits / 2);

        set_siqs_params_to_default(ns[k], &par[k]);
        par[k].nprimes_in_base = (uint32_t) (par[k].nprimes_in_base * fbmul
                                             + 0.5);
        par[k].nprimes_tdiv    = par[k].nprimes_in_base;
        par[k].sieve_half_width = (uint32_t) (par[k].sieve_half_width * hwmul
                                              + 0.5);
        par[k].threshold = (uint32_t) ((long) par[k].threshold + thdelta);
        if (large_primes >= 0) {
            par[k].use_large_primes = (large_primes == 1);
        }
        if (nrelations >= 0) {
            par[k].nrelations = (uint32_t) nrelations;
        }
    }
    mpz_array_t*    const factors = alloc_mpz_array(4);
    uint32_array_t* const multis  = alloc_uint32_array(4);

    unsigned long runs = 0;
    const double t0 = cpu_seconds();
    do {
        for (unsigned long k = 0; k < count; k++) {
            factors->length = 0;
            multis->length  = 0;
            siqs(factors, multis, ns[k], &par[k], SINGLE_RUN);
            runs++;
        }
    } while (cpu_seconds() - t0 < secs);

    const double us = 1e6 * (cpu_seconds() - t0) / (double) runs;
    if (lo == hi) {
        printf("%lu bits: %.2f us per run\n", lo, us);
    } else {
        printf("%lu-%lu bits: %.2f us per run\n", lo, hi, us);
    }
    for (unsigned long k = 0; k < count; k++) {
        mpz_clear(ns[k]);
    }
    free(ns);
    free(par);
    free_mpz_array(factors);
    free_uint32_array(multis);
    gmp_randclear(state);

    return EXIT_SUCCESS;
}
