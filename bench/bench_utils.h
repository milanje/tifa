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
 * \file    bench_utils.h
 * \brief   Shared by the benchmark programs: CPU time, random semiprimes and
 *          argument parsing.
 */

#if !defined(_TIFA_BENCH_UTILS_H_)
#define _TIFA_BENCH_UTILS_H_

//
// For clock_gettime and CLOCK_PROCESS_CPUTIME_ID in C99: included first.
//
#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

static inline double cpu_seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

static inline void random_prime(mpz_t p, gmp_randstate_t state,
                                const unsigned long bits) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);
    } while (!mpz_probab_prime_p(p, 30));
}

//
// n = p.q of exactly bits bits, p of pbits bits (at least 2 on each side).
//
static inline void random_semiprime(mpz_t n, gmp_randstate_t state,
                                    const unsigned long bits,
                                    const unsigned long pbits) {
    mpz_t p, q;
    mpz_inits(p, q, NULL);
    do {
        random_prime(p, state, pbits);
        random_prime(q, state, bits - pbits);
        mpz_mul(n, p, q);
    } while (mpz_sizeinbase(n, 2) != bits);
    mpz_clears(p, q, NULL);
}

static inline bool parse_ulong(const char* const s, unsigned long* const v) {
    char* end = NULL;
    if ((s[0] < '0') || (s[0] > '9')) {
        return false;
    }
    *v = strtoul(s, &end, 10);
    return (*end == '\0');
}

static inline bool parse_long(const char* const s, long* const v) {
    char* end = NULL;
    *v = strtol(s, &end, 10);
    return (end != s) && (*end == '\0');
}

static inline bool parse_double(const char* const s, double* const v) {
    char* end = NULL;
    *v = strtod(s, &end);
    return (end != s) && (*end == '\0');
}

static inline bool is_help(const char* const s) {
    return (0 == strcmp(s, "-h")) || (0 == strcmp(s, "--help"));
}

#endif
