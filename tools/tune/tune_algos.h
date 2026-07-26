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
 * \file    tune_algos.h
 * \brief   Registry describing how to tune each factoring algorithm.
 *
 * Each tunable algorithm is described by one \c tune_algo_desc_t entry giving
 * its tuning-file schema, its parameter grid, and how to run a single timed
 * attempt. Adding an algorithm to the tuner is a matter of adding one entry to
 * the \c TUNE_ALGOS array; the tuner's driver code is algorithm-agnostic.
 */

#if !defined(_TIFA_TUNE_ALGOS_H_)
#define _TIFA_TUNE_ALGOS_H_

#include <inttypes.h>
#include <stdbool.h>

#include <gmp.h>

#include "tune_grid.h"

typedef struct {
        /** Algorithm name; also the tuning file base name and lookup key. */
    const char* name;
        /** Number of tunable columns (matches TIFA_TUNING_NCOLS_*). */
    uint32_t    ncols;
        /** Space-separated column names (matches TIFA_TUNING_COLS_*). */
    const char* colnames;
        /**
         * true for factor-finding methods whose success is probabilistic
         * (scored by expected time); false for the sieve methods.
         */
    bool        probabilistic;
        /** Smallest sensible composite bit-size to tune this algorithm at. */
    uint32_t    min_bits;
        /** Largest sensible composite bit-size to tune this algorithm at. */
    uint32_t    max_bits;
        /**
         * Number of similar-size prime factors in the test composites at or
         * above \c multi_factor_min_bits. Factor-finding methods use 3 (their
         * cost tracks the smallest factor, so a 3-prime composite better models
         * "find a smaller factor of n"); sieve methods use 2 (their cost tracks
         * the size of n itself).
         */
    uint32_t    max_factors;
        /**
         * Bit-size at or above which \c max_factors primes are used; below it a
         * balanced 2-prime composite is used instead, so tiny numbers are not
         * tuned on degenerate, trivially-small factors. Ignored (no effect)
         * when \c max_factors == 2.
         */
    uint32_t    multi_factor_min_bits;
        /** Builds the parameter grid to sweep for a given bit-size. */
    void (*build_grid)(uint32_t bitsize, grid_t* out);
        /**
         * Runs one timed attempt with the given column values on the composite
         * \c n (whose smallest prime factor is \c p), storing the elapsed CPU
         * seconds in \c cpu_seconds_out. Returns true iff a non-trivial factor
         * was found.
         */
    bool (*run_once)(const mpz_t n, const mpz_t p, const uint32_t cols[],
                     double* cpu_seconds_out);
        /**
         * Column index of a "sticky-on" boolean parameter (e.g. the large-prime
         * variation): once its value 1 wins at some bit-size, it is locked to 1
         * for all larger sizes, so the tuner stops re-testing the "off" case
         * where it can no longer win. -1 if the algorithm has no such column.
         */
    int sticky_col;
} tune_algo_desc_t;

extern const tune_algo_desc_t TUNE_ALGOS[];
extern const unsigned int     TUNE_NALGOS;

   /**
    * \brief Returns the descriptor for \c name, or NULL if unknown.
    */
const tune_algo_desc_t* tune_find_algo(const char* name);

#endif
