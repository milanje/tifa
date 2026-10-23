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
 * \file    tune_algos.c
 *
 * The grids here are deliberately modest multiplicative sweeps around a size-
 * dependent center. They are meant to be a sane starting point, not the last
 * word: widen a grid (or the bit-size range and trial count) for a more
 * thorough search. Non-tuned parameters are taken from each algorithm's own
 * set_<algo>_params_to_default so every attempt runs with a fully valid
 * parameter set.
 */

#include <string.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "array.h"
#include "factoring_machine.h"
#include "stopwatch.h"
#include "tuning.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "rho.h"
#include "ecm.h"
#include "cfrac.h"
#include "siqs.h"

#include "tune_algos.h"

// Below this bit-size the factor-finding methods are tuned on 2 balanced primes
// (a 3-prime composite would have degenerate, trivially-small factors); at or
// above it they use 3 similar-size primes, which better models finding a
// smaller factor of a larger n. Used both by the descriptors and by the rho
// grid (whose cost depends on the smallest factor's size).
#define TUNE_MULTI_FACTOR_MIN_BITS 64

// Common multiplicative ratio sets.
static const double R3[]     = { 0.5, 1.0, 2.0 };
static const double R3_B2[]  = { 5.0, 10.0, 20.0 };
// Ratios for the p-1/p+1 stage-1 bound b1. Tuning showed the optimum lives in
// ~50-800*bitsize (small-factor floor ~50-100, rising to ~800 only where p-1
// meets larger factors), so these are anchored on a 200*bitsize center to
// bracket exactly that range in 5 points -- {50, 100, 200, 400, 800} * bitsize
// -- without the wasted high point the earlier exploratory grid carried.
#define TUNE_B1_CENTER_MULT 200.0
static const double R_B1[] = { 0.25, 0.5, 1.0, 2.0, 4.0 };
// Upward-only ratios for CFRAC's factor base: too small a base starves
// relation collection and CFRAC, which has no stall guard, churns indefinitely
// on some inputs.
static const double R_UP[] = { 1.0, 1.5, 2.0 };
// Both ways for SIQS, which gives up on a stalled collection and retries wider.
// The default comes first, so the tuner prunes against it and keeps it on ties.
static const double R_SIQS[] = { 1.0, 0.5, 0.7, 1.4, 2.0 };
// Boolean axis (off, then on): "off" first so it wins exact ties, i.e. the
// large-prime variation is only preferred when it is strictly better.
static const uint32_t BOOL01[] = { 0, 1 };

//
// Builds a representative odd number of the given bit-size, used only to query
// an algorithm's own set_<algo>_params_to_default for a sane grid center.
//
static void representative_n(mpz_t n, uint32_t bitsize) {
    mpz_set_ui(n, 0UL);
    mpz_setbit(n, bitsize - 1);   // most-significant bit -> exact bit-size
    mpz_setbit(n, 0);             // make it odd
}

//
// Success oracle shared by every run_once: a non-trivial divisor of n was
// returned (1 < f < n).
//
static bool found_factor(const mpz_array_t* factors, const mpz_t n) {
    for (uint32_t i = 0; i < factors->length; i++) {
        if ((mpz_cmp_ui(factors->data[i], 1UL) > 0)
            && (mpz_cmp(factors->data[i], n) < 0)) {
            return true;
        }
    }
    return false;
}

//
// Times a single run of any algorithm (already given its filled params) and
// reports whether a factor was found. Centralizes the array/stopwatch dance.
//
#define TIME_RUN(runner_call, secs_out, ok_out)               \
    do {                                                      \
        mpz_array_t*    _f = alloc_mpz_array(4);              \
        uint32_array_t* _m = alloc_uint32_array(4);           \
        stopwatch_t _sw;                                      \
        init_stopwatch(&_sw);                                 \
        reset_stopwatch(&_sw);                                \
        start_stopwatch(&_sw);                                \
        (void) (runner_call);                                 \
        stop_stopwatch(&_sw);                                 \
        *(secs_out) = get_stopwatch_elapsed(&_sw);            \
        *(ok_out)   = found_factor(_f, n);                    \
        free_mpz_array(_f);                                   \
        free_uint32_array(_m);                                \
    } while (0)

//------------------------------------------------------------------------------
//                              Pollard p-1
//------------------------------------------------------------------------------

static void build_grid_p_minus_1(uint32_t bitsize, grid_t* g) {
    grid_clear(g);
    const double b1_center = TUNE_B1_CENTER_MULT * (double) bitsize;
    const double b2_center = 1000.0 * (double) bitsize;
    grid_add_axis_mul(g, b1_center, R_B1,  5, 100, LARGEST_PRIME);   // b1
    grid_add_axis_mul(g, b2_center, R3_B2, 3, 100, LARGEST_PRIME);   // b2
    grid_add_axis_fixed(g, 2);                                       // base
}

static bool run_p_minus_1(const mpz_t n, const mpz_t p, const uint32_t cols[],
                          double* secs) {
    (void) p;
    p_minus_1_params_t params;
    params.b1   = (cols[0] > LARGEST_PRIME) ? LARGEST_PRIME : cols[0];
    params.b2   = (cols[1] > LARGEST_PRIME) ? LARGEST_PRIME : cols[1];
    params.base = cols[2];

    bool ok;
    TIME_RUN(p_minus_1(_f, _m, n, &params, SINGLE_RUN), secs, &ok);
    return ok;
}

//------------------------------------------------------------------------------
//                              Williams p+1
//------------------------------------------------------------------------------

static void build_grid_p_plus_1(uint32_t bitsize, grid_t* g) {
    grid_clear(g);
    const double b1_center = TUNE_B1_CENTER_MULT * (double) bitsize;
    const double b2_center = 1000.0 * (double) bitsize;
    static const uint32_t seeds[] = { 6, 5 };
    grid_add_axis_mul(g, b1_center, R_B1,  5, 100, LARGEST_PRIME);   // b1
    grid_add_axis_mul(g, b2_center, R3_B2, 3, 100, LARGEST_PRIME);   // b2
    grid_add_axis_list(g, seeds, 2);                                // seed
}

static bool run_p_plus_1(const mpz_t n, const mpz_t p, const uint32_t cols[],
                         double* secs) {
    (void) p;
    p_plus_1_params_t params;
    params.b1   = (cols[0] > LARGEST_PRIME) ? LARGEST_PRIME : cols[0];
    params.b2   = (cols[1] > LARGEST_PRIME) ? LARGEST_PRIME : cols[1];
    params.seed = cols[2];

    bool ok;
    TIME_RUN(p_plus_1(_f, _m, n, &params, SINGLE_RUN), secs, &ok);
    return ok;
}

//------------------------------------------------------------------------------
//                              Pollard rho
//------------------------------------------------------------------------------

static void build_grid_rho(uint32_t bitsize, grid_t* g) {
    grid_clear(g);
    // The dominant knob is max_niters, which must scale with sqrt of the
    // smallest factor. That factor is ~bitsize/3 for the 3-prime composites
    // used at/above the crossover, ~bitsize/2 below it -- so size the grid to
    // the actual factor, not to bitsize (which previously pushed the center
    // past the ceiling and collapsed the sweep to a single value).
    unsigned int factor_bits = (bitsize >= TUNE_MULTI_FACTOR_MIN_BITS)
                             ? (bitsize / 3) : (bitsize / 2);
    unsigned int shift = factor_bits / 2;      // ~ log2(sqrt(factor))
    if (shift > 28) {
        shift = 28;
    }
    const double center_max = 4.0 * (double) (1UL << shift);   // ~4*sqrt(factor)

    static const uint32_t a_vals[] = { 1 };
    static const uint32_t b_vals[] = { 1, 2 };
    static const uint32_t acc_vals[] = { 100 };

    grid_add_axis_list(g, a_vals,   1);                             // a
    grid_add_axis_list(g, b_vals,   2);                             // b
    grid_add_axis_list(g, acc_vals, 1);                             // acc_niters
    grid_add_axis_mul(g, center_max, R3, 3, 1000, 200000000);       // max_niters
}

static bool run_rho(const mpz_t n, const mpz_t p, const uint32_t cols[],
                    double* secs) {
    (void) p;
    rho_params_t params;
    set_rho_params_to_default(n, &params);   // fills method (not tuned)
    params.a          = cols[0];
    params.b          = cols[1];
    params.acc_niters = cols[2];
    params.max_niters = cols[3];

    bool ok;
    TIME_RUN(rho(_f, _m, n, &params, SINGLE_RUN), secs, &ok);
    return ok;
}

//------------------------------------------------------------------------------
//                                  ECM
//------------------------------------------------------------------------------

static void build_grid_ecm(uint32_t bitsize, grid_t* g) {
    grid_clear(g);
    // Center on ECM's bounds for the smallest factor of the test composites:
    // ~bitsize/3 for the 3-prime ones used at/above the crossover, ~bitsize/2
    // below it.
    unsigned int factor_bits = (bitsize >= TUNE_MULTI_FACTOR_MIN_BITS)
                             ? (bitsize / 3) : (bitsize / 2);
    ecm_params_t d;
    memset(&d, 0, sizeof(d));
    set_ecm_params_for_factor_size(factor_bits, &d);

    grid_add_axis_mul(g, (double) d.b1,      R3, 3, 100, LARGEST_PRIME);  // b1
    grid_add_axis_mul(g, (double) d.b2,      R3, 3, 200, LARGEST_PRIME);  // b2
    grid_add_axis_mul(g, (double) d.ncurves, R3, 3, 1,   100000);         // ncurves
}

static bool run_ecm(const mpz_t n, const mpz_t p, const uint32_t cols[],
                    double* secs) {
    (void) p;
    ecm_params_t params;
    memset(&params, 0, sizeof(params));       // seed is set internally by ECM
    params.b1      = (cols[0] > LARGEST_PRIME) ? LARGEST_PRIME : cols[0];
    params.b2      = (cols[1] > LARGEST_PRIME) ? LARGEST_PRIME : cols[1];
    params.ncurves = cols[2];
    // ECM requires b2 == 0 or b2 > 105.
    if ((params.b2 != 0) && (params.b2 <= 105)) {
        params.b2 = 0;
    }

    bool ok;
    TIME_RUN(ecm(_f, _m, n, &params, SINGLE_RUN), secs, &ok);
    return ok;
}

//------------------------------------------------------------------------------
//                                 CFRAC
//------------------------------------------------------------------------------

static void build_grid_cfrac(uint32_t bitsize, grid_t* g) {
    grid_clear(g);
    // Center on CFRAC's own size-dependent factor-base default; sweep only the
    // factor-base size and keep the rest at their (sane) defaults.
    mpz_t n;
    mpz_init(n);
    representative_n(n, bitsize);
    cfrac_params_t d;
    set_cfrac_params_to_default(n, &d);
    mpz_clear(n);

    grid_add_axis_mul(g, (double) d.nprimes_in_base, R_UP, 3, 32, NFIRST_PRIMES); // nprimes_in_base
    grid_add_axis_fixed(g, d.nprimes_tdiv);                                       // nprimes_tdiv
    grid_add_axis_fixed(g, d.nrelations);                                         // nrelations
    grid_add_axis_list(g, BOOL01, 2);                                             // use_large_primes
}

static bool run_cfrac(const mpz_t n, const mpz_t p, const uint32_t cols[],
                      double* secs) {
    (void) p;
    cfrac_params_t params;
    set_cfrac_params_to_default(n, &params);   // fills linalg_method, filter_method
    params.nprimes_in_base    = cols[0];
    params.nprimes_tdiv       = cols[1];
    params.nrelations         = cols[2];
    params.use_large_primes   = (cols[3] != 0);

    bool ok;
    TIME_RUN(cfrac(_f, _m, n, &params, SINGLE_RUN), secs, &ok);
    return ok;
}

//------------------------------------------------------------------------------
//                                  SIQS
//------------------------------------------------------------------------------

static void build_grid_siqs(uint32_t bitsize, grid_t* g) {
    grid_clear(g);
    // Center on SIQS's own size-dependent defaults; crucially this uses the
    // algorithm's sane sieve threshold rather than a hardcoded one (a too-low
    // threshold makes sieving pathologically slow on small inputs).
    mpz_t n;
    mpz_init(n);
    representative_n(n, bitsize);
    siqs_params_t d;
    set_siqs_params_to_default(n, &d);
    mpz_clear(n);

    grid_add_axis_mul(g, (double) d.nprimes_in_base,  R_SIQS, 5, 32, NFIRST_PRIMES); // nprimes_in_base
    grid_add_axis_mul(g, (double) d.sieve_half_width, R_SIQS, 5, 1000, 2000000);     // sieve_half_width
    grid_add_axis_fixed(g, d.threshold);                                           // threshold
    grid_add_axis_fixed(g, d.nprimes_tdiv);                                         // nprimes_tdiv
    grid_add_axis_fixed(g, d.nrelations);                                           // nrelations
    grid_add_axis_list(g, BOOL01, 2);                                               // use_large_primes
}

static bool run_siqs(const mpz_t n, const mpz_t p, const uint32_t cols[],
                     double* secs) {
    (void) p;
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);    // fills linalg_method
    params.nprimes_in_base  = cols[0];
    params.sieve_half_width = cols[1];
    params.threshold        = cols[2];
    params.nprimes_tdiv     = cols[3];
    params.nrelations       = cols[4];
    params.use_large_primes = (cols[5] != 0);

    bool ok;
    TIME_RUN(siqs(_f, _m, n, &params, SINGLE_RUN), secs, &ok);
    return ok;
}

//------------------------------------------------------------------------------
//                              Registry
//------------------------------------------------------------------------------

// Trailing field is sticky_col: the use_large_primes column index for the sieve
// methods (cfrac=3, siqs=5), -1 for the rest.
const tune_algo_desc_t TUNE_ALGOS[] = {
    { "p_minus_1", TIFA_TUNING_NCOLS_P_MINUS_1, TIFA_TUNING_COLS_P_MINUS_1,
      true,  40, 140, 3, TUNE_MULTI_FACTOR_MIN_BITS, build_grid_p_minus_1, run_p_minus_1, -1 },
    { "p_plus_1",  TIFA_TUNING_NCOLS_P_PLUS_1,  TIFA_TUNING_COLS_P_PLUS_1,
      true,  40, 140, 3, TUNE_MULTI_FACTOR_MIN_BITS, build_grid_p_plus_1,  run_p_plus_1, -1 },
    { "rho",       TIFA_TUNING_NCOLS_RHO,       TIFA_TUNING_COLS_RHO,
      true,  40, 100, 3, TUNE_MULTI_FACTOR_MIN_BITS, build_grid_rho,       run_rho, -1 },
    { "ecm",       TIFA_TUNING_NCOLS_ECM,       TIFA_TUNING_COLS_ECM,
      true,  40, 140, 3, TUNE_MULTI_FACTOR_MIN_BITS, build_grid_ecm,       run_ecm, -1 },
    { "cfrac",     TIFA_TUNING_NCOLS_CFRAC,     TIFA_TUNING_COLS_CFRAC,
      false, 56, 140, 2, 0, build_grid_cfrac,     run_cfrac, 3 },
    { "siqs",      TIFA_TUNING_NCOLS_SIQS,      TIFA_TUNING_COLS_SIQS,
      false, 65, 200, 2, 0, build_grid_siqs,      run_siqs, 5 },
};

const unsigned int TUNE_NALGOS = sizeof(TUNE_ALGOS) / sizeof(TUNE_ALGOS[0]);

//------------------------------------------------------------------------------
const tune_algo_desc_t* tune_find_algo(const char* name) {
    for (unsigned int i = 0; i < TUNE_NALGOS; i++) {
        if (strcmp(TUNE_ALGOS[i].name, name) == 0) {
            return &TUNE_ALGOS[i];
        }
    }
    return NULL;
}
