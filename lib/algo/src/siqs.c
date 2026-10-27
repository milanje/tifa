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
 * \file    siqs.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>

#include "tifa_config.h"

#define __PREFIX__  "siqs: "
#define __VERBOSE__ TIFA_VERBOSE_SIQS
#define __TIMING__  TIFA_TIMING_SIQS

#if TIFA_VERBOSE_SIQS || TIFA_TIMING_SIQS || TIFA_PRINT_ERROR
    #include <stdio.h>
#endif

#if TIFA_USE_CALLOC_MEMSET
    #include <string.h>
#endif

#include "first_primes.h"
#include "gmp_utils.h"
#include "funcs.h"
#include "macros.h"
#include "x_tree.h"
#include "hashtable.h"
#include "bernsteinisms.h"
#include "factoring_machine.h"
#include "siqs.h"
#include "tuning.h"
#include "tifa_factor.h"
#include "siqs_poly.h"
#include "siqs_sieve.h"
#include "messages.h"
#include "print_error.h"
#include "tifa_alloc.h"

//-----------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------

//
// In order to not degrade too much the performance of the hashtable used
// for the large prime variation (due to too many collisions), we set its size
// to the number of rows of the matrix multiplied by HTABLE_SIZE_MULTIPLIER.
//
#define HTABLE_SIZE_MULTIPLIER  16
//
// Minimum number of candidates for smoothness test
//
#define MIN_BATCH_SIZE      256
//
// A run has stalled after this many batches of candidates in a row without
// a new relation, or this many sieve fills in a row without a candidate
// (healthy runs with the default parameters: none, and at most 6, from 40 to
// 200 bits). It typically happens with a factor base too small to provide
// enough distinct 'a' coefficients: the polynomials repeat and only bring
// duplicate relations.
//
#define MAX_IDLE_BATCHES    16
#define MAX_IDLE_FILLS      4096
//
// Attempts with ever larger factor base and sieve before giving up.
//
#define MAX_ATTEMPTS        3
//
// Largest factor base: about half the tabulated primes have kn as a square
// modulo them, so a quarter is always available.
//
#define MAX_NPRIMES_IN_BASE (NFIRST_PRIMES / 4)
//
// Largest sieve half width, far above any useful one, so that the sieve
// positions and the retries stay within 32-bit integers.
//
#define MAX_HALF_WIDTH      (1U << 26)
//
// Largest number of extra relations, far above any useful one, which keeps
// the matrix size within bounds.
//
#define MAX_NRELATIONS      1024
//
// When no dependency splits n, nrelations more relations are collected, up
// to NRELATIONS_ROUNDS times nrelations in all.
//
#define NRELATIONS_ROUNDS   3
//
// Largest sieve threshold: the sieve bytes start at it, and a position is a
// candidate once the logs subtracted from it set the top bit.
//
#define MAX_THRESHOLD 127
//
// Largest number of primes added to 'a': each one doubles the polynomials
// per 'a', and 2 already tested worse than none at every size.
//
#define MAX_EXTRA_PRIMES_IN_A 2
//
// Number of ranges for the size of the number to factor. Each range is
// define by the "optimal" (read good enough) values of the parameters to use.
//
#define NRANGES             35
//
// Set to non zero to perform additional timing measurements. Has no effect
// if __TIMING__ if set to 0.
//
#define __EXTRA_TIMING__    TIFA_EXTRA_TIMING
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                         NON PUBLIC MACRO(S)
//-----------------------------------------------------------------------------

#if __TIMING__ && __EXTRA_TIMING__
   #define GET_TIMER(NAME)         &(context->NAME ##_timer)
   #define DECL_EXTRA_TIMER(NAME)  stopwatch_t NAME ##_timer;
   #define INIT_EXTRA_TIMER(NAME)  init_stopwatch(GET_TIMER(NAME));
   #define START_EXTRA_TIMER(NAME) start_stopwatch(GET_TIMER(NAME));
   #define STOP_EXTRA_TIMER(NAME)  stop_stopwatch(GET_TIMER(NAME));
   #define RESET_EXTRA_TIMER(NAME) reset_stopwatch(GET_TIMER(NAME));
   #define GET_EXTRA_TIMING(NAME)  get_stopwatch_elapsed(GET_TIMER(NAME))

   #define PRINT_FILL_TIMING_EXTRA_MSG                     \
       print_fill_timing(context->sieve);
   #define PRINT_SCAN_TIMING_EXTRA_MSG                     \
       print_scan_timing(context->sieve);
   #define PRINT_INIT_POLY_TIMING_EXTRA_MSG                \
       print_init_poly_timing(context->sieve);
   #define PRINT_BATCH_TIMING_EXTRA_MSG                                 \
       printf(__PREFIX__ "        Batch:     done in %7.4f seconds\n",  \
              GET_EXTRA_TIMING(batch));
#else
   #define GET_TIMER(NAME)             /* intentionally left empty */
   #define DECL_EXTRA_TIMER(NAME)      /* intentionally left empty */
   #define INIT_EXTRA_TIMER(NAME)      /* intentionally left empty */
   #define START_EXTRA_TIMER(NAME)     /* intentionally left empty */
   #define STOP_EXTRA_TIMER(NAME)      /* intentionally left empty */
   #define RESET_EXTRA_TIMER(NAME)     /* intentionally left empty */
   #define GET_EXTRA_TIMING(NAME)      /* intentionally left empty */
   #define PRINT_FILL_TIMING_EXTRA_MSG       /* intentionally left empty */
   #define PRINT_SCAN_TIMING_EXTRA_MSG       /* intentionally left empty */
   #define PRINT_INIT_POLY_TIMING_EXTRA_MSG  /* intentionally left empty */
   #define PRINT_BATCH_TIMING_EXTRA_MSG      /* intentionally left empty */
#endif

//-----------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//
// An ad-hoc structure holding the data variables used by the
// SIQS implementation.
//
struct struct_siqs_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // A pointer to the SIQS parameters used.
    //
    siqs_params_t* params;
    //
    // Multiplier to use.
    //
    uint32_t multiplier;
    //
    // Number to factor x multiplier.
    //
    mpz_t kn;
    //
    // Floor of the logarithm in base 2 of sqrt(kn).
    //
    uint32_t log_sqrtkn;
    //
    // Factor base used. Only primes pi such that t^2 = N mod (pi) has a
    // solution for t are kept.
    //
    uint32_array_t *factor_base;
    //
    // Logarithms of the elements in the factor base.
    //
    byte_array_t *log_factor_base;
    //
    // For each prime pi in the base, n has a square root mod pi. Keep
    // these modular square roots in the array sqrtm_pi such that
    // sqrtm_pi->data[i]^2 = n mod factor_base->data[i].
    //
    uint32_array_t *sqrtm_pi;
    //
    // Product tree of primes in factor base.
    //
    mpz_tree_t *ptree;
    //
    // Number of columns of the constructed matrix.
    //
    uint32_t matrix_ncols;
    //
    // Maximum number of rows of the constructed matrix.
    //
    uint32_t matrix_nrows;
    //
    // Matrix generated after collecting relations.
    //
    binary_matrix_t *matrix;
    //
    // Number of relations (i.e. of pairs (a.xi+b, g_{a,b}(xi)) with
    // g_{a,b}(xi) smooth) we want to collect.
    //
    uint32_t npairs_wanted;
    //
    // Number of pairs (a.xi+b, g_{a,b}(xi)) to collect before proceeding to
    // a smoothness test.
    //
    uint32_t to_collect;
    //
    // Pool keeping the xi values that has survived the sieving stage, i.e.
    // the xi such that g_{a,b}(xi) _may_ be smooth.
    //
    int32_array_t *xpool;
    //
    // The SIQS sieve
    //
    siqs_sieve_t* sieve;
    //
    // Relations to find: (a.xi+b)^2 = g_{a,b}(xi) (mod kn) with
    // g_{a,b}(xi) smooth. These are the g_{a, b}(xi)/a to check for
    // smoothness (cand_redgx) together with their corresponding ui = a.xi+b
    // and the index of their 'a' parameter in a_values...
    //
    mpz_array_t    *cand_redgx;
    mpz_array_t    *cand_u;
    uint32_array_t *cand_a_idx;
    //
    // All accepted relations, i.e. the g_{a,b}(xi)/a and their
    // corresponding ui = a.xi+b, such that ui^2 = g_{a, b}(xi) (mod n) with
    // g_{a,b}(xi) smooth.
    //
    mpz_array_t *smooth_redgx;
    mpz_array_t *u;
    //
    // Index in a_values of the 'a' parameter associated to each smooth
    // g_{a,b}(xi)/a in the smooth_redgx array.
    //
    uint32_array_t *a_idx_for_smooth;
    //
    // The 'a' parameters used so far, and the polynomial initialization that
    // set the last one (see siqs_poly_t's nfullpolyinit).
    //
    mpz_array_t *a_values;
    uint32_t     a_values_npolyinit;
    //
    // The indexes in the factor base of the primes of each 'a' of a_values,
    // nprimes_in_a (of the polynomial) per 'a'.
    //
    uint32_array_t *a_primes;
    //
    // Hashtable used in large primes variation. It will hold mpz_pair_t's
    // such as ((a.xi + b), g_{a,b}(xi)) with an index in first_primes_array
    // as key.
    //
    hashtable_t *htable;
    //
    // Set once a run has produced only trivial dependencies: from then on,
    // duplicate relations are dropped as they are collected.
    //
    bool drop_duplicates;
    //
    // Set when relation collection stalls (see MAX_IDLE_BATCHES).
    //
    bool stalled;
    //
    // Open-addressing set of the relations kept so far (index in u plus one,
    // zero for an empty slot), used to drop duplicate relations.
    //
    uint32_t* relset;
    uint32_t  relset_mask;
    //
    // Use extra timers to get finer information.
    //
    DECL_EXTRA_TIMER(fill);
    DECL_EXTRA_TIMER(scan);
    DECL_EXTRA_TIMER(init_poly);
    DECL_EXTRA_TIMER(batch);

    uint32_t na_used;
};
//-----------------------------------------------------------------------------
typedef struct struct_siqs_context_t siqs_context_t;
//-----------------------------------------------------------------------------
struct struct_u32triplet_t {
    //
    // If size of number to factor greater than or equal to size_n...
    //
    uint32_t size_n;
    //
    // use size_base as size of factor base.
    //
    uint32_t size_base;
    //
    // use hwidth as sieve half width.
    //
    uint32_t hwidth;
    //
    // use thres as sieve threshold.
    //
    uint32_t thres;
    //
    // use extra as extra_primes_in_a (not interpolated: the row below's).
    //
    uint32_t extra;
};
//------------------------------------------------------------------------------
typedef struct struct_u32triplet_t u32triplet_t;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                         NON PUBLIC DATA
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static const u32triplet_t best_params[NRANGES] = {
    //
    // These "optimal" values were computed by experimental determination.
    // Actually, they should be regarded more as "good enough" values rather
    // than "optimal" values.
    //
    {.size_n =   0, .size_base =   10, .hwidth =   1000, .thres =  1},
    {.size_n =  50, .size_base =   50, .hwidth =  10000, .thres = 23},
    {.size_n =  60, .size_base =   60, .hwidth =  10000, .thres = 26},
    {.size_n =  70, .size_base =   92, .hwidth =  16800, .thres = 26},
    {.size_n =  75, .size_base =  109, .hwidth =  10200, .thres = 28},
    {.size_n =  80, .size_base =  126, .hwidth =  12000, .thres = 30},
    {.size_n =  85, .size_base =  144, .hwidth =   9750, .thres = 31},
    {.size_n =  90, .size_base =  163, .hwidth =  15000, .thres = 33},
    {.size_n = 100, .size_base =  200, .hwidth =  10500, .thres = 35},
    {.size_n = 105, .size_base =  244, .hwidth =  18000, .thres = 38},
    {.size_n = 110, .size_base =  288, .hwidth =  14700, .thres = 40},
    {.size_n = 120, .size_base =  400, .hwidth =  30000, .thres = 44},
    {.size_n = 130, .size_base =  400, .hwidth =  30000, .thres = 44},
    {.size_n = 140, .size_base =  700, .hwidth =  49152, .thres = 48},
    {.size_n = 150, .size_base =  800, .hwidth =  49152, .thres = 51},
    {.size_n = 160, .size_base = 1300, .hwidth =  65536, .thres = 56},
    {.size_n = 170, .size_base = 1900, .hwidth =  65536, .thres = 60},
    {.size_n = 180, .size_base = 2100, .hwidth =  65536, .thres = 62},
    {.size_n = 190, .size_base = 3200, .hwidth =  65536, .thres = 67},
    {.size_n = 200, .size_base = 4600, .hwidth = 196608, .thres = 71},
    //
    // The following parameters were optained by crude extrapolation... and
    // should thus be looked at with a very critical eye...
    //
    {.size_n = 210, .size_base =  6500, .hwidth =  270000, .thres = 76},
    {.size_n = 215, .size_base =  7400, .hwidth =  320000, .thres = 79},
    {.size_n = 220, .size_base =  8800, .hwidth =  375000, .thres = 81},
    {.size_n = 225, .size_base = 10400, .hwidth =  450000, .thres = 84},
    {.size_n = 230, .size_base = 12300, .hwidth =  550000, .thres = 86},
    {.size_n = 235, .size_base = 14500, .hwidth =  650000, .thres = 89},
    {.size_n = 240, .size_base = 17250, .hwidth =  765000, .thres = 91},
    {.size_n = 245, .size_base = 20300, .hwidth =  925000, .thres = 94},
    {.size_n = 250, .size_base = 24000, .hwidth = 1100000, .thres = 96},
    {.size_n = 255, .size_base = 28400, .hwidth = 1300000, .thres = 99},
    {.size_n = 260, .size_base = 33400, .hwidth = 1575000, .thres = 101},
    {.size_n = 265, .size_base = 39500, .hwidth = 1875000, .thres = 104},
    {.size_n = 270, .size_base = 46750, .hwidth = 2250000, .thres = 107},
    {.size_n = 275, .size_base = 55250, .hwidth = 2750000, .thres = 110},
    {.size_n = 280, .size_base = 65250, .hwidth = 3250000, .thres = 113}
};
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_siqs_context(factoring_machine_t* const machine);
static ecode_t clear_siqs_context(factoring_machine_t* const machine);
static ecode_t update_siqs_context(factoring_machine_t* const machine);
static ecode_t perform_siqs(factoring_machine_t* const machine);
static ecode_t collect_and_solve(factoring_machine_t* const machine);
static bool    want_more_relations(siqs_context_t* const context);
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
static void init_factor_base_data(siqs_context_t* const context,
                                  const uint32_t* const nmodp_known,
                                  const int8_t* const legendre_known,
                                  const uint32_t nknown);
static void compute_reduced_polynomial_values(siqs_context_t* const context);
static ecode_t collect_relations(siqs_context_t* const context);
inline static void keep_relations_with_smooth_gx(siqs_context_t* const context);
static void add_primes_of_a(siqs_context_t* const context,
                            byte_matrix_t* const decomp_matrix,
                            const uint32_t first_row);
static void drop_duplicate_relations(
    siqs_context_t* const context,
    const uint32_t first_new
);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t siqs(mpz_array_t* const factors, uint32_array_t* const multis,
             const mpz_t n, const siqs_params_t* const params,
             const factoring_mode_t mode) {

    PRINT_INTRO_MSG("self-initializing quadratic sieve");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    //
    // A copy: the algorithm adjusts it, and the caller's may be shared.
    //
    siqs_params_t own_params = *params;
    //
    // 'a' is the product of determine_nprimes(...) + extra_primes_in_a
    // odd-indexed primes of the base, for a target below sqrt(2kn), with
    // k <= LARGEST_MULTIPLIER < 2^7.
    //
    own_params.extra_primes_in_a = MIN(own_params.extra_primes_in_a,
                                       MAX_EXTRA_PRIMES_IN_A);
    const uint32_t min_base = 2 * (determine_nprimes(
                                       (mpz_sizeinbase(n, 2) + 10) / 2
                                   ) + own_params.extra_primes_in_a) + 1;
    own_params.nprimes_in_base  = MAX(own_params.nprimes_in_base, min_base);
    own_params.nprimes_in_base  = MIN(own_params.nprimes_in_base,
                                      MAX_NPRIMES_IN_BASE);
    own_params.nprimes_tdiv     = MAX(own_params.nprimes_tdiv, 1);
    own_params.nprimes_tdiv     = MIN(own_params.nprimes_tdiv,
                                      own_params.nprimes_in_base);
    own_params.sieve_half_width = MAX(own_params.sieve_half_width, 1);
    own_params.sieve_half_width = MIN(own_params.sieve_half_width,
                                      MAX_HALF_WIDTH);
    own_params.nrelations       = MIN(own_params.nrelations, MAX_NRELATIONS);
    own_params.threshold        = MIN(own_params.threshold, MAX_THRESHOLD);

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) &own_params;
    machine.init_context_func   = init_siqs_context;
    machine.perform_algo_func   = perform_siqs;
    machine.update_context_func = update_siqs_context;
    machine.clear_context_func  = clear_siqs_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    const uint32_t nfactors = factors->length;
    const uint32_t nmultis  = (multis != NULL) ? multis->length : 0;

    for (uint32_t attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
        ecode = run_machine(&machine);
        if (ecode != GIVING_UP) {
            break;
        }
        //
        // Stalled, or only trivial dependencies even without duplicates:
        // more primes give more relations and more distinct polynomials.
        //
        factors->length = nfactors;
        if (multis != NULL) {
            multis->length = nmultis;
        }
        own_params.nprimes_in_base = MIN(own_params.nprimes_in_base * 3 / 2,
                                         MAX_NPRIMES_IN_BASE);
        own_params.nprimes_tdiv    = MIN(own_params.nprimes_tdiv * 3 / 2,
                                         own_params.nprimes_in_base);
        own_params.sieve_half_width = MIN(own_params.sieve_half_width * 3 / 2,
                                          MAX_HALF_WIDTH);
        if (own_params.threshold > 10) {
            own_params.threshold -= 2;
        }
    }

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------
void set_siqs_params_to_default(const mpz_t n, siqs_params_t* const params) {
    //
    // Size of number to factor = size,
    //
    // If:    best_params[i].size_n <= size
    // And:   size < best_params[i+1].size_n
    // Then:  use best_params[i].size_base as size of factor base
    //        and best_params[i].hwidth as sieve half width
    //        and best_params[i].thres as sieve threshold
    // And:   Linear interpolation with best_params[i + 1].*
    //        (but best_params[i].extra as extra_primes_in_a)
    //
    uint32_t size = mpz_sizeinbase(n, 2);
    uint32_t i    = NRANGES - 1;
    while (size < best_params[i].size_n) {
        i--;
    }
    if (i == NRANGES - 1) {
        size = best_params[i].size_n;
    }
    if (size != best_params[i].size_n) {

        //
        // Linear interpolation
        //
        float d = best_params[i+1].size_n - best_params[i].size_n;
        float r = (size - best_params[i].size_n) / d;

        params->nprimes_in_base = (uint32_t) (best_params[i].size_base
            + r * ((float) best_params[i+1].size_base - best_params[i].size_base));

        params->sieve_half_width = (uint32_t) (best_params[i].hwidth
            + r * ((float) best_params[i+1].hwidth - best_params[i].hwidth));

        params->threshold = (uint32_t) (best_params[i].thres
            + r * ((float) best_params[i+1].thres - best_params[i].thres));

    } else {
        params->sieve_half_width = best_params[i].hwidth;
        params->nprimes_in_base  = best_params[i].size_base;
        params->threshold        = best_params[i].thres;
    }
    params->nprimes_tdiv     = params->nprimes_in_base;
    params->nrelations       = SIQS_DFLT_NRELATIONS;
    params->extra_primes_in_a = best_params[i].extra;
    params->linalg_method    = SIQS_DFLT_LINALG_METHOD;
    params->use_large_primes = (mpz_sizeinbase(n, 2)
                                >= SIQS_DFLT_LARGE_PRIMES_MIN_BITS);

    //
    // Override the tabulated/interpolated defaults with machine-tuned values if
    // a tuning file supplies a row for this bit-size. The linear-algebra method
    // is not tuned and keeps its default. Note that the local `size` was
    // clamped to the last row above, so the true bit-size is recomputed.
    //
    tuning_row_t row;
    if (tifa_tuning_get_row(
        "siqs", TIFA_TUNING_NCOLS_SIQS, (uint32_t) mpz_sizeinbase(n, 2),
        TIFA_TUNING_INTERP_SIQS, &row)) {
        params->nprimes_in_base  = row.cols[0];
        params->sieve_half_width = row.cols[1];
        params->threshold        = row.cols[2];
        params->nprimes_tdiv     = row.cols[3];
        params->nrelations       = row.cols[4];
        params->use_large_primes  = (row.cols[5] != 0);
        params->extra_primes_in_a = row.cols[6];
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static ecode_t init_siqs_context(factoring_machine_t* const machine) {
    //
    // Initialize the SIQS implementation specific variables...
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context        = (void*) tifa_malloc(sizeof(siqs_context_t));
    siqs_context_t* context = (siqs_context_t*) machine->context;
    siqs_params_t* params   = (siqs_params_t*)  machine->params;

    context->params = params;

    mpz_init_set(context->n, machine->n);
    mpz_init(context->kn);

    context->factor_base     = alloc_uint32_array(params->nprimes_in_base);
    context->log_factor_base = alloc_byte_array(context->factor_base->alloced);
    context->sqrtm_pi        = alloc_uint32_array(
                                   context->factor_base->alloced
                               );

    uint32_t nmodp_known[KS_NPRIMES];
    int8_t   legendre_known[KS_NPRIMES];
    uint32_t nknown = 0;

    context->multiplier = ks_multiplier_with_residues(
                              context->n,
                              context->factor_base->alloced,
                              nmodp_known, legendre_known, &nknown
                          );

    mpz_mul_ui(context->kn, context->n, context->multiplier);
    //
    // Initialize all the data common to all the polynomials we are likely to
    // use for the sieving.
    //
    init_factor_base_data(context, nmodp_known, legendre_known, nknown);

    context->ptree = prod_tree_ui(context->factor_base);

    //
    // The number of columns is indeed factor_base->length + 1 since the first
    // column is used to keep track of the sign of the residue.
    //
    context->matrix_ncols  = params->nprimes_in_base + 1;
    context->matrix_nrows  = context->matrix_ncols;
    context->matrix_nrows += params->nrelations;

    context->matrix = alloc_binary_matrix(
                          context->matrix_nrows,
                          context->matrix_ncols
                      );
    context->matrix->nrows = 0;

    context->npairs_wanted = context->matrix_nrows;
    uint32_t e             = most_significant_bit(context->npairs_wanted);
    //
    // We'll collect more than the number of relations needed before attempting
    // a batch smoothness test a la Bernstein.
    //
    context->to_collect = 1 << (e + 1);

    context->xpool = alloc_int32_array(context->to_collect);

    context->cand_redgx    = alloc_mpz_array(context->to_collect);
    context->cand_u        = alloc_mpz_array(context->to_collect);
    context->cand_a_idx    = alloc_uint32_array(context->to_collect);

    context->smooth_redgx       = alloc_mpz_array(context->matrix_nrows);
    context->a_idx_for_smooth   = alloc_uint32_array(context->matrix_nrows);
    context->a_values           = alloc_mpz_array(ELONGATION);
    context->a_values_npolyinit = 0;
    context->a_primes           = alloc_uint32_array(ELONGATION);
    context->u                  = alloc_mpz_array(context->matrix_nrows);

    context->drop_duplicates = false;
    context->stalled         = false;
    context->relset          = NULL;
    context->relset_mask     = 0;

    context->htable = NULL;

    if (params->use_large_primes) {
        uint32_t size = HTABLE_SIZE_MULTIPLIER * context->smooth_redgx->alloced;
        context->htable = alloc_init_hashtable(
                              size, uint32_cmp_func, hash_rj_32
                          );
    }

    context->sieve = alloc_siqs_sieve(
                             context->kn,
                             context->factor_base,
                             context->log_factor_base,
                             context->sqrtm_pi,
                             params->sieve_half_width,
                             params->extra_primes_in_a
                         );

    params->sieve_half_width = context->sieve->nchunks
                                    * context->sieve->chunk_size;

    set_siqs_sieve_threshold(context->sieve, params->threshold);

    context->log_sqrtkn = mpz_sizeinbase(context->kn, 2) / 2;

    INIT_EXTRA_TIMER(fill);
    INIT_EXTRA_TIMER(scan);
    INIT_EXTRA_TIMER(init_poly);
    INIT_EXTRA_TIMER(batch);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_siqs_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    siqs_context_t* context = (siqs_context_t*) machine->context;

    if (context != NULL) {

        free_siqs_sieve(context->sieve);

        free_uint32_array(context->cand_a_idx);
        free_mpz_array(context->cand_u);
        free_mpz_array(context->u);
        free_mpz_array(context->cand_redgx);
        free_mpz_array(context->smooth_redgx);
        free_uint32_array(context->a_idx_for_smooth);
        free_mpz_array(context->a_values);
        free_uint32_array(context->a_primes);

        free_mpz_tree(context->ptree);

        free_uint32_array(context->factor_base);
        free_uint32_array(context->sqrtm_pi);

        free_int32_array(context->xpool);

        free_byte_array(context->log_factor_base);

        free_binary_matrix(context->matrix);

        mpz_clear(context->n);
        mpz_clear(context->kn);

        if (context->htable != NULL) {
            free_mpzpair_htable(context->htable);
        }
        free(context->relset);
        free(context);
    }

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_siqs_context(factoring_machine_t* const machine) {
    //
    // All the dependencies were trivial. The usual culprit is duplicate
    // relations (the same u found by several polynomials): drop them, then
    // let the next run collect replacements, dropping duplicates as it goes.
    // Deduplicating only now keeps the common case free of the extra batches
    // that replacing a dropped relation costs.
    //
    siqs_context_t* const context = (siqs_context_t*) machine->context;

    if (context->drop_duplicates || context->stalled) {
        return GIVING_UP;
    }
    context->drop_duplicates = true;

    uint32_t relset_size = 1U << (most_significant_bit(context->matrix_nrows)
                                  + 2);
    context->relset      = tifa_calloc(relset_size, sizeof(uint32_t));
    context->relset_mask = relset_size - 1;

    const uint32_t nrels = context->u->length;

    drop_duplicate_relations(context, 0);

    if (context->u->length == nrels) {
        return GIVING_UP;
    }
    reset_binary_matrix(context->matrix);
    context->matrix->nrows = 0;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t perform_siqs(factoring_machine_t* const machine) {
    //
    // With k extra relations, no dependency splits n = pq with probability
    // 2^-k: rare enough to start with few and collect more only then.
    //
    siqs_context_t* const context = (siqs_context_t*) machine->context;
    for (;;) {
        const ecode_t ecode = collect_and_solve(machine);
        if ((ecode != NO_FACTOR_FOUND) || !want_more_relations(context)) {
            return ecode;
        }
    }
}
//------------------------------------------------------------------------------
static bool want_more_relations(siqs_context_t* const context) {
    //
    // params->nrelations more, up to NRELATIONS_ROUNDS times as many in all.
    // The relations kept are reused: their arrays grow, and the matrix,
    // refilled at each run, is reallocated.
    //
    const uint32_t step  = context->params->nrelations;
    const uint32_t extra = context->matrix_nrows - context->matrix_ncols;
    if ((step == 0) || (extra + step > NRELATIONS_ROUNDS * step)) {
        return false;
    }
    context->matrix_nrows += step;
    resize_mpz_array(context->smooth_redgx, context->matrix_nrows);
    resize_mpz_array(context->u, context->matrix_nrows);
    resize_uint32_array(context->a_idx_for_smooth, context->matrix_nrows);
    //
    // A batch of candidates can be as large as the relations wanted.
    //
    if (context->xpool->alloced < context->matrix_nrows) {
        resize_int32_array(context->xpool, context->matrix_nrows);
        resize_mpz_array(context->cand_redgx, context->matrix_nrows);
        resize_mpz_array(context->cand_u, context->matrix_nrows);
        resize_uint32_array(context->cand_a_idx, context->matrix_nrows);
    }

    free_binary_matrix(context->matrix);
    context->matrix = alloc_binary_matrix(context->matrix_nrows,
                                          context->matrix_ncols);
    context->matrix->nrows = 0;

    return true;
}
//------------------------------------------------------------------------------
static ecode_t collect_and_solve(factoring_machine_t* const machine) {

    INIT_TIMER;

    siqs_context_t* context = (siqs_context_t*) machine->context;
    siqs_params_t*  params  = context->params;

    PRINT_COLLECT_RELS_MSG;
    START_TIMER;

    //
    // Collect congruences relations
    //
    ecode_t ecode = collect_relations(context);

    PRINT_FILL_TIMING_EXTRA_MSG;
    PRINT_SCAN_TIMING_EXTRA_MSG;
    PRINT_INIT_POLY_TIMING_EXTRA_MSG;
    PRINT_BATCH_TIMING_EXTRA_MSG;

    if (ecode != SUCCESS) {
        return ecode;
    }

    STOP_TIMER;
    PRINT_COLLECT_RELS_DONE_MSG;
    PRINT_TIMING;

    //
    // We have now all the relations we need to fill the matrix and solve the
    // resulting linear system. The rest of the algorithm is now completely
    // standard to all congruence of squares methods.
    //
    PRINT_FACTOR_RES_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Cofactors of generated smooth g_{a,b}(xi) after trial divisions.
    //
    mpz_array_t* partial_gx_array = alloc_mpz_array(context->matrix_nrows);
    //
    // Restrict the trial division to only params->nprimes_tdiv primes
    //
    uint32_array_t primes_array;
    primes_array.alloced = 0;
    primes_array.length  = params->nprimes_tdiv;
    primes_array.data    = (uint32_t*)(context->factor_base->data);

    byte_matrix_t* decomp_matrix = alloc_byte_matrix(
                                      context->matrix_nrows,
                                      context->factor_base->length
                                   );

    //
    // Partly fill the matrix via trial divisions of the g_{a,b}/a, then add
    // the primes of their 'a', known: smaller numbers to divide (one limb
    // instead of two at 70 bits), and no prime of 'a' to find again.
    //
    const uint32_t first_row = context->matrix->nrows;
    fill_trial_div_decomp(
        context->matrix,
        decomp_matrix,
        partial_gx_array,
        context->smooth_redgx,
        &primes_array
    );
    add_primes_of_a(context, decomp_matrix, first_row);

    if (params->nprimes_tdiv < params->nprimes_in_base) {
        //
        // The g_{a,b} are not completely factored on our factor base...
        //
        uint32_array_list_t*
        decomp_list = alloc_uint32_array_list(context->smooth_redgx->length);
        //
        // The primes in our factor base that was _not_ used in the
        // previous trial divisions are now pointed by primes_array
        //
        primes_array.length = params->nprimes_in_base - params->nprimes_tdiv;
        primes_array.data = (uint32_t*)
                            &(context->factor_base->data[params->nprimes_tdiv]);
        //
        // Use Bernstein's batch algo for the complete factorization of the
        // residues...
        //
        bern_71(decomp_list, partial_gx_array, &primes_array);
        //
        // ... and finish to fill the matrix...
        //
        fill_matrix_from_list_decomp(
            context->matrix,
            decomp_matrix,
            partial_gx_array,
            decomp_list,
            context->factor_base
        );
        free_uint32_array_list(decomp_list);
    }
    STOP_TIMER;
    PRINT_TIMING;
    PRINT_LIN_ALG_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // The matrix is now completely filled. Solve the linear system to find
    // relations of the form: A^2 = Q^2 mod n where  A^2 = prod(ui)^2 and
    // Q^2  = prod(g_{a,b}(xi)).
    //
    uint32_array_list_t *relations;
    relations = find_dependencies(context->matrix, params->linalg_method);

    STOP_TIMER;
    PRINT_TIMING;
    PRINT_DED_FACTORS_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Deduce factors of n using the found relations.
    //
    ecode = find_factors_decomp(
                machine->factors,
                context->n,
                context->u,
                decomp_matrix,
                relations,
                context->factor_base
            );

    STOP_TIMER;
    PRINT_TIMING;

    free_byte_matrix(decomp_matrix);
    free_uint32_array_list(relations);
    free_mpz_array(partial_gx_array);

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor(factors, multis, n, mode);

    //
    // The following code should be used to factor a number using _only_ SIQS.
    //
    //siqs_params_t params;
    //set_siqs_params_to_default(n, &params);
    //return siqs(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                Determination of factor base and misc. data
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void init_factor_base_data(siqs_context_t* const context,
                                  const uint32_t* const nmodp_known,
                                  const int8_t* const legendre_known,
                                  const uint32_t nknown) {
    //
    // Initializes the factor base factor_base by keeping only those primes
    // pi such that the number to factor n has a square root mod pi, and keeps
    // one of these square roots in the array sqrtm_pi.
    //
    // Also precomputes the logarithms of the primes in the factor base and
    // stores them in the array log_factor_base.
    //
    uint32_array_t* const factor_base     = context->factor_base;
    byte_array_t*   const log_factor_base = context->log_factor_base;
    uint32_array_t* const sqrtm_pi        = context->sqrtm_pi;

    uint32_t fb_alloced = factor_base->alloced;
    const mpz_srcptr kn = context->kn;

    //
    // Always put 2 in the factor base...
    //
    factor_base->data[0] = 2;
    factor_base->length  = 1;

    log_factor_base->data[0] = 1;
    log_factor_base->length  = 1;

    if (0 == mpz_tstbit(kn, 0)) {
        sqrtm_pi->data[0] = 0;
    } else {
        sqrtm_pi->data[0] = 1;
    }
    sqrtm_pi->length  = 1;

    uint32_t t     = 0;
    uint32_t ibase = 1;
    uint32_t nmodp = 0;

    uint32_t* fbptr       = &(factor_base->data[ibase]);
    uint32_t* sqrtptr     = &(sqrtm_pi->data[ibase]);
    const uint32_t* prime = &(first_primes[0]);
    const uint32_t* last  = &(first_primes[NFIRST_PRIMES - 1]);
    unsigned char* logptr = &(log_factor_base->data[ibase]);

    //
    // Determine the factor base and the solutions to the equations
    // t^2 = N (mod p) where the p's are the primes in the factor base
    //
    //
    // For the first primes, ks_multiplier_with_residues already gave n mod p
    // and (kn/p). With b = floor(log2(p)), round(log2(p)) is b, plus 1 when
    // p^2 >= 2^(2b+1), and round(log2(p)/2) is (b + 1)/2.
    //
    while ((ibase != fb_alloced) && (prime != last)) {
        prime++;
        const uint32_t ip = (uint32_t) (prime - first_primes);
        const uint64_t p  = *prime;
        const uint32_t b  = floor_log2(*prime);

        if (ip < nknown) {
            if (legendre_known[ip] == -1) {
                continue;
            }
            nmodp = (uint32_t) ((context->multiplier * (uint64_t) nmodp_known[ip])
                                % p);
        } else {
            nmodp = mpz_fdiv_ui(kn, *prime);
        }
        if (nmodp == 0) {
            //
            // This prime divides n (extremely unlikely) or the multiplier. There
            // is only one (trivial) root to x^2 = nmodp = 0 (mod prime). To
            // keep the code simple we'll do as if there were two (identical)
            // solutions but we divide the contribution from log(prime) by two.
            //
            *fbptr   = *prime;
            *sqrtptr = 0;
            *logptr  = (unsigned char) ((b + 1) / 2);
            fbptr++;
            sqrtptr++;
            logptr++;
            ibase++;

            continue;
        }
        t = (ip < nknown) ? sqrtm_qr(nmodp, *prime) : sqrtm(nmodp, *prime);

        if (t != NO_SQRT_MOD_P) {
            //
            // i.e. if n has a square root mod first_primes[i]...
            //
            *fbptr   = *prime;
            *sqrtptr = t;
            *logptr  = (unsigned char) (b + (p * p >= ((uint64_t) 2 << (2 * b))));
            fbptr++;
            sqrtptr++;
            logptr++;
            ibase++;
        }
    }
    sqrtm_pi->length        = ibase;
    factor_base->length     = ibase;
    log_factor_base->length = ibase;

    prime++;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Miscellaneous utilitarian functions.
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void compute_reduced_polynomial_values(siqs_context_t* const context) {
    //
    // Computes the "reduced" values g_{a,b}(xi)/a = a.xi^2 + 2.b.xi + c for
    // all xi values in the array x (with i >= cand_redgx->length) and stores
    // them in the array cand_redgx. Also keeps the (a.xi + b) in the array
    // cand_u.
    //
    const mpz_ptr a = context->sieve->poly->a;
    const mpz_ptr b = context->sieve->poly->b;
    const mpz_ptr c = context->sieve->poly->c;

    mpz_array_t*    const cand_u     = context->cand_u;
    mpz_array_t*    const cand_redgx = context->cand_redgx;
    uint32_array_t* const cand_a_idx = context->cand_a_idx;

    mpz_t* const cand_u_data     = cand_u->data;
    mpz_t* const cand_redgx_data = cand_redgx->data;

    const uint32_t npolyinit = context->sieve->poly->nfullpolyinit;
    if (npolyinit != context->a_values_npolyinit) {
        const siqs_poly_t* const poly = context->sieve->poly;
        append_mpz_to_array(context->a_values, a);
        for (uint32_t j = 0; j < poly->nprimes_in_a; j++) {
            append_uint32_to_array(context->a_primes, poly->idx_of_a[j]);
        }
        context->a_values_npolyinit = npolyinit;
    }
    const uint32_t a_idx = context->a_values->length - 1;

    const uint32_t x_length = context->xpool->length;
    int32_t* const x_data = context->xpool->data;

    mpz_t rgx;
    mpz_init(rgx);

    for (uint32_t i = cand_redgx->length; i < x_length; i++) {

        int32_t xval = x_data[i];

        if (xval < 0) {
            mpz_mul_ui(rgx, a, (long unsigned int)(-xval));
            mpz_neg(rgx, rgx);
        } else {
            mpz_mul_ui(rgx, a, xval);
        }
        mpz_add(rgx, rgx, b);
        mpz_set(cand_u_data[i], rgx);

        mpz_add(rgx, rgx, b);

        if (xval < 0) {
            mpz_mul_ui(rgx, rgx, (long unsigned int)(-xval));
            mpz_neg(rgx, rgx);
        } else {
            mpz_mul_ui(rgx, rgx, xval);
        }
        mpz_add(rgx, rgx, c);

        mpz_set(cand_redgx_data[i], rgx);
        cand_a_idx->data[i] = a_idx;
    }
    cand_redgx->length = x_length;
    cand_u->length     = x_length;
    cand_a_idx->length = x_length;

    mpz_clear(rgx);
}
//-----------------------------------------------------------------------------
static void add_primes_of_a(siqs_context_t* const context,
                            byte_matrix_t* const decomp_matrix,
                            const uint32_t first_row) {
    //
    // Each prime of a relation's 'a' divides it once more: one more in its
    // decomposition, and its bit of the matrix flipped.
    //
    const uint32_t        npa      = context->sieve->poly->nprimes_in_a;
    const uint32_t* const a_primes = context->a_primes->data;
    const uint32_t* const a_idx    = context->a_idx_for_smooth->data;
    const uint32_t        nrels    = context->matrix->nrows - first_row;

    for (uint32_t i = 0; i < nrels; i++) {
        const uint32_t* const q = &a_primes[a_idx[i] * npa];
        for (uint32_t k = 0; k < npa; k++) {
            decomp_matrix->data[i][q[k]]++;
            flip_matrix_bit(first_row + i, q[k] + 1, context->matrix);
        }
    }
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *          Collection of relations (a.xi + b)^2 = g_{a,b}(xi) mod (kn)
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t collect_relations(siqs_context_t* const context) {
    //
    // Collect the relations (a.xi + b)^2 = g_{a,b}(xi) mod (kn). This is
    // of course the 'core' of the SIQS algorithm...
    //

    uint32_t nrels_to_go  = context->smooth_redgx->alloced;
    uint32_t max_to_go    = nrels_to_go;
    uint32_t idle_batches = 0;
    uint32_t idle_fills   = 0;

    //
    // Collect relations that will ultimately lead to congruences of squares,
    // i.e. keep the (a.x + b) and their associated g_{a,b}(x) values iif
    // g_{a,b}(x) is smooth on the factor base.
    //
    while (nrels_to_go != 0) {
        //
        // Note that we may collect more than the remaining number of relations
        // needed before attempting a batch smoothness test a la Bernstein.
        //
        while (context->xpool->length != context->to_collect) {
            const uint32_t ncands = context->xpool->length;

            fill_sieve(context->sieve);
            scan_sieve(context->sieve, context->xpool, context->to_collect);

            compute_reduced_polynomial_values(context);

            idle_fills = (context->xpool->length == ncands) ? idle_fills + 1 : 0;
            if (idle_fills == MAX_IDLE_FILLS) {
                context->stalled = true;
                return GIVING_UP;
            }
        }
        const uint32_t nrels = context->smooth_redgx->length;
        //
        // After computing the "reduced" g_{a,b}(xi)/a, keep the ones that are
        // smooth in the smooth_redgx array together with the associated
        // u = (a.x + b) values...
        //
        keep_relations_with_smooth_gx(context);

        idle_batches = (context->smooth_redgx->length == nrels)
                       ? idle_batches + 1 : 0;
        if (idle_batches == MAX_IDLE_BATCHES) {
            context->stalled = true;
            return GIVING_UP;
        }
        context->xpool->length        = 0;
        context->cand_redgx->length   = 0;
        context->cand_u->length       = 0;
        context->cand_a_idx->length   = 0;

        PRINT_NRELS_FOUND(
            context->smooth_redgx->length,
            context->smooth_redgx->alloced
        );

        nrels_to_go  = context->smooth_redgx->alloced;
        nrels_to_go -= context->smooth_redgx->length;

        context->to_collect = MAX(nrels_to_go, MIN_BATCH_SIZE);
        context->to_collect = MIN(context->to_collect, max_to_go);
    }
    return SUCCESS;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                        Syntaxic sugar functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void drop_duplicate_relations(siqs_context_t* const context,
                                     const uint32_t first_new) {
    //
    // The same u = a.x + b is found by every polynomial whose 'a' divides
    // u^2 - kn, which is common for small numbers. Duplicate relations (u or
    // -u) only yield trivial dependencies, so keep the first copy of each.
    //
    mpz_t*    const u     = context->u->data;
    mpz_t*    const redgx = context->smooth_redgx->data;
    uint32_t* const a     = context->a_idx_for_smooth->data;

    uint32_t* const relset = context->relset;
    const uint32_t  mask   = context->relset_mask;

    uint32_t nkept = first_new;

    for (uint32_t i = first_new; i < context->u->length; i++) {
        uint32_t slot = (uint32_t) mpz_getlimbn(u[i], 0) & mask;
        bool     dup  = false;

        while (relset[slot] != 0) {
            if (0 == mpz_cmpabs(u[relset[slot] - 1], u[i])) {
                dup = true;
                break;
            }
            slot = (slot + 1) & mask;
        }
        if (dup) {
            continue;
        }
        if (nkept != i) {
            mpz_swap(u[nkept], u[i]);
            mpz_swap(redgx[nkept], redgx[i]);
            a[nkept] = a[i];
        }
        relset[slot] = nkept + 1;
        nkept++;
    }
    context->u->length                = nkept;
    context->smooth_redgx->length     = nkept;
    context->a_idx_for_smooth->length = nkept;
}
//-----------------------------------------------------------------------------
inline static void
keep_relations_with_smooth_gx(siqs_context_t* const context) {

    START_EXTRA_TIMER(batch);

    const uint32_t first_new = context->u->length;

    if (context->params->use_large_primes) {
        bern_21_rt_pairs_lp_siqs(
            context->n,
            context->htable,
            context->u,
            context->smooth_redgx,
            context->a_idx_for_smooth,
            context->cand_u,
            context->cand_redgx,
            context->cand_a_idx,
            context->a_values,
            context->ptree->data[0]
        );
    } else {
        bern_21_rt_pairs_siqs(
            context->u,
            context->smooth_redgx,
            context->a_idx_for_smooth,
            context->cand_u,
            context->cand_redgx,
            context->cand_a_idx,
            context->ptree->data[0]
        );
    }
    if (context->drop_duplicates) {
        drop_duplicate_relations(context, first_new);
    }

    STOP_EXTRA_TIMER(batch);
}
//-----------------------------------------------------------------------------
