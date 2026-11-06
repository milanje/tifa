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
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>

#include "tifa_config.h"

#define __PREFIX__  "siqs: "
#define __VERBOSE__ TIFA_VERBOSE_SIQS
#define __TIMING__  TIFA_TIMING_SIQS

#if TIFA_VERBOSE_SIQS || TIFA_TIMING_SIQS || TIFA_PRINT_ERROR
    #include <stdio.h>
#endif

#include <string.h>

#include "first_primes.h"
#include "gmp_utils.h"
#include "funcs.h"
#include "macros.h"
#include "x_tree.h"
#include "lp_table.h"
#include "bernsteinisms.h"
#include "factoring_machine.h"
#include "siqs.h"
#include "lindep.h"
#include "block_lanczos.h"
#include "tuning.h"
#include "tifa_factor.h"
#include "siqs_poly.h"
#include "siqs_sieve.h"
#include "messages.h"
#include "print_error.h"
#include "tifa_alloc.h"
#include "stopping.h"

//-----------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------

//
// The table of the partial relations starts with this many slots per row of
// the matrix (it grows as needed).
//
#define LP_TABLE_SLOTS_PER_ROW  4
//
// The sieve chunk size the half widths of best_params were measured with.
//
#define SIQS_TABLE_CHUNK_SIZE 65536
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
// Largest factor base: the primes after the tabulated ones are sieved.
//
#define MAX_NPRIMES_IN_BASE (1U << 16)
//
// Odd numbers sieved at a time for the primes after the table.
//
#define FB_SEGMENT          (1U << 15)
//
// Up to DENSE_DECOMP_MAX_BASE primes in the base (about 240 bits), the
// decompositions of the relations are dense, a byte per prime: 0.94 of the
// time of sparse ones at 80 bits, 0.97 at 100. Beyond, they keep only the
// exponents that are not 0 (see decomp_t): a dense one would take 1.5 GB at
// 38280 primes.
//
#define DENSE_DECOMP_MAX_BASE (1U << 14)
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
// From PRECHECK_MIN_BITS, the candidates whose cofactor over the factor base
// is estimated above 2^PRECHECK_MARGIN_BITS times its largest prime leave
// before the batch smoothness test (see drop_hopeless_candidates), and the
// sieve threshold is lowered by PRECHECK_THRESHOLD_DROP (twice that from 200
// bits) for the partial relations the cheaper candidates bring. Measured best
// from 140 to 220 bits on an Apple M5 and x86-64 (below, the check costs more
// than it saves).
//
#define PRECHECK_MIN_BITS       140
#define PRECHECK_MARGIN_BITS    8
#define PRECHECK_THRESHOLD_DROP 2
//
// The smallest primes of the base whose constants the pre-check keeps: at
// least those the sieve may skip.
//
#define PRECHECK_SMALL_PRIMES   64

#if (NPRIMES_TO_SKIP_LARGE >= PRECHECK_SMALL_PRIMES)
    #error "NPRIMES_TO_SKIP_LARGE must be below PRECHECK_SMALL_PRIMES"
#endif
//
// From FILTERED_GAUSS_MIN_BASE primes in the base (about 188 bits), the
// matrix is filtered before the elimination (see FILTERED_GAUSS_ELIM): 0.96
// of the time at 210 and 220 bits on an Apple M5, but up to 1.05 below 160
// bits, where the elimination costs less than the filtering saves.
//
#define FILTERED_GAUSS_MIN_BASE 3000
//
// From BLOCK_LANCZOS_MIN_BASE primes in the base (about 230 bits), block
// Lanczos on a sparse matrix (see BLOCK_LANCZOS) rather than the dense one.
// The linear algebra took the same time at 220 bits on an Apple M5 Pro,
// 0.61 of it at 240, 0.26 at 264, 0.15 at 280 (1.1 s rather than 7.5), and
// the run's peak memory 0.09 to 0.20 of it from 240 to 280 bits.
//
#define BLOCK_LANCZOS_MIN_BASE 12000
//
// With SKIP_LARGE_MIN_BASE primes or more, the sieve skips NPRIMES_TO_SKIP_LARGE
// of them (see siqs_sieve.h): the threshold is SKIP_LARGE_THRESHOLD_DROP lower
// for the logs they no longer add, 0.95 to 0.99 of the time from 170 to 220
// bits on an Apple M5.
//
#define SKIP_LARGE_THRESHOLD_DROP 4
//
// Largest number of primes added to 'a': each one doubles the polynomials
// per 'a', and 2 already tested worse than none at every size.
//
#define MAX_EXTRA_PRIMES_IN_A 2
//
// Number of ranges for the size of the number to factor. Each range is
// define by the "optimal" (read good enough) values of the parameters to use.
//
#define NRANGES             37
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
    // The sieve byte of each xi of xpool.
    //
    byte_array_t *xbytes;
    //
    // Largest cofactor, in bits, a candidate is estimated to keep, or 0 to
    // keep all of them (see drop_hopeless_candidates).
    //
    double max_cofactor_bits;
    //
    // For the first primes p of the factor base: 1/p, log2(p), and p^2 with
    // (2^64 - 1) / p^2 and 2^64 mod p^2 (see drop_hopeless_candidates).
    //
    double   small_inv[PRECHECK_SMALL_PRIMES];
    double   small_log2[PRECHECK_SMALL_PRIMES];
    uint32_t small_p2[PRECHECK_SMALL_PRIMES];
    uint64_t small_m2[PRECHECK_SMALL_PRIMES];
    uint32_t small_r2[PRECHECK_SMALL_PRIMES];
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
    // The partial relations of the large prime variation, by large prime.
    //
    lp_table_t* lpt;
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
    {.size_n =  70, .size_base =   92, .hwidth =  16800, .thres = 27},
    {.size_n =  75, .size_base =  109, .hwidth =  10200, .thres = 29},
    {.size_n =  80, .size_base =  126, .hwidth =  12000, .thres = 30},
    {.size_n =  85, .size_base =  144, .hwidth =   9750, .thres = 32},
    {.size_n =  90, .size_base =  163, .hwidth =  15000, .thres = 33},
    {.size_n = 100, .size_base =  200, .hwidth =  10500, .thres = 36},
    {.size_n = 105, .size_base =  244, .hwidth =  18000, .thres = 39},
    {.size_n = 110, .size_base =  288, .hwidth =  14700, .thres = 41},
    {.size_n = 120, .size_base =  400, .hwidth =  30000, .thres = 45},
    {.size_n = 130, .size_base =  400, .hwidth =  30000, .thres = 45},
    {.size_n = 140, .size_base =  700, .hwidth =  49152, .thres = 48},
    {.size_n = 150, .size_base =  800, .hwidth =  49152, .thres = 51},
    {.size_n = 160, .size_base = 1300, .hwidth =  32768, .thres = 56},
    {.size_n = 170, .size_base = 1900, .hwidth =  32768, .thres = 60},
    {.size_n = 180, .size_base = 2520, .hwidth =  39322, .thres = 62},
    {.size_n = 190, .size_base = 3200, .hwidth =  65536, .thres = 67},
    {.size_n = 195, .size_base = 3900, .hwidth = 131072, .thres = 69},
    {.size_n = 200, .size_base = 4600, .hwidth = 196608, .thres = 73},
    {.size_n = 210, .size_base =  6500, .hwidth =  196608, .thres =  77},
    {.size_n = 215, .size_base =  7400, .hwidth =  196608, .thres =  79},
    {.size_n = 220, .size_base =  8800, .hwidth =  196608, .thres =  81},
    {.size_n = 225, .size_base = 10400, .hwidth =  450000, .thres =  84},
    //
    // From 232 bits, least-squares quadratics in the size of n, fitted on the
    // rows from 140 to 200 bits and on those tuned at 272 and 288 bits on an
    // Apple M5 Pro: log2 of the base, the threshold, and log2 of the half
    // width (from 195 bits, below which a sieve chunk bounds it), which peaks
    // at 268 bits and is kept from there. The base stops at
    // MAX_NPRIMES_IN_BASE, past 288 bits. Below, the fit's half widths and
    // thresholds took 1.06 to 1.18 of the time of the rows above, from 212 to
    // 222 bits. The half widths are half the fit's from 248 to 280 bits, which
    // took 0.89 to 0.96 of its time on the M5 and 0.74 to 0.91 on a Xeon
    // 8573C, but up to 1.09 on an EPYC 7763 from 264 bits; and 0.7 of the
    // fit's from 288 bits, where half took 1.01 to 1.03 on the M5.
    //
    {.size_n = 232, .size_base = 12950, .hwidth =  579000, .thres =  84},
    {.size_n = 240, .size_base = 16600, .hwidth =  693000, .thres =  87},
    {.size_n = 248, .size_base = 21250, .hwidth =  396500, .thres =  90},
    {.size_n = 256, .size_base = 27200, .hwidth =  434000, .thres =  94},
    {.size_n = 264, .size_base = 34700, .hwidth =  454000, .thres =  97},
    {.size_n = 272, .size_base = 44250, .hwidth =  457000, .thres = 100},
    {.size_n = 280, .size_base = 56300, .hwidth =  457000, .thres = 103},
    {.size_n = 288, .size_base = 65536, .hwidth =  640000, .thres = 106},
    {.size_n = 296, .size_base = 65536, .hwidth =  640000, .thres = 110},
    {.size_n = 304, .size_base = 65536, .hwidth =  640000, .thres = 113},
    {.size_n = 312, .size_base = 65536, .hwidth =  640000, .thres = 116},
    {.size_n = 320, .size_base = 65536, .hwidth =  640000, .thres = 119}
};
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
static mpz_array_t* take_array(const int which, const uint32_t length,
                               const bool exact);
static void give_array(const int which, mpz_array_t* const array);
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
static void drop_hopeless_candidates(siqs_context_t* const context,
                                     const uint32_t first);
static ecode_t collect_relations(siqs_context_t* const context);
inline static void keep_relations_with_smooth_gx(siqs_context_t* const context);
static void add_primes_of_a(siqs_context_t* const context,
                            decomp_t* const decomp,
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
        if ((ecode != GIVING_UP) || tifa_should_stop()) {
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
    //
    // The half widths of best_params up to 190 bits fit in one sieve chunk of
    // 64 KiB, the size they were measured with. On smaller chunks (32 KiB on
    // most x86-64 cores), one chunk measured 4% to 12% faster than two from
    // 140 to 190 bits.
    //
    if (   (params->sieve_half_width <= SIQS_TABLE_CHUNK_SIZE)
        && (params->sieve_half_width > siqs_sieve_chunk_size())) {
        params->sieve_half_width = siqs_sieve_chunk_size();
    }
    params->nprimes_tdiv     = params->nprimes_in_base;
    params->nrelations       = SIQS_DFLT_NRELATIONS;
    params->extra_primes_in_a = best_params[i].extra;
    params->linalg_method    = (params->nprimes_in_base >= BLOCK_LANCZOS_MIN_BASE)
                             ? BLOCK_LANCZOS
                             : (params->nprimes_in_base >= FILTERED_GAUSS_MIN_BASE)
                             ? FILTERED_GAUSS_ELIM : SIQS_DFLT_LINALG_METHOD;
    params->use_large_primes = (mpz_sizeinbase(n, 2)
                                >= SIQS_DFLT_LARGE_PRIMES_MIN_BITS);
    if (   params->use_large_primes
        && (mpz_sizeinbase(n, 2) >= PRECHECK_MIN_BITS)) {
        params->threshold -= (mpz_sizeinbase(n, 2) >= 200)
                           ? 2 * PRECHECK_THRESHOLD_DROP
                           : PRECHECK_THRESHOLD_DROP;
    }
    if (params->nprimes_in_base >= SKIP_LARGE_MIN_BASE) {
        params->threshold -= SKIP_LARGE_THRESHOLD_DROP;
    }

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
//
// The integers of a run's largest arrays keep their limbs from one run to the
// next in the same thread: allocating them took 7 to 9% of a run from 66 to 80
// bits. Arrays above POOL_MAX_LENGTH (runs above about 150 bits, too long for
// it to matter) are freed instead, so that a thread keeps 170 KiB at most. A
// thread's arrays are freed when it exits.
//
#define POOL_MAX_LENGTH 1024

enum {
    POOL_CAND_REDGX,
    POOL_CAND_U,
    POOL_REDGX,
    POOL_U,
    POOL_PARTIAL_GX,
    POOL_NARRAYS
};
typedef struct {
    mpz_array_t* arrays[POOL_NARRAYS];
} siqs_pool_t;

static pthread_key_t  pool_key;
static bool           pool_key_ok = false;
static pthread_once_t pool_once   = PTHREAD_ONCE_INIT;

static void free_pool(void* const p) {
    siqs_pool_t* const pool = (siqs_pool_t*) p;
    for (int i = 0; i < POOL_NARRAYS; i++) {
        if (pool->arrays[i] != NULL) {
            free_mpz_array(pool->arrays[i]);
        }
    }
    free(pool);
}

static void create_pool_key(void) {
    pool_key_ok = (0 == pthread_key_create(&pool_key, free_pool));
}

static siqs_pool_t* thread_pool(void) {
    pthread_once(&pool_once, create_pool_key);
    if (!pool_key_ok) {
        return NULL;
    }
    siqs_pool_t* pool = (siqs_pool_t*) pthread_getspecific(pool_key);
    if (pool == NULL) {
        pool = tifa_calloc(1, sizeof(siqs_pool_t));
        if (0 != pthread_setspecific(pool_key, pool)) {
            free(pool);
            return NULL;
        }
    }
    return pool;
}
//------------------------------------------------------------------------------
static mpz_array_t* take_array(const int which, const uint32_t length,
                               const bool exact) {
    //
    // An empty array for 'length' integers, from the thread's pool if it has
    // one. With 'exact', it has room for exactly 'length': the relations
    // arrays stop filling up when full.
    //
    siqs_pool_t* const pool = thread_pool();
    mpz_array_t* array = (pool != NULL) ? pool->arrays[which] : NULL;
    if (array == NULL) {
        return alloc_mpz_array(length);
    }
    pool->arrays[which] = NULL;
    if ((array->alloced < length) || (exact && (array->alloced != length))) {
        resize_mpz_array(array, length);
    }
    array->length = 0;
    return array;
}
//------------------------------------------------------------------------------
static void give_array(const int which, mpz_array_t* const array) {
    //
    // A run nested in another finds the place taken: its arrays are freed.
    //
    siqs_pool_t* const pool = (array->alloced <= POOL_MAX_LENGTH)
                            ? thread_pool() : NULL;
    if ((pool == NULL) || (pool->arrays[which] != NULL)) {
        free_mpz_array(array);
        return;
    }
    pool->arrays[which] = array;
}
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

    const uint32_t* const fb = context->factor_base->data;
    context->max_cofactor_bits =
        (   params->use_large_primes
         && (mpz_sizeinbase(context->n, 2) >= PRECHECK_MIN_BITS))
        ? log2((double) fb[context->factor_base->length - 1])
          + PRECHECK_MARGIN_BITS
        : 0;
    for (uint32_t j = 0; j < MIN(context->factor_base->length,
                                 PRECHECK_SMALL_PRIMES); j++) {
        const uint32_t p  = fb[j];
        const uint32_t p2 = p * p;
        context->small_inv[j]  = 1.0 / p;
        context->small_log2[j] = log2((double) p);
        context->small_p2[j]   = p2;
        context->small_m2[j]   = UINT64_MAX / p2;
        context->small_r2[j]   = (uint32_t) ((UINT64_MAX % p2 + 1) % p2);
    }

    context->ptree = prod_tree_ui(context->factor_base);

    //
    // The number of columns is indeed factor_base->length + 1 since the first
    // column is used to keep track of the sign of the residue.
    //
    context->matrix_ncols  = params->nprimes_in_base + 1;
    context->matrix_nrows  = context->matrix_ncols;
    context->matrix_nrows += params->nrelations;

    //
    // Block Lanczos builds its own sparse matrix from the decompositions.
    //
    context->matrix = NULL;
    if (params->linalg_method != BLOCK_LANCZOS) {
        context->matrix = alloc_binary_matrix(
                              context->matrix_nrows,
                              context->matrix_ncols
                          );
        context->matrix->nrows = 0;
    }

    context->npairs_wanted = context->matrix_nrows;
    uint32_t e             = most_significant_bit(context->npairs_wanted);
    //
    // We'll collect more than the number of relations needed before attempting
    // a batch smoothness test a la Bernstein.
    //
    context->to_collect = 1 << (e + 1);

    context->xpool  = alloc_int32_array(context->to_collect);
    context->xbytes = alloc_byte_array(context->to_collect);

    context->cand_redgx    = take_array(POOL_CAND_REDGX, context->to_collect,
                                        false);
    context->cand_u        = take_array(POOL_CAND_U, context->to_collect,
                                        false);
    context->cand_a_idx    = alloc_uint32_array(context->to_collect);

    context->smooth_redgx       = take_array(POOL_REDGX,
                                             context->matrix_nrows, true);
    context->a_idx_for_smooth   = alloc_uint32_array(context->matrix_nrows);
    context->a_values           = alloc_mpz_array(ELONGATION);
    context->a_values_npolyinit = 0;
    context->a_primes           = alloc_uint32_array(ELONGATION);
    context->u                  = take_array(POOL_U, context->matrix_nrows,
                                             true);

    context->drop_duplicates = false;
    context->stalled         = false;
    context->relset          = NULL;
    context->relset_mask     = 0;

    context->lpt = NULL;

    if (params->use_large_primes) {
        context->lpt = alloc_lp_table(
                           LP_TABLE_SLOTS_PER_ROW * context->smooth_redgx->alloced
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
        give_array(POOL_CAND_U, context->cand_u);
        give_array(POOL_U, context->u);
        give_array(POOL_CAND_REDGX, context->cand_redgx);
        give_array(POOL_REDGX, context->smooth_redgx);
        free_uint32_array(context->a_idx_for_smooth);
        free_mpz_array(context->a_values);
        free_uint32_array(context->a_primes);

        free_mpz_tree(context->ptree);

        free_uint32_array(context->factor_base);
        free_uint32_array(context->sqrtm_pi);

        free_int32_array(context->xpool);
        free_byte_array(context->xbytes);

        free_byte_array(context->log_factor_base);

        if (context->matrix != NULL) {
            free_binary_matrix(context->matrix);
        }

        mpz_clear(context->n);
        mpz_clear(context->kn);

        if (context->lpt != NULL) {
            free_lp_table(context->lpt);
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
    if (context->matrix != NULL) {
        reset_binary_matrix(context->matrix);
        context->matrix->nrows = 0;
    }

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
    // A batch of candidates can be as large as the relations wanted. The
    // candidates' integers, from the thread's pool, may have room already.
    //
    const uint32_t nrows = context->matrix_nrows;
    if (context->xpool->alloced < nrows) {
        resize_int32_array(context->xpool, nrows);
        resize_byte_array(context->xbytes, nrows);
        resize_uint32_array(context->cand_a_idx, nrows);
    }
    if (context->cand_redgx->alloced < nrows) {
        resize_mpz_array(context->cand_redgx, nrows);
    }
    if (context->cand_u->alloced < nrows) {
        resize_mpz_array(context->cand_u, nrows);
    }

    if (context->matrix != NULL) {
        free_binary_matrix(context->matrix);
        context->matrix = alloc_binary_matrix(context->matrix_nrows,
                                              context->matrix_ncols);
        context->matrix->nrows = 0;
    }

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
    mpz_array_t* partial_gx_array = take_array(POOL_PARTIAL_GX,
                                               context->matrix_nrows, false);
    //
    // Restrict the trial division to only params->nprimes_tdiv primes
    //
    uint32_array_t primes_array;
    primes_array.alloced = 0;
    primes_array.length  = params->nprimes_tdiv;
    primes_array.data    = (uint32_t*)(context->factor_base->data);

    //
    // Without a dense matrix (block Lanczos), the decompositions are sparse
    // whatever the size of the base: fill_trial_div_decomp wants them so.
    //
    decomp_t* const decomp = alloc_decomp(
                                 context->matrix_nrows,
                                 context->factor_base->length,
                                    (context->matrix == NULL)
                                 || (context->factor_base->length
                                     > DENSE_DECOMP_MAX_BASE)
                             );

    //
    // Partly fill the matrix via trial divisions of the g_{a,b}/a, then add
    // the primes of their 'a', known: smaller numbers to divide (one limb
    // instead of two at 70 bits), and no prime of 'a' to find again.
    //
    const uint32_t first_row = (context->matrix != NULL)
                             ? context->matrix->nrows : 0;
    fill_trial_div_decomp(
        context->matrix,
        decomp,
        partial_gx_array,
        context->smooth_redgx,
        &primes_array
    );
    add_primes_of_a(context, decomp, first_row);

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
            decomp,
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

    sort_decomp(decomp);
    //
    // The matrix is now completely filled. Solve the linear system to find
    // relations of the form: A^2 = Q^2 mod n where  A^2 = prod(ui)^2 and
    // Q^2  = prod(g_{a,b}(xi)).
    //
    uint32_array_list_t *relations;
    if (context->matrix != NULL) {
        relations = find_dependencies(context->matrix, params->linalg_method);
    } else {
        sparse_gf2_t* const sparse = sparse_gf2_from_decomp(
                                         decomp, context->smooth_redgx
                                     );
        relations = alloc_uint32_array_list(64);
        block_lanczos(sparse, relations);
        free_sparse_gf2(sparse);
    }

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
                decomp,
                relations,
                context->factor_base
            );

    STOP_TIMER;
    PRINT_TIMING;

    free_decomp(decomp);
    free_uint32_array_list(relations);
    give_array(POOL_PARTIAL_GX, partial_gx_array);

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
    //
    // Then the primes after the table, from a sieve of odd numbers by the
    // tabulated primes (enough below their largest squared). On the heap: on
    // the stack, its pages would be probed at each call, 1% of the time of
    // tifa_factor at 80 bits.
    //
    unsigned char* const composite = (ibase != fb_alloced)
                                   ? tifa_malloc(FB_SEGMENT) : NULL;
    uint64_t lo = first_primes[NFIRST_PRIMES - 1] + 2;

    while (ibase != fb_alloced) {
        const uint64_t hi = lo + 2 * (uint64_t) FB_SEGMENT;
        memset(composite, 0, FB_SEGMENT);
        for (uint32_t j = 1; j < NFIRST_PRIMES; j++) {
            const uint64_t q = first_primes[j];
            if (q * q >= hi) {
                break;
            }
            uint64_t m = MAX(q * q, ((lo + q - 1) / q) * q);
            if ((m & 1) == 0) {
                m += q;
            }
            for (; m < hi; m += 2 * q) {
                composite[(m - lo) / 2] = 1;
            }
        }
        for (uint32_t i = 0; (i < FB_SEGMENT) && (ibase != fb_alloced); i++) {
            if (composite[i]) {
                continue;
            }
            const uint32_t p = (uint32_t) (lo + 2 * i);
            const uint32_t b = floor_log2(p);
            nmodp = mpz_fdiv_ui(kn, p);
            t = (nmodp == 0) ? 0 : sqrtm(nmodp, p);
            if (t == NO_SQRT_MOD_P) {
                continue;
            }
            *fbptr   = p;
            *sqrtptr = t;
            *logptr  = (nmodp == 0)
                     ? (unsigned char) ((b + 1) / 2)
                     : (unsigned char) (b + ((uint64_t) p * p
                                             >= ((uint64_t) 2 << (2 * b))));
            fbptr++;
            sqrtptr++;
            logptr++;
            ibase++;
        }
        lo = hi;
    }
    free(composite);
    sqrtm_pi->length        = ibase;
    factor_base->length     = ibase;
    log_factor_base->length = ibase;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Miscellaneous utilitarian functions.
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
static inline void mpz_set_i128(mpz_t z, const __int128 v) {
    const unsigned __int128 m = (v < 0) ? -(unsigned __int128) v
                                        : (unsigned __int128) v;
    mp_limb_t* const d = mpz_limbs_write(z, 2);
    d[0] = (mp_limb_t) m;
    d[1] = (mp_limb_t) (m >> 64);
    mpz_limbs_finish(z, (v < 0) ? -2 : 2);
}

static inline __int128 mpz_get_i128(const mpz_t z) {
    const unsigned __int128 m = ((unsigned __int128) mpz_getlimbn(z, 1) << 64)
                              | mpz_getlimbn(z, 0);
    return (mpz_sgn(z) < 0) ? -(__int128) m : (__int128) m;
}
#endif
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
    //
    // The polynomials of q2 (see siqs_poly_t) are taken at 2x.
    //
    const int32_t xmul = context->sieve->poly->q2 ? 2 : 1;

#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
    //
    // With 'a' and 'b' of one limb, 'c' below 2^125 and |x| at most 2^26
    // (n up to about 128 bits), u = a.x + b and (u + b).x + c stay below
    // 2^127: no GMP arithmetic.
    //
    if (   (mpz_size(a) == 1) && (mpz_size(b) <= 1)
        && (mpz_sizeinbase(c, 2) <= 125)) {
        const __int128 ai = mpz_get_i128(a);
        const __int128 bi = mpz_get_i128(b);
        const __int128 ci = mpz_get_i128(c);
        for (uint32_t i = cand_redgx->length; i < x_length; i++) {
            const __int128 xi = (__int128) x_data[i] * xmul;
            const __int128 u  = ai * xi + bi;
            mpz_set_i128(cand_u_data[i], u);
            mpz_set_i128(cand_redgx_data[i], (u + bi) * xi + ci);
            cand_a_idx->data[i] = a_idx;
        }
        cand_redgx->length = x_length;
        cand_u->length     = x_length;
        cand_a_idx->length = x_length;
        return;
    }
#endif
    mpz_t rgx;
    mpz_init(rgx);

    for (uint32_t i = cand_redgx->length; i < x_length; i++) {

        const int64_t xval = (int64_t) x_data[i] * xmul;

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
//
// log2(m.2^e) for m in [1, 2) given t = m - 1, within 0.01: log2(1 + t) is
// about t.(1.3466 - 0.3466.t).
//
static inline double log2_approx(const int32_t e, const double t) {
    return (double) e + t * (1.3466 - 0.3466 * t);
}

static inline double log2_mpz(const mpz_t y) {
    const size_t   n   = mpz_size(y);
    const uint64_t hi  = mpz_getlimbn(y, n - 1);
    const int      lz  = __builtin_clzll(hi);
    const uint64_t lo  = (n > 1) ? mpz_getlimbn(y, n - 2) : 0;
    const uint64_t top = (lz == 0) ? hi : ((hi << lz) | (lo >> (64 - lz)));

    return log2_approx((int32_t) (64 * n) - lz - 1,
                       (double) (top << 1) * 0x1p-64);
}
//
// e.log2(p) for p^e the largest power of p dividing y, given that p does.
//
static inline double prime_power_bits(const mpz_t y, const uint32_t p) {
    uint32_t      e  = 1;
    unsigned long pk = (unsigned long) p * p;
    while (mpz_divisible_ui_p(y, pk)) {
        e++;
        if (pk > ULONG_MAX / p) {
            break;
        }
        pk *= p;
    }
    const int lz = __builtin_clz(p);
    return e * log2_approx(31 - lz, (double) (p << (lz + 1)) * 0x1p-32);
}
//
// As prime_power_bits for the j-th prime p of the base, below
// PRECHECK_SMALL_PRIMES: p^2 not dividing y (nearly always) is found by a
// Barrett reduction when y has two limbs.
//
static inline double small_prime_power_bits(
    const siqs_context_t* const context, const mpz_t y, const uint32_t j) {
    if (   (mpz_size(y) <= 2)
        && (0 != mod_2limbs(mpz_getlimbn(y, 0), mpz_getlimbn(y, 1),
                            context->small_p2[j], context->small_m2[j],
                            context->small_r2[j]))) {
        return context->small_log2[j];
    }
    return prime_power_bits(y, context->factor_base->data[j]);
}
//
// Whether the i-th prime p of the factor base divides y.
//
static inline bool divides(const siqs_poly_t* const poly, const uint32_t i,
                           const uint32_t p, const mpz_t y) {
    if (mpz_size(y) > 2) {
        return mpz_divisible_ui_p(y, p);
    }
    return 0 == mod_2limbs(mpz_getlimbn(y, 0), mpz_getlimbn(y, 1), p,
                           poly->pinv[i], poly->r64[i]);
}
//-----------------------------------------------------------------------------
static void drop_hopeless_candidates(siqs_context_t* const context,
                                     const uint32_t first) {
    //
    // The cofactor of a candidate y = g(x)/a over the factor base is about y
    // over the primes whose logs the sieve added, and over those it does not
    // sieve with: 2, the smallest primes and the primes q of 'a' (which divide
    // y again for one x mod q). Those are divided out exactly. The candidates
    // whose cofactor is then estimated above max_cofactor_bits (most of them
    // from 180 bits) leave before the batch smoothness test.
    //
    const siqs_sieve_t* const sieve = context->sieve;
    const siqs_poly_t*  const poly  = sieve->poly;

    const uint32_t* const fb     = context->factor_base->data;
    const int32_t*  const sol1   = poly->sol1->data;
    const int32_t*  const sol2   = poly->sol2->data;
    const uint32_t* const idx_a  = poly->idx_of_a;
    const uint32_t        npa    = poly->nprimes_in_a;
    const uint32_t        nsmall = first_sieved_prime(sieve);
    const int32_t         thres  = (int32_t) sieve->threshold;
    const double* const   small_inv = context->small_inv;
    //
    // The primes of 'a' among the small ones are counted with the others.
    //
    uint64_t amask = 0;
    for (uint32_t k = 0; k < npa; k++) {
        if (idx_a[k] < nsmall) {
            amask |= UINT64_C(1) << idx_a[k];
        }
    }

    int32_t*       const x     = context->xpool->data;
    unsigned char* const bytes = context->xbytes->data;
    mpz_t*         const y     = context->cand_redgx->data;
    mpz_t*         const u     = context->cand_u->data;
    uint32_t*      const a_idx = context->cand_a_idx->data;

    const uint32_t length = context->xpool->length;
    uint32_t       kept   = first;

    for (uint32_t i = first; i < length; i++) {
        if (mpz_sgn(y[i]) == 0) {
            continue;
        }
        double bits = log2_mpz(y[i])
                    - (thres - (int32_t) (signed char) bytes[i])
                    - (double) mpz_scan1(y[i], 0);
        if (bits <= context->max_cofactor_bits) {
            goto keep;
        }
        //
        // p divides y iff x mod p is one of its roots: the primes that do,
        // without branches, by bit of hits. |x| / p from 1/p in doubles is
        // the quotient or one less.
        //
        const uint32_t ax   = (uint32_t) ((x[i] < 0) ? -(int64_t) x[i] : x[i]);
        const double   dx   = (double) ax;
        const uint32_t neg  = 0U - (uint32_t) (x[i] < 0);
        uint64_t       hits = 0;

        for (uint32_t j = 1; j < nsmall; j++) {
            const uint32_t p = fb[j];
            const uint32_t q = (uint32_t) (dx * small_inv[j]);
            uint32_t       r = ax - q * p;
            r -= p & (0U - (uint32_t) (r >= p));

            const uint32_t rn = (p - r) & (0U - (uint32_t) (r != 0));
            const uint32_t rx = (rn & neg) | (r & ~neg);
            hits |= (uint64_t) (  (rx == (uint32_t) sol1[j])
                                | (rx == (uint32_t) sol2[j])) << j;
        }
        hits &= ~amask;
        while (hits != 0) {
            const uint32_t j = (uint32_t) __builtin_ctzll(hits);
            hits &= hits - 1;
            bits -= small_prime_power_bits(context, y[i], j);
        }
        for (uint32_t k = 0; k < npa; k++) {
            const uint32_t q = fb[idx_a[k]];
            if (divides(poly, idx_a[k], q, y[i])) {
                bits -= prime_power_bits(y[i], q);
            }
        }
        if (bits > context->max_cofactor_bits) {
            continue;
        }
      keep:
        if (kept != i) {
            x[kept]     = x[i];
            a_idx[kept] = a_idx[i];
            mpz_swap(y[kept], y[i]);
            mpz_swap(u[kept], u[i]);
        }
        kept++;
    }
    context->xpool->length      = kept;
    context->cand_redgx->length = kept;
    context->cand_u->length     = kept;
    context->cand_a_idx->length = kept;
}
//-----------------------------------------------------------------------------
static void add_primes_of_a(siqs_context_t* const context,
                            decomp_t* const decomp,
                            const uint32_t first_row) {
    //
    // Each prime of a relation's 'a' divides it once more: one more in its
    // decomposition, and its bit of the matrix flipped.
    //
    const uint32_t        npa      = context->sieve->poly->nprimes_in_a;
    const uint32_t* const a_primes = context->a_primes->data;
    const uint32_t* const a_idx    = context->a_idx_for_smooth->data;
    if (context->matrix == NULL) {
        const uint32_t nrels = context->smooth_redgx->length;
        for (uint32_t i = 0; i < nrels; i++) {
            const uint32_t* const q = &a_primes[a_idx[i] * npa];
            for (uint32_t k = 0; k < npa; k++) {
                append_to_decomp(decomp, i, q[k], 1);
            }
        }
        return;
    }
    const uint32_t        nrels    = context->matrix->nrows - first_row;
    unsigned char** const dense    = (decomp->dense != NULL)
                                   ? decomp->dense->data : NULL;

    for (uint32_t i = 0; i < nrels; i++) {
        const uint32_t* const q = &a_primes[a_idx[i] * npa];
        for (uint32_t k = 0; k < npa; k++) {
            if (dense != NULL) {
                dense[i][q[k]]++;
            } else {
                append_to_decomp(decomp, i, q[k], 1);
            }
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
    uint32_t fills        = 0;

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
            scan_sieve(context->sieve, context->xpool, context->xbytes->data,
                       context->to_collect);

            compute_reduced_polynomial_values(context);
            if (context->max_cofactor_bits != 0) {
                drop_hopeless_candidates(context, ncands);
            }
            if (((++fills & 15) == 0) && tifa_should_stop()) {
                return GIVING_UP;
            }
            //
            // After the drop, which can leave none of the new candidates.
            //
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
            context->lpt,
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
