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
 * \file    approx.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#include "approx.h"
#include "gmp_utils.h"
#include "macros.h"
#include "funcs.h"
#include "tifa_config.h"

//------------------------------------------------------------------------------
//
// Several log tolerance values depending on how much "room" we have in
// finding an approximation. For example, if the number of possible
// combinations is small, we have to raise the tolerance.
//
#define TINY_DLOG10_TOLERANCE       0.015
#define SMALL_DLOG10_TOLERANCE      0.030
#define MEDIUM_DLOG10_TOLERANCE     0.080
#define LARGE_DLOG10_TOLERANCE      0.200
#define HUGE_DLOG10_TOLERANCE       0.500
#define DESPERATE_DLOG10_TOLERANCE  1.000
//
// Factors used in the approximation are taken between facpool[mid-delta]
// and facpool[mid+delta] for mid such as facpool[mid] is roughly equal to
// the nfactors-th root of the target number to approximate.
//
#define TINY_DELTA     12
#define SMALL_DELTA    22
#define MEDIUM_DELTA   32
#define LARGE_DELTA    48
//
// Increase the tolerance by these increments according to 'delta'.
//
#define TINY_DELTA_TOLERANCE_INCR     0.040
#define SMALL_DELTA_TOLERANCE_INCR    0.020
#define MEDIUM_DELTA_TOLERANCE_INCR   0.010
#define FAILURE_DELTA_TOLERANCE_INCR  0.100
//
// Thresholds used to set 'delta' depending on the value of 'mid', e.g. if
// 'mid' is less that MEDIUM_DELTA_THRESHOLD, use SMALL_DELTA.
//
#define SMALL_DELTA_THRESHOLD      30
#define MEDIUM_DELTA_THRESHOLD     50
#define LARGE_DELTA_THRESHOLD      75
//
// Raise the log tolerance by FAILURE_DELTA_TOLERANCE_INCR after
// NTRIES_BEFORE_RAISING_THRESHOLD consecutive unsuccesful attempts at
// finding a satisfying approximation.
//
#define NTRIES_BEFORE_RAISING_THRESHOLD  10
//
// Increment the subset rank by RANK_INCR.
//
#define RANK_INCR 17

#define MY_DEBUG 0
#if MY_DEBUG
    static unsigned int itries = 0;
    static unsigned int ifails = 0;

    #define FOPEN(...)   fopen(__VA_ARGS__);
    #define FCLOSE(...)  fclose(__VA_ARGS__);
    #define MSG(...)     printf(__VA_ARGS__);fflush(stdout);
    #define FMSG(...)    fprintf(__VA_ARGS__);
    #define PAUSE()      do {                   \
        printf("Press any key to continue..."); \
        char g;                                 \
        g = getchar()                           \
      } while (0)
#else
    #define FOPEN(...)   NULL
    #define FCLOSE(...)  /* voluntarily left empty */
    #define MSG(...)     /* voluntarily left empty */
    #define FMSG(...)    /* voluntarily left empty */
    #define PAUSE(...)   /* voluntarily left empty */
#endif

//------------------------------------------------------------------------------
uint32_tuple_t* find_best_tuple(float, uint32_tuple_t* const, uint32_t);
int uint32_triplet_cmp_func(const void* const a, const void* const b);
//-----------------------------------------------------------------------------
approximer_t* alloc_approximer(mpz_t                 target,
                               uint32_array_t* const facpool,
                               uint32_t              nfactors)
{
    long unsigned int seed = time(NULL);
    srand(seed);

    approximer_t* const aximer = malloc(sizeof(approximer_t));

    mpz_init_set(aximer->target, target);

    aximer->targetlog = mpz_log10(target);
    aximer->facpool   = facpool;
    aximer->nfactors  = nfactors;

    //
    // Find 'mid' : the center of the interval from which the factors will be
    // draw. We want facpool[mid] close to the nfactors-th root of the target
    //
    mpz_t tmp;
    mpz_init_set(tmp, target);
    mpz_root(tmp, tmp, nfactors);

    uint32_t root = mpz_get_ui(tmp);

    mpz_clear(tmp);

    uint32_t mid = 0;
    uint32_t min = 0;
    uint32_t max = aximer->facpool->length;
    uint32_t i;

    uint32_t* fdata = aximer->facpool->data;

    while ((max - min) > 1) {
        i = (max + min) / 2;
        if (root <= fdata[i]) {
            max = i;
        } else {
            min = i;
        }
    }
    uint32_t dmin = ABS((int32_t)root - (int32_t)(fdata[min]));
    uint32_t dmax = ABS((int32_t)root - (int32_t)(fdata[max]));

    if (dmin < dmax) {
        mid = min;
    } else {
        mid = max;
    }
    //
    // Set the delta, the log tolerance and the distribution between odd and
    // even indexes depending on 'mid' and 'nfactors'
    //
    switch (nfactors) {
        //
        // The less the number of factors we can pick, the worse the obtained
        // approximation is.
        //
        case 5:
            aximer->keven = 2;
            aximer->dlog_tolerance = TINY_DLOG10_TOLERANCE;
            break;
        case 4:
            aximer->keven = 2;
            aximer->dlog_tolerance = SMALL_DLOG10_TOLERANCE;
            break;
        case 3:
            aximer->keven = 1;
            aximer->dlog_tolerance = MEDIUM_DLOG10_TOLERANCE;
            break;
        case 2:
            aximer->keven = 1;
            aximer->dlog_tolerance = LARGE_DLOG10_TOLERANCE;
            break;
        case 1:
            aximer->keven = 0;
            aximer->dlog_tolerance = HUGE_DLOG10_TOLERANCE;
            break;
        default:
            aximer->keven = MAX_NPRIMES_IN_TUPLE;
            aximer->dlog_tolerance = TINY_DLOG10_TOLERANCE;
    }
    uint32_t delta = LARGE_DELTA;

    if (mid < SMALL_DELTA_THRESHOLD) {
        delta = TINY_DELTA;
        aximer->dlog_tolerance += TINY_DELTA_TOLERANCE_INCR;

    } else if (mid < MEDIUM_DELTA_THRESHOLD) {
        delta = SMALL_DELTA;
        aximer->dlog_tolerance += SMALL_DELTA_TOLERANCE_INCR;

    } else if (mid < LARGE_DELTA_THRESHOLD) {
        delta = MEDIUM_DELTA;
        aximer->dlog_tolerance += MEDIUM_DELTA_TOLERANCE_INCR;
    }
    aximer->imin = mid - delta;
    aximer->imax = mid + delta;

    aximer->kodd  = nfactors - aximer->keven;
    aximer->neven = delta + 1;
    aximer->nodd  = delta;

    //
    // Switch to the "desperate" mode if needed. Things here are certainly
    // less than optimal but this only happens if the factor base is too small
    // (in other words for tiny factorizations, say n < 2^70).
    //
    bool overflow  = (aximer->imax >= (int32_t)facpool->length);
    bool underflow = (aximer->imin < 0);

    if (overflow) {
        aximer->imax = facpool->length - 1;
    }
    if (underflow) {
        aximer->imin = 1;
    }
    if (overflow || underflow) {
        uint32_t diff = (aximer->imax - aximer->imin + 1);

        aximer->nodd  = diff / 2;
        aximer->neven = diff / 2 + (diff & 1);

        aximer->kodd  = nfactors;
        aximer->keven = 0;

        aximer->dlog_tolerance = DESPERATE_DLOG10_TOLERANCE;
    }

    //
    // Precomputes all the combinations of 'keven' factors from the
    // even-indexed pool and sort them according to the product of their
    // factors (here, by looking at its logarithm in base 10).
    //
    if (aximer->keven > 0) {
        uint32_t *subset = malloc(aximer->keven * sizeof(uint32_t));

        for (uint32_t i = 0; i < aximer->keven; i++) {
            //
            // The first combination: (1, 2, ..., 'keven')
            //
            subset[i] = i + 1;
        }
        aximer->ntuples = n_choose_k(aximer->neven, aximer->keven);
        aximer->tuples  = malloc(aximer->ntuples * sizeof(uint32_tuple_t));

        bool end = false;

        for (uint32_t i = 0; i < aximer->ntuples; i++) {

            aximer->tuples[i].tlog = 0.0;

            for (uint32_t j = 0; j < aximer->keven; j++) {
                uint32_t ifac = aximer->imin + 2 * (subset[j] - 1);

                aximer->tuples[i].tlog     += log10(fdata[ifac]);
                aximer->tuples[i].tuple[j]  = ifac;
            }
            //
            // Generate next combination in lexicographic order
            //
            next_subset_lex(aximer->neven, aximer->keven, subset, &end);
        }
        qsort(
            aximer->tuples,
            aximer->ntuples,
            sizeof(uint32_tuple_t),
            uint32_triplet_cmp_func
        );
        free(subset);
    }
    aximer->nsubsets_odd = n_choose_k(aximer->nodd, aximer->kodd);
    aximer->subset_odd   = malloc(aximer->kodd  * sizeof(uint32_t));

    aximer->rank  = aximer->nsubsets_odd / 4;
    aximer->rank += (rand() % (aximer->nsubsets_odd / 4));

    return aximer;
}
//-----------------------------------------------------------------------------
void free_approximer(approximer_t* aximer)
{
    if (aximer != NULL) {
        mpz_clear(aximer->target);
        free(aximer->subset_odd);

        if (aximer->keven > 0) {
            free(aximer->tuples);
        }
        free(aximer);
    }
}
//-----------------------------------------------------------------------------
void
random_approximation(approximer_t* const aximer,
                     mpz_t               approxed,
                     uint32_t*           indexes)
{

    const uint32_t imin  = aximer->imin;
    const uint32_t kodd  = aximer->kodd;
    const uint32_t keven = aximer->keven;

    const uint32_t* const pool = aximer->facpool->data;
    uint32_t nfails = 0;

  CHOOSE_COMBINATION:
    //
    // Get another rank and generate its associated combination.
    //
    // Note that taking random ranks may not work: because of the birthday
    // paradox, the odds to get twice the same subset is _not_ negligible.
    // For example, with 3 elements out of 48, we can expect a collision
    // with probability 0.5 (clearly not acceptable in our case) after
    // sqrt(2*cnk(48, 3)*ln(1/(1-0.5)) = 155 sets drawn.
    //
    aximer->rank += RANK_INCR;
    while (aximer->rank >= aximer->nsubsets_odd) {
        aximer->rank -= aximer->nsubsets_odd;
    }

    unrank_subset_lex(aximer->nodd, kodd, aximer->rank, aximer->subset_odd);

    mpz_set_ui(approxed, 1);

    for (uint32_t j = 0; j < kodd; j++) {
        uint32_t ip = imin + 2 * aximer->subset_odd[j] - 1;
        mpz_mul_ui(approxed, approxed, pool[ip]);
        indexes[j] = ip;
    }
    float dlog = aximer->targetlog - mpz_log10(approxed);

    if (aximer->keven > 0) {
        //
        // The remaining factors should bring a log10 contribution as close
        // to dlog as possible
        //
        uint32_tuple_t*
            best = find_best_tuple(dlog, aximer->tuples, aximer->ntuples);

        for (uint32_t m = 0; m < keven; m++) {
            uint32_t ip = best->tuple[m];
            mpz_mul_ui(approxed, approxed, pool[ip]);
            indexes[m + kodd] = ip;
        }

        dlog = ABS(aximer->targetlog - mpz_log10(approxed));

        if (nfails == NTRIES_BEFORE_RAISING_THRESHOLD) {
            aximer->dlog_tolerance += FAILURE_DELTA_TOLERANCE_INCR;
            nfails = 0;
        }

        if (dlog > aximer->dlog_tolerance) {
            nfails++;
            goto CHOOSE_COMBINATION;
        }
    }
    qsort(indexes, aximer->nfactors, sizeof(uint32_t), uint32_cmp_func);
}
//-----------------------------------------------------------------------------
int
uint32_triplet_cmp_func(const void* const a, const void* const b)
{
    //
    // Compare a uint32_tuple_t using this function...
    //
    uint32_tuple_t* const ta = (uint32_tuple_t* const) a;
    uint32_tuple_t* const tb = (uint32_tuple_t* const) b;

    if (ta->tlog >= tb->tlog) {
        return 1;
    }
    if (ta->tlog < tb->tlog) {
        return -1;
    }
    if (ta->tuple[0] > tb->tuple[0]) {
        return 1;
    }
    if (ta->tuple[0] < tb->tuple[0]) {
        return -1;
    }
    if (ta->tuple[1] > tb->tuple[1]) {
        return 1;
    }
    if (ta->tuple[1] < tb->tuple[1]) {
        return -1;
    }
    if (ta->tuple[2] > tb->tuple[2]) {
        return 1;
    }
    if (ta->tuple[2] < tb->tuple[2]) {
        return -1;
    }
    return 0;
}
//-----------------------------------------------------------------------------
uint32_tuple_t*
find_best_tuple(float log, uint32_tuple_t* tuples, uint32_t ntuples)
{
    //
    // Returns a pointer to the "best" tuple, i.e. the tuple with its 'tlog'
    // field closer to 'log'.
    //
    //     log    : the target logarithm in base 10
    //     tuples : the list of precomputed tuples
    //     ntuples: the number of tuples in the list
    //
    uint32_t min = 0;
    uint32_t max = ntuples - 1;
    uint32_t best_idx = 0;
    uint32_t index;

    //
    // Simple binary search...
    //
    while ((max - min) > 1) {
        index = (max + min) / 2;
        if (log <= tuples[index].tlog) {
            max = index;
        } else {
            min = index;
        }
    }
    float dmin = ABS(log - tuples[min].tlog);
    float dmax = ABS(log - tuples[max].tlog);

    if (dmin < dmax) {
        best_idx = min;
    } else {
        best_idx = max;
    }
    return &(tuples[best_idx]);
}
//-----------------------------------------------------------------------------
