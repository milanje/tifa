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
 * \file    siqs_sieve.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Structure and functions related to the sieve used in
 *        the SIQS algorithm.
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include "tifa_config.h"
#include "macros.h"
#include "funcs.h"
#include "print_error.h"
#include "siqs_sieve.h"
#include "siqs_poly.h"
#include "gmp_utils.h"

#if !defined(__TIMING__)
    #define __TIMING__ 1
#endif

#if defined(__PREFIX__)
    #undef  __PREFIX__
#endif
#define __PREFIX__ "siqs: "

#if defined(__VERBOSE__)
    #undef  __VERBOSE__
#endif
#define __VERBOSE__ TIFA_VERBOSE_SIQS

#include "timer.h"
#include "messages.h"

//-----------------------------------------------------------------------------
#define __EXTRA_TIMING__ TIFA_EXTRA_TIMING
//
// Macros used to simplify handlings of the extra stopwatches used (if any)
//
#if __EXTRA_TIMING__
   #define GET_TIMER(NAME)         &(sieve->NAME ##_timer)
   #define DECL_EXTRA_TIMER(NAME)  stopwatch_t NAME ##_timer;
   #define INIT_EXTRA_TIMER(NAME)  init_stopwatch(GET_TIMER(NAME));
   #define START_EXTRA_TIMER(NAME) start_stopwatch(GET_TIMER(NAME));
   #define STOP_EXTRA_TIMER(NAME)  stop_stopwatch(GET_TIMER(NAME));
   #define RESET_EXTRA_TIMER(NAME) reset_stopwatch(GET_TIMER(NAME));
   #define GET_EXTRA_TIMING(NAME)  get_stopwatch_elapsed(GET_TIMER(NAME))
#else
   #define GET_TIMER(NAME)             /* intentionally left empty */
   #define DECL_EXTRA_TIMER(NAME)      /* intentionally left empty */
   #define INIT_EXTRA_TIMER(NAME)      /* intentionally left empty */
   #define START_EXTRA_TIMER(NAME)     /* intentionally left empty */
   #define STOP_EXTRA_TIMER(NAME)      /* intentionally left empty */
   #define RESET_EXTRA_TIMER(NAME)     /* intentionally left empty */
   #define GET_EXTRA_TIMING(NAME)      /* intentionally left empty */
#endif
//
// Resets the sieve DEST of length LEN to the threshold value VAL
//
#if TIFA_USE_CALLOC_MEMSET
    #include <string.h>
    #define TIFA_MEMSET(DEST, VAL, LEN) do {      \
        memset((DEST), (VAL), (LEN));             \
    } while (0)
#else
    #define TIFA_MEMSET(DEST, VAL, LEN) do {      \
        for (uint32_t i_ = 0; i_ < LEN; i_++) {   \
            DEST[i_] = (VAL);                     \
        }                                         \
    } while (0)
#endif
//
// Divides X by the length of a chunk sieve.
// WARNING: a chunk sieve must have a length which is a power of two!
//
#define DIV_CHUNK_LEN(X) ((X) >> LOG2_SIEVE_CHUNK_MAX_SIZE)
//
// X modulo the length of a chunk sieve.
// WARNING: a chunk sieve must have a length which is a power of two!
//
#define MOD_CHUNK_LEN(X) ((X) & (SIEVE_CHUNK_MAX_SIZE - 1))
//-----------------------------------------------------------------------------
int compute_threshold_correction(siqs_sieve_t* const sieve);
ecode_t fill_sieve_chunk(siqs_sieve_t* const sieve);
ecode_t scan_sieve_positive(
    siqs_sieve_t* const sieve,
    int32_array_t* const survivors,
    const uint32_t nsurvivors
);
ecode_t scan_sieve_negative(
    siqs_sieve_t* const sieve,
    int32_array_t* const survivors,
    const uint32_t nsurvivors
);
ecode_t update_sieve_params(siqs_sieve_t* const sieve);
void    update_sieve_threshold(siqs_sieve_t* const sieve);
ecode_t reset_poly_sols(siqs_sieve_t* const sieve);
void    switch_polynomial(siqs_sieve_t* const sieve);
void    fill_buckets(siqs_sieve_t* const sieve);
static inline ecode_t empty_bucket(siqs_sieve_t* const sieve);
static inline ecode_t empty_bucket_posi(siqs_sieve_t* const sieve);
static inline ecode_t empty_bucket_nega(siqs_sieve_t* const sieve);
//-----------------------------------------------------------------------------
siqs_sieve_t* alloc_siqs_sieve(mpz_t n,
                               uint32_array_t* const factor_base,
                               byte_array_t*   const log_primes,
                               uint32_array_t* const sqrtm_pi,
                               uint32_t half_width)
{
    siqs_sieve_t* const sieve = malloc(sizeof(siqs_sieve_t));

    INIT_EXTRA_TIMER(init_poly);
    INIT_EXTRA_TIMER(fill);
    INIT_EXTRA_TIMER(scan);

    sieve->log_primes = log_primes;

    if (half_width <= SIEVE_CHUNK_MAX_SIZE) {
        sieve->sieve         = alloc_byte_array(half_width);
        sieve->sieve->length = half_width;
        sieve->chunk_size    = half_width;
        sieve->nchunks       = 1;

#if !ROUND_HALF_WIDTH
        sieve->endlast       = half_width - 1;
#endif

    } else {
        sieve->sieve         = alloc_byte_array(SIEVE_CHUNK_MAX_SIZE);
        sieve->sieve->length = SIEVE_CHUNK_MAX_SIZE;

        uint32_t hwmod = half_width % SIEVE_CHUNK_MAX_SIZE;
        sieve->nchunks = half_width / SIEVE_CHUNK_MAX_SIZE;

#if ROUND_HALF_WIDTH
        //
        // Rounds the half width to the closest multiple of
        // SIEVE_CHUNK_MAX_SIZE.
        //
        if (hwmod > (SIEVE_CHUNK_MAX_SIZE / 2)) {
            sieve->nchunks++;
        }
#else
        //
        // Use the exact half width provided as argument
        //
        if (hwmod != 0) {
            sieve->nchunks++;
            sieve->endlast = hwmod - 1;
        } else {
            sieve->endlast = SIEVE_CHUNK_MAX_SIZE - 1;
        }
#endif

        sieve->chunk_size = SIEVE_CHUNK_MAX_SIZE;
    }
    //
    // Setting scan_begin = chunk_size means that a new chunk should be filled
    //
    sieve->scan_begin = sieve->chunk_size;

    sieve->next_chunkno_to_fill = 1;

    mpz_t target_a;
    mpz_init(target_a);

    mpz_mul_ui(target_a, n, 2);
    mpz_sqrt(target_a, target_a);

#if ROUND_HALF_WIDTH
    mpz_fdiv_q_ui(target_a, target_a, sieve->nchunks * sieve->chunk_size);
#else
    mpz_fdiv_q_ui(target_a, target_a, half_width);
#endif

    sieve->poly = alloc_siqs_poly(target_a, n, factor_base, sqrtm_pi);

    mpz_clear(target_a);

    //
    // Bucket stuff if needed
    //
    uint32_t fblen = factor_base->length;
    uint32_t lastp = factor_base->data[fblen - 1];

    sieve->nprimes_no_buckets = fblen;
    sieve->use_buckets = false;

    if (    (lastp > sieve->chunk_size)
         && (((int32_t)fblen - BUCKETS_SMALLEST_PRIME) > 200)
         && (sieve->nchunks != 1)) {

        sieve->use_buckets = true;
        sieve->nprimes_no_buckets = BUCKETS_SMALLEST_PRIME;
    }

    uint32_t sizecp = sieve->nprimes_no_buckets * sizeof(int32_t);

    sieve->sol1 = alloc_int32_array(sieve->nprimes_no_buckets);
    sieve->sol2 = alloc_int32_array(sieve->nprimes_no_buckets);

    //
    // Keep a copy of the next hit positions
    //
    memcpy(sieve->sol1->data, sieve->poly->sol1->data, sizecp);
    memcpy(sieve->sol2->data, sieve->poly->sol2->data, sizecp);

    sieve->sol1->length = sieve->nprimes_no_buckets;
    sieve->sol2->length = sieve->nprimes_no_buckets;

    if (sieve->use_buckets) {
        //
        // _TO_DO_: The number of primes used for bucket sieving could be
        //          fine-tuned
        //
        sieve->buckets_first_prime = BUCKETS_SMALLEST_PRIME;
        sieve->buckets_positive    = alloc_buckets(2048, sieve->nchunks);
        sieve->buckets_negative    = alloc_buckets(2048, sieve->nchunks);

        fill_buckets(sieve);
    }

    return sieve;
}
//-----------------------------------------------------------------------------
void free_siqs_sieve(siqs_sieve_t* sieve)
{
    if (sieve != NULL) {
        free_siqs_poly(sieve->poly);

        free_byte_array(sieve->sieve);
        free_int32_array(sieve->sol1);
        free_int32_array(sieve->sol2);

        if (sieve->use_buckets) {
            free_buckets(sieve->buckets_positive);
            free_buckets(sieve->buckets_negative);
        }
        free(sieve);
    }
}
//-----------------------------------------------------------------------------
ecode_t fill_sieve(siqs_sieve_t* const sieve)
{
    //
    // The sieve is filled by "chunks" whose size should fit in L1 cache
    // memory for better performance. The sieving process is thus made of
    // entangled calls to fill_sieve and scan_sieve where scan_sieve scans
    // the chunks that was just filled before (see the collect_relations
    // function in lib/algo/src/siqs.c).
    //
    // Since positions that survive the sieve are tested for smoothness in
    // batches, we may have to stop scanning a chunk in the middle to perform
    // such a batch. In such a case sieve->scan_begin != sieve->chunk_size
    // (or != sieve->endlast + 1) and we should _not_ fill a new sieve chunk
    // before resuming the scan.
    //

#if ROUND_HALF_WIDTH
    if (sieve->scan_begin != sieve->chunk_size) {
        return SUCCESS;
    }
#else
    bool islast = (ABS(sieve->next_chunkno_to_fill) - 1 == (int)sieve->nchunks);

    if (islast) {
        //
        // The last chunk (on either side of 0) is not completely sieved
        //
        if (sieve->scan_begin != (sieve->endlast + 1)) {
            return SUCCESS;
        }
    } else {
        if (sieve->scan_begin != sieve->chunk_size) {
            return SUCCESS;
        }
    }
#endif

    //
    // On the positive x-side the sieve array represents positions [from, to]
    //
    //     from = (next_chunkno_to_fill - 1) * chunk_size
    //       to = (next_chunkno_to_fill) * chunk_size - 1
    //
    // with next_chunkno_to_fill > 0, in other words sieve[x] is the value
    // of the sieve at position (x + from).
    //
    // However, on the negative x-side the sieve array represents positions
    // [from, to]
    //
    //     from = (next_chunkno_to_fill + 1) * chunk_size + 1
    //       to = (next_chunkno_to_fill) * chunk_size
    //
    // with next_chunkno_to_fill < 0, in other words sieve[x] is the value
    // of the sieve at position -(x + from).
    //
    // Note that we do not sieve position 0 twice since for
    // next_chunkno_to_fill = -1 (the first chunk on the negative side)
    // from == 1 and sieve[0] is thus the value of the sieve at
    // position x = -1.
    //
    sieve->scan_begin = 0;

    ecode_t ecode;
    update_sieve_params(sieve);
    ecode = fill_sieve_chunk(sieve);

    if (sieve->next_chunkno_to_fill > 0) {
        sieve->next_chunkno_to_fill++;
    } else {
        sieve->next_chunkno_to_fill--;
    }
    return ecode;
}
//-----------------------------------------------------------------------------
int compute_threshold_correction(siqs_sieve_t* const sieve)
{
    //
    // Returns the average contribution of the NPRIMES_TO_SKIP first primes
    // since we won't be sieving with them.
    //
    uint32_t logp     = 0;
    uint32_t curprime = 0;
    uint32_t sindex   = 0;

    siqs_poly_t* const poly = sieve->poly;

    const uint32_array_t* const primes     = poly->factor_base;
    const byte_array_t*   const log_primes = sieve->log_primes;

    const uint32_t sieve_length = sieve->sieve->length;

    const uint32_t* const fbdata = primes->data;

    int32_t* const sol1data = sieve->sol1->data;
    int32_t* const sol2data = sieve->sol2->data;

    unsigned char* const sieve_data = sieve->sieve->data;
    unsigned char* const endptr     = &(sieve_data[sieve_length - 1]);
    unsigned char* sieveptr;

    TIFA_MEMSET(sieve_data, 0, sieve_length);

    //
    // Sieve with the prime 2 only for sol1
    //
    sindex   = sol1data[0];
    sieveptr = &(sieve_data[sindex]);
    while (sieveptr <= endptr) {
        *sieveptr += 1;
        sieveptr  += 2;
        sindex    += 2;
    }

    int32_t* sol1ptr = &(sol1data[1]);
    int32_t* sol2ptr = &(sol2data[1]);

    const uint32_t* fbptr = &(fbdata[1]);
    unsigned char* logptr = &(log_primes->data[1]);

    //
    // Sieve with primes > 2 with sol1 and sol2
    //
    for (uint32_t p = 1; p < NPRIMES_TO_SKIP; p++) {
        curprime = *fbptr;
        logp     = *logptr;
        //
        // Sieve with sol1
        //
        sindex   = *sol1ptr;
        sieveptr = &(sieve_data[sindex]);

        while (sieveptr <= endptr) {
            *sieveptr += logp;
            sieveptr  += curprime;
            sindex    += curprime;
        }
        //
        // Sieve with sol2
        //
        sindex   = *sol2ptr;
        sieveptr = &(sieve_data[sindex]);

        while (sieveptr <= endptr) {
            *sieveptr += logp;
            sieveptr  += curprime;
            sindex    += curprime;
        }
        sol1ptr++;
        sol2ptr++;
        fbptr++;
        logptr++;
    }
    //
    // Return average contribution
    //
    int count = 0;
    for (uint32_t i = 0; i < sieve->chunk_size; i++) {
        count += sieve_data[i];
    }
    count /= sieve->chunk_size;

    return count;
}
//-----------------------------------------------------------------------------
inline ecode_t scan_sieve(siqs_sieve_t* const sieve,
                   int32_array_t* const survivors,
                   const uint32_t nsurvivors)
{
    START_EXTRA_TIMER(scan);

    ecode_t ecode;
    if (sieve->next_chunkno_to_fill > 0) {
        ecode = scan_sieve_positive(sieve, survivors, nsurvivors);
    } else {
        ecode = scan_sieve_negative(sieve, survivors, nsurvivors);
    }

    STOP_EXTRA_TIMER(scan);

    return ecode;
}
//-----------------------------------------------------------------------------
ecode_t update_sieve_params(siqs_sieve_t* const sieve)
{
    int32_t nchunks = (int32_t)sieve->nchunks;
    int32_t next    = sieve->next_chunkno_to_fill;

    if (next > nchunks) {
        reset_poly_sols(sieve);
        sieve->next_chunkno_to_fill = -1;
        return SUCCESS;
    }
    if ((-next) > nchunks) {
        switch_polynomial(sieve);
        sieve->next_chunkno_to_fill = 1;
        return SUCCESS;
    }
    return SUCCESS;
}
//-----------------------------------------------------------------------------
void update_sieve_threshold(siqs_sieve_t* const sieve)
{
    int32_t chunkno = sieve->next_chunkno_to_fill;
    int32_t xval    = 0;

    const mpz_ptr a = sieve->poly->a;
    const mpz_ptr b = sieve->poly->b;
    const mpz_ptr c = sieve->poly->c;

    mpz_t rgx;
    mpz_init(rgx);

    xval = (ABS(chunkno) - 1) * sieve->chunk_size + (sieve->chunk_size / 2);

    if (chunkno < 0) {
        xval = -xval;
    }

    if (xval < 0) {
        mpz_mul_ui(rgx, a, (long unsigned int)(-xval));
        mpz_neg(rgx, rgx);
    } else {
        mpz_mul_ui(rgx, a, xval);
    }
    mpz_add(rgx, rgx, b);
    mpz_add(rgx, rgx, b);

    if (xval < 0) {
        mpz_mul_ui(rgx, rgx, (long unsigned int)(-xval));
        mpz_neg(rgx, rgx);
    } else {
        mpz_mul_ui(rgx, rgx, xval);
    }
    mpz_add(rgx, rgx, c);

    sieve->threshold  = mpz_sizeinbase(rgx, 2);

    mpz_clear(rgx);
}
//-----------------------------------------------------------------------------
ecode_t reset_poly_sols(siqs_sieve_t* const sieve)
{
    //
    // Sieving on the positive side is done. Now computes the first hits
    // on the negative side.
    //
    START_EXTRA_TIMER(init_poly);

    siqs_poly_t* const poly = sieve->poly;

    uint32_t* fb    = poly->factor_base->data;

    int32_t* sol1 = sieve->sol1->data;
    int32_t* sol2 = sieve->sol2->data;

    uint32_t  curprime = 0;
    uint32_t* fbptr    = &(fb[0]);

    uint32_t cpsize = sieve->nprimes_no_buckets * sizeof(int32_t);

    memcpy(sol1, poly->sol1->data, cpsize);
    memcpy(sol2, poly->sol2->data, cpsize);

    for (uint32_t i = 0; i < sieve->nprimes_no_buckets; i++) {
        curprime = *fbptr;
        //
        // Remember that for the first chunk on the negative side, sieve[0]
        // is actually the sieve value at position -1. Hence the "- 1" below...
        //
        sol1[i]  = curprime - sol1[i] - 1;
        sol2[i]  = curprime - sol2[i] - 1;

        fbptr++;
    }

    STOP_EXTRA_TIMER(init_poly);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
void switch_polynomial(siqs_sieve_t* const sieve)
{
    //
    // Switch to a new polynomial, keep a copy of the next hit positions
    // and fill the buckets for this new polynomial if needed.
    //
    START_EXTRA_TIMER(init_poly);

    siqs_poly_t* const poly = sieve->poly;

    int32_t* sol1 = sieve->sol1->data;
    int32_t* sol2 = sieve->sol2->data;

    update_polynomial(poly);

    uint32_t cpsize = sieve->nprimes_no_buckets * sizeof(int32_t);

    memcpy(sol1, poly->sol1->data, cpsize);
    memcpy(sol2, poly->sol2->data, cpsize);

    if (sieve->use_buckets) {
        fill_buckets(sieve);
    }

    STOP_EXTRA_TIMER(init_poly);
}
//-----------------------------------------------------------------------------
ecode_t
fill_sieve_chunk(siqs_sieve_t* const sieve)
{
    START_EXTRA_TIMER(fill);

    uint32_t logp     = 0;
    uint32_t curprime = 0;
    uint32_t sindex   = 0;

    siqs_poly_t* const poly = sieve->poly;

    const uint32_array_t* const primes     = poly->factor_base;
    const byte_array_t*   const log_primes = sieve->log_primes;

    const uint32_t sieve_length = sieve->sieve->length;

    const approximer_t* const aximer = poly->approximer;

    const uint32_t imin = aximer->imin;
    const uint32_t imax = aximer->imax;
    const uint32_t nprimes_in_a = poly->nprimes_in_a;

    const uint32_t* const fbdata = primes->data;

    int32_t* const sol1data = sieve->sol1->data;
    int32_t* const sol2data = sieve->sol2->data;

    unsigned char* const sieve_data = sieve->sieve->data;

    uint32_t iend = sieve_length - 1;

#if !ROUND_HALF_WIDTH
    if (ABS(sieve->next_chunkno_to_fill) == sieve->nchunks) {
        iend = sieve->endlast;
    }
#endif

    unsigned char* const endptr = &(sieve_data[iend]);
    unsigned char* sieveptr;
    unsigned char* endunroll;

    //
    // (Re-)initialize the sieve array with the full log approximation...
    //
    TIFA_MEMSET(sieve_data, sieve->threshold, sieve_length);

    if (sieve->use_buckets) {
        //
        // Reports the hits obtained from bucket sieving in the newly
        // initialized sieve.
        //
        empty_bucket(sieve);
    }

    //
    // Sieve with the prime 2 only for sol1 and only if SIEVE_WITH_2 != 0.
    //
#if SIEVE_WITH_2
    sindex   = sol1data[0];
    sieveptr = &(sieve_data[sindex]);
    while (sieveptr <= endptr) {
        *sieveptr -= 1;
        sieveptr  += 2;
        sindex    += 2;
    }
    sol1data[0] = sindex - sieve_length;
#endif

    uint32_t* avoid_primes = poly->idx_of_a;

    uint32_t rimin = MAX(imin, avoid_primes[0]);
    uint32_t rimax = MIN(imax, avoid_primes[nprimes_in_a - 1]);
    uint32_t start = MIN(NPRIMES_TO_SKIP, rimin);

    int32_t* sol1ptr = &(sol1data[start]);
    int32_t* sol2ptr = &(sol2data[start]);

    const uint32_t* fbptr = &(fbdata[start]);
    unsigned char* logptr = &(log_primes->data[start]);

    unsigned char* bckptr = NULL;

    //
    // Sieve with primes > 2 with sol1 and sol2 but skip the primes that
    // divide the coefficient 'a'
    //
    for (uint32_t p = NPRIMES_TO_SKIP; p < rimin; p++) {

        curprime = *fbptr;
        logp     = *logptr;
        //
        // Sieve with sol1
        //
        sindex   = *sol1ptr;
        sieveptr = &(sieve_data[sindex]);
        bckptr   = sieveptr;

        endunroll = endptr - curprime;
        while (sieveptr <= endunroll) {
            *sieveptr -= logp;
            *(sieveptr + curprime) -= logp;
            sieveptr  += (curprime << 1);
        }
        while (sieveptr <= endptr) {
            *sieveptr -= logp;
            sieveptr  += curprime;
        }
        *sol1ptr = sindex + sieveptr - bckptr - sieve_length;

        //
        // Sieve with sol2
        //
        sindex   = *sol2ptr;
        sieveptr = &(sieve_data[sindex]);
        bckptr   = sieveptr;

        endunroll = endptr - curprime;
        while (sieveptr <= endunroll) {
            *sieveptr -= logp;
            *(sieveptr + curprime) -= logp;
            sieveptr  += (curprime << 1);
        }
        while (sieveptr <= endptr) {
            *sieveptr -= logp;
            sieveptr  += curprime;
        }
        *sol2ptr = sindex + sieveptr - bckptr - sieve_length;

        sol1ptr++;
        sol2ptr++;
        fbptr++;
        logptr++;
    }
    avoid_primes++;
    sol1ptr++;
    sol2ptr++;
    fbptr++;
    logptr++;

    for (uint32_t p = rimin + 1; p < rimax; p++) {
        //
        // Skip the primes that divide the coefficient 'a'
        //
        if (p == *avoid_primes) {
            avoid_primes++;
            sol1ptr++;
            sol2ptr++;
            fbptr++;
            logptr++;
            continue;
        }
        curprime = *fbptr;
        logp     = *logptr;
        //
        // Sieve with sol1
        //
        sindex   = *sol1ptr;
        sieveptr = &(sieve_data[sindex]);
        bckptr   = sieveptr;

        while (sieveptr <= endptr) {
            *sieveptr -= logp;
            sieveptr  += curprime;
        }
        *sol1ptr = sindex + sieveptr - bckptr - sieve_length;

        //
        // Sieve with sol2
        //
        sindex   = *sol2ptr;
        sieveptr = &(sieve_data[sindex]);
        bckptr   = sieveptr;

        while (sieveptr <= endptr) {
            *sieveptr -= logp;
            sieveptr  += curprime;
        }
        *sol2ptr = sindex + sieveptr - bckptr - sieve_length;

        sol1ptr++;
        sol2ptr++;
        fbptr++;
        logptr++;
    }
    sol1ptr++;
    sol2ptr++;
    fbptr++;
    logptr++;

    uint32_t npnb = sieve->nprimes_no_buckets;

    for (uint32_t p = rimax + 1; p < npnb; p++) {
        curprime = *fbptr;
        logp     = *logptr;

        //
        // Sieve with sol1
        //
        sindex    = *sol1ptr;
        sieveptr  = &(sieve_data[*sol1ptr]);
        endunroll = endptr - curprime;

        while (sieveptr <= endunroll) {
            *sieveptr -= logp;
            *(sieveptr + curprime) -= logp;
            sieveptr  += (curprime << 1);
            sindex  += (curprime << 1);
        }

        while (sieveptr <= endptr) {
            *sieveptr -= logp;
            sieveptr  += curprime;
            sindex    += curprime;
        }
        *sol1ptr = sindex - sieve_length;

        //
        // Sieve with sol2
        //
        sindex   = *sol2ptr;
        sieveptr = &(sieve_data[*sol2ptr]);

        while (sieveptr <= endunroll) {
            *sieveptr -= logp;
            *(sieveptr + curprime) -= logp;
            sieveptr  += (curprime << 1);
            sindex    += (curprime << 1);
        }

        while (sieveptr <= endptr) {
            *sieveptr -= logp;
            sieveptr  += curprime;
            sindex    += curprime;
        }
        *sol2ptr = sindex - sieve_length;

        sol1ptr++;
        sol2ptr++;
        fbptr++;
        logptr++;
    }

    STOP_EXTRA_TIMER(fill);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
ecode_t scan_sieve_positive(siqs_sieve_t* const siqs_sieve,
                            int32_array_t* const survivors,
                            const uint32_t nsurvivors)
{
    const byte_array_t*  const sieve      = siqs_sieve->sieve;
    const unsigned char* const sieve_data = sieve->data;

    int32_t*  const xpool = survivors->data;
    uint32_t* const scan_begin = &siqs_sieve->scan_begin;

    int32_t chunkno = siqs_sieve->next_chunkno_to_fill - 1;

    const uint32_t sieve_length = sieve->length;
    const int32_t  sieve_begin  = (chunkno - 1) * sieve_length;

    uint32_t iend = sieve_length - 1;

#if !ROUND_HALF_WIDTH
    if (ABS(chunkno) == siqs_sieve->nchunks) {
        iend = siqs_sieve->endlast;
    }
#endif

    const int32_t  sieve_end = sieve_begin + iend;

    uint32_t tmp_length = survivors->length;
    uint32_t xsieve     = *scan_begin;
    int32_t  xindex     = *scan_begin + sieve_begin;

    int32_t* xpool_to_write = &(xpool[tmp_length]);
    int32_t* xpool_to_stop = &(xpool[nsurvivors]);

    //
    // Make sure &sieve_data[xsieve] points to the beginning of a 64 bit word
    // (because we may resume a previous, incomplete scan)
    //
    while ((xsieve & 7) && (xsieve <= iend)) {
        if (sieve_data[xsieve] & 0x80) {
            *xpool_to_write = (int32_t)xindex;
            xpool_to_write++;
            if (xpool_to_write == xpool_to_stop) {
                xindex++;
                goto STOP_AND_RETURN;
            }
        }
        xindex++;
        xsieve++;
    }
    //
    // The usual trick: casting the sieve array as a 64 bit integers array
    // to read and compare 8 sieve entries in a row.
    //
    uint64_t* sieve64 = (uint64_t*)&sieve_data[xsieve];

    const unsigned char* sdatptr = &sieve_data[xsieve];

    uint32_t* sieve32;

    const int32_t sieve_end2 = sieve_end - 8;

    while (xindex <= sieve_end2) {
        //
        // A byte b is < 0 iif (b & 0x80) != 0
        //
        if ((*sieve64 & 0x8080808080808080) == 0) {
            xindex  += 8;
            sieve64++;
            continue;
        }
        sieve32 = (uint32_t*)sieve64;

        if ((*sieve32 & 0x80808080) != 0) {
            sdatptr = (const unsigned char*)sieve32;
            for (uint8_t i = 0; i < 4; i++) {
                if ((*sdatptr & 0x80) != 0) {
                    *xpool_to_write = (int32_t)xindex;
                    xpool_to_write++;
                    if (xpool_to_write == xpool_to_stop) {
                        xindex++;
                        goto STOP_AND_RETURN;
                    }
                }
                sdatptr++;
                xindex++;
            }
        } else {
            xindex += 4;
        }
        sieve32++;

        if ((*sieve32 & 0x80808080) != 0) {
            sdatptr = (const unsigned char*)sieve32;
            for (uint8_t i = 0; i < 4; i++) {
                if ((*sdatptr & 0x80) != 0) {
                    *xpool_to_write = (int32_t)xindex;
                    xpool_to_write++;
                    if (xpool_to_write == xpool_to_stop) {
                        xindex++;
                        goto STOP_AND_RETURN;
                    }
                }
                sdatptr++;
                xindex++;
            }
        } else {
            xindex += 4;
        }
        sieve64++;
    }

    while (xindex <= sieve_end) {
        if (sieve_data[xindex - sieve_begin] & 0x80) {
            *xpool_to_write = (int32_t)xindex;
            xpool_to_write++;
            if (xpool_to_write == xpool_to_stop) {
                xindex++;
                break;
            }
        }
        xindex++;
    }

  STOP_AND_RETURN:

    *scan_begin = xindex - sieve_begin;
    survivors->length = xpool_to_write - &(xpool[0]);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
ecode_t scan_sieve_negative(siqs_sieve_t* const siqs_sieve,
                            int32_array_t* const survivors,
                            uint32_t nsurvivors)
{
    const byte_array_t*  const sieve      = siqs_sieve->sieve;
    const unsigned char* const sieve_data = sieve->data;

    int32_t*  const xpool = survivors->data;
    uint32_t* const scan_begin = &siqs_sieve->scan_begin;

    int32_t chunkno = siqs_sieve->next_chunkno_to_fill + 1;

    const uint32_t sieve_length = sieve->length;
    const int32_t  sieve_begin  = -((chunkno + 1) * sieve_length - 1);

    uint32_t iend = sieve_length - 1;

#if !ROUND_HALF_WIDTH
    if (ABS(chunkno) == siqs_sieve->nchunks) {
       iend = siqs_sieve->endlast;
    }
#endif

    const int32_t  sieve_end = sieve_begin + iend;

    uint32_t tmp_length = survivors->length;
    uint32_t xsieve     = *scan_begin;
    int32_t  xindex     = *scan_begin + sieve_begin;

    int32_t* xpool_to_write = &(xpool[tmp_length]);
    int32_t* xpool_to_stop = &(xpool[nsurvivors]);

    //
    // Make sure &sieve_data[xsieve] points to the beginning of a 64 bit word
    // (because we may resume a previous, incomplete scan)
    //
    while ((xsieve & 7) && (xsieve <= iend)) {
        if (sieve_data[xsieve] & 0x80) {
            *xpool_to_write = -xindex;
            xpool_to_write++;
            if (xpool_to_write == xpool_to_stop) {
                xindex++;
                goto STOP_AND_RETURN;
            }
        }
        xindex++;
        xsieve++;
    }
    //
    // The usual trick: casting the sieve array as a 64 bit integers array
    // to read and compare 8 sieve entries in a row.
    //
    uint64_t* sieve64 = (uint64_t*)&sieve_data[xsieve];

    const unsigned char* sdatptr = &sieve_data[xsieve];

    uint32_t* sieve32;

    const int32_t  sieve_end2 = sieve_end - 8;

    while (xindex <= sieve_end2) {
        if ((*sieve64 & 0x8080808080808080) == 0) {
            xindex  += 8;
            sieve64++;
            continue;
        }
        sieve32 = (uint32_t*)sieve64;

        if ((*sieve32 & 0x80808080) != 0) {
            sdatptr = (const unsigned char*)sieve32;
            for (uint8_t i = 0; i < 4; i++) {
                if ((*sdatptr & 0x80) != 0) {
                    *xpool_to_write = -xindex;
                    xpool_to_write++;
                    if (xpool_to_write == xpool_to_stop) {
                        xindex++;
                        goto STOP_AND_RETURN;
                    }
                }
                sdatptr++;
                xindex++;
            }
        } else {
            xindex  += 4;
        }
        sieve32++;

        if ((*sieve32 & 0x80808080) != 0) {
            sdatptr = (const unsigned char*)sieve32;
            for (uint8_t i = 0; i < 4; i++) {
                if ((*sdatptr & 0x80) != 0) {
                    *xpool_to_write = -xindex;
                    xpool_to_write++;
                    if (xpool_to_write == xpool_to_stop) {
                        xindex++;
                        goto STOP_AND_RETURN;
                    }
                }
                sdatptr++;
                xindex++;
            }
        } else {
            xindex  += 4;
        }
        sieve64++;
    }

    while (xindex <= sieve_end) {
        if (sieve_data[sieve_end - xindex] & 0x80) {
            *xpool_to_write = -xindex;
            xpool_to_write++;
            if (xpool_to_write == xpool_to_stop) {
                xindex++;
                goto STOP_AND_RETURN;
            }
        }
        xindex++;
    }

  STOP_AND_RETURN:

    *scan_begin = xindex - sieve_begin;
    survivors->length = xpool_to_write - &(xpool[0]);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
void fill_buckets(siqs_sieve_t* const sieve) {
    //
    // Fills all the buckets for the whole sieving interval.
    //
    // Hits on the positive side are stored in sieve->buckets_positive
    // such that buckets_positive->bin[i] covers the sieve interval
    // [i*chunk_size, (i+1)*chunk_size - 1]
    //
    // Hits on the negative side are stored in sieve->buckets_negative
    // such that buckets_negative->bin[i] covers the sieve interval
    // [-(i+1)*chunk_size, -i*chunk_size + 1]
    //
    reset_buckets(sieve->buckets_positive);
    reset_buckets(sieve->buckets_negative);

    const siqs_poly_t*  const poly   = sieve->poly;
    const uint32_t*     const fb     = poly->factor_base->data;
    const byte_array_t* const log_fb = sieve->log_primes;

    const uint32_t fb_length = poly->factor_base->length;
    const uint32_t ipmin     = sieve->nprimes_no_buckets;

    const int32_t* const sol1 = poly->sol1->data;
    const int32_t* const sol2 = poly->sol2->data;

    uint32_t curprime = 0;
    uint8_t  logp;

    int32_t  hw = sieve->nchunks * sieve->chunk_size;
    int32_t  spos;

    buckets_t* bposi = sieve->buckets_positive;
    buckets_t* bnega = sieve->buckets_negative;

    int32_t ibin = 0;
    int32_t rpos = 0;

    const uint32_t* fbptr = &(fb[ipmin]);
    const int32_t*  s1ptr = &(sol1[ipmin]);
    const int32_t*  s2ptr = &(sol2[ipmin]);
    const uint8_t*  lgptr = &(log_fb->data[ipmin]);

    int32_t i = fb_length - ipmin;

    do {
        i--;
        //
        // Sieve with solution 1
        //
        curprime = *fbptr;
        spos = *s1ptr;
        logp = *lgptr;

        while (spos < hw) {
            ibin = DIV_CHUNK_LEN(spos);
            rpos = MOD_CHUNK_LEN(spos);
            add_to_buckets(rpos, logp, ibin, bposi);
            spos += curprime;
        }
        spos = curprime - *s1ptr - 1;

        while (spos < hw) {
            ibin = DIV_CHUNK_LEN(spos);
            rpos = MOD_CHUNK_LEN(spos);
            add_to_buckets(rpos, logp, ibin, bnega);
            spos += curprime;
        }

        //
        // Sieve with solution 2
        //
        spos = *s2ptr;
        while (spos < hw) {
            ibin = DIV_CHUNK_LEN(spos);
            rpos = MOD_CHUNK_LEN(spos);

            add_to_buckets(rpos, logp, ibin, bposi);
            spos += curprime;
        }
        spos = curprime - *s2ptr - 1;

        while (spos < hw) {
            ibin = DIV_CHUNK_LEN(spos);
            rpos = MOD_CHUNK_LEN(spos);

            add_to_buckets(rpos, logp, ibin, bnega);
            spos += curprime;
        }
        fbptr++;
        s1ptr++;
        s2ptr++;
        lgptr++;
    } while (i > 0);
}
//-----------------------------------------------------------------------------
static ecode_t empty_bucket(siqs_sieve_t* const sieve) {
    //
    // Empty the buckets on the current sieve chunk...
    //
    if (sieve->next_chunkno_to_fill > 0) {
        return empty_bucket_posi(sieve);
    } else {
        return empty_bucket_nega(sieve);
    }
}
//------------------------------------------------------------------------------
static inline ecode_t empty_bucket_posi(siqs_sieve_t* const sieve) {
    uint32_t ibin     = sieve->next_chunkno_to_fill - 1;

    unsigned char* sieve_data = sieve->sieve->data;

    ecode_t ecode = SUCCESS;

    buckets_t*      bposi = sieve->buckets_positive;
    uint32_array_t* bin   = bposi->bins[ibin];
    uint32_t* const bdata = bin->data;
    uint32_t        blen  = bin->length;

    uint32_t imod = blen & 7;
    uint32_t idiv = blen >> 3;
    uint32_t j = 0;

    for (uint32_t i = 0; i < idiv; i++) {
        sieve_data[GET_A(bdata[j])]     -= GET_B(bdata[j]);
        sieve_data[GET_A(bdata[j + 1])] -= GET_B(bdata[j + 1]);
        sieve_data[GET_A(bdata[j + 2])] -= GET_B(bdata[j + 2]);
        sieve_data[GET_A(bdata[j + 3])] -= GET_B(bdata[j + 3]);
        sieve_data[GET_A(bdata[j + 4])] -= GET_B(bdata[j + 4]);
        sieve_data[GET_A(bdata[j + 5])] -= GET_B(bdata[j + 5]);
        sieve_data[GET_A(bdata[j + 6])] -= GET_B(bdata[j + 6]);
        sieve_data[GET_A(bdata[j + 7])] -= GET_B(bdata[j + 7]);
        j += 8;
    }
    for (uint32_t i = 0; i < imod; i++) {
        sieve_data[GET_A(bdata[j])] -= GET_B(bdata[j]);
        j++;
    }
    bin->length = 0;

    return ecode;
}
//------------------------------------------------------------------------------
static inline ecode_t empty_bucket_nega(siqs_sieve_t* const sieve) {
    uint32_t ibin   = -(sieve->next_chunkno_to_fill + 1);

    unsigned char* const sieve_data = sieve->sieve->data;

    ecode_t ecode = SUCCESS;

    buckets_t*      bnega = sieve->buckets_negative;
    uint32_array_t* bin   = bnega->bins[ibin];
    uint32_t* const bdata = bin->data;
    uint32_t        blen  = bin->length;

    uint32_t imod = blen & 7;
    uint32_t idiv = blen >> 3;
    uint32_t j = 0;

    for (uint32_t i = 0; i < idiv; i++) {
        sieve_data[GET_A(bdata[j])]     -= GET_B(bdata[j]);
        sieve_data[GET_A(bdata[j + 1])] -= GET_B(bdata[j + 1]);
        sieve_data[GET_A(bdata[j + 2])] -= GET_B(bdata[j + 2]);
        sieve_data[GET_A(bdata[j + 3])] -= GET_B(bdata[j + 3]);
        sieve_data[GET_A(bdata[j + 4])] -= GET_B(bdata[j + 4]);
        sieve_data[GET_A(bdata[j + 5])] -= GET_B(bdata[j + 5]);
        sieve_data[GET_A(bdata[j + 6])] -= GET_B(bdata[j + 6]);
        sieve_data[GET_A(bdata[j + 7])] -= GET_B(bdata[j + 7]);
        j += 8;
    }
    for (uint32_t i = 0; i < imod; i++) {
        sieve_data[GET_A(bdata[j])] -= GET_B(bdata[j]);
        j++;
    }
    bin->length = 0;

    return ecode;
}
//------------------------------------------------------------------------------
void print_init_poly_timing(siqs_sieve_t* const sieve)
{
    PRINT_INIT_POLY_TIMING_MSG;
}
//-----------------------------------------------------------------------------
void print_fill_timing(siqs_sieve_t* const sieve)
{
    PRINT_FILL_TIMING_MSG;
}
//-----------------------------------------------------------------------------
void print_scan_timing(siqs_sieve_t* const sieve)
{
    PRINT_SCAN_TIMING_MSG;
}
//-----------------------------------------------------------------------------
void set_siqs_sieve_threshold(siqs_sieve_t* const sieve, uint32_t threshold)
{
    //
    // For the time being, the threshold is part of the algorithm's
    // parameters. Once more tests have been performed, the threshold will
    // be determined automatically.
    //
    sieve->threshold = threshold;
}
//------------------------------------------------------------------------------
