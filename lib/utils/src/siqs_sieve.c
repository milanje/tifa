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

//
// sysconf's L1 cache query is a glibc extension, hidden in strict C99 mode.
//
#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif
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
#include "tifa_alloc.h"

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
//-----------------------------------------------------------------------------
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
ecode_t reset_poly_sols(siqs_sieve_t* const sieve);
void    switch_polynomial(siqs_sieve_t* const sieve);
void    fill_buckets(siqs_sieve_t* const sieve);
static void size_buckets(siqs_sieve_t* const sieve);
static inline ecode_t empty_bucket(siqs_sieve_t* const sieve);
static inline ecode_t empty_bucket_posi(siqs_sieve_t* const sieve);
static inline ecode_t empty_bucket_nega(siqs_sieve_t* const sieve);
//-----------------------------------------------------------------------------
#if defined(LOG2_SIEVE_CHUNK_MAX_SIZE)
static uint32_t log2_sieve_chunk_size(void)
{
    return LOG2_SIEVE_CHUNK_MAX_SIZE;
}
#else
static uint32_t       log2_chunk_size      = LOG2_SIEVE_CHUNK_MIN;
static pthread_once_t log2_chunk_size_once = PTHREAD_ONCE_INIT;

static void find_log2_sieve_chunk_size(void)
{
    //
    // On CPUs mixing core types the system reports the smallest L1 (64 KB on
    // Apple's efficiency cores against 128 KB on its performance cores), so
    // a chunk stays in cache whichever core the thread runs on.
    //
    unsigned long l1d = 0;
#if defined(__APPLE__)
    uint64_t value  = 0;
    size_t   length = sizeof(value);
    if (0 == sysctlbyname("hw.l1dcachesize", &value, &length, NULL, 0)) {
        l1d = (unsigned long) value;
    }
#elif defined(_SC_LEVEL1_DCACHE_SIZE)
    long value = sysconf(_SC_LEVEL1_DCACHE_SIZE);
    if (value > 0) {
        l1d = (unsigned long) value;
    }
#endif
    uint32_t log2 = LOG2_SIEVE_CHUNK_MIN;
    while ((log2 < LOG2_SIEVE_CHUNK_MAX) && ((1UL << (log2 + 1)) <= l1d)) {
        log2++;
    }
    log2_chunk_size = log2;
}

static uint32_t log2_sieve_chunk_size(void)
{
    //
    // Once per process: sysctlbyname took 1% of a 70-bit SIQS.
    //
    pthread_once(&log2_chunk_size_once, find_log2_sieve_chunk_size);
    return log2_chunk_size;
}
#endif
//
// Bytes past the end of the sieve array, for sieve_primes.
//
#define SIEVE_SLACK 64
//-----------------------------------------------------------------------------
siqs_sieve_t* alloc_siqs_sieve(mpz_t n,
                               uint32_array_t* const factor_base,
                               byte_array_t*   const log_primes,
                               uint32_array_t* const sqrtm_pi,
                               uint32_t half_width,
                               uint32_t extra_primes_in_a)
{
    siqs_sieve_t* const sieve = tifa_malloc(sizeof(siqs_sieve_t));

    INIT_EXTRA_TIMER(init_poly);
    INIT_EXTRA_TIMER(fill);
    INIT_EXTRA_TIMER(scan);

    sieve->log_primes = log_primes;

    sieve->log2_chunk_size   = log2_sieve_chunk_size();
    const uint32_t chunk_max = 1U << sieve->log2_chunk_size;

    if (half_width <= chunk_max) {
        sieve->sieve         = alloc_byte_array(half_width + SIEVE_SLACK);
        sieve->sieve->length = half_width;
        sieve->chunk_size    = half_width;
        sieve->nchunks       = 1;

        sieve->endlast       = half_width - 1;

    } else {
        sieve->sieve         = alloc_byte_array(chunk_max + SIEVE_SLACK);
        sieve->sieve->length = chunk_max;

        uint32_t hwmod = half_width % chunk_max;
        sieve->nchunks = half_width / chunk_max;

        //
        // Use the exact half width provided as argument
        //
        if (hwmod != 0) {
            sieve->nchunks++;
            sieve->endlast = hwmod - 1;
        } else {
            sieve->endlast = chunk_max - 1;
        }

        sieve->chunk_size = chunk_max;
    }
    //
    // Setting scan_begin = chunk_size means that a new chunk should be filled
    //
    sieve->scan_begin = sieve->chunk_size;
    sieve->clean      = false;

    sieve->next_chunkno_to_fill = 1;

    mpz_t target_a;
    mpz_init(target_a);

    mpz_mul_ui(target_a, n, 2);
    mpz_sqrt(target_a, target_a);

    mpz_fdiv_q_ui(target_a, target_a, half_width);

    sieve->poly = alloc_siqs_poly(target_a, n, factor_base, sqrtm_pi,
                                  extra_primes_in_a);

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

        size_buckets(sieve);
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
    // [-to, -from]
    //
    //     from = -(next_chunkno_to_fill + 1) * chunk_size + 1
    //       to = -(next_chunkno_to_fill) * chunk_size
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
//
// Aligned to 64 bytes, as fill_sieve_chunk: otherwise their loops move with
// the code linked before them, for 1% to 2% on SIQS.
//
__attribute__((aligned(64))) inline ecode_t scan_sieve(siqs_sieve_t* const sieve,
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
//
// Sieves positions 0 to iend with the primes of indexes from to to - 1, both
// roots in each pass, and leaves their next hits relative to the following
// chunk of sieve_length positions. Each root of a prime p below iend + 1 hits
// m = (iend + 1) / p times, or m + 1: consecutive primes share m, so the
// loop on the first m hits keeps its count from one prime to the next, and
// the last hit is subtracted without branching, from a byte of the slack
// past the end of the sieve when it falls beyond.
//
static inline void sieve_primes(unsigned char* const sieve_data,
                                const uint32_t iend,
                                const uint32_t sieve_length,
                                const uint32_t* const fbdata,
                                const uint8_t* const logs,
                                int32_t* const sol1data,
                                int32_t* const sol2data,
                                const uint32_t from, const uint32_t to) {
    const uint32_t length = iend + 1;

    uint32_t p = from;
    while ((p < to) && (fbdata[p] <= iend)) {
        const uint32_t m    = length / fbdata[p];
        const uint32_t pmax = length / m;

        for (; (p < to) && (fbdata[p] <= pmax); p++) {
            const uint32_t prime = fbdata[p];
            const uint8_t  lp    = logs[p];

            uint32_t r1 = (uint32_t) sol1data[p];
            uint32_t r2 = (uint32_t) sol2data[p];
            uint32_t j  = m;
            for (; j >= 2; j -= 2) {
                sieve_data[r1]         -= lp;
                sieve_data[r2]         -= lp;
                sieve_data[r1 + prime] -= lp;
                sieve_data[r2 + prime] -= lp;
                r1 += 2 * prime;
                r2 += 2 * prime;
            }
            if (j != 0) {
                sieve_data[r1] -= lp;
                sieve_data[r2] -= lp;
                r1 += prime;
                r2 += prime;
            }
            const uint32_t spill = length + (prime & (SIEVE_SLACK / 2 - 1));
            sieve_data[(r1 <= iend) ? r1 : spill] -= lp;
            sieve_data[(r2 <= iend) ? r2 : spill + SIEVE_SLACK / 2] -= lp;
            r1 += (r1 <= iend) ? prime : 0;
            r2 += (r2 <= iend) ? prime : 0;

            sol1data[p] = (int32_t) r1 - (int32_t) sieve_length;
            sol2data[p] = (int32_t) r2 - (int32_t) sieve_length;
        }
    }
    for (; p < to; p++) {
        const uint32_t prime = fbdata[p];
        const uint8_t  lp    = logs[p];

        uint32_t r1 = (uint32_t) sol1data[p];
        uint32_t r2 = (uint32_t) sol2data[p];
        const uint32_t spill = length + (prime & (SIEVE_SLACK / 2 - 1));
        sieve_data[(r1 <= iend) ? r1 : spill] -= lp;
        sieve_data[(r2 <= iend) ? r2 : spill + SIEVE_SLACK / 2] -= lp;
        r1 += (r1 <= iend) ? prime : 0;
        r2 += (r2 <= iend) ? prime : 0;

        sol1data[p] = (int32_t) r1 - (int32_t) sieve_length;
        sol2data[p] = (int32_t) r2 - (int32_t) sieve_length;
    }
}
//-----------------------------------------------------------------------------
//
// The smallest primes are not sieved, but fewer of them for small factor
// bases: 10 below 150 primes (up to about 85 bits) and 14 below 450 (up to
// about 130 bits), rather than 20. Their hits cost less than the false
// candidates the missing logs let through, the thresholds raised by one
// where it was measured best.
//
static inline uint32_t nprimes_to_skip(const uint32_t fblength) {
    if (fblength < 150) {
        return 10;
    }
    return (fblength < 450) ? 14 : NPRIMES_TO_SKIP;
}
//-----------------------------------------------------------------------------
__attribute__((aligned(64))) ecode_t
fill_sieve_chunk(siqs_sieve_t* const sieve)
{
    START_EXTRA_TIMER(fill);

    siqs_poly_t* const poly = sieve->poly;

    const uint32_array_t* const primes     = poly->factor_base;
    const byte_array_t*   const log_primes = sieve->log_primes;

    const uint32_t sieve_length = sieve->sieve->length;

    const uint32_t imin = poly->approximer->imin;
    const uint32_t nprimes_in_a = poly->nprimes_in_a;

    const uint32_t* const fbdata = primes->data;

    int32_t* const sol1data = sieve->sol1->data;
    int32_t* const sol2data = sieve->sol2->data;

    unsigned char* const sieve_data = sieve->sieve->data;

    uint32_t iend = sieve_length - 1;

    if ((uint32_t)ABS(sieve->next_chunkno_to_fill) == sieve->nchunks) {
        iend = sieve->endlast;
    }


    //
    // (Re-)initialize every byte of the sieve array to the threshold, unless
    // the scan of the previous chunk did.
    //
    if (!sieve->clean) {
        TIFA_MEMSET(sieve_data, sieve->threshold, sieve_length);
    }
    sieve->clean = false;

    if (sieve->use_buckets) {
        //
        // Reports the hits obtained from bucket sieving in the newly
        // initialized sieve.
        //
        empty_bucket(sieve);
    }

    //
    // Sieve with the primes from the first one not skipped, but for those
    // dividing the coefficient 'a': on the ranges between them.
    //
    const uint32_t* const avoid = poly->idx_of_a;

    uint32_t rimin = MAX(imin, avoid[0]);
    const uint32_t skip = nprimes_to_skip(primes->length);
    uint32_t from  = (rimin < skip) ? rimin + 1 : skip;

    for (uint32_t k = 0; k <= nprimes_in_a; k++) {
        const uint32_t to = (k < nprimes_in_a)
                          ? MIN(avoid[k], sieve->nprimes_no_buckets)
                          : sieve->nprimes_no_buckets;
        if (to > from) {
            sieve_primes(sieve_data, iend, sieve_length, fbdata,
                         log_primes->data, sol1data, sol2data, from, to);
        }
        if (to + 1 > from) {
            from = to + 1;
        }
    }
    STOP_EXTRA_TIMER(fill);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
#if defined(__ARM_NEON) && defined(__BYTE_ORDER__) \
    && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    #define SCAN_WITH_NEON 1
    #include <arm_neon.h>
//
// The top bits of the 64 bytes of a, b, c and d, byte i giving bit i.
//
static inline uint64_t top_bits(const uint8x16_t a, const uint8x16_t b,
                                const uint8x16_t c, const uint8x16_t d) {
    static const uint8_t weights[16] = {
        1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128
    };
    const uint8x16_t w  = vld1q_u8(weights);
    const uint8x16_t ta = vandq_u8(vcltzq_s8(vreinterpretq_s8_u8(a)), w);
    const uint8x16_t tb = vandq_u8(vcltzq_s8(vreinterpretq_s8_u8(b)), w);
    const uint8x16_t tc = vandq_u8(vcltzq_s8(vreinterpretq_s8_u8(c)), w);
    const uint8x16_t td = vandq_u8(vcltzq_s8(vreinterpretq_s8_u8(d)), w);
    uint8x16_t sum = vpaddq_u8(vpaddq_u8(ta, tb), vpaddq_u8(tc, td));
    sum = vpaddq_u8(sum, sum);
    return vgetq_lane_u64(vreinterpretq_u64_u8(sum), 0);
}
#else
    #define SCAN_WITH_NEON 0
#endif
//
// Appends sign.(sieve_begin + i) to the survivors for the positions i of
// [*scan_begin, iend] whose sieve byte has its top bit set, until there are
// nsurvivors of them, and leaves *scan_begin past the last position scanned.
// Survivors are rare: 64 bytes are tested at once with NEON, 32 otherwise.
// The bytes scanned are reset to the threshold: those before *scan_begin
// always are, so that a scan stopped in the middle of a block loses nothing.
//
static inline void scan_chunk(unsigned char* const data,
                              const uint32_t iend, const int32_t sieve_begin,
                              const int32_t sign, uint32_t* const scan_begin,
                              int32_array_t* const survivors,
                              const uint32_t nsurvivors,
                              const unsigned char threshold) {
    int32_t*       out  = &(survivors->data[survivors->length]);
    int32_t* const stop = &(survivors->data[nsurvivors]);
    uint32_t       i    = *scan_begin;

#if SCAN_WITH_NEON
    const uint8x16_t t = vdupq_n_u8(threshold);
    for (; i + 64 <= iend + 1; i += 64) {
        const uint8x16_t a = vld1q_u8(data + i);
        const uint8x16_t b = vld1q_u8(data + i + 16);
        const uint8x16_t c = vld1q_u8(data + i + 32);
        const uint8x16_t d = vld1q_u8(data + i + 48);
        if (vmaxvq_u8(vorrq_u8(vorrq_u8(a, b), vorrq_u8(c, d))) >= 0x80) {
            uint64_t mask = top_bits(a, b, c, d);
            while (mask != 0) {
                const uint32_t j = i + (uint32_t) __builtin_ctzll(mask);
                *out = sign * (sieve_begin + (int32_t) j);
                out++;
                if (out == stop) {
                    TIFA_MEMSET((data + i), threshold, (j + 1 - i));
                    i = j + 1;
                    goto STOP;
                }
                mask &= mask - 1;
            }
        }
        vst1q_u8(data + i,      t);
        vst1q_u8(data + i + 16, t);
        vst1q_u8(data + i + 32, t);
        vst1q_u8(data + i + 48, t);
    }
#else
    const uint64_t tops = UINT64_C(0x8080808080808080);
    const uint64_t t8   = UINT64_C(0x0101010101010101) * threshold;
    const uint64_t t32[4] = { t8, t8, t8, t8 };
    while (i + 32 <= iend + 1) {
        uint64_t w[4];
        memcpy(w, data + i, sizeof(w));
        if (((w[0] | w[1] | w[2] | w[3]) & tops) == 0) {
            memcpy(data + i, t32, sizeof(t32));
            i += 32;
            continue;
        }
        if ((w[0] & tops) != 0) {
            for (uint32_t j = i; j < i + 8; j++) {
                if (data[j] & 0x80) {
                    *out = sign * (sieve_begin + (int32_t) j);
                    out++;
                    if (out == stop) {
                        TIFA_MEMSET((data + i), threshold, (j + 1 - i));
                        i = j + 1;
                        goto STOP;
                    }
                }
            }
        }
        memcpy(data + i, t32, sizeof(uint64_t));
        i += 8;
    }
#endif
    for (; i <= iend; i++) {
        const unsigned char byte = data[i];
        data[i] = threshold;
        if (byte & 0x80) {
            *out = sign * (sieve_begin + (int32_t) i);
            out++;
            if (out == stop) {
                i++;
                break;
            }
        }
    }
  STOP:
    *scan_begin = i;
    survivors->length = (uint32_t) (out - &(survivors->data[0]));
}
//-----------------------------------------------------------------------------
//
// Once a chunk is scanned to iend, it holds the threshold everywhere but past
// iend (the last chunk's tail, which buckets may have hit): fill_sieve_chunk
// need not reset it.
//
static inline void end_of_scan(siqs_sieve_t* const sieve, const uint32_t iend) {
    if (sieve->scan_begin != iend + 1) {
        return;
    }
    const uint32_t length = sieve->sieve->length;
    if (iend + 1 < length) {
        unsigned char* const tail = sieve->sieve->data + iend + 1;
        TIFA_MEMSET(tail, sieve->threshold, (length - iend - 1));
    }
    sieve->clean = true;
}
//-----------------------------------------------------------------------------
ecode_t scan_sieve_positive(siqs_sieve_t* const siqs_sieve,
                            int32_array_t* const survivors,
                            const uint32_t nsurvivors)
{
    const int32_t  chunkno      = siqs_sieve->next_chunkno_to_fill - 1;
    const uint32_t sieve_length = siqs_sieve->sieve->length;
    const int32_t  sieve_begin  = (chunkno - 1) * (int32_t) sieve_length;

    uint32_t iend = sieve_length - 1;
    if ((uint32_t)ABS(chunkno) == siqs_sieve->nchunks) {
        iend = siqs_sieve->endlast;
    }
    scan_chunk(siqs_sieve->sieve->data, iend, sieve_begin, 1,
               &siqs_sieve->scan_begin, survivors, nsurvivors,
               (unsigned char) siqs_sieve->threshold);
    end_of_scan(siqs_sieve, iend);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
ecode_t scan_sieve_negative(siqs_sieve_t* const siqs_sieve,
                            int32_array_t* const survivors,
                            uint32_t nsurvivors)
{
    const int32_t  chunkno      = siqs_sieve->next_chunkno_to_fill + 1;
    const uint32_t sieve_length = siqs_sieve->sieve->length;
    const int32_t  sieve_begin  = -((chunkno + 1) * (int32_t) sieve_length - 1);

    uint32_t iend = sieve_length - 1;
    if ((uint32_t)ABS(chunkno) == siqs_sieve->nchunks) {
        iend = siqs_sieve->endlast;
    }
    scan_chunk(siqs_sieve->sieve->data, iend, sieve_begin, -1,
               &siqs_sieve->scan_begin, survivors, nsurvivors,
               (unsigned char) siqs_sieve->threshold);
    end_of_scan(siqs_sieve, iend);

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static void size_buckets(siqs_sieve_t* const sieve) {
    //
    // Makes room in every bin for the most hits of a polynomial, and one more
    // for the last hit fill_buckets writes without branching: no bin grows
    // while it fills them.
    //
    const uint32_t* const fb    = sieve->poly->factor_base->data;
    const uint32_t        fblen = sieve->poly->factor_base->length;

    uint32_t bound = 1;
    for (uint32_t i = sieve->nprimes_no_buckets; i < fblen; i++) {
        bound += 2 * ((sieve->chunk_size - 1) / fb[i] + 1);
    }
    buckets_t* const sides[2] = {
        sieve->buckets_positive, sieve->buckets_negative
    };
    for (int s = 0; s < 2; s++) {
        for (uint32_t i = 0; i < sides[s]->nbins; i++) {
            uint32_array_t* const bin = sides[s]->bins[i];
            if (bin->alloced < bound) {
                bin->data    = tifa_realloc(bin->data, bound * sizeof(uint32_t));
                bin->alloced = bound;
            }
        }
    }
}
//-----------------------------------------------------------------------------
//
// Bins whose write positions fill_buckets keeps on the stack.
//
#define FILL_LOCAL_BINS 64

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
    // [-(i+1)*chunk_size, -i*chunk_size - 1]
    //
    // As in sieve_primes, each root of a prime p hits the hw positions of a
    // side m = hw/p or m + 1 times: m in a loop whose count consecutive
    // primes share, the last one written without branching (past the end of
    // a bin when it falls beyond, see size_buckets).
    //
    const siqs_poly_t*  const poly   = sieve->poly;
    const uint32_t*     const fb     = poly->factor_base->data;
    const byte_array_t* const log_fb = sieve->log_primes;

    const uint32_t fb_length = poly->factor_base->length;
    const uint32_t nbins     = sieve->nchunks;

    const int32_t* const sol1 = poly->sol1->data;
    const int32_t* const sol2 = poly->sol2->data;

    const uint32_t hw    = sieve->nchunks * sieve->chunk_size;
    const uint32_t log2c = sieve->log2_chunk_size;
    const uint32_t maskc = sieve->chunk_size - 1;

    buckets_t* const bposi = sieve->buckets_positive;
    buckets_t* const bnega = sieve->buckets_negative;

    uint32_t*  local[2 * FILL_LOCAL_BINS];
    uint32_t** const cur = (nbins <= FILL_LOCAL_BINS)
                         ? local
                         : tifa_malloc(2 * nbins * sizeof(uint32_t*));
    uint32_t** const curp = cur;
    uint32_t** const curn = cur + nbins;

    for (uint32_t b = 0; b < nbins; b++) {
        curp[b] = bposi->bins[b]->data;
        curn[b] = bnega->bins[b]->data;
    }
    uint32_t i = sieve->nprimes_no_buckets;
    while (i < fb_length) {
        const uint32_t m    = hw / fb[i];
        const uint32_t pmax = (m != 0) ? hw / m : UINT32_MAX;

        for (; (i < fb_length) && (fb[i] <= pmax); i++) {
            const uint32_t p  = fb[i];
            const uint32_t lp = log_fb->data[i];

            const uint32_t starts[4] = {
                (uint32_t) sol1[i], p - (uint32_t) sol1[i] - 1,
                (uint32_t) sol2[i], p - (uint32_t) sol2[i] - 1
            };
            for (int s = 0; s < 4; s++) {
                uint32_t** const c = (s & 1) ? curn : curp;
                uint32_t spos = starts[s];

                for (uint32_t j = 0; j < m; j++) {
                    const uint32_t b = spos >> log2c;
                    *(c[b]) = ((spos & maskc) << 8) | lp;
                    c[b]++;
                    spos += p;
                }
                const bool     in = (spos < hw);
                const uint32_t b  = in ? (spos >> log2c) : 0;
                *(c[b]) = ((spos & maskc) << 8) | lp;
                c[b] += in;
            }
        }
    }
    for (uint32_t b = 0; b < nbins; b++) {
        bposi->bins[b]->length = curp[b] - bposi->bins[b]->data;
        bnega->bins[b]->length = curn[b] - bnega->bins[b]->data;
    }
    if (cur != local) {
        free(cur);
    }
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
void print_init_poly_timing(siqs_sieve_t* const sieve MAYBE_UNUSED)
{
    PRINT_INIT_POLY_TIMING_MSG;
}
//-----------------------------------------------------------------------------
void print_fill_timing(siqs_sieve_t* const sieve MAYBE_UNUSED)
{
    PRINT_FILL_TIMING_MSG;
}
//-----------------------------------------------------------------------------
void print_scan_timing(siqs_sieve_t* const sieve MAYBE_UNUSED)
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
    sieve->clean     = false;
}
//------------------------------------------------------------------------------
uint32_t siqs_sieve_chunk_size(void)
{
    return 1U << log2_sieve_chunk_size();
}
