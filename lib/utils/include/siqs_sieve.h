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
 * \file    siqs_sieve.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Structure and functions related to the sieve used in
 *        the SIQS algorithm.
 */
 
#if !defined(_TIFA_SIQS_SIEVE_H_)
#define _TIFA_SIQS_SIEVE_H_

#include <stdint.h>
#include <stdbool.h>

#include <gmp.h>

#include "exit_codes.h"
#include "array.h"
#include "approx.h"
#include "buckets.h"
#include "siqs_poly.h"
#include "stopwatch.h"

#ifdef __cplusplus
extern "C" {
#endif

//
// Logarithm (in base 2) of the largest size of a sieve chunk. The optimal
// value is architecture dependant, so this value should be tweaked according
// to the target processor.
//
// WARNING: Should be an integer in the [10..20] range.
//
#define LOG2_SIEVE_CHUNK_MAX_SIZE 15

//
// Largest size of a sieve chunk.
//
// SIEVE_CHUNK_MAX_SIZE = 2^LOG2_SIEVE_CHUNK_MAX_SIZE
//
#if   (LOG2_SIEVE_CHUNK_MAX_SIZE == 10)
    #define SIEVE_CHUNK_MAX_SIZE 1024
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 11)
    #define SIEVE_CHUNK_MAX_SIZE 2048
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 12)
    #define SIEVE_CHUNK_MAX_SIZE 4096
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 13)
    #define SIEVE_CHUNK_MAX_SIZE 8192
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 14)
    #define SIEVE_CHUNK_MAX_SIZE 16384
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 15)
    #define SIEVE_CHUNK_MAX_SIZE 32768
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 16)
    #define SIEVE_CHUNK_MAX_SIZE 65536
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 17)
    #define SIEVE_CHUNK_MAX_SIZE 131072
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 18)
    #define SIEVE_CHUNK_MAX_SIZE 262144
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 19)
    #define SIEVE_CHUNK_MAX_SIZE 524288
#elif (LOG2_SIEVE_CHUNK_MAX_SIZE == 20)
    #define SIEVE_CHUNK_MAX_SIZE 1048576
#else
  #error "LOG2_SIEVE_CHUNK_MAX_SIZE should be in [10..20]"
#endif

//
// Do not sieve for the NPRIMES_TO_SKIP smallest primes in the factor base
//
#define NPRIMES_TO_SKIP        20
//
// Smallest prime index for which a bucket sieve is performed
//
#define BUCKETS_SMALLEST_PRIME 2400

#if (NPRIMES_TO_SKIP == 0)
    //
    // Sieving for the prime 2 is done separately...
    //
    #define SIEVE_WITH_2      1
    #undef  NPRIMES_TO_SKIP
    #define NPRIMES_TO_SKIP   1
#else
    #define SIEVE_WITH_2      0
#endif

//
// If ROUND_HALF_WIDTH is not zero, round the half-width to the nearest
// multiple of SIEVE_CHUNK_MAX_SIZE (unless if the half-width is less than 
// SIEVE_CHUNK_MAX_SIZE).
//
#define ROUND_HALF_WIDTH  0

    /**
     * \struct struct_siqs_sieve_t siqs_sieve.h lib/utils/include/siqs_sieve.h
     * \brief  Defines the sieve used by SIQS.
     *
     * This structure defines the sieve used by SIQS together with all
     * its associated data.
     */
struct struct_siqs_sieve_t {
    /**
     * Number of blocks to sieve on each side of zero (\tt{nchunks} for
     * positive \b{x} and \tt{nchunks} for negative \b{x})
     */
    uint32_t nchunks;
    /**
     * Size of a sieve chunk. The total sieving interval is thus given by
     * \tt{2 * chunk_size * nchunks}
     */
    uint32_t chunk_size;
    /**
     * The number of the next sieve chunk to fill (in [\tt{nchunks}, 0[ U
     * ]0, \tt{nchunks}])
     */   
    int32_t  next_chunkno_to_fill;
    /**
     * The position to start scanning in the next sieve chunk to scan
     */   
    uint32_t scan_begin;
    /**
     * The sieve threshold. An \b{xi} will be tested for smoothness if
     * \tt{sieve[}\b{xi}\tt{]} < \tt{threshold}
     */
    uint32_t threshold;
    /**
     * Approximated logarithm (in base 2) of the primes in the factor base
     */
    byte_array_t* log_primes;
    /**
     * The actual sieve array, of size \tt{chunk_size}
     */
    byte_array_t* sieve;
    /**
     * The SIQS polynomial
     */
    siqs_poly_t*  poly;
     /**
      * The first solution to the equation \b{P}(x) = 0 mod \b{pi} for
      * each prime \b{pi} in the factor base
      */
    int32_array_t *sol1;
     /**
      * The first solution to the equation \b{P}(x) = 0 mod \b{pi} for
      * each prime \b{pi} in the factor base
      */
    int32_array_t *sol2;
     /**
      * The index of the last position to sieve in the last sieve chunk (on
      * either side)
      */
#if !ROUND_HALF_WIDTH
    uint32_t endlast;
#endif

     /**
      * Should we use a bucket sieving for the largest primes in the base?
      */
    bool       use_buckets;
     /**
      * Index (in the factor base) of the largest prime for which a standard
      * sieving procedure is use
      */
    uint32_t   nprimes_no_buckets;
     /**
      * Index (in the factor base) of the smallest prime for which the bucket 
      * sieving procedure is use
      */
    uint32_t   buckets_first_prime;
     /**
      * Buckets for the positive \b{x}'s
      */
    buckets_t* buckets_positive;
     /**
      * Buckets for the negative \b{x}'s
      */
    buckets_t* buckets_negative;
     /**
      * Additional stopwatch (to get timing for the polynomial initializations)
      */
    stopwatch_t init_poly_timer;
     /**
      * Additional stopwatch (to get timing for the sieve filling)
      */
    stopwatch_t fill_timer;
     /**
      * Additional stopwatch (to get timing for the sieve scanning)
      */
    stopwatch_t scan_timer;
};

   /**
    * \typedef siqs_sieve_t
    * \brief Equivalent to \tt{struct_siqs_sieve_t}.
    */
typedef struct struct_siqs_sieve_t siqs_sieve_t;

   /**
    * \brief Allocates and returns a new <tt>siqs_sieve_t</tt>.
    *
    * \param n           the number to factor (or a small multiple)
    * \param factor_base the factor base
    * \param log_primes  logarithms (in base 2) of the primes in the base
    * \param sqrtm_pi    modular square roots of \tt{n} for each prime in
    *                    the base
    * \param half_width  the (approximate) half_width of the sieving interval
    *                    (the real half_width will be adjusted to be a multiple
    *                    of \tt{chunk_size} if \tt{ROUND_HALF_WIDTH} is
    *                    defined as non zero)
    *
    * \return A pointer to the newly allocated \c siqs_sieve_t.
    */
siqs_sieve_t* alloc_siqs_sieve(
    mpz_t n,
    uint32_array_t* const factor_base,
    byte_array_t*   const log_primes,
    uint32_array_t* const sqrtm_pi,
    uint32_t half_width
);

   /**
    * \brief Frees a previously allocated <tt>siqs_sieve_t</tt>.
    *
    * Frees all memory used by the pointed <tt>siqs_sieve_t</tt> and then
    * frees the \tt{sieve} pointer.
    *
    * \warning Do not call \tt{free(sieve)} in client code after a call to
    * \tt{free_siqs_sieve(sieve)}: it would result in an error.
    *
    * \param sieve the \tt{siqs_sieve_t} to free.
    */
void free_siqs_sieve(siqs_sieve_t* sieve);

   /**
    * \brief Fills the next chunk of an <tt>siqs_sieve_t</tt>.
    *
    * Fills the next chunk of \tt{sieve}, transparently updating (if needed)
    * the polynomial used.
    *
    * \param sieve the \tt{siqs_sieve_t} to fill.
    * 
    * \return An exit code.
    */
ecode_t fill_sieve(siqs_sieve_t* const sieve);

   /**
    * \brief Scans a chunk of an <tt>siqs_sieve_t</tt>.
    * 
    * Scans the last filled chunk of \tt{sieve}.
    *
    * \param sieve the \tt{siqs_sieve_t} to scan.
    * 
    * \return An exit code.
    */
ecode_t scan_sieve(
    siqs_sieve_t* const sieve,
    int32_array_t* const survivors,
    uint32_t nsurvivors
);

   /**
    * \brief Sets the <tt>siqs_sieve_t</tt>'s threshold.
    * 
    * Sets \tt{sieve}'s threshold (all positions \b{xi} with
    * \tt{sieve[}\b{xi}\tt{]} < threshold will be tested for smoothness).
    *
    * \param sieve     the \tt{siqs_sieve_t} to update.
    * \param threshold the new threshold's value.
    */
void set_siqs_sieve_threshold(siqs_sieve_t* const sieve, uint32_t threshold);

   /**
    * \brief Prints an <tt>siqs_sieve_t</tt>'s poly init timing.
    * 
    * Prints the time taken by \tt{sieve} to initialize its polynomials.
    *
    * \param sieve the \tt{siqs_sieve_t} to read.
    */
void print_init_poly_timing(siqs_sieve_t* const sieve);

   /**
    * \brief Prints an <tt>siqs_sieve_t</tt>'s fill timing.
    * 
    * Prints the time taken by \tt{sieve} to fill its sieve chunks.
    *
    * \param sieve the \tt{siqs_sieve_t} to read.
    */
void print_fill_timing(siqs_sieve_t* const sieve);

   /**
    * \brief Prints an <tt>siqs_sieve_t</tt>'s scan timing.
    * 
    * Prints the time taken by \tt{sieve} to scan its sieve chunks.
    *
    * \param sieve the \tt{siqs_sieve_t} to read.
    */
void print_scan_timing(siqs_sieve_t* const sieve);

#ifdef __cplusplus
}
#endif

#endif


