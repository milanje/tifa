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
 * \file    cofactor.h
 *
 * \brief Cofactorization: does an integer split into primes below a bound?
 *
 * The question asked of the cofactors left by the sieve of the number field
 * sieve (or of any sieve with large primes): given a bound L = 2^lpb, is
 * every prime factor of n below L, and which are they? Unlike
 * \c tifa_factor, the search stops as soon as the answer is known: a prime
 * factor of at least L, or more large primes than allowed.
 *
 * The work that depends only on the bounds is done once, in a strategy that
 * may be shared by several threads.
 */

#if !defined(_TIFA_COFACTOR_H_)
   /**
    * \def _TIFA_COFACTOR_H_
    * Standard include guard.
    */
#define _TIFA_COFACTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <gmp.h>

#include "array.h"

   /**
    * \enum cofactor_status_enum
    * \brief The outcome of \c tifa_cofactor.
    */
enum cofactor_status_enum {
        /**
         * Every prime factor of \c n is below 2^lpb, and there are not more
         * large primes than allowed: the factors are listed.
         */
    COFACTOR_SMOOTH = 0,
        /**
         * \c n has a prime factor of at least 2^lpb, or more large primes
         * than allowed (proven).
         */
    COFACTOR_NOT_SMOOTH,
        /**
         * Undecided: no factorization was completed within the limits set
         * by \c max_siqs_bits.
         */
    COFACTOR_GAVE_UP
};
   /**
    * \typedef cofactor_status_t
    * \brief Equivalent to <tt>enum cofactor_status_enum</tt>.
    */
typedef enum cofactor_status_enum cofactor_status_t;

   /**
    * \struct struct_cofactor_params_t cofactor.h lib/algo/include/cofactor.h
    * \brief  The bounds of a cofactorization.
    */
struct struct_cofactor_params_t {
        /**
         * Prime factors below 2^lpb are accepted (1 <= lpb <= 63).
         */
    uint32_t lpb;
        /**
         * The numbers to cofactor have no prime factor below this bound (the
         * factor base bound of a sieve). Below 1000 (0 if nothing is known),
         * the primes below 1000 are divided out first. A wrong bound gives
         * unspecified results.
         */
    uint64_t fb_bound;
        /**
         * Largest number of prime factors of at least \c fb_bound, counted
         * with multiplicity, or 0 for no limit.
         */
    uint32_t max_large_primes;
        /**
         * Above this many bits, a composite is never handed to SIQS: once
         * ECM has searched it for as long as SIQS would take, without
         * finding a factor, it is given up (\c COFACTOR_GAVE_UP). With
         * three large primes and <tt>max_siqs_bits = 2 lpb</tt>, this takes
         * 60% to 65% of the time of exact answers and misses 0.25% to 0.30%
         * of the smooth numbers (readme/news.md). 0 for no limit: every
         * answer is then exact. Ignored without 128-bit integers (the ECM
         * used here needs them).
         */
    uint32_t max_siqs_bits;
};
   /**
    * \typedef cofactor_params_t
    * \brief Equivalent to <tt>struct struct_cofactor_params_t</tt>.
    */
typedef struct struct_cofactor_params_t cofactor_params_t;

   /**
    * \typedef cofactor_strategy_t
    * \brief The bounds and the precomputations shared by all the calls.
    */
typedef struct struct_cofactor_strategy_t cofactor_strategy_t;

   /**
    * \brief Allocates the strategy for the bounds in \c params.
    *
    * \return The strategy, or NULL if \c lpb is out of range.
    */
cofactor_strategy_t* alloc_cofactor_strategy(const cofactor_params_t* params);

   /**
    * \brief Frees a strategy allocated by \c alloc_cofactor_strategy.
    */
void free_cofactor_strategy(cofactor_strategy_t* strategy);

   /**
    * \brief Finds whether the prime factors of \c n are all below 2^lpb.
    *
    * \param[out] factors The distinct prime factors of \c n, on
    *             \c COFACTOR_SMOOTH (unspecified otherwise).
    * \param[out] multis Their multiplicities.
    * \param[in]  n The positive integer to cofactor.
    * \param[in]  strategy The bounds (several threads may share it).
    * \return The outcome.
    */
cofactor_status_t tifa_cofactor(mpz_array_t* factors, uint32_array_t* multis,
                                const mpz_t n,
                                const cofactor_strategy_t* strategy);

#ifdef __cplusplus
}
#endif

#endif
