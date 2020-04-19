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
 * \file    approx.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Approximate a value by multiplying some numbers from a pool.
 *
 * This provides a structure \tt{approximer_t} and associated functions that
 * can be used to approximate a target value by multiplying a given number of
 * factors from a given base. Each factor is allowed to appear only once in
 * the decomposition of the approximation on the given base.
 *
 * This is used in TIFA's SIQS implementation where we need to find a
 * polynomial coefficient of a given order from the product of some
 * prime numbers.
 *
 * The strategy used to reach a good approximation is adapted from the
 * Carrier-Wagstaff method.
 *
 * \li Let \b{t} be the target number to be obtained by multiplying \b{n}
 * distinct numbers from a given base \b{B}=[p_1, p_2, p_3, ...]
 * (with p_i < p_{i+1}).
 *
 * \li Set \b{m} so that \b{B[m]} is roughly equal to the
 * \b{n}<sup>th</sup>-root of \b{t}. Factors will be chosen from the set
 * \b{B[m-d]}, ..., \b{B[m+d]} with \b{d} suitably chosen.
 *
 * \li The approximation \b{a} of the target \b{t} is obtained by choosing
 * a random combination of (\b{n-i}) factors completed by \b{i} other factors
 * so as to obtain the best approximation as possible.
 *
 * \li Since all numbers must be distinct the (\b{n-i}) randomly chosen
 * factors are picked from the set of \b{B[j]} with \b{j} odd and the remaining
 * \b{i} factors from the set of \b{B[k]} with \b{k} even.
 *
 * \li The best remaining \b{i} factors are obtained by precomputing and
 * sorting all combinations of \b{i} factors, then picking the most suitable
 * one.
 */

/**
 * \def _TIFA_APPROX_H_
 * Standard include guard.
 */
#if !defined(_TIFA_APPROX_H_)
#define _TIFA_APPROX_H_

#include <stdint.h>
#include <stdbool.h>

#include "exit_codes.h"
#include "array.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * \def MAX_NPRIMES_IN_TUPLE
     * Maximum number of factors in the sorted factor combinations.
     */
#define MAX_NPRIMES_IN_TUPLE 3

    /**
     * \struct struct_uint32_tuple_t approx.h lib/utils/include/approx.h
     * \brief  Defines a tuple of integers together with a sorting key.
     *
     * This structure defines a tuple of up to \tt{MAX_NPRIMES_IN_TUPLE}
     * integers together with a sorting key given as a float.
     */
struct struct_uint32_tuple_t {
      /**
       * The value of the tuple.
       */
    uint32_t tuple[MAX_NPRIMES_IN_TUPLE];
       /**
        * Used as a key to sort tuples.
        */
    float    tlog;
};

    /**
     * \typedef uint32_tuple_t
     * \brief Equivalent to \tt{struct_uint32_tuple_t}.
     */
typedef struct struct_uint32_tuple_t uint32_tuple_t;

    /**
     * \struct struct_approximer_t approx.h lib/utils/include/approx.h
     * \brief  Structure used to find number approximation.
     *
     * The structure \tt{approximer_t} (and its associated functions)
     * is used to approximate a target value by multiplying a given number of
     * factors from a given base. Each factor is allowed to appear only once in
     * the decomposition of the approximation on the given base.
     */
struct struct_approximer_t {
       /**
        * The target number to approximate.
        */
    mpz_t target;
       /**
        * Logarithm in base 10 of the target number.
        */
    float targetlog;
       /**
        * The logarithm in base 10 of the approximation should be within
        * \tt{dlog_tolerance} of \tt{targetlog}.
        */
    float dlog_tolerance;
       /**
        * Pool of potential factors.
        */
    uint32_array_t* facpool;
       /**
        * Numbers of factors from \tt{facpool} in the approximation.
        */
    uint32_t  nfactors;
       /**
        * \tt{facpool[imin]} is the smallest factor allowed.
        */
    int32_t imin;
       /**
        * \tt{facpool[imax]} is the largest factor allowed.
        */
    int32_t imax;
       /**
        * Number of factors to choose with even indexes.
        */
    uint32_t keven;
       /**
        * Number of available factors with even indexes.
        */
    uint32_t neven;
       /**
        * Number of factors to choose with odd indexes.
        */
    uint32_t kodd;
       /**
        * Number of available factors with odd indexes.
        */
    uint32_t nodd;
       /**
        * Number of distinct combination of \tt{kodd} factors from
        * the pool of odd indexed factors (in other
        * words C(\tt{nodd}, \tt{kodd}) ).
        */
    uint32_t  nsubsets_odd;
       /**
        * A subset of \tt{kodd} odd-indexes.
        */
    uint32_t* subset_odd;
       /**
        * The rank of the subset of \tt{kodd} odd-indexes in lexicographic
        * order.
        */
    uint32_t rank;
       /**
        * List of tuples obtained by precomputing all combinations of
        * \tt{keven} factors.
        */
    uint32_tuple_t* tuples;
       /**
        * Number of all possible combinations of \tt{keven} factors.
        */
    uint32_t        ntuples;
};

    /**
     * \typedef struct_approximer_t
     * \brief Equivalent to \tt{approximer_t}.
     */
typedef struct struct_approximer_t approximer_t;

   /**
    * \brief Allocates and returns a new <tt>approximer_t</tt>.
    *
    * \param target the target number to approximate.
    * \param facpool the pool of available factors.
    * \param nfactors the number of factors from \tt{facpool} to use.
    *
    * \return A pointer to the newly allocated \c approximer_t.
    */
approximer_t* alloc_approximer(
    mpz_t                 target,
    uint32_array_t* const facpool,
    uint32_t              nfactors
);

   /**
    * \brief Frees a previously allocated <tt>approximer_t</tt>.
    *
    * Frees all memory used by the pointed <tt>approximer_t</tt> and then
    * frees the \tt{aximer} pointer.
    *
    * \warning Do not call \tt{free(aximer)} in client code after a call to
    * \tt{free_approximer(aximer)}: it would result in an error.
    *
    * \param aximer the \tt{approximer_t} to free.
    */
void free_approximer(approximer_t* aximer);

   /**
    * \brief Generates a "random" approximation.
    *
    * \param[in]  aximer   the \tt{approximer_t} to use.
    * \param[out] approxed the approximation obtained.
    * \param[out] indexes  the (sorted) indexes of the factors making up the
    *                      approximation.
    */
void random_approximation(
  approximer_t* const aximer,
  mpz_t               approxed,
  uint32_t*           indexes
);

#ifdef __cplusplus
}
#endif

#endif
