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
 * \file    siqs.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief The Self-Initializing Quadratic Sieve factorization algorithm.
 */

#if !defined(_TIFA_SIQS_H_)
   /**
    * \def _TIFA_SIQS_H_
    * Standard include guard.
    */
#define _TIFA_SIQS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <inttypes.h>
#include <gmp.h>

#include "array.h"
#include "methods.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \def SIQS_DFLT_NRELATIONS
    * Default number of relations to collect beyond the size of the factor
    * base (see \c nrelations).
    */
#define SIQS_DFLT_NRELATIONS       8
   /**
    * \def SIQS_DFLT_LINALG_METHOD
    * Default linear system resolution method to use.
    */
#define SIQS_DFLT_LINALG_METHOD    SMART_GAUSS_ELIM
   /**
    * \def SIQS_DFLT_LARGE_PRIMES_MIN_BITS
    * Use the large prime variation by default from this many bits.
    */
#define SIQS_DFLT_LARGE_PRIMES_MIN_BITS 103
   /**
    * \def SIQS_FILTERED_GAUSS_MIN_BASE
    * By default, the matrix is filtered before the elimination
    * (\c FILTERED_GAUSS_ELIM) from this many primes in the factor base (about
    * 188 bits): 0.96 of the time at 210 and 220 bits on an Apple M5, but up
    * to 1.05 below 160 bits, where the elimination costs less than the
    * filtering saves.
    */
#define SIQS_FILTERED_GAUSS_MIN_BASE 3000
   /**
    * \def SIQS_BLOCK_LANCZOS_MIN_BASE
    * By default, block Lanczos (\c BLOCK_LANCZOS) on a sparse matrix, without
    * any dense one, from this many primes in the factor base (about 230
    * bits). The linear algebra took the same time at 220 bits on an Apple M5
    * Pro, 0.61 of it at 240, 0.26 at 264, 0.15 at 280 (1.1 s rather than
    * 7.5), and the run's peak memory 0.09 to 0.20 of it from 240 to 280 bits.
    */
#define SIQS_BLOCK_LANCZOS_MIN_BASE 12000

   /**
    * \struct struct_siqs_params_t siqs.h lib/algo/include/siqs.h
    * \brief  Defines the variable parameters used in the SIQS algorithm.
    *
    * This structure defines the set of the variable parameters used in the
    * SIQS algorithm.
    */
struct struct_siqs_params_t {
       /**
        * Sieve's half-width, i.e. the SIQS will sieve in the interval
        * <tt>[-sieve_half_width, sieve_half_width]</tt> (from 1 to 2^26).
        */
    uint32_t sieve_half_width;
       /**
        * Number of prime numbers composing the factor base on which to factor
        * the residues, at most 65536, and enough to choose the coefficient
        * \c a.
        */
    uint32_t nprimes_in_base;
      /**
       * Sieve threshold (at most 127).
       */
   uint32_t threshold;
       /**
        * Number of the first primes to use in the trial division
        * of the residues (from 1 to \c nprimes_in_base).
        */
    uint32_t nprimes_tdiv;
       /**
        * Number of relations to collect beyond the size of the factor base
        * (plus one, for the sign), before looking for linear dependences (at
        * most 1024). When none of them splits the number, as many more are
        * collected, up to three times as many in all.
        */
    uint32_t nrelations;
       /**
        * Linear system resolution method to use: by default the one of
        * \c siqs_dflt_linalg_method() for \c nprimes_in_base.
        */
    linalg_method_t linalg_method;
       /**
        * True if we use the large prime variation.
        * False otherwise.
        */
    bool use_large_primes;
       /**
        * Number of primes added to the polynomials' coefficient \c a, beyond
        * the usual number for its size (at most 2). Each added prime doubles
        * the number of polynomials sharing an \c a, whose setup is costly,
        * but makes these primes smaller. The defaults add none.
        */
    uint32_t extra_primes_in_a;
};

   /**
    * \typedef siqs_params_t
    * \brief Equivalent to <tt>struct struct_siqs_params_t</tt>.
    */
typedef struct struct_siqs_params_t siqs_params_t;

   /**
    * \brief The default linear algebra of SIQS for a factor base of
    * \c nprimes_in_base primes.
    *
    * The Gaussian elimination, filtered from \c SIQS_FILTERED_GAUSS_MIN_BASE
    * primes, and block Lanczos from \c SIQS_BLOCK_LANCZOS_MIN_BASE.
    * \c set_siqs_params_to_default() takes it for the factor base it sets,
    * tuning files included; take it again after changing \c nprimes_in_base.
    *
    * \param[in] nprimes_in_base The number of primes in the factor base.
    * \return The linear algebra method.
    */
static inline linalg_method_t siqs_dflt_linalg_method(
    const uint32_t nprimes_in_base) {
    if (nprimes_in_base >= SIQS_BLOCK_LANCZOS_MIN_BASE) {
        return BLOCK_LANCZOS;
    }
    if (nprimes_in_base >= SIQS_FILTERED_GAUSS_MIN_BASE) {
        return FILTERED_GAUSS_ELIM;
    }
    return SIQS_DFLT_LINALG_METHOD;
}

   /**
    * \brief Fills a \c siqs_params_t with default values.
    *
    * Fills a \c siqs_params_t with default values.
    *
    * \param[in]  n      The number to factor.
    * \param[out] params A pointer to the \c siqs_params_t structure to fill.
    */
TIFA_API void set_siqs_params_to_default(const mpz_t n,
                                         siqs_params_t* const params);

   /**
    * \brief Integer factorization via the self-initializing quadratic sieve
    * (SIQS) algorithm.
    *
    * Attempts to factor the non perfect square integer \c n with the
    * SIQS algorithm, using the set of parameters given by <tt>params</tt>
    * and the factoring mode given by <tt>mode</tt>. Found factors are then
    * stored in \c factors. Additionally, if the factoring mode used is set
    * to FIND_COMPLETE_FACTORIZATION, factors' multiplicities are stored in
    * the array <tt>multis</tt>.
    *
    * \note If the factoring mode used is different from
    * FIND_COMPLETE_FACTORIZATION, \c multis is allowed to be a NULL pointer.
    * Otherwise, using a NULL pointer will lead to a fatal error.
    *
    * \warning If the \c factors and \c multis arrays have not enough room
    * to store the found factors (and the multiplicities, if any), they will
    * be automatically resized to accommodate the data. This has to be kept
    * in mind when trying to do ingenious stuff with memory management (hint:
    * don't try to be clever here).
    *
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis Pointer to the multiplicities of the found factors
    *                    (only computed if \c mode is set to
    *                     FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The non perfect square integer to factor.
    * \param[in] params Pointer to the values of the parameters used in the
    *                   SIQS algorithm.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
TIFA_API ecode_t siqs(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const siqs_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
