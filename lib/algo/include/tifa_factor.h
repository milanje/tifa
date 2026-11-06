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
 * \file    tifa_factor.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief TIFA's generic factorization function.
 *
 * This is the TIFA library's generic factorization function: it picks the most
 * suitable factoring algorithm depending on the size of the number to factor.
 */

#if !defined(_TIFA_TIFA_FACTOR_H_)
   /**
    * \def _TIFA_TIFA_FACTOR_H_
    * Standard include guard.
    */
#define _TIFA_TIFA_FACTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \brief Generic Integer factorization.
    *
    * Attempts to factor the integer \c n with the most suitable algorithms
    * (chosen according to the size of \c n) and with the factoring mode given
    * by <tt>mode</tt>. Found factors are then stored in \c factors.
    * Additionally, if the factoring mode used is set to
    * \c FIND_COMPLETE_FACTORIZATION, factors' multiplicities are stored in
    * the array <tt>multis</tt>. Numbers below 2^64 are factored with
    * machine-word arithmetic where the compiler supports it (see
    * factor_u64.h). Otherwise, prime factors below 1000 (below 8192 for
    * numbers of at most 26 bits) are found by trial division before any
    * other algorithm runs. The composites left go to SQUFOF up to 58 bits;
    * above, a short ECM search (see ecm_schedule.h, up to 512 bits) looks
    * for a factor of up to 60 bits (85 bits above 256), or, where it can not
    * run, single p-1 and p+1 runs from 90 bits, before SIQS. Above 512 bits,
    * those p-1 and p+1 runs, then \c ecm up to factors of 85 bits, without a
    * budget (SIQS would take years). \c tifa_factor_with_hint changes this
    * when more is known of <tt>n</tt>.
    *
    * \note If the factoring mode used is different from
    * \c FIND_COMPLETE_FACTORIZATION, \c multis is allowed to be a NULL
    * pointer. Otherwise, using a NULL pointer will lead to a fatal error.
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
    *                     \c FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The integer to factor.
    * \param[in] mode The factoring mode to use.
    * \return An exit code: \c NO_FACTOR_FOUND if \c n is prime.
    */
TIFA_API ecode_t tifa_factor(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const factoring_mode_t mode
);

   /**
    * \enum tifa_factor_hint_enum
    * \brief What is known of the number given to \c tifa_factor_with_hint.
    *
    * A hint never changes the factors found, only the time taken to find
    * them. It applies to the prime factors of n above 1000, and is ignored up
    * to 64 bits (58 without 128-bit integers), where machine words (or
    * SQUFOF) factor n.
    */
enum tifa_factor_hint_enum {
        /**
         * Nothing is known: ECM searches for a factor for 20% of the time
         * SIQS would take, then SIQS runs. This is what \c tifa_factor does.
         */
    TIFA_FACTOR_RANDOM = 0,
        /**
         * n has two prime factors of about the same size (an RSA modulus):
         * SIQS runs right away.
         */
    TIFA_FACTOR_BALANCED,
        /**
         * n has a prime factor of less than about a third of its size: ECM
         * searches for as long as SIQS would take (with 128-bit integers,
         * up to 512 bits; otherwise, and above, as \c TIFA_FACTOR_RANDOM).
         */
    TIFA_FACTOR_UNBALANCED
};
   /**
    * \typedef tifa_factor_hint_t
    * \brief Equivalent to <tt>enum tifa_factor_hint_enum</tt>.
    */
typedef enum tifa_factor_hint_enum tifa_factor_hint_t;

   /**
    * \brief Integer factorization, given what is known of the number.
    *
    * Same as \c tifa_factor, with the methods tried before SIQS chosen
    * according to <tt>hint</tt>. On numbers of the shape it names, a hint
    * saves time; on others, it can cost much more (readme/news.md). The
    * hint carries over to what is left of n after trial division, and to
    * the root of a perfect power, but the parts of n split by the other
    * methods are factored as by \c tifa_factor.
    *
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis Pointer to the multiplicities of the found factors
    *                    (only computed if \c mode is set to
    *                     \c FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The integer to factor.
    * \param[in] mode The factoring mode to use.
    * \param[in] hint What is known of <tt>n</tt>.
    * \return An exit code: \c NO_FACTOR_FOUND if \c n is prime,
    *         \c FAILURE if \c hint is not a \c tifa_factor_hint_t.
    */
TIFA_API ecode_t tifa_factor_with_hint(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const factoring_mode_t mode,
    const tifa_factor_hint_t hint
);

#ifdef __cplusplus
}
#endif

#endif

