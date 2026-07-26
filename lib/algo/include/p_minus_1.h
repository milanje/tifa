// Copyright (C) 2026 Jerome Milan.
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
 * \file    p_minus_1.h
 * \author  Jerome Milan
 * \date    2026
 * \version 1.0
 *
 * \brief Pollard's p-1 factorization algorithm.
 *
 * This is the TIFA library's implementation of the (first stage of) Pollard's
 * p-1 factorization algorithm. It recovers a prime factor \c p of the number
 * to factor whenever <tt>p - 1</tt> is B1-smooth, by computing
 * <tt>gcd(a^E - 1, n)</tt> where \c E is the product of the prime powers up
 * to the bound \c B1.
 */

#if !defined(_TIFA_P_MINUS_1_H_)
#define _TIFA_P_MINUS_1_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <gmp.h>

#include "array.h"
#include "exit_codes.h"

   /**
    * \struct struct_p_minus_1_params_t p_minus_1.h lib/algo/include/p_minus_1.h
    * \brief  Defines the variable parameters used in Pollard's p-1 algorithm.
    *
    * This structure defines the set of the variable parameters used in the
    * first stage of Pollard's p-1 algorithm.
    */
struct struct_p_minus_1_params_t {
        /**
         * Stage 1 smoothness bound: all prime powers up to \c b1 are
         * accumulated in the exponent. A prime factor \c p is found when
         * <tt>p - 1</tt> is <tt>b1</tt>-smooth. Must be no greater than the
         * largest tabulated prime (see \c LARGEST_PRIME in first_primes.h).
         */
    uint32_t b1;
        /**
         * Stage 2 bound. Stage 2 recovers a factor \c p when <tt>p - 1</tt> is
         * a \c b1-smooth number times a single prime in <tt>(b1, b2]</tt>. If
         * <tt>b2 <= b1</tt>, stage 2 is skipped. Must be no greater than the
         * largest tabulated prime (see \c LARGEST_PRIME in first_primes.h).
         */
    uint32_t b2;
        /**
         * The base \c a of the exponentiation <tt>a^E mod n</tt>. Usually 2.
         */
    uint32_t base;
};

   /**
    * \typedef p_minus_1_params_t
    * \brief Equivalent to <tt>struct struct_p_minus_1_params_t</tt>.
    */
typedef struct struct_p_minus_1_params_t p_minus_1_params_t;

   /**
    * \brief Fills a \c p_minus_1_params_t with default values.
    *
    * This function is intended to fill a previously allocated
    * \c p_minus_1_params_t with default values.
    *
    * \param[in]  n The number to factor. Its bit-size selects the row from a
    *               tuning file when one is available (see tuning.h); otherwise
    *               size-independent built-in defaults are used.
    * \param[out] params A pointer to the \c p_minus_1_params_t structure to fill.
    */
void set_p_minus_1_params_to_default(const mpz_t n, p_minus_1_params_t* const params);

   /**
    * \brief Integer factorization via Pollard's p-1 factorization algorithm.
    *
    * Attempts to factor the non perfect square integer \c n with the first
    * stage of Pollard's p-1 algorithm, using the factoring mode given by
    * <tt>mode</tt>. Found factors are then stored in \c factors. Additionally,
    * if the factoring mode used is set to FIND_COMPLETE_FACTORIZATION,
    * factors' multiplicities are stored in the array <tt>multis</tt>.
    *
    * \note If the factoring mode used is different from
    * FIND_COMPLETE_FACTORIZATION, \c multis is allowed to be a NULL pointer.
    * Otherwise, using a NULL pointer will lead to a fatal error.
    *
    * \warning If the \c factors and \c multis arrays have not enough room
    * to store the found factors (and the multiplicities, if any), they will
    * be automatically resized to accommodate the data.
    *
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis Pointer to the multiplicities of the found factors
    *                    (only computed if \c mode is set to
    *                     FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The non perfect square integer to factor.
    * \param[in] params Pollard's p-1 parameters.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t p_minus_1(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const p_minus_1_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
