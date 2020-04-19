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
 * \file    rho.h
 * \author  Jerome Milan
 * \date    Mon Apr 28 2014
 * \version 2014-04-28
 *
 * \brief Pollard's rho factorization algorithm.
 *
 * This is the TIFA library's naive implementation of Pollard's rho
 * factorization algorithm, using Floyd's cycle finding algorithm.
 */

#if !defined(_TIFA_RHO_H_)
#define _TIFA_RHO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <gmp.h>

#include "array.h"
#include "exit_codes.h"

   /**
    * \enum rho_method_enum
    *
    * Enumeration listing the different cycle-finding method implemented.
    */
enum rho_method_enum {
        /**
         * Floyd's algorithm, as used in Pollard's original paper "A Monte Carlo
         * method for factorization", BIT Numerical Mathematics 1975, Volume 15,
         * pp 331-334.
         */
    FLOYD_CYCLE_FINDING = 0,
        /**
         * Brent's improvement, as described in "An improved Monte Carlo
         * factorization algoritm", BIT Numerical Mathematics 1980, Volume 20,
         * Issue 2, pp 176-184.
         */
    BRENT_CYCLE_FINDING
};

   /**
    * \typedef rho_method_t
    * \brief Equivalent to <tt>enum rho_method_enum</tt>.
    */
typedef enum rho_method_enum rho_method_t;

static const char* const rho_method_to_str[2] = {
    "Floyd's cycle detection",
    "Brent's cycle detection"
};

   /**
    * \struct struct_rho_params_t rho.h lib/algo/include/rho.h
    * \brief  Defines the variable parameters used in Pollard's rho algorithm.
    *
    * This structure is intended to define the set of the variable parameters
    * used in Pollard's rho algorithm.
    */
struct struct_rho_params_t {
        /**
         * The 'a' coefficient in the iteration function given by
         * (a * x^2 + b) mod n, where n is the number to factor.
         */
    uint32_t a;
        /**
         * The 'b' coefficient in the iteration function given by
         * (a * x^2 + b) mod n, where n is the number to factor.
         */
    uint32_t b;
        /**
         * Number of (x_{j} - x_{i}) to multiply together before computing
         * gcd(prod (x_{j} - x_{i})), n).
         */
    uint32_t acc_niters;
        /**
         * Maximum number of evaluations of the iteration function to compute
         * before giving up the factorization.
         */
    uint32_t max_niters;
        /**
         * The cycle finding method to use.
         */
    rho_method_t method;
};

   /**
    * \typedef rho_params_t
    * \brief Equivalent to <tt>struct struct_rho_params_t</tt>.
    */
typedef struct struct_rho_params_t rho_params_t;

   /**
    * \brief Fills a \c rho_params_t with default values.
    *
    * This function is intended to fill a previously allocated \c rho_params_t
    * with default values.
    *
    * \param params A pointer to the \c rho_params_t structure to fill.
    */
void set_rho_params_to_default(const mpz_t n, rho_params_t* const params);

   /**
    * \brief Integer factorization via Pollard's rho factorization algorithm.
    *
    * Attempts to factor the non perfect square integer \c n with
    * Pollard's rho algorithm, using the factoring mode given by <tt>mode</tt>.
    * Found factors are then stored in \c factors. Additionally, if the
    * factoring mode used is set to FIND_COMPLETE_FACTORIZATION,
    * factors' multiplicities are stored in the array <tt>multis</tt>.
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
    * \param[in] params Pollard's rho parameters.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t rho(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const rho_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif

