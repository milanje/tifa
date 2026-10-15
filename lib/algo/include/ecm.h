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
 * \file    ecm.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief The elliptic curve method of integer factorization (ECM).
 *
 * This is the TIFA library's implementation of the 'ECM' factorization
 * algorithm. The second phase of the algorithm follows the standard
 * continuation and is implemented in a way reminiscent of the description
 * given in the article "Implementing the Elliptic Curve Method of Factoring
 * in Reconfigurable Hardware" by Kris Gaj et al.
 *
 * \see "Implementing the Elliptic Curve Method of Factoring in Reconfigurable
 * Hardware", K. Gaj et al., <i>Cryptographic Hardware and Embedded Systems -
 * CHES 2006</i>.
 *
 * \warning This is merely a toy-implementation of ECM without any smart
 * optimizations. More work is certainly needed to make it competitive
 * for small numbers. Large numbers are, of course, out of the scope of this
 * library.
 */

#if !defined(_TIFA_ECM_H_)
   /**
    * \def _TIFA_ECM_H_
    * Standard include guard.
    */
#define _TIFA_ECM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \struct struct_ecm_params_t ecm.h lib/algo/include/ecm.h
    * \brief  Defines the variable parameters used in ECM.
    *
    * This structure defines the set of the variable parameters used in the
    * elliptic curve method (ECM).
    */
struct struct_ecm_params_t {
        /**
         * Bound used in the first phase of ECM. Values above the largest
         * tabulated prime (see \c LARGEST_PRIME in first_primes.h) are
         * lowered to it.
         */
    uint32_t b1;
        /**
         * Bound used in the second phase of ECM. If set to 0, no second phase
         * is performed.
         *
         * \warning Due to a current limitation in the code, it is \c required
         * than \c b2 == 0 (no second phase) or \c b2 > 105. Failure to assess
         * such a condition will lead to unpredictable behaviour.
         */
    uint32_t b2;
        /**
         * Number of curves per round: the only round in the SINGLE_RUN
         * factoring mode, the first one otherwise (see
         * \c set_ecm_params_to_default).
         */
    uint32_t ncurves;
        /**
         * Seed of the curve choice: two runs on the same number with the
         * same nonzero seed use the same curves. 0 derives the seed from the
         * number to factor.
         */
    uint32_t seed;
};
   /**
    * \typedef ecm_params_t
    * \brief Equivalent to <tt>struct struct_ecm_params_t</tt>.
    */
typedef struct struct_ecm_params_t ecm_params_t;

   /**
    * \brief Fills an \c ecm_params_t with default values.
    *
    * Sets the bounds suited to the smallest factors (about 20 bits). In
    * modes other than SINGLE_RUN, \c ecm then moves up to the bounds for
    * larger factors, 5 bits at a time, as long as no factor is found, up to
    * factors of half the size of \c n, and tries that last level a few more
    * times before giving up. A tuning file row for the size of \c n (see
    * tuning.h) overrides the starting values.
    *
    * \param[in]  n The \c mpz_t integer to factor.
    * \param[out] params A pointer to the \c ecm_params_t structure to fill.
    */
void set_ecm_params_to_default(const mpz_t n, ecm_params_t* const params);

   /**
    * \brief Fills an \c ecm_params_t for factors of a given size.
    *
    * Sets \c b1 and \c b2 to the bounds minimizing the expected time to
    * find a prime factor of \c fbits bits, and \c ncurves to the mean number
    * of curves that takes.
    *
    * \param[in]  fbits The size of the factor sought, in bits.
    * \param[out] params A pointer to the \c ecm_params_t structure to fill.
    */
void set_ecm_params_for_factor_size(uint32_t fbits,
                                    ecm_params_t* const params);

   /**
    * \brief Integer factorization with the elliptic curve method (ECM).
    *
    * Attempts to factor the non perfect square integer \c n with the
    * ECM, using the set of parameters given by <tt>params</tt> and the
    * factoring mode given by <tt>mode</tt>. Found factors are then
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
    * \param[in] params Pointer to the values of the parameters used in the ECM.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t ecm(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const ecm_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
