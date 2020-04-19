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
 * \file    siqs_poly.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Structure and functions related to the polynomials used in
 *        the SIQS algorithm.
 */
 
#if !defined(_TIFA_SIQS_POLY_H_)
#define _TIFA_SIQS_POLY_H_

#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

#include "exit_codes.h"
#include "array.h"
#include "approx.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \struct struct_siqs_poly_t siqs_poly.h lib/utils/include/siqs_poly.h
    * \brief  Defines polynomials used by SIQS.
    *
    * This structure defines the polynomials used by SIQS together with all
    * its associated data.
    *
    * A polynomial \b{P} is given by \b{a}X^2 + \b{b}X + \b{c}.
    */
struct struct_siqs_poly_t {
    /**
     * The \b{a} coefficient
     */
  mpz_t a;
    /**
     * The \b{b} coefficient
     */
  mpz_t b;
    /**
     * The \b{c} coefficient
     */
  mpz_t c;
    /**
     * The number to factor (or a small multiple if a multiplier is used)
     */
  mpz_t n;
    /**
     * Logarithm of \b{a} in base 2
     */
  uint32_t loga;
    /**
     * Logarithm of \b{b} in base 2
     */
  uint32_t logb;
    /**
     * Logarithm of \b{c} in base 2
     */
  uint32_t logc;
    /**
     * An \b{approximer_t} used to determine suitable values of the \b{a}
     * coefficient
     */
  approximer_t* approximer;
    /**
     * The factor base
     */
  uint32_array_t *factor_base;
    /**
     * The modular square roots of \tt{n} modulo each primes in the factor base
     */
  uint32_array_t *sqrtm_pi;
    /**
     * The first solution to the equation \b{P}(x) = 0 mod \b{pi} for
     * each prime \b{pi} in the factor base
     */
  int32_array_t *sol1;
    /**
     * The second solution to the equation \b{P}(x) = 0 mod \b{pi} for
     * each prime \b{pi} in the factor base
     */
  int32_array_t *sol2;
    /**
     * For all \b{l} in [0, \tt{npolys}-1] solutions to
     *
     * \li \tt{Bl}^2 = \tt{n} mod \b{q_l}
     * \li \tt{Bl} = 0 mod \b{q_j} for all \b{j} != \b{l}
     *
     * where \b{q_i} are the primes from the factor base such that
     * \b{a} = \b{q_0} x \b{q_1} x ... x \b{q_}\tt{npolys}
     */
  mpz_array_t* Bl;
    /**
     * \tt{Bainv2[i][j]} = 2 x \tt{Bl[i]} x inv(\tt{a}) mod \b{pj} for
     * \b{i} in [0, \tt{npolys}-1] and \b{pj} in the factor base 
     */
  uint32_t**   Bainv2;
    /**
     * Number of different polynomials having the same \tt{a} coefficient
     */
  uint32_t npolys;
    /**
     * Current polynomial number (from 1 to \tt{npolys})
     */
  uint32_t polyno;
    /**
     * Number of primes in the prime decomposition of \tt{a}
     */  
  uint32_t  nprimes_in_a;
    /**
     * Indexes (in the factor base) of the (prime) factors of \tt{a}
     */
  uint32_t* idx_of_a;
    /**
     * Number of "full" polynomial initializations performed
     */
  uint32_t nfullpolyinit;
};

   /**
    * \typedef siqs_poly_t
    * \brief Equivalent to \tt{struct_siqs_poly_t}.
    */
typedef struct struct_siqs_poly_t siqs_poly_t;

   /**
    * \brief Allocates and returns a new <tt>siqs_poly_t</tt>.
    *
    * \param target_a    the target leading coefficient to approximate.
    * \param n           the number to factor (or a small multiple).
    * \param factor_base the factor base.
    * \param sqrtm_pi    the modular square roots of n.
    *
    * \return A pointer to the newly allocated \c siqs_poly_t.
    */
siqs_poly_t* alloc_siqs_poly(
  mpz_t target_a,
  mpz_t n,
  uint32_array_t* const factor_base,
  uint32_array_t* const sqrtm_pi
);

   /**
    * \brief Frees a previously allocated <tt>siqs_poly_t</tt>.
    *
    * Frees all memory used by the pointed <tt>siqs_poly_t</tt> and then
    * frees the \tt{poly} pointer.
    *
    * \warning Do not call \tt{free(poly)} in client code after a call to
    * \tt{free_siqs_poly(poly)}: it would result in an error.
    *
    * \param poly the \tt{siqs_poly_t} to free.
    */
void free_siqs_poly(siqs_poly_t* poly);

   /**
    * \brief Updates a polynomial.
    *
    * Updates the polynomial \tt{poly} by either, deriving a new \tt{b} value
    * (the so-called "fast" initialization) or by computing a new leading
    * coefficient (the "full" or "slow" initilization).
    * 
    * \param poly the polynomial to update.
    *
    * \return An error code (either \tt{SUCCESS} or \tt{FATAL_INTERNAL_ERROR})
    */
ecode_t update_polynomial(siqs_poly_t* const poly);

   /**
    * \brief Returns the number of "full" initialization performed.
    *
    * This is also the number of distinct \tt{a} used.
    * 
    * \param poly the polynomial used.
    *
    * \return The number of "full" initialization performed.
    */
int na_used(siqs_poly_t* const poly);

#ifdef __cplusplus
}
#endif

#endif
