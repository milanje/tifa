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
 * \file    lindep.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Functions used in the resolution of the linear systems.
 *
 * Functions used in the resolution of the linear systems over GF(2) found
 * in factorization problems and in the very last stage of the factorization
 * process (determination of factors after linear algebra phase).
 */

#if !defined(_TIFA_LINDEP_H_)
   /**
    * \def _TIFA_LINDEP_H_
    * Standard include guard.
    */
#define _TIFA_LINDEP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>

#include "array.h"
#include "x_array_list.h"
#include "matrix.h"
#include "exit_codes.h"

   /**
    * \enum linalg_method_enum
    *
    * Enumeration listing the different linear system resolution method
    * implemented.
    *
    * For the time being, only one method is available.
    */
enum linalg_method_enum {
       /**
        * "Smart" gaussian elimination described in:
        * "A compact algorithm for Gaussian elimination over GF(2) implemented
        * on highly parallel computers", by D. Parkinson and M. Wunderlich
        * (Parallel Computing 1 (1984) 65-73).
        */
    SMART_GAUSS_ELIM = 0
};

   /**
    * \typedef linalg_method_t
    * \brief Equivalent to <tt>enum linalg_method_enum</tt>.
    */
typedef enum linalg_method_enum linalg_method_t;

   /**
    * Global constant array mapping linalg methods to their string 
    * representations.
    */
static const char* const linalg_method_to_str[1] = {
    "smart gaussian elimination"
};

   /**
    * \brief Fills a binary matrix via trial divisions.
    *
    * Fills the binary matrix \c matrix by trial divisions of the integers
    * listed in \c to_factor by the integers listed in <tt>factor_base</tt>.
    * After these trial divisions, each partially factored integer
    * from \c to_factor are stored in <tt>partially_factored</tt>.
    *
    * \note The binary matrix is filled so that:
    *
    * \li There is a '1' in the (i-th row, 1st col) position in the matrix
    * if <tt>to_factor->data[i]</tt> is negative.
    * \li There is a '1' in the (i-th row, j-th col) position in the matrix
    * if <tt>to_factor->data[i]</tt> is divisible by an odd power of
    * <tt>factor_base->data[j]</tt>.
    * \li In all other cases, the (i-th row, j-th col) position in the matrix
    * contains a 0.
    *
    * \param[out] matrix A pointer to the binary matrix to fill.
    * \param[out] partially_factored A pointer to the partially factored
    *                                integers.
    * \param[in] to_factor A pointer to the array listing the integers to
    *                      factor.
    * \param[in] factor_base  A pointer to the array listing the integers
    *                         to trial divide by.
    */
void fill_matrix_trial_div(
          binary_matrix_t* const matrix,
          mpz_array_t*     const partially_factored,
    const mpz_array_t*     const to_factor,
    const uint32_array_t*  const factor_base
);
   
   /**
    * \brief Similar to \c fill_matrix_trial_div but also stores valuations.
    *
    * Fills the binary matrix \c matrix by trial divisions of the integers
    * listed in \c to_factor by the integers listed in <tt>factor_base</tt>.
    * After these trial divisions, each partially factored integer
    * from \c to_factor are stored in <tt>partially_factored</tt>.
    *
    * Also stores in \c decomp_matrix the valuation of each integer
    * from \c to_factor for each prime in the factor base. For example the
    * valuation of <tt>to_factor->data[i]</tt> for the prime given by
    * <tt>factor_base->data[j]</tt> will be stored in
    * <tt>decomp_matrix->data[i][j]</tt>.
    *
    * \note The binary matrix is filled so that:
    *
    * \li There is a '1' in the (i-th row, 1st col) position in the matrix
    * if <tt>to_factor->data[i]</tt> is negative.
    * \li There is a '1' in the (i-th row, j-th col) position in the matrix
    * if <tt>to_factor->data[i]</tt> is divisible by an odd power of
    * <tt>factor_base->data[j]</tt>.
    * \li In all other cases, the (i-th row, j-th col) position in the matrix
    * contains a 0.
    *
    * \param[out] matrix        A pointer to the binary matrix to fill.
    * \param[out] decomp_matrix A pointer to the byte matrix to fill.
    * \param[out] partially_factored A pointer to the partially factored
    *                                integers.
    * \param[in] to_factor A pointer to the array listing the integers to
    *                      factor.
    * \param[in] factor_base  A pointer to the array listing the integers
    *                         to trial divide by.
    */
void fill_trial_div_decomp(
          binary_matrix_t* const matrix,
          byte_matrix_t*   const decomp_matrix,
          mpz_array_t*     const partially_factored,
    const mpz_array_t*     const to_factor,
    const uint32_array_t*  const factor_base
);

   /**
    * \brief Fills a binary matrix from a list of factors.
    *
    * Fills the binary matrix \c matrix from a previously computed list
    * giving all known factors.
    *
    * \note \c list countains the previously computed factors of each
    * integers in <tt>smooth_array</tt>, in other words, we know that
    * <tt>smooth_array->data[i]</tt> is divisible by all the integers
    * of the \c uint32_array_t given by <tt>list->data[i]</tt>.
    *
    * The binary matrix is filled so that:
    *
    * \li There is a '1' in the (i-th row, 1st col) position in the matrix
    * if <tt>smooth_factor->data[i]</tt> is negative.
    * \li There is a '1' in the (i-th row, j-th col) position in the matrix
    * if <tt>smooth_factor->data[i]</tt> is divisible by an odd power of
    * <tt>list->data[i]->data[k]</tt>. Here j is found so that
    * <tt>factor_base->data[j] = list->data[i]->data[k]</tt>.
    * \li In all other cases, the (i-th row, j-th col) position in the matrix
    * contains a 0.
    *
    * \param[out] matrix A pointer to the binary matrix to fill.
    * \param[out] smooth_array A pointer to the array giving the integers
    *                          to factor.
    * \param[in] list A pointer to the factor list for each integer to factor.
    * \param[in] factor_base  A pointer to the array listing the integers
    *                         to trial divide by.
    */
void fill_matrix_from_list(
          binary_matrix_t*     const matrix,
    const mpz_array_t*         const smooth_array,
    const uint32_array_list_t* const list,
    const uint32_array_t*      const factor_base
);

   /**
    * \brief Similar to \c fill_matrix_from_list but also stores valuations.
    *
    * Fills the binary matrix \c matrix from a previously computed list
    * giving all known factors. 
    *
    * Also stores in \c decomp_matrix the valuation of each integer
    * from \c smooth_array for each prime in the factor base. For example the
    * valuation of <tt>smooth_array->data[i]</tt> for the prime given by
    * <tt>factor_base->data[j]</tt> will be stored in
    * <tt>decomp_matrix->data[i][j]</tt>.
    *
    * \note \c list countains the previously computed factors of each
    * integers in <tt>smooth_array</tt>, in other words, we know that
    * <tt>smooth_array->data[i]</tt> is divisible by all the integers
    * of the \c uint32_array_t given by <tt>list->data[i]</tt>.
    *
    * The binary matrix is filled so that:
    *
    * \li There is a '1' in the (i-th row, 1st col) position in the matrix
    * if <tt>smooth_factor->data[i]</tt> is negative.
    * \li There is a '1' in the (i-th row, j-th col) position in the matrix
    * if <tt>smooth_factor->data[i]</tt> is divisible by an odd power of
    * <tt>list->data[i]->data[k]</tt>. Here j is found so that
    * <tt>factor_base->data[j] = list->data[i]->data[k]</tt>.
    * \li In all other cases, the (i-th row, j-th col) position in the matrix
    * contains a 0.
    *
    * \param[out] matrix        A pointer to the binary matrix to fill.
    * \param[out] decomp_matrix A pointer to the byte matrix to fill.
    * \param[out] smooth_array  A pointer to the array giving the integers
    *                           to factor.
    * \param[in] list A pointer to the factor list for each integer to factor.
    * \param[in] factor_base  A pointer to the array listing the integers
    *                         to trial divide by.
    */
void fill_matrix_from_list_decomp(
          binary_matrix_t*     const matrix,
          byte_matrix_t*       const decomp_matrix,
    const mpz_array_t*         const smooth_array,
    const uint32_array_list_t* const list,
    const uint32_array_t*      const factor_base
);

   /**
    * \brief Solves a linear system over GF(2).
    *
    * Solves the linear system over GF(2) given by the binary matrix
    * <tt>matrix</tt>, using the resolution method <tt>method</tt>.
    *
    * \note For the time being, the only implemented method is
    * <tt>SMART_GAUSS_ELIM</tt>.
    *
    * \param[in,out] matrix A pointer to the binary matrix giving the system
    *                       to solve.
    * \param[in] method The linear system resolution method to use.
    */
uint32_array_list_t* find_dependencies(
    binary_matrix_t* const matrix,
    linalg_method_t        method
);

   /**
    * \brief Find factors of an integer from congruence relations.
    *
    * Find factors of the integer \c n from congruence relations of the form
    * <tt>{(x_i)^2}_i = {y_i}_i (mod n)</tt> where <tt>{y_i}_i</tt> is
    * a perfect square.
    *
    * Each entry in \c dependencies gives the list of the aforementioned indexes
    * i so that such a previous relation will hold.
    *
    * Upon termination, returns the following <tt>ecode_t</tt>:
    *
    * \li \c SOME_FACTORS_FOUND   if some factors were found
    * \li \c NO_FACTOR FOUND      is no factor was found
    * \li \c FATAL_INTERNAL_ERROR is case of... a really ugly error!
    *
    * \param[out] factors The factors of \c n found.
    * \param[in] n The integer to factor.
    * \param[in] xi_array A pointer to an array giving the avalaible
    *                     <tt>x_i</tt> values.
    * \param[in] yi_array A pointer to an array giving the avalaible
    *                     <tt>y_i</tt> values.
    * \param[in] dependencies A pointer to an array list giving the sets
    *                         of indexes from which a congruence relation
    *                         can be computed.
    * \return An exit code.
    */
ecode_t find_factors(
          mpz_array_t*         const factors,
    const mpz_t                      n,
    const mpz_array_t*         const xi_array,
    const mpz_array_t*         const yi_array,
    const uint32_array_list_t* const dependencies
);

   /**
    * \brief Similar to \c find_factors but uses the factorization of
    * each \c y_i.
    *
    * Find factors of the integer \c n from congruence relations of the form
    * <tt>{(x_i)^2}_i = {y_i}_i (mod n)</tt> where <tt>{y_i}_i</tt> is
    * a perfect square.
    *
    * Each entry in \c dependencies gives the list of the aforementioned indexes
    * i so that such a previous relation will hold.
    *
    * The difference with the \c find_factors function is that here the \c y_i
    * are not given directly but rather by their factorizations on the factor 
    * base by \c yi_decomp_matrix. For example the valuation of \c y_i for the 
    * prime given by <tt>factor_base->data[j]</tt> is stored in
    * <tt>yi_decomp_matrix->data[i][j]</tt>.
    *
    * Upon termination, returns the following <tt>ecode_t</tt>:
    *
    * \li \c SOME_FACTORS_FOUND   if some factors were found
    * \li \c NO_FACTOR FOUND      is no factor was found
    * \li \c FATAL_INTERNAL_ERROR is case of... a really ugly error!
    *
    * \param[out] factors The factors of \c n found.
    * \param[in] n The integer to factor.
    * \param[in] xi_array A pointer to an array giving the avalaible
    *                     <tt>x_i</tt> values.
    * \param[in] yi_decomp_matrix A pointer to a byte matrix giving the
    *                             factorization of the \c y_i.
    * \param[in] dependencies A pointer to an array list giving the sets
    *                         of indexes from which a congruence relation
    *                         can be computed.
    * \param[in] factor_base  A pointer to the array listing the primes
    *                         in the factor base.
    * \return An exit code.
    */
ecode_t find_factors_decomp(
          mpz_array_t*         const factors,
    const mpz_t                      n,
    const mpz_array_t*         const xi_array,
    const byte_matrix_t*       const yi_decomp_matrix,
    const uint32_array_list_t* const dependencies,
    const uint32_array_t*      const factor_base
);

#ifdef __cplusplus
}
#endif

#endif
