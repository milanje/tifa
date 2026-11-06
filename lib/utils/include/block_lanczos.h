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
 * \file    block_lanczos.h
 * \author  Jerome Milan
 *
 * \brief Montgomery's block Lanczos over GF(2), on sparse matrices.
 *
 * The dependencies among the relations of a factorization: sets of
 * relations whose exponent vectors sum to zero modulo 2. The relations are
 * the columns of a sparse matrix whose rows are the primes (and the sign).
 * Unlike the Gaussian eliminations of gauss_elim.h, nothing dense is
 * stored: memory in the number of nonzero entries, time in the number of
 * relations times that number.
 */

#if !defined(_TIFA_BLOCK_LANCZOS_H_)
   /**
    * \def _TIFA_BLOCK_LANCZOS_H_
    * Standard include guard.
    */
#define _TIFA_BLOCK_LANCZOS_H_

#include <inttypes.h>

#include <gmp.h>

#include "array.h"
#include "lindep.h"
#include "matrix.h"
#include "x_array_list.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \struct struct_sparse_gf2_t block_lanczos.h
    * \brief A matrix over GF(2) by columns: the rows of the ones of column
    * \c c are <tt>rows[start[c]]</tt> to <tt>rows[start[c + 1] - 1]</tt>.
    */
struct struct_sparse_gf2_t {
    uint32_t  nrows;
    uint32_t  ncols;
    uint32_t* start;
    uint32_t* rows;
    uint32_t  ncols_alloced;
    uint32_t  nnz_alloced;
};
   /**
    * \typedef sparse_gf2_t
    * \brief Equivalent to <tt>struct struct_sparse_gf2_t</tt>.
    */
typedef struct struct_sparse_gf2_t sparse_gf2_t;

   /**
    * \brief An empty matrix of \c nrows rows, its columns to be appended.
    */
sparse_gf2_t* alloc_sparse_gf2(uint32_t nrows, uint32_t ncols_alloced);

   /**
    * \brief Frees a matrix allocated by \c alloc_sparse_gf2.
    */
void free_sparse_gf2(sparse_gf2_t* m);

   /**
    * \brief Appends a column with ones at the \c n rows of \c rows.
    */
void append_sparse_gf2_column(sparse_gf2_t* m, const uint32_t* rows,
                              uint32_t n);

   /**
    * \brief The transpose of a dense matrix whose rows are relations.
    *
    * The rows of \c dense (its first \c dense->nrows) become the columns.
    */
sparse_gf2_t* sparse_gf2_from_relations(const binary_matrix_t* dense);

   /**
    * \brief The matrix of the relations of a factorization, from their
    * exponents.
    *
    * Column \c i is the i-th relation, for the first
    * <tt>values->length</tt> of \c decomp: a one at row 0 if
    * <tt>values->data[i]</tt> is negative, and at row <tt>j + 1</tt> if the
    * exponent of the j-th prime of \c decomp in it is odd. \c decomp is
    * sorted (see \c sort_decomp) if sparse.
    */
sparse_gf2_t* sparse_gf2_from_decomp(const decomp_t* decomp,
                                     const mpz_array_t* values);

   /**
    * \brief Dependencies among the columns of \c m, by block Lanczos.
    *
    * Appends to \c deps up to its free room (at most 64) dependencies, each
    * the sorted indices of columns of \c m whose sum is zero, every one
    * checked. Columns in no dependency (with a row of their own) are left
    * out first, then the heaviest ones beyond an excess of a hundred columns
    * over the rows they use. The time limit and the interrupt check of the
    * calling thread stop it, with what it has found then (nothing).
    *
    * \param[in] m The matrix.
    * \param[in,out] deps The list of dependencies to append to.
    * \return The number of dependencies appended.
    */
uint32_t block_lanczos(const sparse_gf2_t* m, uint32_array_list_t* deps);

#ifdef __cplusplus
}
#endif

#endif
