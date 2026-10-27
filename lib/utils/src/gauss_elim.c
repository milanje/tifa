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
 * \file    gauss_elim.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdio.h>
#include <stdlib.h>

#include "gauss_elim.h"
#include "array.h"
#include "tifa_alloc.h"

//-----------------------------------------------------------------------------
// Ref. "A compact algorithm for Gaussian elimination over GF(2)
//       implemented on highly parallel computers",
//       Dennis Parkinson, Marvin Wunderlich,
//       Parallel Computing 1 (1984) 65-73
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
#if TIFA_SIZEOF_BITSTRING_T <= 8
//
// The rightmost column before c (the largest index below c) with a one in
// the row r, or NO_SUCH_ROW. Column j is bit BITSTRING_T_BITSIZE - 1 -
// j % BITSTRING_T_BITSIZE of word j / BITSTRING_T_BITSIZE.
//
static inline uint32_t last_col_before(const TIFA_BITSTRING_T* const r,
                                       const uint32_t c) {
    if (c == 0) {
        return NO_SUCH_ROW;
    }
    uint32_t w = (c - 1) / BITSTRING_T_BITSIZE;
    const uint32_t pos = BITSTRING_T_BITSIZE - 1 - (c - 1) % BITSTRING_T_BITSIZE;
    unsigned long long bits = r[w] & ~((((TIFA_BITSTRING_T) 1) << pos) - 1);
    while (bits == 0) {
        if (w == 0) {
            return NO_SUCH_ROW;
        }
        bits = r[--w];
    }
    return (w + 1) * BITSTRING_T_BITSIZE - 1 - (uint32_t) __builtin_ctzll(bits);
}
#endif
//-----------------------------------------------------------------------------
void gaussian_elim(uint32_array_list_t* relations,
                   binary_matrix_t* const matrix) {
    //
    // Parkinson and Wunderlich's algorithm is straitforwardly implemented
    // without any kind of enhancement.
    //
    // _WARNING_: The aformentionned paper contains an omission in the
    //            description of the algorithm. See comments for futher
    //            explanation...
    //
    const uint32_t nrows = matrix->nrows;
    const uint32_t ncols = matrix->ncols;

    uint32_t row = 0;

    binary_array_t* used  = alloc_binary_array(nrows);
    uint32_array_t* pivot = alloc_uint32_array(ncols);
    used->length  = nrows;
    
    //
    // Stage 1: Process the matrix
    //
    // _NOTE_: We process the matrix from right to left since the rightmost
    //         part of the matrix is sparser than the leftmost part
    //
    TIFA_BITSTRING_T** const data = matrix->data;
    const uint32_t nwords = matrix->ncols_alloced;
#if TIFA_SIZEOF_BITSTRING_T <= 8
    //
    // Columns are processed from right to left, so the one column a row can
    // be a pivot for next is its rightmost unprocessed one: rows are kept in
    // a list per such column, and a column only visits the rows with a one
    // on it. Any of them is as good a pivot as P & W's first: the pivot just
    // leaves the lists, instead of being cleared so that no column picks it.
    //
    uint32_t* const head = tifa_malloc(ncols * sizeof(uint32_t));
    uint32_t* const next = tifa_malloc(nrows * sizeof(uint32_t));
    for (uint32_t icol = 0U; icol < ncols; icol++) {
        head[icol] = NO_SUCH_ROW;
    }
    for (uint32_t irow = nrows - 1; irow != UINT32_MAX; irow--) {
        const uint32_t c = last_col_before(data[irow], ncols);
        if (c != NO_SUCH_ROW) {
            next[irow] = head[c];
            head[c]    = irow;
        }
    }
    for (uint32_t icol = ncols - 1; icol != UINT32_MAX; icol--) {
        row = head[icol];
        if (row == NO_SUCH_ROW) {
            continue;
        }
        pivot->data[icol] = row;
        set_array_bit_to_one(row, used);

        const uint32_t         w    = icol / BITSTRING_T_BITSIZE;
        const TIFA_BITSTRING_T mask = ((TIFA_BITSTRING_T) 1)
                                      << (BITSTRING_T_BITSIZE - 1
                                          - icol % BITSTRING_T_BITSIZE);
        const TIFA_BITSTRING_T* const prow = data[row];
        uint32_t irow = next[row];
        while (irow != NO_SUCH_ROW) {
            const uint32_t inext = next[irow];
            TIFA_BITSTRING_T* const r = data[irow];
            for (uint32_t jw = 0U; jw < nwords; jw++) {
                r[jw] ^= prow[jw];
            }
            r[w] |= mask;
            const uint32_t c = last_col_before(r, icol);
            if (c != NO_SUCH_ROW) {
                next[irow] = head[c];
                head[c]    = irow;
            }
            irow = inext;
        }
    }
    free(head);
    free(next);
#else
    for (uint32_t icol = ncols - 1; icol != UINT32_MAX; icol--) {
        row = first_row_with_one_on_col(icol, matrix);
        if (row == NO_SUCH_ROW) {
            continue;
        }
        pivot->data[icol] = row;

        for (uint32_t irow = row + 1; irow < nrows; irow++) {
            if (1 == get_matrix_bit(irow, icol, matrix)) {
                for (uint32_t jw = 0U; jw < nwords; jw++) {
                    data[irow][jw] ^= data[row][jw];
                }
                set_matrix_bit_to_one(irow, icol, matrix);
            }
        }
        //
        // _WARNING_: The following is completely omitted in the algorithm
        //            description of P & W's paper: clear the pivot row so
        //            that no later column picks it.
        //
        for (uint32_t jw = 0U; jw < nwords; jw++) {
            data[row][jw] = 0U;
        }
        set_array_bit_to_one(row, used);
    }
#endif
    pivot->length = ncols;
    
    //
    // Stage 2: Find dependencies
    //
    for (uint32_t irow = 0U; irow < nrows; irow++) {
        if (relations->length == relations->alloced) {
            break;
        }
        if (0U == get_array_bit(irow, used)) {
            //
            // The ones of the row: one pivot row each, plus the row itself.
            //
            uint32_t nrows_in_rel = 0U;
#if TIFA_SIZEOF_BITSTRING_T <= 8
            //
            // A word at a time (bits past ncols are zero): column j is bit
            // BITSTRING_T_BITSIZE - 1 - j % BITSTRING_T_BITSIZE of its word.
            //
            const TIFA_BITSTRING_T* const r = data[irow];
            for (uint32_t jw = 0U; jw < nwords; jw++) {
                nrows_in_rel += (uint32_t) __builtin_popcountll(r[jw]);
            }
            uint32_array_t* rel = alloc_uint32_array(nrows_in_rel + 1);
            rel->length = nrows_in_rel + 1;

            nrows_in_rel = 0U;
            for (uint32_t jw = 0U; jw < nwords; jw++) {
                unsigned long long bits = r[jw];
                while (bits != 0) {
                    const int top = 63 - __builtin_clzll(bits);
                    rel->data[nrows_in_rel++]
                        = pivot->data[(jw + 1) * BITSTRING_T_BITSIZE - 1 - top];
                    bits &= ~(1ULL << top);
                }
            }
#else
            for (uint32_t icol = 0U; icol < ncols; icol++) {
                if (1 == get_matrix_bit(irow, icol, matrix)) {
                    nrows_in_rel++;
                }
            }
            uint32_array_t* rel = alloc_uint32_array(nrows_in_rel + 1);
            rel->length = nrows_in_rel + 1;

            nrows_in_rel = 0U;
            for (uint32_t icol = 0U; icol < ncols; icol++) {
                if (1 == get_matrix_bit(irow, icol, matrix)) {
                    rel->data[nrows_in_rel++] = pivot->data[icol];
                }
            }
#endif
            rel->data[nrows_in_rel] = irow;
            add_entry_in_uint32_array_list(rel, relations);
        }
    }
    free_binary_array(used);
    free_uint32_array(pivot);
}
//-----------------------------------------------------------------------------
