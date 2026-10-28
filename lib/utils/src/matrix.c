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
 * \file    matrix.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "matrix.h"

#include "bitstring_t.h"
#include "tifa_alloc.h"

/*
 *-----------------------------------------------------------------------------
 *              binary_matrix_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
void reset_binary_matrix(binary_matrix_t* const matrix) {
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
#if TIFA_USE_CALLOC_MEMSET
        memset(matrix->data[i], 0x00,
               matrix->ncols_alloced * SIZEOF_BITSTRING_T);
#else
        for (uint32_t j = 0U; j < matrix->ncols_alloced; j++) {
            matrix->data[i][j] = 0;
            //
            // Use an xor instead of simply 0 to make sure we have a bitstring
            // made of 0 bits only. Granted, TIFA is not likely to compile on
            // the exotic architectures where the value 0 is not made of 0 bits
            // only anyways, but still...
            //
            matrix->data[i][j] ^= 0;
        }
#endif
    }
}
//-----------------------------------------------------------------------------
binary_matrix_t* alloc_binary_matrix(uint32_t nrows, uint32_t ncols) {

    binary_matrix_t* matrix = tifa_malloc(sizeof(binary_matrix_t));

    matrix->nrows_alloced = nrows;

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    matrix->ncols_alloced = ((ncols - 1) >> POW_TWO_BITSTRING_T_SIZE) + 1;
#else
    matrix->ncols_alloced = (ncols - 1)/BITSTRING_T_BITSIZE + 1;
#endif

    matrix->nrows = nrows;
    matrix->ncols = ncols;
    //
    // The rows are one block, kept past the last row: one allocation, not
    // one per row (2.6% of a 70-bit SIQS run with the byte matrices).
    //
    matrix->data  = tifa_malloc((matrix->nrows_alloced + 1)
                                * sizeof(TIFA_BITSTRING_T*));
    const size_t nwords = (size_t) matrix->nrows_alloced
                        * matrix->ncols_alloced;
#if TIFA_USE_CALLOC_MEMSET
    TIFA_BITSTRING_T* const block = tifa_calloc(nwords, SIZEOF_BITSTRING_T);
#else
    TIFA_BITSTRING_T* const block = tifa_malloc(nwords * SIZEOF_BITSTRING_T);
#endif
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
        matrix->data[i] = block + (size_t) i * matrix->ncols_alloced;
    }
    matrix->data[matrix->nrows_alloced] = block;

#if ! TIFA_USE_CALLOC_MEMSET
    reset_binary_matrix(matrix);
#endif

    return matrix;
}
//-----------------------------------------------------------------------------
binary_matrix_t* clone_binary_matrix(const binary_matrix_t * const matrix) {

    binary_matrix_t* clone = tifa_malloc(sizeof(binary_matrix_t));

    clone->nrows_alloced = matrix->nrows_alloced;
    clone->ncols_alloced = matrix->ncols_alloced;
    clone->nrows = matrix->nrows;
    clone->ncols = matrix->ncols;
    clone->data  = tifa_malloc((clone->nrows_alloced + 1)
                               * sizeof(TIFA_BITSTRING_T*));

    size_t rowsize = clone->ncols_alloced * SIZEOF_BITSTRING_T;
    TIFA_BITSTRING_T* const block = tifa_malloc(clone->nrows_alloced
                                                * rowsize);

    for (uint32_t i = 0U; i < clone->nrows_alloced; i++) {
        clone->data[i] = block + (size_t) i * clone->ncols_alloced;
        memcpy(clone->data[i], matrix->data[i], rowsize);
    }
    clone->data[clone->nrows_alloced] = block;
    return clone;
}
//-----------------------------------------------------------------------------
void free_binary_matrix(binary_matrix_t* const matrix) {
    free(matrix->data[matrix->nrows_alloced]);
    free(matrix->data);
    free(matrix);
}
//-----------------------------------------------------------------------------
void print_binary_matrix(const binary_matrix_t* const matrix) {
    //
    // Mostly for debugging purposes...
    //
    for (uint32_t i = 0U; i < matrix->nrows; i++) {
        printf("row %4u : ", i);
        for (uint32_t j = 0U; j < matrix->ncols; j++) {
            if (0 != get_matrix_bit(i, j, matrix)) {
                printf("1");
            } else {
                printf(".");
            }
        }
        printf("\n");
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              byte_matrix_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
void reset_byte_matrix(byte_matrix_t* const matrix) {
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
#if TIFA_USE_CALLOC_MEMSET
        memset(
            matrix->data[i],
            0x00,
            matrix->ncols_alloced * sizeof(unsigned char)
        );
#else
        for (uint32_t j = 0U; j < matrix->ncols_alloced; j++) {
            matrix->data[i][j] = 0;
        }
#endif
    }
}
//-----------------------------------------------------------------------------
byte_matrix_t* alloc_byte_matrix(uint32_t nrows, uint32_t ncols) {

    byte_matrix_t* matrix = tifa_malloc(sizeof(byte_matrix_t));

    matrix->nrows_alloced = nrows;
    matrix->ncols_alloced = ncols;

    matrix->nrows = nrows;
    matrix->ncols = ncols;

    //
    // The rows are one block, kept past the last row (see
    // alloc_binary_matrix).
    //
    matrix->data  = tifa_malloc((matrix->nrows_alloced + 1)
                                * sizeof(unsigned char*));
    const size_t nbytes = (size_t) matrix->nrows_alloced
                        * matrix->ncols_alloced;
#if TIFA_USE_CALLOC_MEMSET
    unsigned char* const block = tifa_calloc(nbytes, sizeof(unsigned char));
#else
    unsigned char* const block = tifa_malloc(nbytes);
#endif
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
        matrix->data[i] = block + (size_t) i * matrix->ncols_alloced;
    }
    matrix->data[matrix->nrows_alloced] = block;

#if ! TIFA_USE_CALLOC_MEMSET
    reset_byte_matrix(matrix);
#endif

    return matrix;
}
//-----------------------------------------------------------------------------
byte_matrix_t* clone_byte_matrix(const byte_matrix_t * const matrix) {

    byte_matrix_t* clone = tifa_malloc(sizeof(byte_matrix_t));

    clone->nrows_alloced = matrix->nrows_alloced;
    clone->ncols_alloced = matrix->ncols_alloced;
    clone->nrows = matrix->nrows;
    clone->ncols = matrix->ncols;
    clone->data  = tifa_malloc((clone->nrows_alloced + 1)
                               * sizeof(unsigned char*));

    size_t rowsize = clone->ncols_alloced * sizeof(unsigned char);
    unsigned char* const block = tifa_malloc(clone->nrows_alloced * rowsize);

    for (uint32_t i = 0U; i < clone->nrows_alloced; i++) {
        clone->data[i] = block + (size_t) i * clone->ncols_alloced;
        memcpy(clone->data[i], matrix->data[i], rowsize);
    }
    clone->data[clone->nrows_alloced] = block;
    return clone;
}
//-----------------------------------------------------------------------------
void free_byte_matrix(byte_matrix_t* const matrix) {
    free(matrix->data[matrix->nrows_alloced]);
    free(matrix->data);
    free(matrix);
}
//-----------------------------------------------------------------------------
void print_byte_matrix(const byte_matrix_t* const matrix) {
    //
    // Mostly for debugging purposes...
    //
    for (uint32_t i = 0U; i < matrix->nrows; i++) {
        printf("row %4u : ", i);
        for (uint32_t j = 0U; j < matrix->ncols; j++) {
            printf("%4u", (unsigned int)matrix->data[i][j]);
        }
        printf("\n");
    }
}
//-----------------------------------------------------------------------------
