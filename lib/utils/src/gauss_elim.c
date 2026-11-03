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
#include <stdbool.h>
#include <string.h>

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
#if TIFA_SIZEOF_BITSTRING_T <= 8
//
// The ones of the row r on the columns before c (still to process).
//
static inline __attribute__((always_inline))
uint32_t weight_before(const TIFA_BITSTRING_T* const r, const uint32_t c) {
    const uint32_t w   = c / BITSTRING_T_BITSIZE;
    const uint32_t pos = BITSTRING_T_BITSIZE - 1 - c % BITSTRING_T_BITSIZE;
    uint32_t       k   = 0;
    for (uint32_t jw = 0; jw < w; jw++) {
        k += (uint32_t) __builtin_popcountll(r[jw]);
    }
    const TIFA_BITSTRING_T above = (pos + 1 < BITSTRING_T_BITSIZE)
                                 ? ~((((TIFA_BITSTRING_T) 1) << (pos + 1)) - 1)
                                 : 0;
    return k + (uint32_t) __builtin_popcountll(r[w] & above);
}
//
// Parkinson and Wunderlich's elimination (see gaussian_elim), with the
// lightest row of a column's list as its pivot when lightest is set: its ones
// on the columns still to process are what it adds to the other rows.
//
static inline __attribute__((always_inline))
void pw_elim(uint32_array_list_t* relations, binary_matrix_t* const matrix,
             const bool lightest) {
    const uint32_t nrows = matrix->nrows;
    const uint32_t ncols = matrix->ncols;

    uint32_t row = 0;

    binary_array_t* used  = alloc_binary_array(nrows);
    uint32_array_t* pivot = alloc_uint32_array(ncols);
    used->length  = nrows;

    TIFA_BITSTRING_T** const data = matrix->data;
    const uint32_t nwords = matrix->ncols_alloced;
    //
    // Columns are processed from right to left, so the one column a row can
    // be a pivot for next is its rightmost unprocessed one: rows are kept in
    // a list per such column, and a column only visits the rows with a one
    // on it. Any of them is as good a pivot as P & W's first: the pivot just
    // leaves the lists, instead of being cleared so that no column picks it.
    //
    uint32_t* const head   = tifa_malloc(ncols * sizeof(uint32_t));
    uint32_t* const next   = tifa_malloc(nrows * sizeof(uint32_t));
    uint32_t* const weight = lightest ? tifa_malloc(nrows * sizeof(uint32_t))
                                      : NULL;
    for (uint32_t icol = 0U; icol < ncols; icol++) {
        head[icol] = NO_SUCH_ROW;
    }
    for (uint32_t irow = nrows - 1; irow != UINT32_MAX; irow--) {
        const uint32_t c = last_col_before(data[irow], ncols);
        if (c != NO_SUCH_ROW) {
            next[irow] = head[c];
            head[c]    = irow;
            if (lightest) {
                weight[irow] = weight_before(data[irow], c);
            }
        }
    }
    for (uint32_t icol = ncols - 1; icol != UINT32_MAX; icol--) {
        row = head[icol];
        if (row == NO_SUCH_ROW) {
            continue;
        }
        if (lightest) {
            uint32_t prev  = NO_SUCH_ROW;
            uint32_t best  = row;
            uint32_t bprev = NO_SUCH_ROW;
            for (uint32_t r = row; r != NO_SUCH_ROW; prev = r, r = next[r]) {
                if (weight[r] < weight[best]) {
                    best  = r;
                    bprev = prev;
                }
            }
            if (best != row) {
                next[bprev] = next[best];
                next[best]  = row;
                row         = best;
            }
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
                if (lightest) {
                    weight[irow] = weight_before(r, c);
                }
            }
            irow = inext;
        }
    }
    free(head);
    free(next);
    free(weight);
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
            // A word at a time (bits past ncols are zero): column j is bit
            // BITSTRING_T_BITSIZE - 1 - j % BITSTRING_T_BITSIZE of its word.
            //
            uint32_t nrows_in_rel = 0U;
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
            rel->data[nrows_in_rel] = irow;
            add_entry_in_uint32_array_list(rel, relations);
        }
    }
    free_binary_array(used);
    free_uint32_array(pivot);
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
#if TIFA_SIZEOF_BITSTRING_T <= 8
    pw_elim(relations, matrix, false);
#else
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
            rel->data[nrows_in_rel] = irow;
            add_entry_in_uint32_array_list(rel, relations);
        }
    }
    free_binary_array(used);
    free_uint32_array(pivot);
#endif
}
//-----------------------------------------------------------------------------
#if TIFA_SIZEOF_BITSTRING_T <= 8
#if defined(__x86_64__) && !defined(__POPCNT__) && defined(__GNUC__)
//
// The lightest pivots weigh the rows by popcounts: without POPCNT, which
// the x86-64 baseline lacks, GCC calls a function for each word. A copy
// compiled with it runs where the processor has it.
//
#define PW_ELIM_POPCNT 1

__attribute__((target("popcnt")))
static void pw_elim_lightest_popcnt(uint32_array_list_t* relations,
                                    binary_matrix_t* const matrix) {
    pw_elim(relations, matrix, true);
}
#endif

static void pw_elim_lightest(uint32_array_list_t* relations,
                             binary_matrix_t* const matrix) {
#if defined(PW_ELIM_POPCNT)
    if (__builtin_cpu_supports("popcnt")) {
        pw_elim_lightest_popcnt(relations, matrix);
        return;
    }
#endif
    pw_elim(relations, matrix, true);
}
#endif
#if TIFA_SIZEOF_BITSTRING_T <= 8
//
// The columns of the ones of a row of nwords words, one after the other:
// column j is bit BITSTRING_T_BITSIZE - 1 - j % BITSTRING_T_BITSIZE of word
// j / BITSTRING_T_BITSIZE.
//
#define FOR_EACH_COLUMN(ROW, NWORDS, COL, BODY)                              \
    for (uint32_t jw_ = 0U; jw_ < (NWORDS); jw_++) {                          \
        unsigned long long bits_ = (ROW)[jw_];                                \
        while (bits_ != 0) {                                                  \
            const int top_ = 63 - __builtin_clzll(bits_);                     \
            bits_ &= ~(1ULL << top_);                                         \
            const uint32_t COL = (jw_ + 1) * BITSTRING_T_BITSIZE - 1 - top_;  \
            BODY                                                              \
        }                                                                     \
    }

typedef struct {
    const binary_matrix_t* matrix;
    uint32_t*              colweight;
    uint8_t*               dead;
    const uint32_t*        start;
    const uint32_t*        rows;
    uint32_t*              stack;
    uint32_t               sp;
    uint32_t               nalive;
    uint32_t               nactive;
} filter_t;

static void kill_row(filter_t* const f, const uint32_t row) {
    f->dead[row] = 1;
    f->nalive--;
    FOR_EACH_COLUMN(f->matrix->data[row], f->matrix->ncols_alloced, c, {
        f->colweight[c]--;
        if (f->colweight[c] == 1) {
            f->stack[f->sp++] = c;
        } else if (f->colweight[c] == 0) {
            f->nactive--;
        }
    })
}
//
// A row with a column of its own is in no dependency: those rows go, and so
// do the rows the columns they leave alone are then of.
//
static void kill_singletons(filter_t* const f) {
    while (f->sp != 0) {
        const uint32_t c = f->stack[--f->sp];
        if (f->colweight[c] != 1) {
            continue;
        }
        for (uint32_t k = f->start[c]; k < f->start[c + 1]; k++) {
            if (!f->dead[f->rows[k]]) {
                kill_row(f, f->rows[k]);
                break;
            }
        }
    }
}

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*) a;
    const uint64_t y = *(const uint64_t*) b;
    return (x > y) - (x < y);
}
//
// The rows equal to another one, but for the first of them: their pairs are
// dependencies of their own, trivial ones when a relation was found twice,
// which would take the place of the others (see filtered_gaussian_elim).
//
static void kill_duplicates(filter_t* const f) {
    const binary_matrix_t* const m = f->matrix;
    const uint32_t nwords = m->ncols_alloced;

    uint64_t* const keys = tifa_malloc((m->nrows + 1) * sizeof(uint64_t));
    uint32_t n = 0;
    for (uint32_t r = 0; r < m->nrows; r++) {
        if (f->dead[r]) {
            continue;
        }
        uint64_t h = 0x9E3779B97F4A7C15ULL;
        for (uint32_t jw = 0; jw < nwords; jw++) {
            h = (h ^ (uint64_t) m->data[r][jw]) * 0x100000001B3ULL;
            h ^= h >> 29;
        }
        //
        // The hash in the high bits, the row in the low ones: equal rows end
        // up next to one another, the first first.
        //
        keys[n++] = (h & ~(uint64_t) UINT32_MAX) | r;
    }
    qsort(keys, n, sizeof(uint64_t), compare_u64);
    for (uint32_t k = 1; k < n; k++) {
        const uint32_t r = (uint32_t) keys[k];
        for (uint32_t j = k; j > 0; j--) {
            if ((keys[j - 1] ^ keys[k]) >> 32) {
                break;
            }
            const uint32_t q = (uint32_t) keys[j - 1];
            if (!f->dead[q]
                && (0 == memcmp(m->data[q], m->data[r],
                                nwords * sizeof(TIFA_BITSTRING_T)))) {
                kill_row(f, r);
                break;
            }
        }
    }
    free(keys);
}

//
// The indexes i in [0, n) by decreasing weights[i] (then by increasing i),
// into order.
//
static void sort_heaviest_first(uint32_t* const order,
                                const uint32_t* const indexes,
                                const uint32_t n,
                                const uint32_t* const weights) {
    uint64_t* const keys = tifa_malloc((n + 1) * sizeof(uint64_t));
    for (uint32_t k = 0; k < n; k++) {
        keys[k] = ((uint64_t) (UINT32_MAX - weights[indexes[k]]) << 32)
                | indexes[k];
    }
    qsort(keys, n, sizeof(uint64_t), compare_u64);
    for (uint32_t k = 0; k < n; k++) {
        order[k] = (uint32_t) keys[k];
    }
    free(keys);
}
#endif
//-----------------------------------------------------------------------------
void filtered_gaussian_elim(uint32_array_list_t* relations,
                            binary_matrix_t* const matrix) {
#if TIFA_SIZEOF_BITSTRING_T <= 8
    //
    // The rows in no dependency go first (see kill_singletons), and the
    // copies of others (see kill_duplicates), then the heaviest ones while
    // there are more rows than columns plus the excess of the matrix (nrows -
    // ncols, the dependencies wanted); the columns left are packed, the
    // heaviest first, so that the elimination from the right meets them last.
    // On the SIQS matrices of 180 to 220 bits, with the lightest pivots
    // (see pw_elim): 0.29 to 0.46 of gaussian_elim's time.
    //
    const uint32_t nrows  = matrix->nrows;
    const uint32_t ncols  = matrix->ncols;
    const uint32_t nwords = matrix->ncols_alloced;

    uint32_t* const colweight = tifa_calloc(ncols + 1, sizeof(uint32_t));
    uint32_t* const rowweight = tifa_calloc(nrows + 1, sizeof(uint32_t));
    uint8_t*  const dead      = tifa_calloc(nrows + 1, 1);
    uint32_t* const start     = tifa_calloc(ncols + 1, sizeof(uint32_t));
    uint32_t* const stack     = tifa_malloc((ncols + 1) * sizeof(uint32_t));

    uint32_t nones = 0;
    for (uint32_t r = 0; r < nrows; r++) {
        FOR_EACH_COLUMN(matrix->data[r], nwords, c, {
            colweight[c]++;
            rowweight[r]++;
            nones++;
        })
    }
    //
    // The rows of each column, column after column.
    //
    for (uint32_t c = 0; c < ncols; c++) {
        start[c + 1] = start[c] + colweight[c];
    }
    uint32_t* const fill = tifa_malloc((ncols + 1) * sizeof(uint32_t));
    uint32_t* const rows = tifa_malloc((nones + 1) * sizeof(uint32_t));
    memcpy(fill, start, ncols * sizeof(uint32_t));
    for (uint32_t r = 0; r < nrows; r++) {
        FOR_EACH_COLUMN(matrix->data[r], nwords, c, {
            rows[fill[c]++] = r;
        })
    }
    free(fill);

    filter_t f = {
        .matrix = matrix, .colweight = colweight, .dead = dead,
        .start = start, .rows = rows, .stack = stack, .sp = 0,
        .nalive = nrows, .nactive = 0
    };
    for (uint32_t c = 0; c < ncols; c++) {
        f.nactive += (colweight[c] != 0);
        if (colweight[c] == 1) {
            stack[f.sp++] = c;
        }
    }
    kill_singletons(&f);
    kill_duplicates(&f);
    kill_singletons(&f);

    if (nrows > ncols) {
        const uint32_t excess = nrows - ncols;
        uint32_t* const order = tifa_malloc(nrows * sizeof(uint32_t));
        for (uint32_t r = 0; r < nrows; r++) {
            order[r] = r;
        }
        sort_heaviest_first(order, order, nrows, rowweight);
        for (uint32_t k = 0; (k < nrows) && (f.nalive > f.nactive + excess);
             k++) {
            if (!dead[order[k]]) {
                kill_row(&f, order[k]);
                kill_singletons(&f);
            }
        }
        free(order);
    }
    //
    // The columns left, the heaviest first, and the rows left.
    //
    uint32_t* const order = tifa_malloc((ncols + 1) * sizeof(uint32_t));
    uint32_t* const cmap  = tifa_malloc((ncols + 1) * sizeof(uint32_t));
    uint32_t ncols_left = 0;
    for (uint32_t c = 0; c < ncols; c++) {
        if (colweight[c] != 0) {
            order[ncols_left++] = c;
        }
    }
    sort_heaviest_first(order, order, ncols_left, colweight);
    for (uint32_t k = 0; k < ncols_left; k++) {
        cmap[order[k]] = k;
    }
    uint32_t* const rmap = tifa_malloc((f.nalive + 1) * sizeof(uint32_t));
    uint32_t nrows_left = 0;
    for (uint32_t r = 0; r < nrows; r++) {
        if (!dead[r]) {
            rmap[nrows_left++] = r;
        }
    }
    if (ncols_left == 0) {
        //
        // Every row left is zero: a dependency on its own.
        //
        for (uint32_t k = 0; (k < nrows_left)
                             && (relations->length < relations->alloced); k++) {
            uint32_array_t* rel = alloc_uint32_array(1);
            rel->data[0] = rmap[k];
            rel->length  = 1;
            add_entry_in_uint32_array_list(rel, relations);
        }
    } else if (nrows_left != 0) {
        binary_matrix_t* const packed = alloc_binary_matrix(nrows_left,
                                                            ncols_left);
        packed->nrows = nrows_left;
        packed->ncols = ncols_left;
        for (uint32_t k = 0; k < nrows_left; k++) {
            TIFA_BITSTRING_T* const dst = packed->data[k];
            memset(dst, 0, packed->ncols_alloced * sizeof(TIFA_BITSTRING_T));
            FOR_EACH_COLUMN(matrix->data[rmap[k]], nwords, c, {
                if (colweight[c] != 0) {
                    const uint32_t j = cmap[c];
                    dst[j / BITSTRING_T_BITSIZE] |=
                        ((TIFA_BITSTRING_T) 1)
                        << (BITSTRING_T_BITSIZE - 1 - j % BITSTRING_T_BITSIZE);
                }
            })
        }
        const uint32_t first = relations->length;
        pw_elim_lightest(relations, packed);
        for (uint32_t d = first; d < relations->length; d++) {
            uint32_array_t* const rel = relations->data[d];
            for (uint32_t k = 0; k < rel->length; k++) {
                rel->data[k] = rmap[rel->data[k]];
            }
        }
        free_binary_matrix(packed);
    }
    free(order);
    free(cmap);
    free(rmap);
    free(colweight);
    free(rowweight);
    free(dead);
    free(start);
    free(stack);
    free(rows);
#else
    gaussian_elim(relations, matrix);
#endif
}
//-----------------------------------------------------------------------------
