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
 * \file    block_lanczos.c
 * \author  Jerome Milan
 *
 * Montgomery's block Lanczos ("A block Lanczos algorithm for finding
 * dependencies over GF(2)", Eurocrypt 1995), on blocks of 64 vectors, after
 * the presentation of msieve's lanczos.c (Jason Papadopoulos, public domain).
 *
 * The matrix B has the relations as its columns; A = B^T.B is symmetric,
 * and a vector X with A.X = 0 gives, after a last elimination on at most
 * 128 vectors, ones with B.X = 0: dependencies among the relations.
 */

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "tifa_config.h"
#include "tifa_alloc.h"
#include "block_lanczos.h"
#include "stopping.h"

//------------------------------------------------------------------------------
sparse_gf2_t* alloc_sparse_gf2(const uint32_t nrows,
                               const uint32_t ncols_alloced) {
    sparse_gf2_t* const m = tifa_malloc(sizeof(sparse_gf2_t));
    m->nrows         = nrows;
    m->ncols         = 0;
    m->ncols_alloced = (ncols_alloced != 0) ? ncols_alloced : 16;
    m->nnz_alloced   = 16 * m->ncols_alloced;
    m->start         = tifa_malloc(((size_t) m->ncols_alloced + 1)
                                   * sizeof(uint32_t));
    m->rows          = tifa_malloc((size_t) m->nnz_alloced * sizeof(uint32_t));
    m->start[0]      = 0;
    return m;
}
//------------------------------------------------------------------------------
void free_sparse_gf2(sparse_gf2_t* const m) {
    if (m != NULL) {
        free(m->start);
        free(m->rows);
        free(m);
    }
}
//------------------------------------------------------------------------------
void append_sparse_gf2_column(sparse_gf2_t* const m,
                              const uint32_t* const rows, const uint32_t n) {
    if (m->ncols == m->ncols_alloced) {
        m->ncols_alloced *= 2;
        m->start = tifa_realloc(m->start, ((size_t) m->ncols_alloced + 1)
                                          * sizeof(uint32_t));
    }
    const uint32_t used = m->start[m->ncols];
    if (used + n > m->nnz_alloced) {
        while (used + n > m->nnz_alloced) {
            m->nnz_alloced *= 2;
        }
        m->rows = tifa_realloc(m->rows, (size_t) m->nnz_alloced
                                        * sizeof(uint32_t));
    }
    memcpy(m->rows + used, rows, (size_t) n * sizeof(uint32_t));
    m->ncols++;
    m->start[m->ncols] = used + n;
}
//------------------------------------------------------------------------------
sparse_gf2_t* sparse_gf2_from_relations(const binary_matrix_t* const dense) {
    //
    // Column j of a dense row is bit BITSTRING_T_BITSIZE - 1 - j % size of
    // its word j / size (see matrix.h).
    //
    sparse_gf2_t* const m = alloc_sparse_gf2(dense->ncols, dense->nrows);
    uint32_t* const ones  = tifa_malloc(((size_t) dense->ncols + 1)
                                        * sizeof(uint32_t));
    for (uint32_t i = 0; i < dense->nrows; i++) {
        const TIFA_BITSTRING_T* const r = dense->data[i];
        uint32_t n = 0;
        for (uint32_t w = 0; w < dense->ncols_alloced; w++) {
            const TIFA_BITSTRING_T word = r[w];
            if (word == 0) {
                continue;
            }
            for (uint32_t b = 0; b < BITSTRING_T_BITSIZE; b++) {
                if ((word >> (BITSTRING_T_BITSIZE - 1 - b)) & 1) {
                    const uint32_t col = w * BITSTRING_T_BITSIZE + b;
                    if (col < dense->ncols) {
                        ones[n++] = col;
                    }
                }
            }
        }
        append_sparse_gf2_column(m, ones, n);
    }
    free(ones);
    return m;
}

//------------------------------------------------------------------------------
sparse_gf2_t* sparse_gf2_from_decomp(const decomp_t* const decomp,
                                     const mpz_array_t* const values) {
    //
    // The relations are those with a value: decomp has room for more.
    //
    const uint32_t nprimes = decomp->nprimes;
    const uint32_t nrels   = (values->length < decomp->nrows)
                           ? values->length : decomp->nrows;
    sparse_gf2_t* const m  = alloc_sparse_gf2(nprimes + 1, nrels);
    uint32_t* const ones   = tifa_malloc(((size_t) nprimes + 2)
                                         * sizeof(uint32_t));
    //
    // A prime can come more than once in a sparse relation (its 'a', an
    // exponent above 255): the parities of its entries are summed.
    //
    uint8_t* const odd = tifa_calloc((size_t) nprimes + 1, 1);
    for (uint32_t i = 0; i < nrels; i++) {
        uint32_t n = 0;
        if (mpz_sgn(values->data[i]) < 0) {
            ones[n++] = 0;
        }
        if (decomp->dense != NULL) {
            const unsigned char* const row = decomp->dense->data[i];
            for (uint32_t j = 0; j < nprimes; j++) {
                if (row[j] & 1) {
                    ones[n++] = j + 1;
                }
            }
        } else {
            const uint32_t first = n;
            for (uint32_t k = decomp->starts[i]; k < decomp->starts[i + 1];
                 k++) {
                const uint32_t j = decomp->sorted[k] >> 8;
                if ((decomp->sorted[k] & 1) == 0) {
                    continue;
                }
                if (odd[j] == 0) {
                    ones[n++] = j + 1;
                }
                odd[j] ^= 1;
            }
            uint32_t kept = first;
            for (uint32_t k = first; k < n; k++) {
                const uint32_t j = ones[k] - 1;
                if (odd[j]) {
                    ones[kept++] = ones[k];
                }
                odd[j] = 0;
            }
            n = kept;
        }
        append_sparse_gf2_column(m, ones, n);
    }
    free(odd);
    free(ones);
    return m;
}

//------------------------------------------------------------------------------
//                           64 x 64 BLOCK ARITHMETIC
//------------------------------------------------------------------------------
//
// A 64 x 64 matrix is 64 rows of 64 bits, bit j of row i its entry (i, j);
// an N x 64 block is N such rows.
//
#define BIT(i) (UINT64_C(1) << (i))

static void mul_64x64(const uint64_t* const a, const uint64_t* const b,
                      uint64_t* const c) {
    uint64_t r[64];
    for (uint32_t i = 0; i < 64; i++) {
        uint64_t acc = 0;
        uint64_t bits = a[i];
        while (bits != 0) {
            acc ^= b[__builtin_ctzll(bits)];
            bits &= bits - 1;
        }
        r[i] = acc;
    }
    memcpy(c, r, sizeof(r));
}
//
// c = a^T.b for N x 64 blocks a and b: a byte of a[k] at a time, eight
// tables of 256 sums of rows of b, then the rows of c from the tables.
//
static void mul_trans_Nx64(const uint64_t* const a, const uint64_t* const b,
                           const uint32_t n, uint64_t* const c) {
    uint64_t (*const t)[256] = tifa_calloc(8, sizeof(*t));
    for (uint32_t k = 0; k < n; k++) {
        const uint64_t x = a[k];
        const uint64_t y = b[k];
        t[0][(x >>  0) & 0xFF] ^= y;
        t[1][(x >>  8) & 0xFF] ^= y;
        t[2][(x >> 16) & 0xFF] ^= y;
        t[3][(x >> 24) & 0xFF] ^= y;
        t[4][(x >> 32) & 0xFF] ^= y;
        t[5][(x >> 40) & 0xFF] ^= y;
        t[6][(x >> 48) & 0xFF] ^= y;
        t[7][(x >> 56) & 0xFF] ^= y;
    }
    for (uint32_t i = 0; i < 8; i++) {
        for (uint32_t j = 0; j < 8; j++) {
            uint64_t acc = 0;
            for (uint32_t v = 0; v < 256; v++) {
                if (v & BIT(j)) {
                    acc ^= t[i][v];
                }
            }
            c[8 * i + j] = acc;
        }
    }
    free(t);
}
//
// y ^= a.m for an N x 64 block a and a 64 x 64 matrix m: eight tables of
// the 256 sums of eight rows of m.
//
static void mul_Nx64_64x64_acc(const uint64_t* const a,
                               const uint64_t* const m, uint64_t* const y,
                               const uint32_t n) {
    uint64_t (*const t)[256] = tifa_malloc(8 * sizeof(*t));
    for (uint32_t i = 0; i < 8; i++) {
        t[i][0] = 0;
        for (uint32_t v = 1; v < 256; v++) {
            t[i][v] = t[i][v & (v - 1)] ^ m[8 * i + __builtin_ctz(v)];
        }
    }
    for (uint32_t k = 0; k < n; k++) {
        const uint64_t x = a[k];
        y[k] ^= t[0][(x >>  0) & 0xFF] ^ t[1][(x >>  8) & 0xFF]
              ^ t[2][(x >> 16) & 0xFF] ^ t[3][(x >> 24) & 0xFF]
              ^ t[4][(x >> 32) & 0xFF] ^ t[5][(x >> 40) & 0xFF]
              ^ t[6][(x >> 48) & 0xFF] ^ t[7][(x >> 56) & 0xFF];
    }
    free(t);
}
//
// w = B.v (one word per row), then y = B^T.w (one word per column).
//
static void mul_B(const sparse_gf2_t* const m, const uint64_t* const v,
                  uint64_t* const w) {
    memset(w, 0, (size_t) m->nrows * sizeof(uint64_t));
    for (uint32_t c = 0; c < m->ncols; c++) {
        const uint64_t x = v[c];
        for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
            w[m->rows[k]] ^= x;
        }
    }
}

static void mul_BtB(const sparse_gf2_t* const m, const uint64_t* const v,
                    uint64_t* const w, uint64_t* const y) {
    mul_B(m, v, w);
    for (uint32_t c = 0; c < m->ncols; c++) {
        uint64_t acc = 0;
        for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
            acc ^= w[m->rows[k]];
        }
        y[c] = acc;
    }
}
//
// The columns s of t to keep, and w the inverse of t restricted to them:
// Montgomery's choice, which keeps every column not in the last choice
// (last_s), then as many others as leave the restriction invertible. 0 if
// none fits (the iteration has broken down).
//
static uint32_t find_nonsingular_sub(const uint64_t* const t,
                                     uint32_t* const s,
                                     const uint32_t* const last_s,
                                     const uint32_t last_dim,
                                     uint64_t* const w) {
    uint64_t M[64][2];
    for (uint32_t i = 0; i < 64; i++) {
        M[i][0] = t[i];
        M[i][1] = BIT(i);
    }
    uint64_t mask = 0;
    for (uint32_t i = 0; i < last_dim; i++) {
        mask |= BIT(last_s[i]);
    }
    uint32_t order[64];
    uint32_t n = 0;
    for (uint32_t i = 0; i < 64; i++) {
        if ((mask & BIT(i)) == 0) {
            order[n++] = i;
        }
    }
    for (uint32_t i = 0; i < last_dim; i++) {
        order[n++] = last_s[i];
    }
    uint32_t dim = 0;
    for (uint32_t i = 0; i < 64; i++) {
        const uint64_t bit = BIT(order[i]);
        uint64_t* const row_i = M[order[i]];
        uint32_t j;
        for (j = i; j < 64; j++) {
            uint64_t* const row_j = M[order[j]];
            if (row_j[0] & bit) {
                const uint64_t m0 = row_j[0];
                const uint64_t m1 = row_j[1];
                row_j[0] = row_i[0];
                row_j[1] = row_i[1];
                row_i[0] = m0;
                row_i[1] = m1;
                break;
            }
        }
        if (j < 64) {
            for (j = 0; j < 64; j++) {
                uint64_t* const row_j = M[order[j]];
                if ((row_j != row_i) && (row_j[0] & bit)) {
                    row_j[0] ^= row_i[0];
                    row_j[1] ^= row_i[1];
                }
            }
            s[dim++] = order[i];
            continue;
        }
        for (j = i; j < 64; j++) {
            uint64_t* const row_j = M[order[j]];
            if (row_j[1] & bit) {
                const uint64_t m0 = row_j[0];
                const uint64_t m1 = row_j[1];
                row_j[0] = row_i[0];
                row_j[1] = row_i[1];
                row_i[0] = m0;
                row_i[1] = m1;
                break;
            }
        }
        if (j == 64) {
            return 0;
        }
        for (j = 0; j < 64; j++) {
            uint64_t* const row_j = M[order[j]];
            if ((row_j != row_i) && (row_j[1] & bit)) {
                row_j[0] ^= row_i[0];
                row_j[1] ^= row_i[1];
            }
        }
        row_i[0] = 0;
        row_i[1] = 0;
    }
    for (uint32_t i = 0; i < 64; i++) {
        w[i] = M[i][1];
    }
    //
    // The recurrence needs every column of t not chosen last time.
    //
    mask = 0;
    for (uint32_t i = 0; i < dim; i++) {
        mask |= BIT(s[i]);
    }
    for (uint32_t i = 0; i < last_dim; i++) {
        mask |= BIT(last_s[i]);
    }
    return (mask == UINT64_MAX) ? dim : 0;
}

//------------------------------------------------------------------------------
//                               THE ITERATION
//------------------------------------------------------------------------------
//
// x and v, N x 64 blocks, such that the dependencies are in the span of
// their 128 columns; false if the iteration broke down or was stopped.
//
static bool lanczos_core(const sparse_gf2_t* const m, uint64_t seed,
                         uint64_t* const x, uint64_t* const v_last) {
    const uint32_t n = m->ncols;
    uint64_t* const y       = tifa_malloc((size_t) n * sizeof(uint64_t));
    uint64_t* const v0_init = tifa_malloc((size_t) n * sizeof(uint64_t));
    uint64_t* v[3];
    for (uint32_t i = 0; i < 3; i++) {
        v[i] = tifa_calloc(n, sizeof(uint64_t));
    }
    uint64_t* vnext   = tifa_malloc((size_t) n * sizeof(uint64_t));
    uint64_t* const w = tifa_malloc(((size_t) m->nrows + 1)
                                    * sizeof(uint64_t));
    //
    // xorshift64*, for the random start.
    //
    for (uint32_t k = 0; k < n; k++) {
        seed ^= seed >> 12;
        seed ^= seed << 25;
        seed ^= seed >> 27;
        y[k] = seed * UINT64_C(2685821657736338717);
    }
    mul_BtB(m, y, w, v[0]);
    memcpy(v0_init, v[0], (size_t) n * sizeof(uint64_t));
    memset(x, 0, (size_t) n * sizeof(uint64_t));

    uint64_t vt_a_v[2][64], vt_a2_v[2][64], winv[3][64];
    uint64_t d[64], e[64], f[64], f2[64];
    uint32_t s[2][64];
    uint32_t dim1 = 64;
    uint64_t mask1 = UINT64_MAX;
    memset(vt_a_v[1], 0, sizeof(vt_a_v[1]));
    memset(vt_a2_v[1], 0, sizeof(vt_a2_v[1]));
    memset(winv, 0, sizeof(winv));
    for (uint32_t i = 0; i < 64; i++) {
        s[1][i] = i;
    }
    //
    // Each iteration takes up to 64 dimensions out of n; a few more than
    // n / 63 means a breakdown.
    //
    const uint32_t max_iter = n / 60 + 20;
    bool ok = false;

    for (uint32_t iter = 0; iter < max_iter; iter++) {
        if (tifa_should_stop()) {
            break;
        }
        mul_BtB(m, v[0], w, vnext);
        mul_trans_Nx64(v[0], vnext, n, vt_a_v[0]);
        mul_trans_Nx64(vnext, vnext, n, vt_a2_v[0]);

        uint64_t any = 0;
        for (uint32_t i = 0; i < 64; i++) {
            any |= vt_a_v[0][i];
        }
        if (any == 0) {
            ok = true;
            break;
        }
        const uint32_t dim0 = find_nonsingular_sub(vt_a_v[0], s[0], s[1],
                                                   dim1, winv[0]);
        if (dim0 == 0) {
            break;
        }
        uint64_t mask0 = 0;
        for (uint32_t i = 0; i < dim0; i++) {
            mask0 |= BIT(s[0][i]);
        }
        //
        // D = I + Winv0.(V0^T.A^2.V0.S0.S0^T + V0^T.A.V0)
        //
        for (uint32_t i = 0; i < 64; i++) {
            d[i] = (vt_a2_v[0][i] & mask0) ^ vt_a_v[0][i];
        }
        mul_64x64(winv[0], d, d);
        for (uint32_t i = 0; i < 64; i++) {
            d[i] ^= BIT(i);
        }
        //
        // E = Winv1.V0^T.A.V0.S0.S0^T
        //
        mul_64x64(winv[1], vt_a_v[0], e);
        for (uint32_t i = 0; i < 64; i++) {
            e[i] &= mask0;
        }
        //
        // F = Winv2.(I + V1^T.A.V1.Winv1).(V1^T.A^2.V1.S1.S1^T
        //     + V1^T.A.V1).S0.S0^T
        //
        mul_64x64(vt_a_v[1], winv[1], f);
        for (uint32_t i = 0; i < 64; i++) {
            f[i] ^= BIT(i);
        }
        mul_64x64(winv[2], f, f);
        for (uint32_t i = 0; i < 64; i++) {
            f2[i] = ((vt_a2_v[1][i] & mask1) ^ vt_a_v[1][i]) & mask0;
        }
        mul_64x64(f, f2, f);
        //
        // V_next = A.V0.S0.S0^T + V0.D + V1.E + V2.F
        //
        for (uint32_t k = 0; k < n; k++) {
            vnext[k] &= mask0;
        }
        mul_Nx64_64x64_acc(v[0], d, vnext, n);
        mul_Nx64_64x64_acc(v[1], e, vnext, n);
        mul_Nx64_64x64_acc(v[2], f, vnext, n);
        //
        // X += V0.Winv0.V0^T.V_init
        //
        mul_trans_Nx64(v[0], v0_init, n, d);
        mul_64x64(winv[0], d, d);
        mul_Nx64_64x64_acc(v[0], d, x, n);

        uint64_t* const old = v[2];
        v[2] = v[1];
        v[1] = v[0];
        v[0] = vnext;
        vnext = old;
        memcpy(winv[2], winv[1], sizeof(winv[1]));
        memcpy(winv[1], winv[0], sizeof(winv[0]));
        memcpy(vt_a_v[1], vt_a_v[0], sizeof(vt_a_v[0]));
        memcpy(vt_a2_v[1], vt_a2_v[0], sizeof(vt_a2_v[0]));
        memcpy(s[1], s[0], sizeof(s[0]));
        mask1 = mask0;
        dim1  = dim0;
    }
    if (ok) {
        //
        // A.X = A.Y, so that A.(X + Y) = 0 up to V_last.
        //
        for (uint32_t k = 0; k < n; k++) {
            x[k] ^= y[k];
        }
        memcpy(v_last, v[0], (size_t) n * sizeof(uint64_t));
    }
    free(y);
    free(v0_init);
    for (uint32_t i = 0; i < 3; i++) {
        free(v[i]);
    }
    free(vnext);
    free(w);
    return ok;
}
//
// The combinations of the 128 columns of x and v that B maps to zero, as
// up to 64 vectors deps (N x 64, bit j of deps[k] whether column k is in
// dependency j); returns the mask of the nonzero ones.
//
static uint64_t combine(const sparse_gf2_t* const m, const uint64_t* const x,
                        const uint64_t* const v, uint64_t* const deps) {
    const uint32_t n      = m->ncols;
    const uint32_t nrows  = m->nrows;
    const uint32_t nwords = (nrows + 63) / 64;

    uint64_t* const bx = tifa_malloc(((size_t) nrows + 1) * sizeof(uint64_t));
    uint64_t* const bv = tifa_malloc(((size_t) nrows + 1) * sizeof(uint64_t));
    mul_B(m, x, bx);
    mul_B(m, v, bv);
    //
    // The 128 columns of [B.x | B.v] as rows of nrows bits, each with the
    // combination of x and v it stands for (two words).
    //
    uint64_t* const rows = tifa_calloc((size_t) 128 * nwords,
                                       sizeof(uint64_t));
    uint64_t tags[128][2];
    for (uint32_t r = 0; r < nrows; r++) {
        uint64_t bits = bx[r];
        while (bits != 0) {
            const uint32_t j = (uint32_t) __builtin_ctzll(bits);
            rows[(size_t) j * nwords + r / 64] |= BIT(r % 64);
            bits &= bits - 1;
        }
        bits = bv[r];
        while (bits != 0) {
            const uint32_t j = (uint32_t) __builtin_ctzll(bits);
            rows[(size_t) (64 + j) * nwords + r / 64] |= BIT(r % 64);
            bits &= bits - 1;
        }
    }
    for (uint32_t j = 0; j < 128; j++) {
        tags[j][0] = (j < 64) ? BIT(j) : 0;
        tags[j][1] = (j < 64) ? 0 : BIT(j - 64);
    }
    //
    // Elimination: the rows reduced to zero give the combinations.
    //
    uint32_t pivot_word[128];
    uint64_t pivot_bit[128];
    uint32_t pivots[128];
    uint32_t npivots = 0;
    uint64_t combos[64][2];
    uint32_t ncombos = 0;

    for (uint32_t j = 0; j < 128; j++) {
        uint64_t* const row = rows + (size_t) j * nwords;
        for (uint32_t p = 0; p < npivots; p++) {
            if (row[pivot_word[p]] & pivot_bit[p]) {
                const uint64_t* const prow = rows + (size_t) pivots[p] * nwords;
                for (uint32_t k = 0; k < nwords; k++) {
                    row[k] ^= prow[k];
                }
                tags[j][0] ^= tags[pivots[p]][0];
                tags[j][1] ^= tags[pivots[p]][1];
            }
        }
        uint32_t k = 0;
        while ((k < nwords) && (row[k] == 0)) {
            k++;
        }
        if (k < nwords) {
            pivots[npivots]     = j;
            pivot_word[npivots] = k;
            pivot_bit[npivots]  = row[k] & (~row[k] + 1);
            npivots++;
        } else if (ncombos < 64) {
            combos[ncombos][0] = tags[j][0];
            combos[ncombos][1] = tags[j][1];
            ncombos++;
        }
    }
    free(rows);
    free(bx);
    free(bv);
    //
    // deps = x.Tx + v.Tv, column j of Tx and Tv the j-th combination.
    //
    uint64_t tx[64], tv[64];
    memset(tx, 0, sizeof(tx));
    memset(tv, 0, sizeof(tv));
    for (uint32_t j = 0; j < ncombos; j++) {
        for (uint32_t i = 0; i < 64; i++) {
            if (combos[j][0] & BIT(i)) {
                tx[i] |= BIT(j);
            }
            if (combos[j][1] & BIT(i)) {
                tv[i] |= BIT(j);
            }
        }
    }
    memset(deps, 0, (size_t) n * sizeof(uint64_t));
    mul_Nx64_64x64_acc(x, tx, deps, n);
    mul_Nx64_64x64_acc(v, tv, deps, n);

    uint64_t nonzero = 0;
    for (uint32_t k = 0; k < n; k++) {
        nonzero |= deps[k];
    }
    return nonzero;
}

//------------------------------------------------------------------------------
//                                 FILTERING
//------------------------------------------------------------------------------
//
// Columns beyond this many over the rows they use are left out, the
// heaviest first: block Lanczos needs some excess (more relations than
// primes), and finds at most 64 dependencies anyway.
//
#define BL_EXCESS 100
//
// The columns of m kept (their indices in m), in a new matrix whose rows
// are renumbered to those still used. NULL if none are left.
//
static sparse_gf2_t* filter(const sparse_gf2_t* const m,
                            uint32_t** const kept_out) {
    const uint32_t ncols = m->ncols;
    uint32_t* const weight = tifa_calloc((size_t) m->nrows + 1,
                                         sizeof(uint32_t));
    bool* const dead = tifa_calloc((size_t) ncols + 1, sizeof(bool));
    for (uint32_t k = 0; k < m->start[ncols]; k++) {
        weight[m->rows[k]]++;
    }
    //
    // A column with a row of its own is in no dependency: those go, and so
    // do the columns this leaves with a row of their own, and so on. The
    // rows of weight one, found again after each pass.
    //
    uint32_t nalive = ncols;
    for (;;) {
        bool changed = false;
        for (uint32_t c = 0; c < ncols; c++) {
            if (dead[c]) {
                continue;
            }
            for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
                if (weight[m->rows[k]] == 1) {
                    dead[c] = true;
                    break;
                }
            }
            if (dead[c]) {
                for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
                    weight[m->rows[k]]--;
                }
                nalive--;
                changed = true;
            }
        }
        uint32_t nused = 0;
        for (uint32_t r = 0; r < m->nrows; r++) {
            nused += (weight[r] != 0);
        }
        if (!changed && (nalive <= nused + BL_EXCESS)) {
            break;
        }
        if (nalive > nused + BL_EXCESS) {
            //
            // The heaviest columns go, down to the excess wanted, then the
            // singletons again.
            //
            uint32_t extra = nalive - (nused + BL_EXCESS);
            uint32_t maxw = 0;
            for (uint32_t c = 0; c < ncols; c++) {
                const uint32_t cw = m->start[c + 1] - m->start[c];
                if (!dead[c] && (cw > maxw)) {
                    maxw = cw;
                }
            }
            for (uint32_t target = maxw; (target > 0) && (extra > 0);
                 target--) {
                for (uint32_t c = 0; (c < ncols) && (extra > 0); c++) {
                    if (dead[c] || (m->start[c + 1] - m->start[c] != target)) {
                        continue;
                    }
                    dead[c] = true;
                    for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
                        weight[m->rows[k]]--;
                    }
                    nalive--;
                    extra--;
                }
            }
        }
    }
    sparse_gf2_t* f = NULL;
    uint32_t* kept = NULL;
    if (nalive != 0) {
        uint32_t* const renum = tifa_malloc(((size_t) m->nrows + 1)
                                            * sizeof(uint32_t));
        uint32_t nused = 0;
        for (uint32_t r = 0; r < m->nrows; r++) {
            renum[r] = (weight[r] != 0) ? nused++ : UINT32_MAX;
        }
        f    = alloc_sparse_gf2(nused, nalive);
        kept = tifa_malloc((size_t) nalive * sizeof(uint32_t));
        uint32_t* const col = tifa_malloc(((size_t) m->nrows + 1)
                                          * sizeof(uint32_t));
        for (uint32_t c = 0; c < ncols; c++) {
            if (dead[c]) {
                continue;
            }
            uint32_t len = 0;
            for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
                col[len++] = renum[m->rows[k]];
            }
            kept[f->ncols] = c;
            append_sparse_gf2_column(f, col, len);
        }
        free(col);
        free(renum);
    }
    free(weight);
    free(dead);
    *kept_out = kept;
    return f;
}

//------------------------------------------------------------------------------
//                                THE FUNCTION
//------------------------------------------------------------------------------
//
// Tries with these many random starts before giving up.
//
#define BL_ATTEMPTS 3

uint32_t block_lanczos(const sparse_gf2_t* const m,
                       uint32_array_list_t* const deps) {
    uint32_t* kept = NULL;
    sparse_gf2_t* const f = filter(m, &kept);
    if (f == NULL) {
        return 0;
    }
    const uint32_t n = f->ncols;
    uint64_t* const x    = tifa_malloc((size_t) n * sizeof(uint64_t));
    uint64_t* const v    = tifa_malloc((size_t) n * sizeof(uint64_t));
    uint64_t* const dep  = tifa_malloc((size_t) n * sizeof(uint64_t));
    uint64_t* const w    = tifa_malloc(((size_t) f->nrows + 1)
                                       * sizeof(uint64_t));
    uint32_t nfound = 0;

    for (uint32_t attempt = 0; (attempt < BL_ATTEMPTS) && (nfound == 0);
         attempt++) {
        const uint64_t seed = UINT64_C(0x9E3779B97F4A7C15) * (attempt + 1)
                            ^ ((uint64_t) n << 32) ^ f->nrows;
        if (!lanczos_core(f, seed, x, v)) {
            if (tifa_should_stop()) {
                break;
            }
            continue;
        }
        uint64_t nonzero = combine(f, x, v, dep);
        //
        // Checked: B.dep must be zero for each dependency kept.
        //
        mul_B(f, dep, w);
        for (uint32_t r = 0; r < f->nrows; r++) {
            nonzero &= ~w[r];
        }
        for (uint32_t j = 0; (j < 64) && (deps->length < deps->alloced);
             j++) {
            if ((nonzero & BIT(j)) == 0) {
                continue;
            }
            uint32_t len = 0;
            for (uint32_t k = 0; k < n; k++) {
                len += (uint32_t) ((dep[k] >> j) & 1);
            }
            uint32_array_t* const rel = alloc_uint32_array(len);
            for (uint32_t k = 0; k < n; k++) {
                if ((dep[k] >> j) & 1) {
                    rel->data[rel->length++] = kept[k];
                }
            }
            add_entry_in_uint32_array_list(rel, deps);
            nfound++;
        }
    }
    free(x);
    free(v);
    free(dep);
    free(w);
    free(kept);
    free_sparse_gf2(f);
    return nfound;
}
