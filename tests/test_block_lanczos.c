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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "tifa_test.h"
#include "block_lanczos.h"
#include "lindep.h"

static uint64_t rng = 88172645463325252ULL;

static uint32_t next_rand(void) {
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (uint32_t) (rng >> 32);
}
//
// Columns as a sieve makes them: row r taken with a chance of about
// 1 / (r + 2), so that the first rows are heavy and the last ones light,
// a dozen ones each.
//
static sparse_gf2_t* random_relations(const uint32_t nrows,
                                      const uint32_t ncols) {
    sparse_gf2_t* const m = alloc_sparse_gf2(nrows, ncols);
    uint32_t* const col  = malloc((nrows + 1) * sizeof(uint32_t));
    uint8_t*  const seen = calloc(nrows, 1);
    for (uint32_t c = 0; c < ncols; c++) {
        uint32_t n = 0;
        const uint32_t want = 8 + next_rand() % 10;
        while (n < want) {
            //
            // 1 / (r + 2) by inverse transform: r = 2^u - 2, u uniform.
            //
            const double u = (double) next_rand() / 4294967296.0;
            const double lg = log2((double) nrows + 2);
            uint32_t r = (uint32_t) (exp2(u * lg) - 2);
            if (r >= nrows) {
                r = nrows - 1;
            }
            if (!seen[r]) {
                seen[r] = 1;
                col[n++] = r;
            }
        }
        for (uint32_t k = 0; k < n; k++) {
            seen[col[k]] = 0;
        }
        append_sparse_gf2_column(m, col, n);
    }
    free(col);
    free(seen);
    return m;
}
//
// Whether every dependency sums to zero, uses distinct columns of m, and
// is not empty.
//
static bool all_sum_to_zero(const sparse_gf2_t* const m,
                            const uint32_array_list_t* const deps) {
    uint8_t* const sum  = calloc(m->nrows, 1);
    uint8_t* const used = calloc(m->ncols, 1);
    bool ok = true;
    for (uint32_t d = 0; ok && (d < deps->length); d++) {
        const uint32_array_t* const dep = deps->data[d];
        if (dep->length == 0) {
            ok = false;
        }
        for (uint32_t i = 0; ok && (i < dep->length); i++) {
            const uint32_t c = dep->data[i];
            if ((c >= m->ncols) || used[c]) {
                ok = false;
                break;
            }
            used[c] = 1;
            for (uint32_t k = m->start[c]; k < m->start[c + 1]; k++) {
                sum[m->rows[k]] ^= 1;
            }
        }
        for (uint32_t r = 0; r < m->nrows; r++) {
            ok = ok && (sum[r] == 0);
            sum[r] = 0;
        }
        for (uint32_t i = 0; i < dep->length; i++) {
            if (dep->data[i] < m->ncols) {
                used[dep->data[i]] = 0;
            }
        }
    }
    free(sum);
    free(used);
    return ok;
}

static void check_random(const uint32_t nrows, const uint32_t ncols,
                         const uint32_t min_deps) {
    sparse_gf2_t* const m = random_relations(nrows, ncols);
    uint32_array_list_t* const deps = alloc_uint32_array_list(64);
    const uint32_t n = block_lanczos(m, deps);
    char msg[128];
    snprintf(msg, sizeof(msg), "%u x %u: %u dependencies, at least %u",
             nrows, ncols, n, min_deps);
    CHECK((n == deps->length) && (n >= min_deps), msg);
    snprintf(msg, sizeof(msg), "%u x %u: every dependency sums to zero",
             nrows, ncols);
    CHECK(all_sum_to_zero(m, deps), msg);
    free_uint32_array_list(deps);
    free_sparse_gf2(m);
}

int main(void) {
    TIFA_TEST_BEGIN("block Lanczos");

    check_random(300, 400, 32);
    check_random(2000, 2100, 32);
    check_random(5000, 5040, 20);
    check_random(30000, 30200, 32);
    //
    // Small matrices, where a breakdown near the end lost every dependency
    // of three attempts in a few percent of them.
    //
    for (int i = 0; i < 40; i++) {
        check_random(130, 194, 32);
    }
    //
    // More than 100 empty columns (relations that are squares): the
    // filter left them, and the excess above its bound, forever.
    //
    {
        sparse_gf2_t* const m = alloc_sparse_gf2(10, 106);
        uint32_t rows[2] = { 0, 1 };
        for (uint32_t c = 0; c < 101; c++) {
            append_sparse_gf2_column(m, rows, 0);
        }
        for (uint32_t c = 0; c < 5; c++) {
            append_sparse_gf2_column(m, rows, 2);
        }
        uint32_array_list_t* const deps = alloc_uint32_array_list(64);
        CHECK(block_lanczos(m, deps) >= 32, "101 empty columns");
        CHECK(all_sum_to_zero(m, deps), "their dependencies sum to zero");
        free_uint32_array_list(deps);
        free_sparse_gf2(m);
    }
    //
    // Fewer columns than rows: dependencies only if some rows are unused
    // or the columns happen to be dependent; none must be wrong.
    //
    {
        sparse_gf2_t* const m = random_relations(1000, 600);
        uint32_array_list_t* const deps = alloc_uint32_array_list(64);
        block_lanczos(m, deps);
        CHECK(all_sum_to_zero(m, deps), "600 columns over 1000 rows: right");
        free_uint32_array_list(deps);
        free_sparse_gf2(m);
    }
    //
    // The room left in the list bounds what is appended.
    //
    {
        sparse_gf2_t* const m = random_relations(500, 700);
        uint32_array_list_t* const deps = alloc_uint32_array_list(3);
        CHECK(block_lanczos(m, deps) == 3, "no more than the room left");
        CHECK(all_sum_to_zero(m, deps), "the three sum to zero");
        free_uint32_array_list(deps);
        free_sparse_gf2(m);
    }
    //
    // From a dense matrix whose rows are the relations, as the
    // factorization functions build them.
    //
    {
        const uint32_t nrel = 900, nprimes = 800;
        binary_matrix_t* const dense = alloc_binary_matrix(nrel, nprimes);
        sparse_gf2_t* const ref = random_relations(nprimes, nrel);
        for (uint32_t c = 0; c < ref->ncols; c++) {
            for (uint32_t k = ref->start[c]; k < ref->start[c + 1]; k++) {
                flip_matrix_bit(c, ref->rows[k], dense);
            }
        }
        sparse_gf2_t* const m = sparse_gf2_from_relations(dense);
        bool same = (m->ncols == ref->ncols) && (m->nrows == ref->nrows);
        uint32_array_list_t* const deps = alloc_uint32_array_list(64);
        CHECK(same && (block_lanczos(m, deps) >= 32), "from a dense matrix");
        CHECK(all_sum_to_zero(ref, deps), "its dependencies sum to zero");
        uint32_array_list_t* const found = find_dependencies(dense,
                                                             BLOCK_LANCZOS);
        CHECK((found->length > 0) && all_sum_to_zero(ref, found),
              "find_dependencies with BLOCK_LANCZOS");
        free_uint32_array_list(found);
        free_uint32_array_list(deps);
        free_sparse_gf2(m);
        free_sparse_gf2(ref);
        free_binary_matrix(dense);
    }

    TIFA_TEST_END();
}
