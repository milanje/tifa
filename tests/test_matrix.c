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

// Tests for matrix cloning in lib/utils/src/matrix.c, and for the
// dependencies gaussian_elim finds (lib/utils/src/gauss_elim.c).
//
// test_byte_matrix_clone is a regression test for the row-pointer
// under-allocation in clone_byte_matrix (malloc'd with sizeof(unsigned char)
// instead of sizeof(unsigned char*)). Pre-fix this over-writes the heap; it is
// caught here by comparing every cell (and cleanly under -fsanitize=address).

#include <stdlib.h>

#include "tifa_test.h"
#include "matrix.h"
#include "gauss_elim.h"

static void test_byte_matrix_clone(void) {
    const uint32_t nrows = 5;
    const uint32_t ncols = 3;

    byte_matrix_t* m = alloc_byte_matrix(nrows, ncols);
    for (uint32_t i = 0; i < nrows; i++) {
        for (uint32_t j = 0; j < ncols; j++) {
            m->data[i][j] = (unsigned char)(i * ncols + j + 1);
        }
    }

    byte_matrix_t* clone = clone_byte_matrix(m);

    CHECK_EQ_UI(clone->nrows, nrows, "byte clone nrows");
    CHECK_EQ_UI(clone->ncols, ncols, "byte clone ncols");
    for (uint32_t i = 0; i < nrows; i++) {
        for (uint32_t j = 0; j < ncols; j++) {
            CHECK_EQ_UI(clone->data[i][j], m->data[i][j],
                        "byte clone cell matches original");
        }
    }

    free_byte_matrix(clone);
    free_byte_matrix(m);
}

static void test_binary_matrix_clone(void) {
    const uint32_t nrows = 4;
    const uint32_t ncols = 70;  // spans two 64-bit words

    binary_matrix_t* m = alloc_binary_matrix(nrows, ncols);
    set_matrix_bit_to_one(0, 0,  m);
    set_matrix_bit_to_one(1, 63, m);
    set_matrix_bit_to_one(2, 64, m);
    set_matrix_bit_to_one(3, 69, m);

    binary_matrix_t* clone = clone_binary_matrix(m);

    CHECK_EQ_UI(clone->nrows, nrows, "binary clone nrows");
    CHECK_EQ_UI(clone->ncols, ncols, "binary clone ncols");
    for (uint32_t i = 0; i < nrows; i++) {
        for (uint32_t j = 0; j < ncols; j++) {
            CHECK_EQ_UI(get_matrix_bit(i, j, clone), get_matrix_bit(i, j, m),
                        "binary clone bit matches original");
        }
    }

    free_binary_matrix(clone);
    free_binary_matrix(m);
}

// The rank of a 0/1 matrix, by plain elimination.
static uint32_t rank_of(unsigned char* m, uint32_t nrows, uint32_t ncols) {
    uint32_t rank = 0;
    for (uint32_t c = 0; (c < ncols) && (rank < nrows); c++) {
        uint32_t r = rank;
        while ((r < nrows) && (m[r * ncols + c] == 0)) {
            r++;
        }
        if (r == nrows) {
            continue;
        }
        for (uint32_t k = 0; k < ncols; k++) {
            unsigned char tmp = m[r * ncols + k];
            m[r * ncols + k] = m[rank * ncols + k];
            m[rank * ncols + k] = tmp;
        }
        for (uint32_t i = 0; i < nrows; i++) {
            if ((i != rank) && m[i * ncols + c]) {
                for (uint32_t k = 0; k < ncols; k++) {
                    m[i * ncols + k] ^= m[rank * ncols + k];
                }
            }
        }
        rank++;
    }
    return rank;
}

// Every dependency adds up to zero, and there are as many as the nullity,
// on random matrices shaped like the SIQS and CFRAC ones (a few dense
// columns for the sign and the smallest primes, the others sparse).
static void test_gaussian_elim(void) {
    static const uint32_t shapes[][2] = {
        { 10, 4 }, { 70, 64 }, { 117, 93 }, { 230, 201 }, { 300, 130 },
        { 64, 70 }, { 129, 1 }
    };
    srand(12345);
    for (unsigned s = 0; s < sizeof(shapes) / sizeof(shapes[0]); s++) {
        const uint32_t nrows = shapes[s][0];
        const uint32_t ncols = shapes[s][1];
        unsigned char* dense = calloc((size_t) nrows * ncols, 1);
        binary_matrix_t* m = alloc_binary_matrix(nrows, ncols);
        m->nrows = nrows;
        for (uint32_t i = 0; i < nrows; i++) {
            for (uint32_t c = 0; c < ncols; c++) {
                const int odds = (c < 3) ? 2 : 12;
                if (rand() % odds == 0) {
                    set_matrix_bit_to_one(i, c, m);
                    dense[i * ncols + c] = 1;
                }
            }
        }
        binary_matrix_t* orig = clone_binary_matrix(m);
        uint32_array_list_t* deps = alloc_uint32_array_list(nrows);
        gaussian_elim(deps, m);

        const uint32_t nullity = nrows - rank_of(dense, nrows, ncols);
        CHECK_EQ_UI(deps->length, nullity, "as many dependencies as the nullity");

        unsigned char* seen_last = calloc(nrows, 1);
        for (uint32_t d = 0; d < deps->length; d++) {
            const uint32_array_t* dep = deps->data[d];
            unsigned char* used = calloc(nrows, 1);
            int ok = (dep->length > 0);
            for (uint32_t c = 0; ok && (c < ncols); c++) {
                uint32_t sum = 0;
                for (uint32_t k = 0; k < dep->length; k++) {
                    sum ^= get_matrix_bit(dep->data[k], c, orig);
                }
                ok = (sum == 0);
            }
            for (uint32_t k = 0; ok && (k < dep->length); k++) {
                ok = (dep->data[k] < nrows) && !used[dep->data[k]];
                used[dep->data[k]] = 1;
            }
            CHECK(ok, "a dependency adds up to zero, with distinct rows");
            const uint32_t last = dep->data[dep->length - 1];
            CHECK(!seen_last[last], "dependencies end with distinct free rows");
            seen_last[last] = 1;
            free(used);
        }
        free(seen_last);
        free_uint32_array_list(deps);
        free_binary_matrix(orig);
        free_binary_matrix(m);
        free(dense);
    }
}

int main(void) {
    TIFA_TEST_BEGIN("matrix clone and gaussian elimination");
    test_byte_matrix_clone();
    test_binary_matrix_clone();
    test_gaussian_elim();
    TIFA_TEST_END();
}
