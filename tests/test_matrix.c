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

// Tests for matrix cloning in lib/utils/src/matrix.c.
//
// test_byte_matrix_clone is a regression test for the row-pointer
// under-allocation in clone_byte_matrix (malloc'd with sizeof(unsigned char)
// instead of sizeof(unsigned char*)). Pre-fix this over-writes the heap; it is
// caught here by comparing every cell (and cleanly under -fsanitize=address).

#include "tifa_test.h"
#include "matrix.h"

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

int main(void) {
    TIFA_TEST_BEGIN("matrix clone");
    test_byte_matrix_clone();
    test_binary_matrix_clone();
    TIFA_TEST_END();
}
