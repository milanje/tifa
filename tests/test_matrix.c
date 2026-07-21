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
