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

// Regression tests for the container helpers: clearing one bit of a binary
// array or matrix (once cleared the whole word), sorting a byte array (once
// compared four bytes per element), the hashtable's entry count (never
// counted additions), and adding to a full mpz_array_list_t (wrote past it).

#include <stdlib.h>

#include "tifa_test.h"
#include "array.h"
#include "matrix.h"
#include "hashtable.h"
#include "x_array_list.h"
#include "funcs.h"
#include "tifa_alloc.h"

static void test_binary_array_bits(void) {
    static const unsigned int bits[] = { 1, 0, 1, 1, 0, 0, 1, 0, 1 };
    const uint32_t nbits = sizeof(bits) / sizeof(bits[0]);

    binary_array_t* a = alloc_binary_array(nbits);
    for (uint32_t i = 0; i < nbits; i++) {
        append_bit_to_array(a, bits[i]);
    }
    for (uint32_t i = 0; i < nbits; i++) {
        CHECK_EQ_UI(get_array_bit(i, a), bits[i], "appended bit kept");
    }
    set_array_bit_to_zero(2, a);
    CHECK_EQ_UI(get_array_bit(2, a), 0, "cleared bit");
    CHECK_EQ_UI(get_array_bit(0, a), 1, "bit before the cleared one kept");
    CHECK_EQ_UI(get_array_bit(3, a), 1, "bit after the cleared one kept");
    free_binary_array(a);
}

static void test_binary_matrix_bits(void) {
    binary_matrix_t* m = alloc_binary_matrix(2, 70);
    set_matrix_bit_to_one(1, 3, m);
    set_matrix_bit_to_one(1, 5, m);
    set_matrix_bit_to_one(1, 66, m);
    set_matrix_bit_to_zero(1, 5, m);
    CHECK_EQ_UI(get_matrix_bit(1, 5, m), 0, "cleared matrix bit");
    CHECK_EQ_UI(get_matrix_bit(1, 3, m), 1, "same-word matrix bit kept");
    CHECK_EQ_UI(get_matrix_bit(1, 66, m), 1, "next-word matrix bit kept");
    free_binary_matrix(m);
}

static void test_byte_array_sort(void) {
    static const unsigned char in[] = { 200, 3, 255, 0, 17, 3, 128 };
    static const unsigned char out[] = { 0, 3, 3, 17, 128, 200, 255 };
    const uint32_t n = sizeof(in);

    byte_array_t* a = alloc_byte_array(n);
    for (uint32_t i = 0; i < n; i++) {
        append_byte_to_array(a, in[i]);
    }
    qsort_byte_array(a);
    for (uint32_t i = 0; i < n; i++) {
        CHECK_EQ_UI(a->data[i], out[i], "byte array sorted");
    }
    free_byte_array(a);
}

static void test_hashtable_count(void) {
    hashtable_t* h = alloc_init_hashtable(16, uint32_cmp_func, hash_rj_32);
    static uint32_t values[3] = { 10, 20, 30 };
    for (uint32_t i = 0; i < 3; i++) {
        uint32_t* key = tifa_malloc(sizeof(uint32_t));
        *key = i + 1;
        add_entry_in_hashtable(h, key, &values[i]);
    }
    CHECK_EQ_UI(h->nentries, 3, "entries counted on addition");

    const uint32_t key = 2;
    CHECK(remove_entry_in_hashtable(h, &key) == &values[1], "removed entry");
    CHECK_EQ_UI(h->nentries, 2, "entries counted on removal");

    // free_hashtable leaves the keys to their owner: remove the others.
    for (uint32_t k = 1; k <= 3; k += 2) {
        remove_entry_in_hashtable(h, &k);
    }
    free_hashtable(h);
}

static void test_mpz_array_list_full(void) {
    mpz_array_list_t* list = alloc_mpz_array_list(1);
    mpz_array_t* first  = alloc_mpz_array(1);
    mpz_array_t* second = alloc_mpz_array(1);

    add_entry_in_mpz_array_list(first, list);
    add_entry_in_mpz_array_list(second, list);
    CHECK_EQ_UI(list->length, 1, "a full list takes no more entries");
    CHECK(list->data[0] == first, "the first entry stays");

    free_mpz_array(second);
    free_mpz_array_list(list);
}

int main(void) {
    TIFA_TEST_BEGIN("containers");
    test_binary_array_bits();
    test_binary_matrix_bits();
    test_byte_array_sort();
    test_hashtable_count();
    test_mpz_array_list_full();
    TIFA_TEST_END();
}
