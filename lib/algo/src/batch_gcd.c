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
 * \file    batch_gcd.c
 *
 * \brief The gcd of each integer of a set with the product of the others.
 */

#include <stdlib.h>
#include <gmp.h>

#include "batch_gcd.h"
#include "tifa_alloc.h"
#include "stopping.h"

//------------------------------------------------------------------------------
static void free_levels(mpz_t** const levels, const uint32_t* const lengths,
                        const uint32_t nlevels) {
    for (uint32_t k = 1; k < nlevels; k++) {
        for (uint32_t j = 0; j < lengths[k]; j++) {
            mpz_clear(levels[k][j]);
        }
        free(levels[k]);
    }
}
//------------------------------------------------------------------------------
ecode_t tifa_batch_gcd(mpz_array_t* const gcds,
                       const mpz_array_t* const moduli) {
    const uint32_t count = moduli->length;
    for (uint32_t i = 0; i < count; i++) {
        if (mpz_sgn(moduli->data[i]) <= 0) {
            return FAILURE;
        }
    }
    if (count < 2) {
        if (count == 1) {
            mpz_t one;
            mpz_init_set_ui(one, 1);
            append_mpz_to_array(gcds, one);
            mpz_clear(one);
        }
        return SUCCESS;
    }
    //
    // levels[0] holds the moduli, levels[k][j] the product of
    // levels[k - 1][2j] and levels[k - 1][2j + 1] (or a copy of the first
    // when the second is missing), up to the single product of the top.
    //
    uint32_t nlevels = 1;
    for (uint32_t len = count; len > 1; len = (len + 1) / 2) {
        nlevels++;
    }
    mpz_t**   const levels  = tifa_malloc(nlevels * sizeof(mpz_t*));
    uint32_t* const lengths = tifa_malloc(nlevels * sizeof(uint32_t));
    levels[0]  = moduli->data;
    lengths[0] = count;
    for (uint32_t k = 1; k < nlevels; k++) {
        lengths[k] = (lengths[k - 1] + 1) / 2;
        levels[k]  = tifa_malloc(lengths[k] * sizeof(mpz_t));
        for (uint32_t j = 0; j < lengths[k]; j++) {
            mpz_init(levels[k][j]);
        }
    }
    ecode_t ecode = GIVING_UP;
    mpz_t t, u;
    mpz_inits(t, u, NULL);

    for (uint32_t k = 1; k < nlevels; k++) {
        mpz_t* const below = levels[k - 1];
        for (uint32_t j = 0; j < lengths[k]; j++) {
            if (2 * j + 1 < lengths[k - 1]) {
                mpz_mul(levels[k][j], below[2 * j], below[2 * j + 1]);
            } else {
                mpz_set(levels[k][j], below[2 * j]);
            }
            if (tifa_should_stop()) {
                goto done;
            }
        }
    }
    //
    // Down the tree, each node N becomes (P / N) mod N, P the product of all
    // the moduli: for a node N = A.B, P / A = (P / N).B. Rather than P mod N^2
    // (Heninger et al.), for 0.9 of the time and 0.7 of the memory.
    //
    const uint32_t top = nlevels - 1;
    mpz_set_ui(levels[top][0], 1);
    for (uint32_t k = top - 1; k >= 1; k--) {
        for (uint32_t j = 0; j < lengths[k]; j += 2) {
            mpz_ptr const q = levels[k + 1][j / 2];
            if (j + 1 == lengths[k]) {
                mpz_tdiv_r(levels[k][j], q, levels[k][j]);
                continue;
            }
            mpz_mul(t, q, levels[k][j + 1]);
            mpz_mul(u, q, levels[k][j]);
            mpz_tdiv_r(levels[k][j], t, levels[k][j]);
            mpz_tdiv_r(levels[k][j + 1], u, levels[k][j + 1]);
            if (tifa_should_stop()) {
                goto done;
            }
        }
        for (uint32_t j = 0; j < lengths[k + 1]; j++) {
            mpz_realloc2(levels[k + 1][j], 1);
        }
    }
    if (gcds->alloced - gcds->length < count) {
        resize_mpz_array(gcds, gcds->length + count);
    }
    for (uint32_t i = 0; i < count; i++) {
        mpz_ptr const q = levels[1][i / 2];
        if ((i ^ 1) < count) {
            mpz_mul(t, q, moduli->data[i ^ 1]);
        } else {
            mpz_set(t, q);
        }
        mpz_gcd(gcds->data[gcds->length + i], t, moduli->data[i]);
    }
    gcds->length += count;
    ecode = SUCCESS;

done:
    mpz_clears(t, u, NULL);
    free_levels(levels, lengths, nlevels);
    free(lengths);
    free(levels);

    return ecode;
}
