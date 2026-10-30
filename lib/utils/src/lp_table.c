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
 * \file    lp_table.c
 * \brief   The partial relations of SIQS's large prime variation, by large
 *          prime.
 */

#include <stdlib.h>
#include <string.h>

#include "tifa_config.h"
#include "tifa_alloc.h"
#include "lp_table.h"

//------------------------------------------------------------------------------
static inline uint32_t slot_of(const lp_table_t* const table,
                               const uint32_t prime) {
    return (uint32_t) (prime * UINT32_C(2654435761)) >> table->shift;
}
//------------------------------------------------------------------------------
static void alloc_slots(lp_table_t* const table, uint32_t nslots) {
    uint32_t log2 = 4;
    while ((1U << log2) < nslots) {
        log2++;
    }
    table->nslots = 1U << log2;
    table->shift  = 32 - log2;
    table->primes = tifa_calloc(table->nslots, sizeof(uint32_t));
    table->index  = tifa_malloc(table->nslots * sizeof(uint32_t));
}
//------------------------------------------------------------------------------
lp_table_t* alloc_lp_table(uint32_t nslots) {
    lp_table_t* const table = tifa_malloc(sizeof(lp_table_t));
    alloc_slots(table, nslots);
    table->nrelations    = 0;
    table->alloced       = 0;
    table->offset        = NULL;
    table->xsize         = NULL;
    table->ysize         = NULL;
    table->limbs         = NULL;
    table->nlimbs        = 0;
    table->alloced_limbs = 0;
    return table;
}
//------------------------------------------------------------------------------
void free_lp_table(lp_table_t* const table) {
    if (table != NULL) {
        free(table->primes);
        free(table->index);
        free(table->offset);
        free(table->xsize);
        free(table->ysize);
        free(table->limbs);
        free(table);
    }
}
//------------------------------------------------------------------------------
uint32_t find_in_lp_table(const lp_table_t* const table, const uint32_t prime) {
    const uint32_t mask = table->nslots - 1;
    uint32_t slot = slot_of(table, prime);
    while (table->primes[slot] != 0) {
        if (table->primes[slot] == prime) {
            return table->index[slot];
        }
        slot = (slot + 1) & mask;
    }
    return LP_TABLE_NONE;
}
//------------------------------------------------------------------------------
void get_from_lp_table(const lp_table_t* const table, const uint32_t k,
                       mpz_t x, mpz_t y) {
    const mp_limb_t* const xl = table->limbs + table->offset[k];
    const int32_t          xs = table->xsize[k];
    mpz_roinit_n(x, xl, xs);
    mpz_roinit_n(y, xl + abs(xs), table->ysize[k]);
}
//------------------------------------------------------------------------------
static void put_in_slots(lp_table_t* const table, const uint32_t prime,
                         const uint32_t index) {
    const uint32_t mask = table->nslots - 1;
    uint32_t slot = slot_of(table, prime);
    while (table->primes[slot] != 0) {
        slot = (slot + 1) & mask;
    }
    table->primes[slot] = prime;
    table->index[slot]  = index;
}
//------------------------------------------------------------------------------
void add_to_lp_table(lp_table_t* const table, const uint32_t prime,
                     const mpz_t x, const mpz_t y) {
    const uint32_t k = table->nrelations;
    //
    // At most half full: twice the slots, the relations put back in them.
    //
    if (2 * (k + 1) > table->nslots) {
        uint32_t* const primes = table->primes;
        uint32_t* const index  = table->index;
        const uint32_t  nslots = table->nslots;
        alloc_slots(table, 2 * nslots);
        for (uint32_t s = 0; s < nslots; s++) {
            if (primes[s] != 0) {
                put_in_slots(table, primes[s], index[s]);
            }
        }
        free(primes);
        free(index);
    }
    if (k == table->alloced) {
        table->alloced = (k < 256) ? 512 : 2 * k;
        table->offset  = tifa_realloc(table->offset,
                                      table->alloced * sizeof(size_t));
        table->xsize   = tifa_realloc(table->xsize,
                                      table->alloced * sizeof(int32_t));
        table->ysize   = tifa_realloc(table->ysize,
                                      table->alloced * sizeof(int32_t));
    }
    const size_t nx = mpz_size(x);
    const size_t ny = mpz_size(y);
    if (table->nlimbs + nx + ny > table->alloced_limbs) {
        size_t alloced = (table->alloced_limbs < 2048) ? 4096
                                                       : 2 * table->alloced_limbs;
        while (alloced < table->nlimbs + nx + ny) {
            alloced *= 2;
        }
        table->limbs = tifa_realloc(table->limbs, alloced * sizeof(mp_limb_t));
        table->alloced_limbs = alloced;
    }
    mp_limb_t* const dest = table->limbs + table->nlimbs;
    memcpy(dest, mpz_limbs_read(x), nx * sizeof(mp_limb_t));
    memcpy(dest + nx, mpz_limbs_read(y), ny * sizeof(mp_limb_t));

    table->offset[k] = table->nlimbs;
    table->xsize[k]  = (int32_t) ((mpz_sgn(x) < 0) ? -(long) nx : (long) nx);
    table->ysize[k]  = (int32_t) ((mpz_sgn(y) < 0) ? -(long) ny : (long) ny);
    table->nlimbs   += nx + ny;
    table->nrelations = k + 1;

    put_in_slots(table, prime, k);
}
//------------------------------------------------------------------------------
