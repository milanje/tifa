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
 * \file    buckets.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Structure and inline functions to implement bucket sieving.
 */
 
#include <stdlib.h>

#if !defined(_TIFA_BUCKETS_H_)
#define _TIFA_BUCKETS_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
#define A_MASK      0xFFFFFF00
#define B_MASK      0x000000FF
#define GET_A(X)    (((X) & A_MASK) >> 8)
#define GET_B(X)    ((X) & B_MASK)
//-----------------------------------------------------------------------------
struct struct_buckets_t {
    uint16_t         nbins;
    uint32_array_t** bins;
};
//-----------------------------------------------------------------------------
typedef struct struct_buckets_t buckets_t;
//-----------------------------------------------------------------------------
static inline buckets_t*
alloc_buckets(uint16_t balloced, uint8_t nbins) {
    
    buckets_t* buckets = malloc(sizeof(buckets_t));
    
    buckets->nbins = nbins;
    buckets->bins  = malloc(nbins * sizeof(uint32_array_t*));
    
    for (uint8_t i = 0; i < nbins; i++) {
      buckets->bins[i] = alloc_uint32_array(balloced);
    }
    return buckets;
}
//-----------------------------------------------------------------------------
static inline void
free_buckets(buckets_t* const buckets) {
    if (buckets != NULL) {
        for (uint8_t i = 0; i < buckets->nbins; i++) {
            free_uint32_array(buckets->bins[i]);
        }
        free(buckets->bins);
        free(buckets);
    }
}
//-----------------------------------------------------------------------------
static inline void
reset_buckets(buckets_t* const buckets) {
    for (uint8_t i = 0; i < buckets->nbins; i++) {
        buckets->bins[i]->length = 0;
    }
}
//-----------------------------------------------------------------------------
static inline void
add_to_buckets(const uint32_t a, const uint8_t b, const uint8_t ibin,
               buckets_t* const buckets) {

    uint32_array_t* bin  = buckets->bins[ibin];
    uint32_t        blen = bin->length;
    
    if (blen == bin->alloced) {
        bin->data    = realloc(bin->data, 2 * blen * sizeof(uint32_t));
        bin->alloced = 2 * blen;
    }
    bin->data[blen] = (a << 8) ^ (b & 0xFF);
    bin->length++;
}
//-----------------------------------------------------------------------------
static inline void
print_buckets(buckets_t* const buckets) {
    for (uint32_t i = 0; i < buckets->nbins; i++) {
        printf("bin %2i:\n", i);
        printf("-------\n");
        uint32_array_t* bin = buckets->bins[i]; 
        for (uint32_t j = 0; j < bin->length; j++) {
            printf("  (%5i, %2i)\n", GET_A(bin->data[j]), GET_B(bin->data[j]));
        }
    }
}
//-----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif

