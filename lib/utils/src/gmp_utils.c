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
 * \file    gmp_utils.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

#include "tifa_config.h"
#include "macros.h"
#include "funcs.h"
#include "gmp_utils.h"

//-----------------------------------------------------------------------------
void empty_mpzpair_htable(hashtable_t* const htable) {
    //
    // Empties the hashtable and all of its entries... It's a bit messy
    // so should we...
    //
    // _NOTE_: Consider adding in the hashtable_t structure (and hence
    //         in the linked_list_t structure) a pointer to a function
    //         clearing the memory space used by the entries?
    //
    for (uint32_t i = 0; i < htable->alloced; i++) {
        linked_list_node_t *node = htable->buckets[i].head;
        linked_list_node_t *next = NULL;
        while (NULL != node) {
            hashtable_entry_t* entry = (hashtable_entry_t*)node->data;
            clear_mpz_pair((mpz_pair_t*)entry->data);
            free(entry->data);
            free(entry->key);
            free(entry);
            next = node->next;
            free(node);
            node = next;
        }
        htable->buckets[i].head   = NULL;
        htable->buckets[i].tail   = NULL;
        htable->buckets[i].length = 0;
    }
}
//-----------------------------------------------------------------------------
float mpz_log10(mpz_t n) {
    if (SIZ(n) == 1) {
        return log10(MPZ_LAST_LIMB_VALUE(n));
    }
    //
    // The log approximation is obtained by taking the log10 of n's most   
    // significant limb left-shifted and completed with bits from its second
    // most significant limb for a better approximation.
    //
    mp_limb_t tmp   = MPZ_LAST_LIMB_VALUE(n);
    uint32_t  shift = GMP_NUMB_BITS - ceil_log2_mp_limb(tmp) - 1;
    mp_limb_t tmp2  = (PTR(n)[SIZ(n) - 2] & GMP_NUMB_MASK);

    tmp  <<= shift;
    tmp2 >>= (GMP_NUMB_BITS - shift);    
    tmp   ^= tmp2;
    
    float res = LOG10_GMP_NUMB_MAX * (SIZ(n) - 1);
    
    res += log10((double)(tmp));
    res -= log10((unsigned long int)1 << shift);
    
    return res;
}
//-----------------------------------------------------------------------------
