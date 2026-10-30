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
 * \file    lp_table.h
 * \brief   The partial relations of SIQS's large prime variation, by large
 *          prime.
 *
 * An open-addressing table of the large primes; the relations' integers are
 * copied, limbs and signed sizes, into one buffer, and read back as read-only
 * views (mpz_roinit_n): no allocation per relation.
 */

#if !defined(_TIFA_LP_TABLE_H_)
#define _TIFA_LP_TABLE_H_

#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \def LP_TABLE_NONE
    * Returned by \c find_in_lp_table for a large prime not in the table.
    */
#define LP_TABLE_NONE UINT32_MAX

   /**
    * \struct struct_lp_table_t lp_table.h lib/utils/include/lp_table.h
    * \brief  The partial relations of the large prime variation.
    */
struct struct_lp_table_t {
       /**
        * Number of slots, a power of two.
        */
    uint32_t nslots;
       /**
        * 32 minus the base 2 logarithm of \c nslots, for the hash.
        */
    uint32_t shift;
       /**
        * The large prime of the relation stored at each slot (0 for none).
        */
    uint32_t* primes;
       /**
        * The index of the relation stored at each slot.
        */
    uint32_t* index;
       /**
        * Number of relations stored, and room for them.
        */
    uint32_t nrelations;
    uint32_t alloced;
       /**
        * For each relation, the offset of its x in \c limbs (its y follows),
        * and the signed sizes of x and y (as \c mpz_size, negative for a
        * negative integer).
        */
    size_t*  offset;
    int32_t* xsize;
    int32_t* ysize;
       /**
        * The limbs of the relations, and the number used and allocated.
        */
    mp_limb_t* limbs;
    size_t     nlimbs;
    size_t     alloced_limbs;
};

   /**
    * \typedef lp_table_t
    * \brief Equivalent to <tt>struct struct_lp_table_t</tt>.
    */
typedef struct struct_lp_table_t lp_table_t;

   /**
    * \brief Allocates an empty \c lp_table_t.
    *
    * \param[in] nslots The initial number of slots (rounded up to a power of
    *                   two); the table grows as relations are added.
    * \return The table.
    */
lp_table_t* alloc_lp_table(uint32_t nslots);

   /**
    * \brief Frees an \c lp_table_t.
    *
    * \param[in] table The table to free.
    */
void free_lp_table(lp_table_t* const table);

   /**
    * \brief The index of the relation of a large prime, or \c LP_TABLE_NONE.
    *
    * \param[in] table The table.
    * \param[in] prime The large prime (not 0).
    * \return The index of its relation, or \c LP_TABLE_NONE if there is none.
    */
uint32_t find_in_lp_table(const lp_table_t* const table, const uint32_t prime);

   /**
    * \brief Read-only views of the x and y of a stored relation.
    *
    * \c x and \c y need not be initialized and must be neither modified nor
    * cleared; they stay valid until the next relation is added.
    *
    * \param[in] table The table.
    * \param[in] k The index of the relation.
    * \param[out] x Its x.
    * \param[out] y Its y.
    */
void get_from_lp_table(const lp_table_t* const table, const uint32_t k,
                       mpz_t x, mpz_t y);

   /**
    * \brief Stores the relation (x, y) of a large prime not yet in the table.
    *
    * \param[in] table The table.
    * \param[in] prime The large prime (not 0), not in the table.
    * \param[in] x The x of the relation.
    * \param[in] y The y of the relation.
    */
void add_to_lp_table(lp_table_t* const table, const uint32_t prime,
                     const mpz_t x, const mpz_t y);

#ifdef __cplusplus
}
#endif

#endif
