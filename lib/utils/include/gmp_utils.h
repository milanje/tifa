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
 * \file    gmp_utils.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Various GMP small utilities.
 *
 * GMP small utilities' definitions should go here.
 */

#if !defined(_TIFA_GMP_UTILS_H_)
   /**
    * \def _TIFA_GMP_UTILS_H_
    * Standard include guard.
    */
#define _TIFA_GMP_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>
#include "hashtable.h"

   /**
    * \struct struct_mpz_pair_t array.h lib/utils/include/gmp_utils.h
    * \brief  A pair of \c mpz_t integers.
    *
    * This very simple structure defines a pair of \c mpz_t integers.
    */
struct struct_mpz_pair_t {
        /**
         * The first \c mpz_t integer of the pair.
         */
    mpz_t x;
        /**
         * The second \c mpz_t integer of the pair.
         */
    mpz_t y;
};

   /**
    * \typedef mpz_pair_t
    * \brief Equivalent to <tt>struct struct_mpz_pair_t</tt>.
    */
typedef struct struct_mpz_pair_t mpz_pair_t;

   /**
    * \brief inits a <tt>mpz_pair_t</tt>.
    *
    * Inits a <tt>mpz_pair_t</tt> by initializing each of its \c mpz_t
    * element.
    *
    * \param[in] pair A pointer to the <tt>mpz_pair_t</tt> to init.
    */
inline static void init_mpz_pair(mpz_pair_t* pair) {
    mpz_init(pair->x);
    mpz_init(pair->y);
}

   /**
    * \brief Clears a <tt>mpz_pair_t</tt>.
    *
    * Clears a <tt>mpz_pair_t</tt>.
    *
    * \param[in] pair A pointer to the <tt>mpz_pair_t</tt> to clear.
    */
inline static void clear_mpz_pair(mpz_pair_t* pair) {
    mpz_clear(pair->x);
    mpz_clear(pair->y);
}

   /**
    * \brief Empties a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s.
    *
    * Empties a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s and clears
    * the memory associated to the keys and their associated data.
    *
    * \param[in] htable A pointer to the <tt>hashtable_t</tt> to empty.
    */
void empty_mpzpair_htable(hashtable_t* const htable);

   /**
    * \brief Clears a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s.
    *
    * Clears a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s. It clears
    * the memory associated to the keys, their associated data and the
    * hashtable itself.
    *
    * \param[in] htable A pointer to the <tt>hashtable_t</tt> to clear.
    */
inline static void free_mpzpair_htable(hashtable_t* htable) {
    empty_mpzpair_htable(htable);
    free_hashtable(htable);
}

  /**
    * \brief Logarithm in base 10 of a multi-precision integer.
    *
    * Returns an <em>crude approximation</em> of the logarithm (in base 10) of
    * a positive multi-precision integer \c n. The approximation is usually
    * valid up to the 4th decimal.
    *
    * \param[in] n A positive multi-precision integer.
    * \returns An approximation of log(<tt>n</tt>) in base 10.
    */
float mpz_log10(mpz_t n);

#ifdef __cplusplus
}
#endif

#endif
