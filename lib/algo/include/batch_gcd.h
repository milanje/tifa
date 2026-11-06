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
 * \file    batch_gcd.h
 *
 * \brief The gcd of each integer of a set with the product of the others.
 *
 * Two RSA moduli that share a prime both fall to a gcd, and among millions
 * of moduli those few pairs are found without computing every pairwise gcd:
 * with Bernstein's product and remainder trees, as Heninger, Durumeric,
 * Wustrow and Halderman did on the keys of the Internet ("Mining Your Ps and
 * Qs", USENIX Security 2012). The product P of the n_i is computed by a
 * product tree, then P modulo each n_i^2 by a remainder tree, and the gcd of
 * n_i with the product of the others is that of n_i and
 * (P mod n_i^2) / n_i.
 */

#if !defined(_TIFA_BATCH_GCD_H_)
   /**
    * \def _TIFA_BATCH_GCD_H_
    * Standard include guard.
    */
#define _TIFA_BATCH_GCD_H_

#include "tifa_config.h"
#include "array.h"
#include "exit_codes.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \brief Appends to \c gcds the gcd of each integer of \c moduli with the
    *        product of the others.
    *
    * For each n_i of \c moduli, in order, appends gcd(n_i, the product of
    * the n_j, j != i): 1 when n_i shares no prime with the others, n_i
    * itself when all its primes are shared (as with a duplicate). The time
    * is that of a few multiplications of numbers the size of the product
    * of the n_i at each of the log2(count) levels of the trees, and the
    * memory about log2(count) times the size of that product: for 100,000
    * moduli of 2048 bits, ... s and ... MB on an Apple M5 Pro.
    *
    * Gives up, appending nothing, when the calling thread's time limit
    * passes or its interrupt check returns true (see time_limit.h).
    *
    * \param[out] gcds The array the gcds are appended to.
    * \param[in] moduli The integers, all positive.
    * \return \c SUCCESS, \c FAILURE if some integer is not positive (nothing
    *         appended), or \c GIVING_UP.
    */
TIFA_API ecode_t tifa_batch_gcd(mpz_array_t* gcds, const mpz_array_t* moduli);

#ifdef __cplusplus
}
#endif

#endif
