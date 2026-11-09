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
 * \file    primes.h
 *
 * \brief The primes up to 2^32 - 1, from the table of first_primes.h and
 *        then from a segmented sieve, and the stage-1 exponent of p-1 and
 *        ECM over them.
 */

#if !defined(_TIFA_PRIMES_H_)
   /**
    * \def _TIFA_PRIMES_H_
    * Standard include guard.
    */
#define _TIFA_PRIMES_H_

#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \struct prime_iter_t primes.h lib/utils/include/primes.h
    * \brief The state of an enumeration of the primes in [lo, hi].
    */
typedef struct {
    uint64_t hi;
    uint64_t start;
    uint64_t seg_lo;
    uint32_t ip;
    uint32_t pos;
    uint32_t len;
    bool     done;
    uint8_t* composite;
} prime_iter_t;

   /**
    * \brief Starts the enumeration of the primes in [lo, hi].
    */
void prime_iter_init(prime_iter_t* it, uint32_t lo, uint32_t hi);

   /**
    * \brief The next prime of the enumeration, in increasing order, or 0
    * after the last one.
    */
uint32_t prime_iter_next(prime_iter_t* it);

   /**
    * \brief Frees the memory of an enumeration.
    */
void prime_iter_clear(prime_iter_t* it);

   /**
    * \brief Sets \c e to the product of the largest powers of the primes up
    * to \c b that do not exceed \c b.
    *
    * \return false, with \c e unspecified, if \c stoppable and the calling
    *         thread's call should stop (see stopping.h) before the end.
    */
bool prime_powers_product(mpz_t e, uint32_t b, bool stoppable);

#ifdef __cplusplus
}
#endif

#endif
