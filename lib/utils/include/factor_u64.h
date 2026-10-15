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
 * \file    factor_u64.h
 *
 * \brief Primality testing and factorization of integers below 2^64 with
 *        machine-word arithmetic.
 *
 * Trial division by the primes below 1000, then Hart's one line factoring
 * (up to 40 bits), ECM on Montgomery curves and Brent's rho, all in 64-bit
 * Montgomery arithmetic. Only available where the compiler provides 128-bit
 * integers and \c unsigned \c long is 64 bits wide.
 */

#if !defined(_TIFA_FACTOR_U64_H_)
   /**
    * \def _TIFA_FACTOR_U64_H_
    * Standard include guard.
    */
#define _TIFA_FACTOR_U64_H_

#include <stdint.h>
#include <stdbool.h>

#include "tifa_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SIZEOF_INT128__) && (TIFA_SIZEOF_ULONG_T >= 8)
    #define TIFA_HAVE_FACTOR_U64 1
#else
    #define TIFA_HAVE_FACTOR_U64 0
#endif

   /**
    * \def FACTOR_U64_MAX_PRIMES
    * Largest number of distinct prime factors of an integer below 2^64.
    */
#define FACTOR_U64_MAX_PRIMES 15

#if TIFA_HAVE_FACTOR_U64

   /**
    * \brief Deterministic primality test for 64-bit integers.
    *
    * \param[in] n The integer to test.
    * \return true if \c n is prime, false otherwise.
    */
bool is_prime_u64(uint64_t n);

   /**
    * \brief Complete factorization of a 64-bit integer.
    *
    * Stores the distinct prime factors of \c n in increasing order in
    * \c primes and their multiplicities in \c multis (both arrays must hold
    * at least \c FACTOR_U64_MAX_PRIMES entries).
    *
    * \param[in]  n The integer to factor, at least 2.
    * \param[out] primes The distinct prime factors of \c n.
    * \param[out] multis Their multiplicities.
    * \return The number of distinct prime factors, or 0 if \c n < 2 or a
    *         composite factor could not be split.
    */
uint32_t factor_u64(uint64_t n, uint64_t* primes, uint32_t* multis);

#endif

#ifdef __cplusplus
}
#endif

#endif
