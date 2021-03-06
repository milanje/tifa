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
 * \file    first_primes.h
 * \author  Automatically generated by genprimes.pl
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Precomputed small primes.
 *
 * This is a list of the precomputed small primes together with
 * a \c uint32_array_t wrapper.
 */

#if !defined(_TIFA_FIRST_PRIMES_H_)
   /**
    * \def _TIFA_FIRST_PRIMES_H_
    * Standard include guard.
    */
#define _TIFA_FIRST_PRIMES_H_

#include <inttypes.h>

#include "array.h"
#include "tifa_config.h"

/**
 * \def NFIRST_PRIMES
 * Number of precomputed primes in the \c first_primes array.
 */
#define NFIRST_PRIMES 65536

/**
 * The \c first_primes array is a global array of \c uint32_t elements
 * containing the first \c NFIRST_PRIMES prime numbers (from 2 and beyond).
 */
extern const uint32_t first_primes[NFIRST_PRIMES] MAYBE_UNUSED;

/**
 * The largest prime in the \c first_primes array.
 */
extern const uint32_t LARGEST_PRIME MAYBE_UNUSED;

/**
 * \c first_primes_array is a \c uint32_array_t wrapper to the array
 * \c first_primes.
 *
 * \note \c first_primes_array 's \c alloced field is set to zero. Indeed,
 * \c first_primes_array is merely a \c uint32_array_t wrapper for
 * \c first_primes, and as such, it has no real "alloced" memory. Setting
 * \c first_primes_array.alloced to 0 will prevent errors
 * if \c free_mpz_array is inadvertently called on \c first_primes_array.
 */
extern const uint32_array_t first_primes_array MAYBE_UNUSED;

#endif

