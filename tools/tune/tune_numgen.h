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
 * \file    tune_numgen.h
 * \brief   Test-composite generator for the tuning tool.
 *
 * Generates composites that are the product of \c nfactors primes of similar
 * bit-size (a balanced semiprime for <tt>nfactors == 2</tt>).
 */

#if !defined(_TIFA_TUNE_NUMGEN_H_)
#define _TIFA_TUNE_NUMGEN_H_

#include <gmp.h>

   /**
    * \def TUNE_MAX_NFACTORS
    * Maximum number of prime factors a generated composite may have.
    */
#define TUNE_MAX_NFACTORS 8

   /**
    * \brief Builds a composite that is a product of similar-size primes.
    *
    * On return <tt>n</tt> is the product of \c nfactors distinct primes, each
    * of roughly <tt>n_bits / nfactors</tt> bits, with \c n having exactly
    * \c n_bits bits (retried internally until the size matches). \c smallest
    * receives the least of those prime factors.
    *
    * The caller owns \c n and \c smallest (already \c mpz_init'd) and the random
    * state \c rng. Requires <tt>n_bits >= 2 * nfactors</tt>.
    *
    * \param[out]    n        The generated composite.
    * \param[out]    smallest The smallest prime factor of \c n.
    * \param[in]     nfactors Number of prime factors (2..TUNE_MAX_NFACTORS).
    * \param[in]     n_bits   Target bit-size of \c n.
    * \param[in,out] rng      A seeded GMP random state.
    */
void tune_gen_composite(mpz_t n, mpz_t smallest,
                        unsigned int nfactors, unsigned int n_bits,
                        gmp_randstate_t rng);

#endif
