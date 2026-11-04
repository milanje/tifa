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
 * \file    prime_cert.h
 *
 * \brief Primality certificates: proofs that the factors found are prime.
 *
 * \c tifa_factor tells primes from composites with a Baillie-PSW test,
 * which has no known counterexample but is no proof. \c tifa_prove_prime
 * proves that an integer is prime and returns a certificate that anyone can
 * check with modular exponentiations and GCDs, without trusting TIFA:
 * \c tifa_check_prime_cert is such a check, and \c tifa.check_certificate
 * another, in Python. Nothing here is used by \c tifa_factor, whose speed
 * is the same with or without it.
 *
 * A certificate is a list of steps. Each step proves a prime p from primes
 * q_1, ..., q_k of p - 1, each below 2^64 or proven by an earlier step,
 * and integers a_1, ..., a_k. Let F be the part of p - 1 made of the q_i
 * (the product of the largest powers of the q_i dividing p - 1). If, for
 * each i,
 *
 *     a_i^(p - 1) = 1 (mod p) and gcd(a_i^((p - 1) / q_i) - 1, p) = 1,
 *
 * every prime factor of p is 1 mod F (Pocklington, 1914), and p is prime if
 *
 *   - (F + 1)^2 > p, or
 *   - F^3 > p and c1^2 - 4 c2 is not a square, where
 *     p = c2 F^2 + c1 F + 1 and 0 <= c1 < F (Brillhart, Lehmer and
 *     Selfridge, 1975; Theorem 4.1.5 of Crandall and Pomerance's "Prime
 *     Numbers: A Computational Perspective").
 *
 * The last step proves the integer certified. Primes below 2^64 need no
 * step: they are strong probable primes to the bases 2, 3, 5, ..., 37 (the
 * first 12 primes), which no composite below 2^64 is (Sorenson and Webster,
 * 2017, "Strong pseudoprimes to twelve prime bases").
 *
 * To prove p, \c tifa_prove_prime factors p - 1 with \c tifa_factor, keeps
 * its smallest prime factors until F^3 > p, and proves those above 2^64 in
 * the same way. The time it takes is mostly that of factoring p - 1: on an
 * Apple M5 Pro, 0.4 ms on average for random primes of 128 bits, 15 ms at
 * 192 bits and 0.2 s at 224 bits, but up to 3 s when p - 1 has two large
 * prime factors, which SIQS has to split; about 2 s at 256 bits.
 * There is no time limit. A check takes less than 0.1 ms.
 */

#if !defined(_TIFA_PRIME_CERT_H_)
   /**
    * \def _TIFA_PRIME_CERT_H_
    * Standard include guard.
    */
#define _TIFA_PRIME_CERT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

#include "array.h"

   /**
    * \enum prime_cert_status_enum
    * \brief The outcome of \c tifa_prove_prime.
    */
enum prime_cert_status_enum {
        /**
         * The integer is prime, and the certificate proves it.
         */
    PRIME_PROVEN = 0,
        /**
         * The integer is not prime (proven: it is below 2, or a witness
         * of its compositeness was found).
         */
    NOT_PRIME,
        /**
         * Undecided: n - 1, or the p - 1 of a prime below, could not be
         * factored far enough.
         */
    PRIME_NOT_PROVEN
};
   /**
    * \typedef prime_cert_status_t
    * \brief Equivalent to <tt>enum prime_cert_status_enum</tt>.
    */
typedef enum prime_cert_status_enum prime_cert_status_t;

   /**
    * \struct struct_prime_cert_step_t prime_cert.h lib/algo/include/prime_cert.h
    * \brief  One step of a certificate: the proof that \c p is prime.
    */
struct struct_prime_cert_step_t {
        /**
         * The prime proven.
         */
    mpz_t p;
        /**
         * The primes q_i of p - 1, in increasing order.
         */
    mpz_array_t* primes;
        /**
         * The a_i, one for each q_i.
         */
    uint32_array_t* witnesses;
};
   /**
    * \typedef prime_cert_step_t
    * \brief Equivalent to <tt>struct struct_prime_cert_step_t</tt>.
    */
typedef struct struct_prime_cert_step_t prime_cert_step_t;

   /**
    * \struct struct_prime_cert_t prime_cert.h lib/algo/include/prime_cert.h
    * \brief  A primality certificate: steps, each using primes of earlier
    *         ones.
    */
struct struct_prime_cert_t {
        /**
         * Number of steps, none for a prime below 2^64.
         */
    uint32_t nsteps;
        /**
         * Number of steps allocated.
         */
    uint32_t alloced;
        /**
         * The steps, the last one proving the integer certified.
         */
    prime_cert_step_t* steps;
};
   /**
    * \typedef prime_cert_t
    * \brief Equivalent to <tt>struct struct_prime_cert_t</tt>.
    */
typedef struct struct_prime_cert_t prime_cert_t;

   /**
    * \brief Allocates an empty certificate.
    */
TIFA_API prime_cert_t* alloc_prime_cert(void);

   /**
    * \brief Frees a certificate allocated by \c alloc_prime_cert.
    */
TIFA_API void free_prime_cert(prime_cert_t* cert);

   /**
    * \brief Proves that \c n is prime.
    *
    * \param[out] cert The certificate of \c n on \c PRIME_PROVEN (its
    *             previous steps are dropped).
    * \param[in]  n The integer to prove prime.
    * \return The outcome.
    */
TIFA_API prime_cert_status_t tifa_prove_prime(prime_cert_t* cert,
                                              const mpz_t n);

   /**
    * \brief Checks a certificate.
    *
    * \param[in] cert A certificate.
    * \param[in] n The integer it should prove prime.
    * \return \c true if and only if every step of \c cert holds and the
    *         last one proves \c n prime (or, with no steps, if \c n is a
    *         prime below 2^64).
    */
TIFA_API bool tifa_check_prime_cert(const prime_cert_t* cert, const mpz_t n);

#ifdef __cplusplus
}
#endif

#endif
