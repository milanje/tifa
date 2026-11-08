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
 * A certificate is a list of steps, each proving a prime p from primes of
 * p - 1 and of p + 1, each below 2^64 or proven by an earlier step.
 *
 * The primes q_1, ..., q_k of p - 1 come with integers a_1, ..., a_k. Let
 * F1 be the part of p - 1 made of the q_i (the product of the largest
 * powers of the q_i dividing p - 1). If, for each i,
 *
 *     a_i^(p - 1) = 1 (mod p) and gcd(a_i^((p - 1) / q_i) - 1, p) = 1,
 *
 * every prime factor of p is 1 mod F1 (Pocklington, 1914).
 *
 * The primes r_1, ..., r_l of p + 1 come with an integer D and integers
 * P_1, ..., P_l. Let F2 be the part of p + 1 made of the r_j, Q_j =
 * (P_j^2 - D) / 4 and U(j) the Lucas sequence of parameters P_j and Q_j
 * (U_0 = 0, U_1 = 1, U_(m+1) = P_j U_m - Q_j U_(m-1)). If gcd(p, D) = 1
 * and, for each j, P_j^2 - D is a multiple of 4, gcd(p, Q_j) = 1 and
 *
 *     p divides U(j)_(p + 1) and gcd(U(j)_((p + 1) / r_j), p) = 1,
 *
 * every prime factor s of p is (D / s) mod F2, that is 1 or -1 (Morrison,
 * 1975: the rank of apparition of s divides s - (D / s)).
 *
 * Then p is prime if
 *
 *   - (F1 + 1)^2 > p, or
 *   - F1^3 > p and c1^2 - 4 c2 is not a square, where
 *     p = c2 F1^2 + c1 F1 + 1 and 0 <= c1 < F1 (Brillhart, Lehmer and
 *     Selfridge, 1975; Theorem 4.1.5 of Crandall and Pomerance's "Prime
 *     Numbers: A Computational Perspective"), or
 *   - G^2 > p, G being the lcm of F1 and F2, and c = p mod G is 1, or at
 *     least p, or does not divide p: a prime factor of p below its square
 *     root would be 1 mod G, impossible, or c, the residue that is 1 mod F1
 *     and -1 mod F2.
 *
 * The last step proves the integer certified. Primes below 2^64 need no
 * step: they are strong probable primes to the bases 2, 3, 5, ..., 37 (the
 * first 12 primes), which no composite below 2^64 is (Sorenson and Webster,
 * 2017, "Strong pseudoprimes to twelve prime bases").
 * \c tifa_prime_cert_to_str writes a certificate as text, in JSON, and
 * \c tifa_prime_cert_from_str reads it back.
 *
 * To prove p, \c tifa_prove_prime searches p - 1 for factors by trial
 * division and ECM, then p + 1 if they are not enough, and, if still not,
 * factors the smallest composite left by SIQS. It keeps the smallest primes
 * found until one of the conditions above holds (from p - 1 alone if it
 * can), and proves those above 2^64 in the same way. On an Apple M5 Pro, a
 * proof took 0.5 ms on average for random primes of 128 bits, 1.4 ms at 160
 * bits, 5.5 ms at 192 bits, 53 ms at 224 bits (1.1 s at most, when both
 * p - 1 and p + 1 needed SIQS) and 66 ms at 256 bits; a check takes less
 * than 0.5 ms. Above, the composites left grow with p: a proof can take as
 * long as SIQS on a composite of almost its size (minutes at 280 bits, far
 * longer above), unless a time limit (see tifa_options.h) stops it with
 * \c PRIME_NOT_PROVEN; composites above the options' \c siqs_max_bits (320
 * by default) never go to SIQS.
 */

#if !defined(_TIFA_PRIME_CERT_H_)
   /**
    * \def _TIFA_PRIME_CERT_H_
    * Standard include guard.
    */
#define _TIFA_PRIME_CERT_H_

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#include "array.h"
#include "tifa_options.h"

#ifdef __cplusplus
extern "C" {
#endif

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
    mpz_array_t* pm1_primes;
        /**
         * The a_i, one for each q_i.
         */
    uint32_array_t* pm1_witnesses;
        /**
         * D, or 0 if no prime of p + 1 is used.
         */
    int32_t pp1_d;
        /**
         * The primes r_j of p + 1, in increasing order.
         */
    mpz_array_t* pp1_primes;
        /**
         * The P_j, one for each r_j.
         */
    uint32_array_t* pp1_witnesses;
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
    * \param[in]  options The options of the call (time limit, interrupt
    *             check, threads, sizes: see tifa_options.h), or NULL for the
    *             defaults.
    * \return The outcome: \c PRIME_NOT_PROVEN if the options are invalid.
    */
TIFA_API prime_cert_status_t tifa_prove_prime(
    prime_cert_t* cert,
    const mpz_t n,
    const tifa_options_t* const options
);

   /**
    * \brief Checks a certificate.
    *
    * Each prime has at most one step: the last for \c n, the others for
    * primes below \c n. The time is then at most that of the step for \c n
    * times the number of steps.
    *
    * \param[in] cert A certificate.
    * \param[in] n The integer it should prove prime.
    * \return \c true if and only if every step of \c cert holds and the
    *         last one proves \c n prime (or, with no steps, if \c n is a
    *         prime below 2^64).
    */
TIFA_API bool tifa_check_prime_cert(const prime_cert_t* cert, const mpz_t n);

   /**
    * \brief Writes a certificate as text.
    *
    * The text is JSON, the format of the Python package's
    * \c tifa.check_certificate: a list of steps, one per line, each
    * <tt>[p, [[q_1, a_1], ...]]</tt>, or
    * <tt>[p, [[q_1, a_1], ...], D, [[r_1, P_1], ...]]</tt> with primes of
    * p + 1. A prime below 2^64 has <tt>[]</tt>.
    *
    * As \c snprintf does, it writes at most \c size characters, the
    * terminating null one included. A step whose arrays differ in length,
    * or with primes of p + 1 and D = 0, which \c tifa_check_prime_cert
    * refuses, is written so that \c tifa_prime_cert_from_str refuses it.
    *
    * \param[out] str The text, which can be \c NULL if \c size is 0.
    * \param[in]  size The size of \c str.
    * \param[in]  cert A certificate.
    * \return The length of the whole text, the terminating null excluded:
    *         the text was truncated if it is \c size or more.
    */
TIFA_API size_t tifa_prime_cert_to_str(char* str, size_t size,
                                       const prime_cert_t* cert);

   /**
    * \brief Reads a certificate from text.
    *
    * The text is that of \c tifa_prime_cert_to_str, or the same JSON with
    * other spaces, tabs or newlines between its numbers and brackets (that
    * of \c json.dumps of a Python certificate, for example). Only its form
    * is checked: \c tifa_check_prime_cert checks the certificate.
    *
    * \param[out] cert The certificate read (its previous steps are dropped),
    *             without steps if the text is not one.
    * \param[in]  str The text, null-terminated.
    * \return \c true if \c str is a list of steps and nothing else but
    *         spaces: p, the q_i, a_i, r_j and P_j at least 0, the a_i and
    *         P_j below 2^32, and D, if there is one, not 0 and below 2^31 in
    *         absolute value.
    */
TIFA_API bool tifa_prime_cert_from_str(prime_cert_t* cert, const char* str);

#ifdef __cplusplus
}
#endif

#endif
