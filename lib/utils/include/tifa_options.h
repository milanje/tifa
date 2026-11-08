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
 * \file    tifa_options.h
 * \brief   The options of a factorization: time limit, interrupt, threads.
 */

#if !defined(_TIFA_TIFA_OPTIONS_H_)
   /**
    * \def _TIFA_TIFA_OPTIONS_H_
    * Standard include guard.
    */
#define _TIFA_TIFA_OPTIONS_H_

#include <stdbool.h>
#include <stdint.h>

#include "tifa_config.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \typedef tifa_interrupt_check_t
    * \brief A function telling a factorization whether to stop.
    */
typedef bool (*tifa_interrupt_check_t)(void* arg);

   /**
    * \enum tifa_factor_hint_enum
    * \brief What is known of the number given to \c tifa_factor.
    *
    * A hint never changes the factors found, only the time taken to find
    * them. It applies to the prime factors of n above 1000, and is ignored up
    * to 64 bits (58 without 128-bit integers), where machine words (or
    * SQUFOF) factor n, and above \c siqs_max_bits, where ECM alone searches.
    */
enum tifa_factor_hint_enum {
        /**
         * Nothing is known: ECM searches for a factor for 20% of the time
         * SIQS would take, then SIQS runs. The default.
         */
    TIFA_FACTOR_RANDOM = 0,
        /**
         * n has two prime factors of about the same size (an RSA modulus):
         * SIQS runs right away.
         */
    TIFA_FACTOR_BALANCED,
        /**
         * n has a prime factor of less than about a third of its size: ECM
         * searches for as long as SIQS would take (with 128-bit integers,
         * up to 512 bits; otherwise, and above, as \c TIFA_FACTOR_RANDOM).
         */
    TIFA_FACTOR_UNBALANCED
};
   /**
    * \typedef tifa_factor_hint_t
    * \brief Equivalent to <tt>enum tifa_factor_hint_enum</tt>.
    */
typedef enum tifa_factor_hint_enum tifa_factor_hint_t;

   /**
    * \def TIFA_SIQS_MAX_BITS
    * The default of \c siqs_max_bits: SIQS would take about 40 minutes at
    * 320 bits on an Apple M5 Pro (extrapolated from 288 bits), twice as long
    * every 8 bits above.
    */
#define TIFA_SIQS_MAX_BITS 320

   /**
    * \def TIFA_ECM_MAX_BITS
    * The default of \c ecm_max_bits.
    */
#define TIFA_ECM_MAX_BITS 85

   /**
    * \def TIFA_ECM_TOP_BITS
    * The largest \c ecm_max_bits: ECM's levels stop at factors of this size
    * (stage 1 takes its primes from the table of first_primes.h).
    */
#define TIFA_ECM_TOP_BITS 85

   /**
    * \struct struct_tifa_options_t tifa_options.h
    *         lib/utils/include/tifa_options.h
    * \brief  The options of one call: all zero means the defaults.
    *
    * The functions that can run long take a pointer to them, or NULL for the
    * defaults: \c tifa_factor, \c tifa_prove_prime, \c tifa_cofactor,
    * \c tifa_batch_gcd and the algorithms (\c siqs, \c ecm, ...). They hold
    * for that call only, the calls it makes included. A call made from within
    * another (from its interrupt check, for instance) has its own, or the
    * defaults, and gives the outer call back its own when it returns; the
    * interrupt check is not called again from there. Several threads can
    * share the same options: their interrupt check is then called from each
    * of them, with the same argument. A time limit counts from each call:
    * for one budget over several calls, give each the time left.
    *
    * \code
    * tifa_options_t options = { .time_limit = 10, .threads = 4 };
    * ecode_t e = tifa_factor(factors, multis, n,
    *                         FIND_COMPLETE_FACTORIZATION, &options);
    * \endcode
    */
struct struct_tifa_options_t {
        /**
         * The time limit, in seconds from the call: 0 or less, or not
         * finite, for none. Past it, the call gives up when it next checks
         * the time, and returns what it found: a partial factorization at
         * best. SIQS checks it every 16 chunks of its sieve, ECM after each
         * curve (and within it above 256 bits), CFRAC after each batch of
         * residues, the elimination every 256 columns. The other algorithms
         * check only between the runs of a factoring mode other than
         * \c SINGLE_RUN.
         */
    double time_limit;
        /**
         * The number of threads, the calling one included (0 or 1: the
         * calling thread alone, as fast as without threads). SIQS collects
         * its relations in up to that many, at most one per 100 primes of
         * its factor base (none below about 100 bits) and 64 in all. The ECM
         * search that \c tifa_factor and \c tifa_prove_prime run before SIQS
         * shares its curves among as many, at most one per 0.1 ms of curves:
         * from about 120 bits for \c tifa_factor without a hint, 100 with
         * the unbalanced one, 135 for \c tifa_prove_prime. The calling thread
         * alone checks the time limit and the interrupt check, and stops the
         * others. To factor many numbers, one factorization per thread uses
         * the cores best.
         */
    uint32_t threads;
        /**
         * A function called where the time limit is checked, every few
         * milliseconds in SIQS and ECM, in the calling thread, with
         * \c interrupt_arg: once it returns true, the call gives up as at
         * its time limit. A program can so stop a factorization from
         * another thread, through a flag that it reads, or on a signal. It
         * should take little time. NULL for none.
         */
    tifa_interrupt_check_t interrupt;
        /**
         * The argument of \c interrupt.
         */
    void* interrupt_arg;
        /**
         * What is known of the number given to \c tifa_factor, which refuses
         * an unknown one with \c FAILURE (the other functions ignore it).
         */
    tifa_factor_hint_t hint;
        /**
         * The largest composite, in bits, that \c tifa_factor (and
         * \c tifa_prove_prime) hands to SIQS: above, ECM alone searches it,
         * for factors of up to \c ecm_max_bits bits, and a composite left
         * gives \c FACTORS_TOO_LARGE. 0 for \c TIFA_SIQS_MAX_BITS.
         */
    uint32_t siqs_max_bits;
        /**
         * The size, in bits, of the largest factors ECM searches for, by
         * itself (\c ecm) and on the composites above \c siqs_max_bits: 0
         * for \c TIFA_ECM_MAX_BITS (85 bits, also the largest: 12 minutes at
         * 520 bits on an Apple M5 Pro). Each 5 bits less halve about the
         * time, for an earlier \c FACTORS_TOO_LARGE.
         */
    uint32_t ecm_max_bits;
        /**
         * For the options of later versions: must be zero (the functions
         * refuse other values).
         */
    uint64_t reserved[8];
};
   /**
    * \typedef tifa_options_t
    * \brief Equivalent to <tt>struct struct_tifa_options_t</tt>.
    */
typedef struct struct_tifa_options_t tifa_options_t;

#ifdef __cplusplus
}
#endif

#endif
