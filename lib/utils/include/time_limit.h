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
 * \file    time_limit.h
 * \brief   A time limit on the factorizations of a thread, and interrupts.
 */

#if !defined(_TIFA_TIME_LIMIT_H_)
#define _TIFA_TIME_LIMIT_H_

#include <stdbool.h>

#include "tifa_config.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \brief Limits the time of the calling thread's factorizations.
    *
    * From \c seconds after this call, \c tifa_factor, \c tifa_prove_prime
    * and the algorithms that take long (SIQS, its linear algebra and ECM's
    * curves) give up when they next check the time, and return what they
    * found: a partial factorization at best. The limit is approximate: SIQS
    * checks it every 16 chunks of its sieve, ECM after each curve, the
    * elimination every 256 columns. The other algorithms check only between
    * the runs of a factoring mode other than \c SINGLE_RUN. Each thread has
    * its own limit, which holds until the next call: \c seconds of 0 or
    * less, or not finite, removes it.
    *
    * \param[in] seconds The time left from now, in seconds.
    */
TIFA_API void tifa_set_time_limit(double seconds);

   /**
    * \brief Returns whether the calling thread's time limit has passed.
    *
    * \return true if the thread has a time limit and it has passed.
    */
TIFA_API bool tifa_time_limit_reached(void);

   /**
    * \typedef tifa_interrupt_check_t
    * \brief A function telling a thread's factorizations whether to stop.
    */
typedef bool (*tifa_interrupt_check_t)(void* arg);

   /**
    * \brief Lets a function interrupt the calling thread's factorizations.
    *
    * The factorizations of the calling thread call <tt>check(arg)</tt> where
    * they check their time limit (see \c tifa_set_time_limit): every few
    * milliseconds in SIQS and ECM. Once it returns true, they give up as at
    * the time limit, until the next call. A program can so stop them from
    * another thread, through a flag that \c check reads, or on a signal.
    * \c check runs in the calling thread, and should take little time. NULL
    * removes it.
    *
    * \param[in] check The function, or NULL.
    * \param[in] arg   Its argument.
    */
TIFA_API void tifa_set_interrupt_check(tifa_interrupt_check_t check,
                                       void* arg);

   /**
    * \brief Returns whether the calling thread's factorizations were
    * interrupted.
    *
    * \return true if the function set by \c tifa_set_interrupt_check has
    *         returned true since.
    */
TIFA_API bool tifa_interrupted(void);

#ifdef __cplusplus
}
#endif

#endif
