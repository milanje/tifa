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
 * \file    ecm_schedule.h
 *
 * \brief A bounded search for a factor of up to some size: p-1, then ECM
 *        curves level by level (fixed-limb arithmetic, see mont_ecm.h).
 *
 * Used where a composite would otherwise go to SIQS: the search runs as long
 * as it costs a given fraction of that SIQS run, estimated from the size of
 * the composite.
 */

#if !defined(_TIFA_ECM_SCHEDULE_H_)
   /**
    * \def _TIFA_ECM_SCHEDULE_H_
    * Standard include guard.
    */
#define _TIFA_ECM_SCHEDULE_H_

#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

#include "mont_ecm.h"

#ifdef __cplusplus
extern "C" {
#endif

#if TIFA_HAVE_MONT_ECM

   /**
    * \typedef ecm_schedule_t
    * \brief The plans of the p-1 run and of the ECM levels.
    */
typedef struct struct_ecm_schedule_t ecm_schedule_t;

   /**
    * \brief Allocates the schedule for factors of up to \c max_bits bits.
    *
    * The plans of the levels above 60 bits are built when a search first
    * reaches them; those above 85 bits, freed when the last search using
    * them ends.
    */
ecm_schedule_t* alloc_ecm_schedule(uint32_t max_bits);

   /**
    * \brief Frees a schedule allocated by \c alloc_ecm_schedule.
    */
void free_ecm_schedule(ecm_schedule_t* schedule);

   /**
    * \brief Searches \c n for a factor of at most \c target_bits bits.
    *
    * Runs p-1, then ECM curves level by level (by 5 bits of factor size),
    * skipping the levels below \c min_bits (\c n has no prime factor below
    * 2^min_bits) and stopping after the one that holds \c target_bits, or
    * before the curves' estimated time exceeds \c budget times that of SIQS
    * on \c n.
    *
    * \param[out] g A proper factor of \c n, if one is found.
    * \param[in]  n The odd composite to search (see \c mont_ecm_handles).
    * \return true if a proper factor was found.
    */
bool ecm_schedule_split(mpz_t g, const mpz_t n,
                        ecm_schedule_t* schedule, double min_bits,
                        double target_bits, double budget);

   /**
    * \brief Searches \c n for a factor of at most half its size, and of at
    * most \c max_bits bits, before SIQS (with an infinite budget: instead).
    *
    * As \c ecm_schedule_split, except that the first level runs whole as
    * soon as one of its curves fits in the budget. With a budget of 0.2, the
    * level does not fit from 74 to 83 bits, but at 80 bits its two curves
    * find two thirds of the 22-bit factors that p-1 misses, for a quarter of
    * the SIQS run.
    *
    * The curves run in as many threads as the options allow (see
    * tifa_options.h), at most one per 0.1 ms of them, the first factor
    * found stopping all.
    *
    * \param[out] g A proper factor of \c n, if one is found.
    * \param[in]  n The odd composite to search (see \c mont_ecm_handles).
    * \return true if a proper factor was found.
    */
bool ecm_schedule_prepass(mpz_t g, const mpz_t n,
                          ecm_schedule_t* schedule, double budget,
                          double max_bits);


#endif

#ifdef __cplusplus
}
#endif

#endif
