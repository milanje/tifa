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
 * \file    mont_ecm.h
 *
 * \brief ECM curves, p-1 and p+1 in fixed-size Montgomery arithmetic, for
 *        moduli of one to four 64-bit limbs.
 *
 * Stage 1 raises to the product E of the prime powers up to B1 (Montgomery
 * ladder for ECM, binary powering for p-1, Lucas ladder for p+1); stage 2 is
 * the standard continuation with baby steps and giant steps. The work that
 * depends only on the bounds (E, stage 2 schedule) is done once in a plan
 * shared by all the runs.
 */

#if !defined(_TIFA_MONT_ECM_H_)
   /**
    * \def _TIFA_MONT_ECM_H_
    * Standard include guard.
    */
#define _TIFA_MONT_ECM_H_

#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

#if    defined(__SIZEOF_INT128__) && (GMP_LIMB_BITS == 64) \
    && (GMP_NAIL_BITS == 0)
    #define TIFA_HAVE_MONT_ECM 1
#else
    #define TIFA_HAVE_MONT_ECM 0
#endif

   /**
    * \def MONT_ECM_MAX_LIMBS
    * Largest modulus size, in 64-bit limbs.
    */
#define MONT_ECM_MAX_LIMBS 4

#if TIFA_HAVE_MONT_ECM

   /**
    * \typedef mont_ecm_plan_t
    * \brief The precomputations shared by all the curves for given bounds.
    */
typedef struct struct_mont_ecm_plan_t mont_ecm_plan_t;

   /**
    * \brief Allocates the plan for the bounds \c b1 and \c b2.
    *
    * \c b1 is lowered to the largest tabulated prime (see
    * first_primes.h) if needed. If <tt>b2 > b1</tt>, \c b1 is raised to at
    * least 105 (half the stage 2 step); no stage 2 is planned if \c b2 does
    * not exceed the resulting \c b1.
    */
mont_ecm_plan_t* alloc_mont_ecm_plan(uint32_t b1, uint32_t b2);

   /**
    * \brief Frees a plan allocated by \c alloc_mont_ecm_plan.
    */
void free_mont_ecm_plan(mont_ecm_plan_t* plan);

   /**
    * \brief Whether \c mont_ecm_curve handles \c n (odd, one to four limbs).
    */
bool mont_ecm_handles(const mpz_t n);

   /**
    * \brief Runs one ECM curve on \c n.
    *
    * Uses the curve of Suyama's family with parameter \c sigma (at least 6).
    *
    * \param[out] factor A proper factor of \c n, if one is found.
    * \param[in]  n The number to factor (see \c mont_ecm_handles).
    * \param[in]  sigma The curve parameter.
    * \param[in]  plan The bounds and precomputations.
    * \return true if a proper factor was found.
    */
bool mont_ecm_curve(mpz_t factor, const mpz_t n, uint32_t sigma,
                    const mont_ecm_plan_t* plan);

   /**
    * \brief Runs the k-th curve of TIFA's sequence on \c n.
    *
    * The curves with torsion Z/2 x Z/8 of z2z8_curves.h first (orders
    * divisible by 16), then Suyama's curves from sigma = 6.
    *
    * \param[out] factor A proper factor of \c n, if one is found.
    * \param[in]  n The number to factor (see \c mont_ecm_handles).
    * \param[in]  k The index of the curve in the sequence.
    * \param[in]  plan The bounds and precomputations.
    * \return true if a proper factor was found.
    */
bool mont_ecm_nth_curve(mpz_t factor, const mpz_t n, uint32_t k,
                        const mont_ecm_plan_t* plan);

   /**
    * \brief Runs Pollard's p-1 on \c n.
    *
    * Stage 1 computes a = base^E; stage 2 looks for a prime q of the plan
    * with a^q = 1 modulo a prime factor of \c n.
    *
    * \param[out] factor A proper factor of \c n, if one is found.
    * \param[in]  n The number to factor (see \c mont_ecm_handles).
    * \param[in]  base The base, smaller than \c n.
    * \param[in]  plan The bounds and precomputations.
    * \return 1 if a proper factor was found, -1 if a gcd was \c n (all the
    *         prime factors were caught at once), 0 otherwise.
    */
int mont_pm1(mpz_t factor, const mpz_t n, uint32_t base,
             const mont_ecm_plan_t* plan);

   /**
    * \brief Runs Williams's p+1 on \c n.
    *
    * Stage 1 computes the Lucas sequence term V_E(seed); stage 2 works as
    * for \c mont_pm1.
    *
    * \param[out] factor A proper factor of \c n, if one is found.
    * \param[in]  n The number to factor (see \c mont_ecm_handles).
    * \param[in]  seed The seed P of V_n(P, 1), smaller than \c n.
    * \param[in]  plan The bounds and precomputations.
    * \return As for \c mont_pm1.
    */
int mont_pp1(mpz_t factor, const mpz_t n, uint32_t seed,
             const mont_ecm_plan_t* plan);

#endif

#ifdef __cplusplus
}
#endif

#endif
