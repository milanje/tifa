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
 * \file    bernsteinisms.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Algorithms from two D. J. Bernstein's papers on the factorization
 * of small integers.
 *
 * Algorithms from two D. J. Bernstein's papers on the factorization
 * of small integers:
 * \li "How to find small factors of integers", http://cr.yp.to/papers/sf.pdf
 * \li "How to find smooth parts of integers",
 *   http://cr.yp.to/factorization/smoothparts-20040510.pdf
 */

#if !defined(_TIFA_BERNSTEINISMS_H_)
   /**
    * \def _TIFA_BERNSTEINISMS_H_
    * Standard include guard.
    */
#define _TIFA_BERNSTEINISMS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include "array.h"
#include "x_array_list.h"
#include "hashtable.h"
#include "smooth_filter.h"

//-----------------------------------------------------------------------------
// Ref. "How to find small factors of integers", Daniel J. Bernstein
//       http://cr.yp.to/papers/sf.pdf
//-----------------------------------------------------------------------------

   /**
    * \brief Daniel J. Bernstein's algorithm 5.1.
    *
    * Given a positive integer \c b and an odd positive integer \c u, returns
    * a non negative integer \c v < 2^\c b such that
    * 1 + <tt>u*v</tt> = 0 (mod 2^\c b).
    *
    * This is the algorithm 5.1 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    * http://cr.yp.to/papers/sf.pdf
    *
    * \param[in] b A positive integer.
    * \param[in] u An odd positive \c mpz_t integer.
    * \return A non negative \c mpz_t integer \c v < 2^\c b such that
    *         1 + <tt>u*v</tt> = 0 (mod 2^\c b).
    */
mpz_t* bern_51(uint32_t b, const mpz_t u);

   /**
    * \brief Daniel J. Bernstein's algorithm 5.3.
    *
    * Given an odd positive integer \c u < 2^\c c and a non negative integer
    * \c x < 2^(\c b + \c c), returns a non negative integer
    * \c r < 2^(\c c + 1) such that \c r*2^\c b = \c x (mod \c u).
    *
    * This is the algorithm 5.3 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    *  http://cr.yp.to/papers/sf.pdf
    *
    * \param[in] b A positive integer.
    * \param[in] u An odd positive \c mpz_t integer.
    * \param[in] x An non negative \c mpz_t integer.
    * \return A non negative \c mpz_t integer \c r such that
    *         \c r*2^\c b = \c x (mod \c u).
    */
mpz_t* bern_53(uint32_t b, const mpz_t u, const mpz_t x);

   /**
    * \brief Daniel J. Bernstein's algorithm 6.3.
    *
    * Given a non negative integer \c x and given the product tree \c tree
    * of a sequence of odd positive integers p_i, returns the integers
    * p_i such that: \c x mod p_i = 0.
    *
    * This is the algorithm 6.3 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    *  http://cr.yp.to/papers/sf.pdf
    *
    * \param[in] x A non negative positive integer.
    * \param[in] tree The product tree of a sequence of odd positive
    *                 integers p_i.
    * \return A pointer to an \c uint32_array_t holding the integers
    *         p_i such that: \c x mod p_i = 0.
    */
uint32_array_t* bern_63(const mpz_t x, const mpz_array_t* const tree);

   /**
    * \brief Daniel J. Bernstein's algorithm 7.1.
    *
    * Given a sequence of odd primes p_j given by \c odd_primes and a
    * set of integers n_i given by \c to_be_factored, determines, for each
    * n_i, the list of odd primes p_j such that (n_i mod p_j = 0) and stores
    * them in <tt>decomp_list</tt>. Each entry in \c decomp_list->data[i]
    * is a pointer to a \c mpz_array_t listing the p_j for the integer
    * <tt>to_be_factored->data[i]</tt>.
    *
    * This is the algorithm 7.1 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    * http://cr.yp.to/papers/sf.pdf
    *
    * \param[out] decomp_list A pointer to the list of matching p_j for
    *                         each n_i.
    * \param[in] to_be_factored A pointer to the set of integers n_i.
    * \param[in] odd_primes A pointer to the set of integers p_j.
    */
void bern_71(
    uint32_array_list_t* const decomp_list,
    const mpz_array_t* const to_be_factored,
    const uint32_array_t* const odd_primes
);

//-----------------------------------------------------------------------------
// Ref. "How to find smooth parts of integers", Daniel J. Bernstein
//       http://cr.yp.to/factorization/smoothparts-20040510.pdf
//-----------------------------------------------------------------------------

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified for SIQS (with a
    * remainder tree only when \c z has more than
    * \c BERN_SIQS_DIRECT_MAX_LIMBS (40) limbs; otherwise, each y_i is reduced
    * modulo \c z directly).
    *
    * Given <tt>z</tt>, the product of prime numbers p_j and the positive
    * integers y_i listed by <tt>cand_yi</tt>, determines the y_i that are
    * {p_j}-smooth and stores them in <tt>smooth_yi</tt>, so that
    * \c smooth_yi->data[i] is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem, particularly
    * to our SIQS implementation where we need to keep track of the
    * coefficient \c a of the polynomial each \c y_i comes from, given by its
    * index in the table of the coefficients used: \c cand_a lists these
    * indexes for the \c y_i, and those of the smooth y_i are appended to
    * \c a_for_smooth_gx.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[out] xi A pointer to the \c x_i associated to the {p_j}-smooth
    *                \c y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth \c y_i integer.
    * \param[out] a_for_smooth_gx The indexes of the \c a of the smooth \c y_i.
    * \param[in] cand_xi A pointer to the list of the \c x_i integers.
    * \param[in] cand_yi A pointer to the list of the \c y_i integers.
    * \param[in] cand_a The indexes of the \c a of the \c y_i.
    * \param[in] z The product of the the \c p_j prime numbers in the factor
    *              base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt_pairs_siqs(
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    uint32_array_t* const a_for_smooth_gx,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const uint32_array_t* const cand_a,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified for SIQS, with large
    * primes variation (with a remainder tree only when \c z has more than
    * \c BERN_SIQS_DIRECT_MAX_LIMBS (40) limbs; otherwise, each y_i is reduced
    * modulo \c z directly).
    *
    * Given <tt>z</tt>, the product of prime numbers p_j and the positive
    * integers y_i listed by <tt>cand_yi</tt>, determines the y_i that are
    * {p_j}-smooth and stores them in <tt>smooth_yi</tt>, so that
    * \c smooth_yi->data[i] is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * Moreover, this function implements the so-called large primes variation.
    * If a given y_i is not {p_j}-smooth but is the product of a prime \c p by
    * a {p_j}-smooth number, it is stored in the hashtable <tt>htable</tt>.
    * Subsequently, if another y_j is the product of a {p_j}-smooth number
    * by the same prime number \c p, then <tt>y_i*y_j/(p^2)</tt> is stored in
    * \c smooth_yi and <tt>x_i*x_j*pinv</tt> is stored in \c xi where
    * \c pinv is the inverse of \c p in Z/<tt>c</tt>Z.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem, particularly
    * to our SIQS implementation where we need to keep track of the
    * coefficient \c a of the polynomial each \c y_i comes from, given by its
    * index in the table of the coefficients used: \c cand_a lists these
    * indexes for the \c y_i, and those of the smooth y_i are appended to
    * \c a_for_smooth_gx.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[in] n The integer to factor.
    * \param[in] htable A pointer to the hashtable used for the large prime
    *                   variation.
    * \param[out] xi A pointer to the \c x_i associated to the {p_j}-smooth
    *                \c y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth \c y_i integer.
    * \param[out] a_for_smooth_yi The indexes of the \c a of the smooth \c y_i.
    * \param[in] cand_xi A pointer to the list of the \c x_i integers.
    * \param[in] cand_yi A pointer to the list of the \c y_i integers.
    * \param[in] cand_a The indexes of the \c a of the \c y_i.
    * \param[in] a_values The coefficients \c a the indexes refer to.
    * \param[in] z The product of the the \c p_j prime numbers in the factor
    *              base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt_pairs_lp_siqs(
    const mpz_t n,
    hashtable_t* const htable,
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    uint32_array_t* const a_for_smooth_yi,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const uint32_array_t* const cand_a,
    const mpz_array_t* const a_values,
    const mpz_t z
);

  /**
    * \brief Daniel J. Bernstein's algorithm 2.1 adapted to be used with a
    * \c smooth_filter_t.
    *
    * No early abort strategy is performed (\c complete_filter_init sets
    * \c filter->nsteps to 0 for this method). Appends to
    * \c filter->accepted_xi and \c filter->accepted_yi the relations of
    * \c filter->candidate_xi and \c filter->candidate_yi whose y is smooth
    * over the primes of \c filter->prod_pj[0]. With \c filter->htable, a
    * relation whose y is smooth but for one tabulated prime is kept there
    * until another one with the same prime comes, and the two are combined.
    *
    * \param filter pointer to the \c smooth_filter_t to use
    * \param step unused
    * \return The number of relations used from the "candidate" arrays.
    */
uint32_t djb_batch_rt(
    smooth_filter_t* const filter,
    unsigned long int step
);

#ifdef __cplusplus
}
#endif

#endif
