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
 * \file    bernsteinisms.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "tifa_config.h"
#include "gmp_utils.h"
#include "bernsteinisms.h"
#include "x_tree.h"
#include "funcs.h"

#include "first_primes.h"
#include "hashtable.h"
#include "macros.h"
#include "tifa_alloc.h"

#define DEBUG 0
#if DEBUG
    #include <stdio.h>
    #define PRINTF(...)     printf(__VA_ARGS__); fflush(stdout); 
    #define GMP_PRINTF(...) gmp_printf(__VA_ARGS__); fflush(stdout); 
#else
    #define PRINTF(...)     /* */
    #define GMP_PRINTF(...) /* */
#endif

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
void bern_63_rec(
    uint32_array_t*,
    const mpz_t,
    const mpz_t[],
    uint32_t nb_nodes,
    uint32_t curnode
);
//------------------------------------------------------------------------------
void bern_71_rec(
    uint32_array_list_t* const,
    uint32_t* const,
    uint32_array_t*,
    const mpz_array_t* const,
    uint32_t
);
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_no_ea(smooth_filter_t* const filter);
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Ref. "How to find small factors of integers", Daniel J. Bernstein
//       http://cr.yp.to/papers/sf.pdf
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
mpz_t* bern_51(uint32_t b, const mpz_t u) {
    //
    // _NOTE_: This is the algorithm 5.1 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are the ones used in the paper.
    //
    mpz_t* v = tifa_malloc(sizeof(mpz_t));
    mpz_init(*v);
    //
    // Step 1
    //
    if (1 == b) {
        mpz_set_ui(*v, 1);
        return v;
    }
    //
    // Step 2
    //
    uint32_t c;
    c = (b & 1) + (b >> 1);

    //
    // Step 3
    //
    mpz_t* v0 = bern_51(c, u);
    //
    // Step 4
    //
    mpz_t u0;
    mpz_t u1;
    mpz_init(u0);
    mpz_init(u1);

    mpz_tdiv_r_2exp(u0, u, c);
    mpz_tdiv_q_2exp(u1, u, c);
    mpz_tdiv_r_2exp(u1, u1, c);
    //
    // Step 5
    //
    mpz_t z;
    mpz_init_set_ui(z, 1);

    mpz_addmul(z, u0, *v0);
    mpz_tdiv_q_2exp(z, z, c);
    mpz_addmul(z, u1, *v0);
    mpz_tdiv_r_2exp(z, z, c);
    //
    // Step 6
    //
    mpz_mul(u0, *v0, z);
    mpz_mul_2exp(u0, u0, c);
    mpz_tdiv_r_2exp(u1, u0, b);
    mpz_add(*v, *v0, u1);
    //
    // Clean up
    //
    mpz_clear(*v0);
    free(v0);
    mpz_clear(u0);
    mpz_clear(u1);
    mpz_clear(z);

    return v;
}
//------------------------------------------------------------------------------
mpz_t* bern_53(uint32_t b, const mpz_t u, const mpz_t x) {
    //
    // _NOTE_: This is the algorithm 5.3 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are the ones used in the paper.
    //
    mpz_t *r = tifa_malloc(sizeof(mpz_t));
    mpz_init(*r);
    //
    // Step 1
    //
    mpz_t *v = bern_51(b, u);
    //
    // Step 2
    //
    mpz_t x0;
    mpz_t x1;
    mpz_init(x0);
    mpz_init(x1);

    mpz_tdiv_q_2exp(x1, x, b);
    mpz_tdiv_r_2exp(x0, x, b);
    //
    // Step 3
    //
    mpz_t q;
    mpz_init(q);

    mpz_mul(q, *v, x0);
    mpz_tdiv_r_2exp(q, q, b);
    //
    // Step 4
    //
    mpz_set(*r, x0);
    mpz_addmul(*r, u, q);
    mpz_tdiv_q_2exp(*r, *r, b);
    mpz_add(*r, *r, x1);
    //
    // Clean up
    //
    mpz_clear(*v);
    free(v);
    mpz_clear(x0);
    mpz_clear(x1);
    mpz_clear(q);

    return r;
}
//------------------------------------------------------------------------------
uint32_array_t* bern_63(const mpz_t x, const mpz_tree_t* const tree) {
    //
    // _NOTE_: This is the algorithm 6.3 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //

    //
    // Bootstrap...
    //
    uint32_array_t* result = tifa_malloc(sizeof(mpz_array_t));

    if (0 == tree->length) {
        result->alloced = 0;
        result->length = 0;
        result->data = NULL;
        return result;
    }
    result->alloced = (tree->length + 1) / 2;
    result->length  = 0;
    result->data    = tifa_malloc(result->alloced*sizeof(uint32_t));
    //
    // Recursion: steps 1 to 7
    //
    bern_63_rec(result, x, (const mpz_t*)tree->data, tree->length, 0);

    return result;
}
//------------------------------------------------------------------------------
void bern_63_rec(uint32_array_t* result, const mpz_t x,
                 const mpz_t tree[],
                 uint32_t nb_nodes, uint32_t curnode) {
    //
    // _NOTE_: This is the algorithm 6.3 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //

    //
    // Step 1
    //
    const mpz_t* const u = &(tree[curnode]);
    //
    // Step 2
    //
    uint32_t c = mpz_sizeinbase(*u, 2);
    uint32_t d = mpz_sizeinbase(x, 2);

    mpz_t *r;
    if (d > (c + 1)) {
        //
        // Step 3
        //
        r = bern_53(d - c, *u, x);

    } else {
        //
        // Step 4
        //
        r = tifa_malloc(sizeof(mpz_t));
        mpz_init_set(*r, x);
    }
    //
    // Step 5
    //
    uint8_t flag = 0;

    if ((1 == nb_nodes)) {
        //
        // The root of the tree has no children
        //
        mpz_t uc;
        mpz_init_set(uc, *u);
        //
        // If (u != 1) : check if r is in {0, u, 2*u, 3*u}
        // If (u == 1) : useless trivial case...
        //
        if (mpz_cmp_ui(*u, 1) != 0) {

            if ((mpz_cmp_ui(*r, 0) == 0) || (mpz_cmp(*r, uc) == 0)) {
                flag++;
            } else {
                mpz_add(uc, uc, *u);
                if ((mpz_cmp(*r, uc) == 0)) {
                    flag++;
                } else {
                    mpz_add(uc, uc, *u);
                    if ((mpz_cmp(*r, uc) == 0)) {
                        flag++;
                    }
                }
            }
        }
        if (flag != 0) {
            result->data[result->length] = mpz_get_ui(*u);
            result->length++;
        }
        mpz_clear(uc);

    } else {
        //
        // Step 6 and 7
        //
        bern_63_rec(result, x, tree, (nb_nodes - 1)/2, 2*curnode + 1);
        bern_63_rec(result, x, tree, (nb_nodes - 1)/2, 2*curnode + 2);
    }
    mpz_clear(*r);
    free(r);
}
//------------------------------------------------------------------------------
void bern_71(uint32_array_list_t* const decomp_list,
             const mpz_array_t* const to_be_factored,
             const uint32_array_t* const odd_primes) {
    //
    // _NOTE_: This is the algorithm 7.1 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //

    //
    // Bootstrap. Step 1
    //
    if (0U == to_be_factored->length) {
        for (uint32_t i = 0; i < to_be_factored->length; i++) {
            add_entry_in_uint32_array_list(NULL, decomp_list);
        }
        return;
    }
    //
    // Step 2 - processed only once.
    //
    mpz_tree_t* x = prod_tree(to_be_factored);
    uint32_t nb_int_to_factor = to_be_factored->length;
    //
    // Recursion: steps 3 to 8
    //
    bern_71_rec(decomp_list, &nb_int_to_factor, (uint32_array_t*)odd_primes,
                x, 0);

    free_mpz_tree(x);
}
//------------------------------------------------------------------------------
void bern_71_rec(uint32_array_list_t* const decomp_list,
                 uint32_t * const nb_int_to_factor,
                 uint32_array_t*  odd_primes,
                 const mpz_array_t* const x, uint32_t curnode) {

    //
    // _TO_DO_: Lots of redundant computations here! Fix this!
    //

    //
    // Step 3
    //
    mpz_tree_t *T = prod_tree_ui((const uint32_array_t*)odd_primes);
    //
    // Step 4
    //
    uint32_array_t *P = bern_63(x->data[curnode], (const mpz_array_t*)T);

    if (0 == P->length) {
        //
        // Step 1 (Yes, indeed! Why call bern_71_rec another time since we
        //         can handle it now...)
        //
        // The list of primes such that x->data[curnode] mod p = 0 is empty.
        // Thus, for every node i in the subtree of root curnode, the list
        // such that subtree[i] mod p = 0 will be empty. Consequently, we have
        // to add in our list as many null lists as there is nodes in the
        // subtree of root curnode.
        //
        // Instead of invoking recursively bern_71_rec, compute this number of
        // nodes by finding the depth of curnode is the global x tree...
        //
        uint8_t depth_curnode = 0;
        while ((1U<<depth_curnode) <= (curnode+1)) {
            depth_curnode++;
        }
        depth_curnode--;

        uint32_t nb_leafs = (x->length + 1)/2;
        uint32_t nb_leafs_subtree = nb_leafs/(1<<depth_curnode);
        //
        // ... and add as many null list as required...
        //
        for (uint32_t i = 0; i < nb_leafs_subtree; i++) {
            if (*nb_int_to_factor != 0U) {
                add_entry_in_uint32_array_list(NULL, decomp_list);
                (*nb_int_to_factor)--;
            }
        }
        free_uint32_array(P);
        free_mpz_tree(T);

        return;
    }
    if (curnode >= (x->length/2)) {
        //
        // i.e. curnode is a leaf...
        //
        // Now, check *nb_int_to_factor as there may be useless nodes (i.e.
        // nodes with value 1) in the product tree of all the integers to
        // factor...
        //
        if (*nb_int_to_factor != 0U) {
            //
            // Step 5.
            //
            add_entry_in_uint32_array_list(P, decomp_list);
            (*nb_int_to_factor)--;
        } else {
            free_uint32_array(P);
        }
    } else {
        //
        // Step 6, 7 and 8
        //
        bern_71_rec(decomp_list, nb_int_to_factor, P, x, 2*curnode+1);
        bern_71_rec(decomp_list, nb_int_to_factor, P, x, 2*curnode+2);

        free_uint32_array(P);
    }
    free_mpz_tree(T);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Ref. "How to find smooth part of integers", Daniel J. Bernstein
//       http://cr.yp.to/factorization/smoothparts-20040510.pdf
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//
// For the SIQS batches. While z has at most BERN_SIQS_DIRECT_MAX_LIMBS limbs
// (factor bases of up to about 250 primes), z mod y for each candidate y is
// cheaper computed directly (linear in the size of z) than through the
// product and remainder trees.
//
#if !defined(BERN_SIQS_DIRECT_MAX_LIMBS)
#define BERN_SIQS_DIRECT_MAX_LIMBS 40
#endif

#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
//
// For y = 2^t.yo of one limb, yo odd, and r = z mod y: r^(2^e) mod yo up to
// a unit, by Montgomery squarings of r left unconverted (each multiplies the
// result by 2^-64). It is 0 iff y is smooth (z holds 2), and its gcd with
// yo is the odd part of the smooth part of y.
//
static inline uint64_t pow2e_mod_odd(const uint64_t r, const uint64_t yo,
                                     const uint32_t e) {
    uint64_t ninv = yo;
    for (int i = 0; i < 5; i++) {
        ninv *= 2 - yo * ninv;
    }
    uint64_t x = r % yo;
    for (uint32_t i = 0; i < e; i++) {
        const unsigned __int128 t = (unsigned __int128) x * x;
        const uint64_t lo = (uint64_t) t;
        const uint64_t hi = (uint64_t) (t >> 64);
        const uint64_t qn = (uint64_t) (((unsigned __int128) (lo * ninv) * yo)
                                        >> 64);
        x = (hi < qn) ? (hi - qn + yo) : (hi - qn);
    }
    return x;
}

//
// As pow2e_mod_odd, for yo = (y0, y1) of two limbs (y1 != 0) and r < yo:
// Montgomery squarings with 2^128 (CIOS), into x.
//
static inline void pow2e_mod_odd_2(uint64_t x[2], const uint64_t r[2],
                                   const uint64_t yo[2], const uint32_t e) {
    uint64_t ninv = yo[0];
    for (int i = 0; i < 5; i++) {
        ninv *= 2 - yo[0] * ninv;
    }
    ninv = 0 - ninv;
    x[0] = r[0];
    x[1] = r[1];
    for (uint32_t k = 0; k < e; k++) {
        uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0;
        for (int i = 0; i < 2; i++) {
            unsigned __int128 c = (unsigned __int128) x[0] * x[i] + t0;
            t0 = (uint64_t) c;
            c  = (c >> 64) + (unsigned __int128) x[1] * x[i] + t1;
            t1 = (uint64_t) c;
            c  = (c >> 64) + t2;
            t2 = (uint64_t) c;
            t3 = (uint64_t) (c >> 64);

            const uint64_t m = t0 * ninv;
            c  = ((unsigned __int128) m * yo[0] + t0) >> 64;
            c += (unsigned __int128) m * yo[1] + t1;
            t0 = (uint64_t) c;
            c  = (c >> 64) + t2;
            t1 = (uint64_t) c;
            t2 = t3 + (uint64_t) (c >> 64);
        }
        //
        // t < 2.yo: subtract yo once if needed.
        //
        const unsigned __int128 d0 = (unsigned __int128) t0 - yo[0];
        const uint64_t b0 = (uint64_t) (d0 >> 64) & 1;
        const unsigned __int128 d1 = (unsigned __int128) t1 - yo[1] - b0;
        const uint64_t b1 = (uint64_t) (d1 >> 64) & 1;
        const bool keep = (t2 == 0) && (b1 != 0);
        x[0] = keep ? t0 : (uint64_t) d0;
        x[1] = keep ? t1 : (uint64_t) d1;
    }
}

static inline uint64_t gcd_u64(uint64_t a, uint64_t b) {
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    const int shift = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) {
            const uint64_t t = a;
            a = b;
            b = t;
        }
        b -= a;
    } while (b != 0);
    return a << shift;
}
#endif
//
// Whether the candidate y is smooth over the primes of z, given r = z mod |y|
// (if rleaf is NULL, it is computed from z). If not and nsp is not NULL, nsp
// is set to |y| over its smooth part. yk and yo are scratch.
//
static bool siqs_smooth(const mpz_t y, const mpz_t z, const mpz_t rleaf,
                        mpz_t yk, mpz_t yo, mpz_t nsp) {
    const uint32_t msb = mpz_sizeinbase(y, 2) - 1;
    const uint32_t e   = most_significant_bit(msb) + 1;
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
    if (mpz_size(y) == 1) {
        const uint64_t ya = mpz_getlimbn(y, 0);
        const uint64_t r  = (rleaf != NULL)
                          ? mpz_getlimbn(rleaf, 0)
                          : mpn_mod_1(mpz_limbs_read(z), mpz_size(z), ya);
        const uint64_t yo = ya >> __builtin_ctzll(ya);
        const uint64_t x  = pow2e_mod_odd(r, yo, e);
        if (x == 0) {
            return true;
        }
        if (nsp != NULL) {
            mpz_set_ui(nsp, (unsigned long int) (yo / gcd_u64(x, yo)));
        }
        return false;
    }
    if (mpz_size(y) == 2) {
        if (rleaf != NULL) {
            mpz_set(yk, rleaf);
        } else {
            mpz_tdiv_r(yk, z, y);
            mpz_abs(yk, yk);
        }
        //
        // yk = r mod yo, with yo the odd part of |y|.
        //
        mpz_abs(yo, y);
        mpz_tdiv_q_2exp(yo, yo, mpz_scan1(yo, 0));
        mpz_tdiv_r(yk, yk, yo);
        uint64_t x[2];
        if (mpz_size(yo) == 1) {
            const uint64_t yo1 = mpz_getlimbn(yo, 0);
            x[0] = pow2e_mod_odd(mpz_getlimbn(yk, 0), yo1, e);
            x[1] = 0;
        } else {
            const uint64_t r2[2]  = { mpz_getlimbn(yk, 0), mpz_getlimbn(yk, 1) };
            const uint64_t yo2[2] = { mpz_getlimbn(yo, 0), mpz_getlimbn(yo, 1) };
            pow2e_mod_odd_2(x, r2, yo2, e);
        }
        const bool smooth = (x[0] == 0) && (x[1] == 0);
        if (!smooth && (nsp != NULL)) {
            mpz_import(yk, 2, -1, sizeof(uint64_t), 0, 0, x);
            mpz_gcd(yk, yk, yo);
            mpz_divexact(nsp, yo, yk);
        }
        return smooth;
    }
#endif
    if (rleaf != NULL) {
        mpz_set(yk, rleaf);
    } else {
        mpz_tdiv_r(yk, z, y);
        mpz_abs(yk, yk);
    }
    mpz_powm_ui(yk, yk, 1<<e, y);
    if (0 == mpz_sgn(yk)) {
        return true;
    }
    if (nsp != NULL) {
        //
        // The gcd is the smooth part of y.
        //
        mpz_gcd(nsp, yk, y);
        mpz_divexact(nsp, y, nsp);
        mpz_abs(nsp, nsp);
    }
    return false;
}
//------------------------------------------------------------------------------
uint32_t bern_21_rt_pairs_siqs(mpz_array_t* const xi,
                               mpz_array_t* const smooth_yi,
                               uint32_array_t* const a_for_smooth_gx,
                               const mpz_array_t* const cand_xi,
                               const mpz_array_t* const cand_yi,
                               const uint32_array_t* const cand_a,
                               const mpz_t z) {
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* rtree = NULL;
    if (mpz_size(z) > BERN_SIQS_DIRECT_MAX_LIMBS) {
        mpz_tree_t* pyitree = prod_tree(cand_yi);
        rtree = rem_tree(z, pyitree);
        free_mpz_tree(pyitree);
    }
    //
    // Step 3
    //
    mpz_t yk, yo;
    mpz_init(yk);
    mpz_init(yo);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        const mpz_srcptr rleaf = (rtree != NULL)
                               ? rtree->data[rtree->length/2 + k] : NULL;
        if (siqs_smooth(cand_yi->data[k], z, rleaf, yk, yo, NULL)) {
			//
            // We found a good pair (cand_yi->data[k], cand_xi->data[k])!
            // Add cand_yi->data[k] is the smooth_yi array and add
            // cand_xi->data[k] is the xi array.
            //
            mpz_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_set(xi->data[xi->length], cand_xi->data[k]);
            a_for_smooth_gx->data[a_for_smooth_gx->length] = cand_a->data[k];
            smooth_yi->length++;
            xi->length++;
            a_for_smooth_gx->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        }
    }
    mpz_clear(yk);
    mpz_clear(yo);

    if (rtree != NULL) {
        free_mpz_tree(rtree);
    }
    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_rt_pairs_lp_siqs(const mpz_t n,
                                  hashtable_t* const htable,
                                  mpz_array_t* const xi,
                                  mpz_array_t* const smooth_yi,
                                  uint32_array_t* const a_for_smooth_gx,
                                  const mpz_array_t* const cand_xi,
                                  const mpz_array_t* const cand_yi,
                                  const uint32_array_t* const cand_a,
                                  const mpz_array_t* const a_values,
                                  const mpz_t z) {

    //
    // _NOTE_: The coefficient "a" of the SIQS polynomial of each relation is
    //         stored along with it.
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
        
    uint32_t retval = cand_xi->length;
    
    //printf("In batch: %i candidates\n", (int)retval);
    
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* rtree = NULL;
    if (mpz_size(z) > BERN_SIQS_DIRECT_MAX_LIMBS) {
        mpz_tree_t* pyitree = prod_tree(cand_yi);
        rtree = rem_tree(z, pyitree);
        free_mpz_tree(pyitree);
    }
    //
    // Step 3
    //
    mpz_t yk, yo;
    mpz_init(yk);
    mpz_init(yo);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);
    
    const uint32_t cand_yi_length = cand_yi->length;

    mpz_t* candidate_y = &(cand_yi->data[0]);
    mpz_t* candidate_x = &(cand_xi->data[0]);
    const uint32_t* candidate_a = &(cand_a->data[0]);
    mpz_t* rtreeptr    = (rtree != NULL) ? &(rtree->data[rtree->length/2])
                                         : NULL;
    mpz_t* smooth_y    = &(smooth_yi->data[smooth_yi->length]);
    mpz_t* x           = &(xi->data[xi->length]);
    uint32_t* matching_a = &(a_for_smooth_gx->data[a_for_smooth_gx->length]);
    
    for (uint32_t k = 0U; k < cand_yi_length; k++) {
        //
        // _NOTE_: nsp is set to cand_yi->data[k] over its smooth part when
        //         it is not smooth.
        //
        if (siqs_smooth(*candidate_y, z, (rtreeptr != NULL) ? *rtreeptr : NULL,
                        yk, yo, nsp)) {
			//
            // We found a good pair (cand_yi->data[k], cand_xi->data[k])!
            // Add cand_yi->data[k] is the smooth_yi array and add
            // cand_xi->data[k] is the xi array.
            //
            mpz_set(*smooth_y, *candidate_y);
            mpz_set(*x, *candidate_x);
            *matching_a = *candidate_a;
            
            smooth_yi->length++;
            xi->length++;
            a_for_smooth_gx->length++;
            
            smooth_y++;
            x++;
            matching_a++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation: is cand_yi->data[k] the
            // product of its smooth part by a large prime?
            //
            if (mpz_cmp_ui(nsp, 15000000) <= 0 ) {
                uint32_t nsp_ui = mpz_get_ui(nsp);
                
                //
                // cand_yi->data[k] is the product of a smooth number with
                // the prime number first_primes_array->data[ind]
                //
                uint32_t *key = tifa_malloc(sizeof(uint32_t));
                *key  = nsp_ui;
                
                //
                // _NOTE_: Actually it may be possible to leave the pair
                //         in the hashtable...
                //
                // _TO_DO: Try to keep the entry in the hashtable...
                //
                //mpz_pair_t *found = remove_entry_in_hashtable(htable,
                //key);
                
                mpz_pair_t *found = get_entry_in_hashtable(htable,key);
                
                if (found != NULL) {
                    //
                    // The hashtable did contain another entry for this
                    // particular prime number: we can compute a new
                    // pair to add in our arrays xi and smooth_yi
                    //
                    if (0 == mpz_cmpabs(*candidate_x, found->x)) {
                        //
                        // The same relation, found again with another
                        // polynomial: pairing it with itself would only
                        // give a trivial relation.
                        //
                        free(key);
                        goto next_candidate;
                    }

                    *matching_a = *candidate_a;
                
                    mpz_mul(*smooth_y, *candidate_y, found->y);
                    mpz_mul(*x, *candidate_x, found->x);
                    
                    mpz_divexact_ui(*smooth_y, *smooth_y, nsp_ui);
                    mpz_divexact_ui(*smooth_y, *smooth_y, nsp_ui);
                
                    if (0 != mpz_invert(nsp_inv, nsp, n)) {
                        mpz_mul(*x, *x, nsp_inv);
                    } else {
                        //
                        // This should not happen: since nsp is a prime,
                        // it has an inverse in Z/nZ... unless n is a
                        // multiple of nsp! In such a (rare) case, we have
                        // found a factor of n, but handling it is quite
                        // cumbersome, so we just ignore it...
                        //
                        free(key);
                        //clear_mpz_pair(found);
                        //free(found);
                        goto next_candidate;
                    }
                    smooth_yi->length++;
                    xi->length++;
                    a_for_smooth_gx->length++;
                    
                    smooth_y++;
                    x++;
                    matching_a++;
                
                    free(key);
                    //clear_mpz_pair(found);
                    //free(found);
                
                    if (smooth_yi->length == smooth_yi->alloced) {
                        retval = k + 1;
                        break;
                    }
                } else {
                    //
                    // Add this pair in the hashtable with the position
                    // of the prime in first_primes_array as key...
                    //
                    mpz_pair_t* pair = tifa_malloc(sizeof(mpz_pair_t));
                    mpz_init_set(pair->x, *candidate_x);
                    mpz_init_set(pair->y, *candidate_y);
                    mpz_mul(pair->y, pair->y, a_values->data[*candidate_a]);
                    add_entry_in_hashtable(htable, (void*)key, (void*)pair);
                    //
                    // _WARNING_: Do not free key and pair as they are now
                    //            referenced by the hashtable!
                    //
                }
            }
        }
      next_candidate:
        candidate_y++;
        candidate_x++;
        candidate_a++;
        if (rtreeptr != NULL) {
            rtreeptr++;
        }
    }
    mpz_clear(yk);
    mpz_clear(yo);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    if (rtree != NULL) {
        free_mpz_tree(rtree);
    }
    return retval;
}
//------------------------------------------------------------------------------
uint32_t djb_batch_rt(smooth_filter_t* const filter, unsigned long int step) {
    //
    // complete_filter_init sets nsteps to 0 for this method: no early abort.
    //
    (void) step;
    return djb_batch_rt_no_ea(filter);
}
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
uint32_t djb_batch_rt_no_ea(smooth_filter_t* const filter) {
        
    mpz_array_t* const in_xi  = filter->candidate_xi;
    mpz_array_t* const in_yi  = filter->candidate_yi;
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    uint32_t next;
        
    hashtable_t* const htable = filter->htable;
    mpz_ptr n = filter->n;
    
    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }
    uint32_t retval = in_xi->length;
    //
    // Step 2
    //
    mpz_tree_t* const pyitree = prod_tree(in_yi);
    mpz_tree_t* const rtree   = rem_tree(filter->prod_pj[0], pyitree);
    
    free_mpz_tree(pyitree);
    //
    // Step 3
    //
    uint32_t e   = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = in_yi->length - 1; k != UINT32_MAX; k--) {
        msb = mpz_sizeinbase(in_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;
            
        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, in_yi->data[k]);
        
        next = acc_yi->length;    
        //
        // _NOTE_: in_yi->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether in_yi->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {

            mpz_set(acc_yi->data[next], in_yi->data[k]);
            mpz_set(acc_xi->data[next], in_xi->data[k]);
        
            acc_xi->length++;
            acc_yi->length++;
                    
            if (acc_yi->length == acc_yi->alloced) {
                retval = in_yi->length - k;
                break;
            }
        } else {            
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, in_yi->data[k]);
    	    //
    	    // gcd is now the smooth part of in_yi->data[k].
    	    // Is in_yi->data[k] the product of gcd by a large prime?
    	    //
            mpz_divexact(nsp, in_yi->data[k], gcd);
            
            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0 ) {
            
                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (in_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (in_yi->data[k] / gcd) is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length
                                );
                
                if ((htable != NULL) && (ind != NOT_IN_ARRAY)) {
                    //
                    // in_yi->data[k] is the product of a smooth number with
                    // the prime number first_primes_array->data[ind]
                    //
                    uint32_t *key = tifa_malloc(sizeof(uint32_t));
                    *key  = ind;
            
                    //
                    // _NOTE_: Actually it may be possible to leave the pair
                    //         in the hashtable...
                    //
                    // _TO_DO: Try to keep the entry in the hashtable...
                    //

                //mpz_pair_t *found = remove_entry_in_hashtable(htable, key);
      
                    mpz_pair_t *found = get_entry_in_hashtable(htable,key);
            
                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our accepted_* arrays.
                        //
                        mpz_mul(acc_yi->data[next], in_yi->data[k], found->y);
                        mpz_mul(acc_xi->data[next], in_xi->data[k], found->x);
                        
                        mpz_divexact_ui(
                            acc_yi->data[next],
                            acc_yi->data[next],
                            nsp_ui
                        );
                        mpz_divexact_ui(
                            acc_yi->data[next],
                            acc_yi->data[next],
                            nsp_ui
                        );
            
                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(
                                acc_xi->data[next], 
                                acc_xi->data[next],
                                nsp_inv
                            );
                        
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            
                            //clear_mpz_pair(found);
                            //free(found);
                            
                            continue;
                        }
                        acc_yi->length++;
                        acc_xi->length++;
                        
                        free(key);
                    
                        //clear_mpz_pair(found);
                        //free(found);
            
                        if (acc_yi->length == acc_yi->alloced) {
                            retval = in_yi->length - k;
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = tifa_malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, in_xi->data[k]);
                        mpz_init_set(pair->y, in_yi->data[k]);
                        add_entry_in_hashtable(htable, (void*)key, (void*)pair);
                        //
                        // _WARNING_: Do not free key and pair as they are now
                        //            referenced by the hashtable!
                        //
                    }
                }
            }
        }
    }
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    free_mpz_tree(rtree);

    return retval;   
}
//------------------------------------------------------------------------------
