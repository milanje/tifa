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
 * \file    lindep.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>

#include "tifa_config.h"
#include <gmp.h>
#if TIFA_USE_GMP_INTERNAL_FUNCS
    #include "gmp-impl.h"
#endif

#include "macros.h"
#include "funcs.h"
#include "lindep.h"
#include "gauss_elim.h"
#include "first_primes.h"
#include "print_error.h"
#include "tifa_alloc.h"

//-----------------------------------------------------------------------------
void fill_matrix_trial_div(binary_matrix_t* const matrix,
                           mpz_array_t* const partially_factored,
                           const mpz_array_t* const to_factor,
                           const uint32_array_t* const factor_base) {

    bool flip_bit = false;

    uint32_t nrows        = matrix->nrows;
    uint32_t nrows_to_add = to_factor->length;

    if ((nrows_to_add + nrows) > matrix->nrows_alloced) {
        nrows_to_add = matrix->nrows_alloced - nrows;
    }

    for (uint32_t i = 0; i < nrows_to_add; i++) {

        mpz_set(partially_factored->data[i], to_factor->data[i]);

        if (mpz_sgn(partially_factored->data[i]) == -1) {
            flip_matrix_bit(nrows, 0, matrix);
            mpz_neg(partially_factored->data[i], partially_factored->data[i]);
        }
        for (uint32_t j = 0; j < factor_base->length; j++) {

            while (0 != mpz_divisible_ui_p(partially_factored->data[i],
                                           factor_base->data[j]) ) {

                mpz_divexact_ui(
                    partially_factored->data[i],
                    partially_factored->data[i],
                    factor_base->data[j]
                );
                flip_bit = !flip_bit;
            }
            if (flip_bit) {
                flip_matrix_bit(nrows, j + 1, matrix);
                flip_bit = false;
            }
        }
        nrows++;
    }
    partially_factored->length = to_factor->length;
    matrix->nrows = nrows;
}
//-----------------------------------------------------------------------------
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
//
// The trial divisions of fill_trial_div_decomp on the limbs of the residue,
// least significant first, for residues of at most TDIV_MAX_LIMBS limbs.
// The odd primes are tested and divided out as by GMP's divexact, with
// pinv = 1/p mod 2^64 computed once per prime: from the least significant
// limb, each quotient limb is (limb - c).pinv and leaves a borrow c, and p
// divides x iff it divides (top limb - c).
//
#define TDIV_MAX_LIMBS 16

static inline bool divisible_odd(const uint64_t* const l, const size_t n,
                                 const uint64_t p, const uint64_t pinv,
                                 const uint64_t plim) {
    uint64_t c = 0;
    for (size_t i = 0; i + 1 < n; i++) {
        const uint64_t borrow = (l[i] < c);
        const uint64_t q      = (l[i] - c) * pinv;
        c = (uint64_t) (((unsigned __int128) q * p) >> 64) + borrow;
    }
    const uint64_t d = (l[n - 1] >= c) ? (l[n - 1] - c) : (c - l[n - 1]);
    return d * pinv <= plim;
}

static inline void divexact_odd(uint64_t* const l, size_t* const n,
                                const uint64_t p, const uint64_t pinv) {
    uint64_t c = 0;
    for (size_t i = 0; i < *n; i++) {
        const uint64_t borrow = (l[i] < c);
        const uint64_t q      = (l[i] - c) * pinv;
        l[i] = q;
        c = (uint64_t) (((unsigned __int128) q * p) >> 64) + borrow;
    }
    while ((*n > 1) && (l[*n - 1] == 0)) {
        (*n)--;
    }
}

static inline uint32_t remove_twos(uint64_t* const l, size_t* const n) {
    size_t z = 0;
    while (l[z] == 0) {
        z++;
    }
    const uint32_t b = (uint32_t) __builtin_ctzll(l[z]);
    for (size_t i = 0; i + z < *n; i++) {
        const uint64_t hi = (i + z + 1 < *n) ? l[i + z + 1] : 0;
        l[i] = (b == 0) ? l[i + z] : ((l[i + z] >> b) | (hi << (64 - b)));
    }
    *n -= z;
    while ((*n > 1) && (l[*n - 1] == 0)) {
        (*n)--;
    }
    return (uint32_t) (64 * z) + b;
}
//
// Montgomery arithmetic modulo an odd n of at most two limbs, R = 2^128, for
// the products of find_factors_decomp: each product leaves a factor 1/R, and
// those of a whole product are taken out at once at its end.
//
typedef struct {
    uint64_t n[2];
    uint64_t ninv;      // -1/n mod 2^64
    uint64_t one[2];    // R mod n
    uint64_t r2[2];     // R^2 mod n
} mont2_t;

static bool mont2_init(mont2_t* const m, const mpz_t n) {
    if ((mpz_size(n) > 2) || mpz_even_p(n)) {
        return false;
    }
    m->n[0] = mpz_getlimbn(n, 0);
    m->n[1] = mpz_getlimbn(n, 1);
    uint64_t inv = m->n[0];
    for (int i = 0; i < 5; i++) {
        inv *= 2 - m->n[0] * inv;
    }
    m->ninv = 0 - inv;

    mpz_t t;
    mpz_init(t);
    mpz_setbit(t, 128);
    mpz_mod(t, t, n);
    m->one[0] = mpz_getlimbn(t, 0);
    m->one[1] = mpz_getlimbn(t, 1);
    mpz_set_ui(t, 0);
    mpz_setbit(t, 256);
    mpz_mod(t, t, n);
    m->r2[0] = mpz_getlimbn(t, 0);
    m->r2[1] = mpz_getlimbn(t, 1);
    mpz_clear(t);

    return true;
}
//
// r = a.b/R mod n, for a < n and a.b < n.R (CIOS).
//
static inline void mont2_mul(uint64_t r[2], const uint64_t a[2],
                             const uint64_t b[2], const mont2_t* const m) {
    uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0;
    for (int i = 0; i < 2; i++) {
        unsigned __int128 c = (unsigned __int128) a[0] * b[i] + t0;
        t0 = (uint64_t) c;
        c  = (c >> 64) + (unsigned __int128) a[1] * b[i] + t1;
        t1 = (uint64_t) c;
        c  = (c >> 64) + t2;
        t2 = (uint64_t) c;
        t3 = (uint64_t) (c >> 64);

        const uint64_t q = t0 * m->ninv;
        c  = ((unsigned __int128) q * m->n[0] + t0) >> 64;
        c += (unsigned __int128) q * m->n[1] + t1;
        t0 = (uint64_t) c;
        c  = (c >> 64) + t2;
        t1 = (uint64_t) c;
        t2 = t3 + (uint64_t) (c >> 64);
    }
    //
    // t < 2n: subtract n once if needed.
    //
    const unsigned __int128 d0 = (unsigned __int128) t0 - m->n[0];
    const uint64_t b0 = (uint64_t) (d0 >> 64) & 1;
    const unsigned __int128 d1 = (unsigned __int128) t1 - m->n[1] - b0;
    const uint64_t b1 = (uint64_t) (d1 >> 64) & 1;
    const bool keep = (t2 == 0) && (b1 != 0);
    r[0] = keep ? t0 : (uint64_t) d0;
    r[1] = keep ? t1 : (uint64_t) d1;
}
//
// x = x.R^k mod n, which takes out the factors 1/R of k products: x times
// R^(k+1), the Montgomery form of R^k, divided by R.
//
static void mont2_unscale(uint64_t x[2], const uint32_t k,
                          const mont2_t* const m) {
    uint64_t c[2] = { m->one[0], m->one[1] };
    uint32_t bit = 1U << 31;
    while ((bit != 0) && ((k & bit) == 0)) {
        bit >>= 1;
    }
    for (; bit != 0; bit >>= 1) {
        mont2_mul(c, c, c, m);
        if (k & bit) {
            mont2_mul(c, c, m->r2, m);
        }
    }
    mont2_mul(x, x, c, m);
}
//
// r = u mod n, for the u = a.x + b of the relations.
//
static inline void mont2_residue(uint64_t r[2], const mpz_t u, const mpz_t n,
                                 const mont2_t* const m, mpz_t scratch) {
    if ((mpz_size(u) > 2) || (mpz_cmpabs(u, n) >= 0)) {
        mpz_mod(scratch, u, n);
        r[0] = mpz_getlimbn(scratch, 0);
        r[1] = mpz_getlimbn(scratch, 1);
        return;
    }
    r[0] = mpz_getlimbn(u, 0);
    r[1] = mpz_getlimbn(u, 1);
    if ((mpz_sgn(u) < 0) && ((r[0] | r[1]) != 0)) {
        const unsigned __int128 d0 = (unsigned __int128) m->n[0] - r[0];
        const uint64_t b0 = (uint64_t) (d0 >> 64) & 1;
        r[1] = m->n[1] - r[1] - b0;
        r[0] = (uint64_t) d0;
    }
}

static void mpz_set_mont2(mpz_t z, const uint64_t v[2]) {
    mp_limb_t* const d = mpz_limbs_write(z, 2);
    d[0] = v[0];
    d[1] = v[1];
    mpz_limbs_finish(z, 2);
}
#endif
//-----------------------------------------------------------------------------
decomp_t* alloc_decomp(const uint32_t nrows, const uint32_t nprimes,
                       const bool sparse) {
    decomp_t* const decomp = tifa_malloc(sizeof(decomp_t));
    decomp->nrows    = nrows;
    decomp->nprimes  = nprimes;
    decomp->dense    = sparse ? NULL : alloc_byte_matrix(nrows, nprimes);
    decomp->alloced  = sparse ? 32 * nrows + 32 : 0;
    decomp->entries  = sparse ? tifa_malloc(decomp->alloced * sizeof(uint64_t))
                              : NULL;
    decomp->nentries = 0;
    decomp->starts   = NULL;
    decomp->sorted   = NULL;
    return decomp;
}
//-----------------------------------------------------------------------------
void free_decomp(decomp_t* const decomp) {
    if (decomp->dense != NULL) {
        free_byte_matrix(decomp->dense);
    }
    free(decomp->entries);
    free(decomp->starts);
    free(decomp->sorted);
    free(decomp);
}
//-----------------------------------------------------------------------------
void append_to_decomp(decomp_t* const decomp, const uint32_t i,
                      const uint32_t j, uint32_t e) {
    while (e != 0) {
        const uint32_t part = MIN(e, 255);
        if (decomp->nentries == decomp->alloced) {
            decomp->alloced *= 2;
            decomp->entries  = tifa_realloc(
                                   decomp->entries,
                                   decomp->alloced * sizeof(uint64_t)
                               );
        }
        decomp->entries[decomp->nentries++] =
            ((uint64_t) i << 32) | ((uint64_t) j << 8) | part;
        e -= part;
    }
}
//-----------------------------------------------------------------------------
void sort_decomp(decomp_t* const decomp) {
    if ((decomp->dense != NULL) || (decomp->starts != NULL)) {
        return;
    }
    uint32_t* const starts = tifa_calloc(decomp->nrows + 1, sizeof(uint32_t));
    for (uint32_t k = 0; k < decomp->nentries; k++) {
        starts[(decomp->entries[k] >> 32) + 1]++;
    }
    for (uint32_t i = 0; i < decomp->nrows; i++) {
        starts[i + 1] += starts[i];
    }
    uint32_t* const next   = tifa_malloc((decomp->nrows + 1) * sizeof(uint32_t));
    uint32_t* const sorted = tifa_malloc((decomp->nentries + 1)
                                         * sizeof(uint32_t));
    memcpy(next, starts, (decomp->nrows + 1) * sizeof(uint32_t));
    for (uint32_t k = 0; k < decomp->nentries; k++) {
        const uint64_t entry = decomp->entries[k];
        sorted[next[entry >> 32]++] = (uint32_t) entry;
    }
    free(next);
    free(decomp->entries);
    decomp->entries = NULL;
    decomp->starts  = starts;
    decomp->sorted  = sorted;
}
//-----------------------------------------------------------------------------
//
// add_to_decomp with the rows of a dense decomp_t loaded once.
//
static inline void add_exponent(unsigned char** const dense,
                                decomp_t* const decomp, const bool sparse,
                                const uint32_t i, const uint32_t j,
                                const uint32_t e) {
    if (sparse) {
        append_to_decomp(decomp, i, j, e);
    } else {
        dense[i][j] += (unsigned char) e;
    }
}
//-----------------------------------------------------------------------------
//
// Inlined once per kind of decomp_t: a test of the kind per exponent added
// cost 0.7% of the time of tifa_factor at 80 bits.
//
static inline __attribute__((always_inline))
void fill_trial_div_rows(binary_matrix_t* const matrix,
                         decomp_t* const decomp,
                         mpz_array_t* const partially_factored,
                         const mpz_array_t* const to_factor,
                         const uint32_array_t* const factor_base,
                         const bool sparse) {

    unsigned char** const dense = sparse ? NULL : decomp->dense->data;
    bool flip_bit = false;

    uint32_t nrows        = matrix->nrows;
    uint32_t nrows_to_add = to_factor->length;

    if ((nrows_to_add + nrows) > matrix->nrows_alloced) {
        nrows_to_add = matrix->nrows_alloced - nrows;
    }
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
    uint64_t* const pinv = tifa_malloc(factor_base->length * sizeof(uint64_t));
    uint64_t* const plim = tifa_malloc(factor_base->length * sizeof(uint64_t));
    for (uint32_t j = 0; j < factor_base->length; j++) {
        const uint64_t p = factor_base->data[j];
        uint64_t inv = p;
        for (int k = 0; k < 5; k++) {
            inv *= 2 - p * inv;
        }
        pinv[j] = inv;
        plim[j] = UINT64_MAX / p;
    }
#endif
    for (uint32_t i = 0; i < nrows_to_add; i++) {

        mpz_ptr const v = partially_factored->data[i];
        mpz_set(v, to_factor->data[i]);

        if (mpz_sgn(v) == -1) {
            flip_matrix_bit(nrows, 0, matrix);
            mpz_neg(v, v);
        }
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
        size_t n = mpz_size(v);
        if ((n != 0) && (n <= TDIV_MAX_LIMBS)) {
            uint64_t l[TDIV_MAX_LIMBS];
            const mp_limb_t* const vl = mpz_limbs_read(v);
            for (size_t k = 0; k < n; k++) {
                l[k] = vl[k];
            }
            //
            // The factor base is sorted: only its first prime can be 2.
            //
            uint32_t j = 0;
            if ((factor_base->length != 0) && (factor_base->data[0] == 2)) {
                uint32_t e;
                if (n == 1) {
                    e = (uint32_t) __builtin_ctzll(l[0]);
                    l[0] >>= e;
                } else {
                    e = remove_twos(l, &n);
                }
                add_exponent(dense, decomp, sparse, i, 0, e);
                if (e & 1) {
                    flip_matrix_bit(nrows, 1, matrix);
                }
                j = 1;
            }
            while (j < factor_base->length) {
                //
                // Once the residue is 1, the primes left cannot divide it.
                // Once it is below p^2, with no factor below p left, it is 1
                // or a prime, found in the factor base by a binary search.
                //
                const uint64_t p = factor_base->data[j];
                if ((n == 1) && (l[0] == 1)) {
                    break;
                }
                if ((n == 1) && (l[0] < p * p)) {
                    const uint32_t t = index_in_sorted_uint32_array(
                                           (uint32_t) l[0], factor_base, j,
                                           factor_base->length
                                       );
                    if ((l[0] <= UINT32_MAX) && (t != NOT_IN_ARRAY)) {
                        add_exponent(dense, decomp, sparse, i, t, 1);
                        flip_matrix_bit(nrows, t + 1, matrix);
                        l[0] = 1;
                    }
                    break;
                }
                //
                // Most primes divide nothing, and a branch per prime
                // mispredicts: test up to 64 primes without branches, then
                // divide out the divisors only.
                //
                const uint32_t end = MIN(j + 64, factor_base->length);
                uint64_t divides = 0;
                if (n == 1) {
                    const uint64_t x = l[0];
                    for (uint32_t k = j; k < end; k++) {
                        divides |= (uint64_t) (x * pinv[k] <= plim[k])
                                   << (k - j);
                    }
                } else {
                    for (uint32_t k = j; k < end; k++) {
                        divides |= (uint64_t) divisible_odd(
                                       l, n, factor_base->data[k], pinv[k],
                                       plim[k]
                                   ) << (k - j);
                    }
                }
                while (divides != 0) {
                    const uint32_t k = j + (uint32_t) __builtin_ctzll(divides);
                    const uint64_t q = factor_base->data[k];
                    divides &= divides - 1;
                    uint32_t e = 0;
                    if (n == 1) {
                        uint64_t x = l[0];
                        do {
                            x *= pinv[k];
                            e++;
                        } while (x * pinv[k] <= plim[k]);
                        l[0] = x;
                    } else {
                        do {
                            divexact_odd(l, &n, q, pinv[k]);
                            e++;
                        } while (divisible_odd(l, n, q, pinv[k], plim[k]));
                    }
                    add_exponent(dense, decomp, sparse, i, k, e);
                    if (e & 1) {
                        flip_matrix_bit(nrows, k + 1, matrix);
                    }
                }
                j = end;
            }
            mp_limb_t* const wl = mpz_limbs_write(v, (mp_size_t) n);
            for (size_t k = 0; k < n; k++) {
                wl[k] = l[k];
            }
            mpz_limbs_finish(v, (mp_size_t) n);
            nrows++;
            continue;
        }
#endif
        for (uint32_t j = 0; j < factor_base->length; j++) {

            uint32_t e = 0;
            while (0 != mpz_divisible_ui_p(v, factor_base->data[j])) {
                mpz_divexact_ui(v, v, factor_base->data[j]);
                e++;
                flip_bit = !flip_bit;
            }
            add_exponent(dense, decomp, sparse, i, j, e);
            if (flip_bit) {
                flip_matrix_bit(nrows, j + 1, matrix);
                flip_bit = false;
            }
        }
        nrows++;
    }
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
    free(pinv);
    free(plim);
#endif
    partially_factored->length = to_factor->length;
    matrix->nrows = nrows;
}
//-----------------------------------------------------------------------------
void fill_trial_div_decomp(binary_matrix_t* const matrix,
                           decomp_t* const decomp,
                           mpz_array_t* const partially_factored,
                           const mpz_array_t* const to_factor,
                           const uint32_array_t* const factor_base
                          ) {
    if (decomp->dense != NULL) {
        fill_trial_div_rows(matrix, decomp, partially_factored, to_factor,
                            factor_base, false);
    } else {
        fill_trial_div_rows(matrix, decomp, partially_factored, to_factor,
                            factor_base, true);
    }
}
//-----------------------------------------------------------------------------
void fill_matrix_from_list(binary_matrix_t* const matrix,
                           const mpz_array_t* const smooth_array,
                           const uint32_array_list_t* const list,
                           const uint32_array_t* const factor_base) {
    mpz_t tmp;
    mpz_init(tmp);
    uint32_t ind_p = 0;

    bool flip_bit = false;

    for (uint32_t i = 0; i < smooth_array->length; i++) {

        mpz_set(tmp, smooth_array->data[i]);
        if (list->data[i] != NULL) {
            //
            // Simple trial divisions by the primes given by the
            // uint32_array_t list->data[i].
            //
            for (uint32_t j = 0; j < list->data[i]->length; j++) {

                if (0 != mpz_divisible_ui_p(tmp, list->data[i]->data[j])) {
                    //
                    // We need the index of this prime in the factor base to
                    // properly flip the matrix's bits...
                    //
                    ind_p = index_in_sorted_uint32_array(
                                list->data[i]->data[j], factor_base,
                                0, factor_base->length
                            );

                    mpz_divexact_ui(tmp, tmp, factor_base->data[ind_p]);

                    flip_bit = !flip_bit;
                }
                while (0 != mpz_divisible_ui_p(tmp, list->data[i]->data[j])) {
                    mpz_divexact_ui(tmp, tmp, factor_base->data[ind_p]);
                    flip_bit = !flip_bit;
                }
                if (flip_bit) {
                    flip_matrix_bit(i, ind_p+1, matrix);
                    flip_bit = false;
                }
            }
        }
    }
    mpz_clear(tmp);
}
//-----------------------------------------------------------------------------
void fill_matrix_from_list_decomp(binary_matrix_t* const matrix,
                                  decomp_t* const decomp,
                                  const mpz_array_t* const smooth_array,
                                  const uint32_array_list_t* const list,
                                  const uint32_array_t* const factor_base) {
    mpz_t tmp;
    mpz_init(tmp);
    uint32_t ind_p = 0;

    bool flip_bit = false;

    for (uint32_t i = 0; i < smooth_array->length; i++) {

        mpz_set(tmp, smooth_array->data[i]);
        if (list->data[i] != NULL) {
            //
            // Simple trial divisions by the primes given by the
            // uint32_array_t list->data[i].
            //
            for (uint32_t j = 0; j < list->data[i]->length; j++) {

                uint32_t e = 0;
                if (0 != mpz_divisible_ui_p(tmp, list->data[i]->data[j])) {
                    //
                    // We need the index of this prime in the factor base to
                    // properly flip the matrix's bits...
                    //
                    ind_p = index_in_sorted_uint32_array(
                                list->data[i]->data[j], factor_base,
                                0, factor_base->length
                            );

                    mpz_divexact_ui(tmp, tmp, factor_base->data[ind_p]);

                    flip_bit = !flip_bit;
                    e++;
                }
                while (0 != mpz_divisible_ui_p(tmp, list->data[i]->data[j])) {
                    mpz_divexact_ui(tmp, tmp, factor_base->data[ind_p]);
                    flip_bit = !flip_bit;
                    e++;
                }
                if (e != 0) {
                    add_to_decomp(decomp, i, ind_p, e);
                }
                if (flip_bit) {
                    flip_matrix_bit(i, ind_p+1, matrix);
                    flip_bit = false;
                }
            }
        }
    }
    mpz_clear(tmp);
}
//-----------------------------------------------------------------------------
uint32_array_list_t*
find_dependencies(binary_matrix_t* const matrix, linalg_method_t method) {

    uint32_array_list_t* relations;
    //
    // There are nrows - rank dependencies: at least nrows - ncols, which is
    // plenty when positive; otherwise there may still be some.
    //
    relations = alloc_uint32_array_list(
                    (matrix->nrows > matrix->ncols)
                    ? matrix->nrows - matrix->ncols : matrix->nrows
                );

    switch (method) {
        case SMART_GAUSS_ELIM:
            gaussian_elim(relations, matrix);
            break;
        case FILTERED_GAUSS_ELIM:
            filtered_gaussian_elim(relations, matrix);
            break;
        default:
            gaussian_elim(relations, matrix);
    }
    return relations;
}
//-----------------------------------------------------------------------------
ecode_t find_factors(mpz_array_t*  const factors,
                     const mpz_t n,
                     const mpz_array_t* const xi_array,
                     const mpz_array_t* const yi_array,
                     const uint32_array_list_t* const dependencies) {

    ecode_t ecode = NO_FACTOR_FOUND;

    mpz_t x;
    mpz_init(x);

    mpz_t y;
    mpz_init(y);
    
    mpz_t rem;
    mpz_init(rem);

    mpz_t gcd;
    mpz_init(gcd);

    uint32_t irow = 0;

    for (uint32_t idep = 0; idep < dependencies->length; idep++) {
        if (factors->length == factors->alloced) {
            break;
        }
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 1);
        //
        // Compute y   = Sqrt(Prod(yi))
        // Compute x   = Prod(xi) + y = Prod(xi) + Sqrt(Prod(yi))
        // Compute gcd = gcd(x, n)
        //
        for (uint32_t irel = 0;
             irel < dependencies->data[idep]->length;
             irel++) {

            irow = dependencies->data[idep]->data[irel];
            mpz_mul(x, x, xi_array->data[irow]);
            mpz_mod(x, x, n);
            //
            // _TO_DO_: We should use a product tree here, but since this
            //          function is only used with CFRAC (where the size of
            //          the residues are bounded) this is not critical...
            //
            mpz_mul(y, y, yi_array->data[irow]);
        }
        mpz_sqrtrem(y, rem, y);
        
        if (0 != mpz_sgn(rem)) {
            //
            // At this stage, if we made no mistake, y should be a product
            // of squares. If that's not the case, something went terribly
            // wrong and it's no use to proceed any further...
            //
            PRINTF_STDERR("\n");
            PRINT_ERROR("What should be Y^2 is not a perfect square!\n");
            PRINT_ERROR("This should not happen - Factorization aborted!\n");
            ecode = FATAL_INTERNAL_ERROR;
            goto clean_and_return;
        }
        mpz_mod(y, y, n);
        mpz_add(x, x, y);
        mpz_gcd(gcd, x, n);

        if ( (0 != mpz_cmp(gcd, n)) && (0 != mpz_cmp_ui(gcd, 1)) ) {
            //
            // gcd is a non trivial factor of n. Add it to our list if it's
            // not already there...
            //
            if (NOT_IN_ARRAY == index_in_mpz_array(gcd, factors)) {
                append_mpz_to_array(factors, gcd);
                ecode = SOME_FACTORS_FOUND;
            }
        }
    }

  clean_and_return:

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(rem);
    mpz_clear(gcd);

    return ecode;
}
//-----------------------------------------------------------------------------
ecode_t find_factors_decomp(mpz_array_t*  const factors,
                            const mpz_t n,
                            const mpz_array_t* const xi_array,
                            const decomp_t* const yi_decomp,
                            const uint32_array_list_t* const dependencies,
                            const uint32_array_t* const factor_base) {

    ecode_t ecode = NO_FACTOR_FOUND;

    mpz_t x;
    mpz_t y;
    mpz_t gcd;

    mpz_init(x);
    mpz_init(y);
    mpz_init(gcd);

    uint32_t irow = 0;

    const uint32_t  nprimes   = yi_decomp->nprimes;
    unsigned int*   decomp    = tifa_malloc(nprimes * sizeof(unsigned int));
    uint32_t*       base_data = factor_base->data;

#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
    //
    // Up to two limbs, the products are Montgomery products, without any of
    // GMP's divisions: they took half of this function at 70 bits.
    //
    mont2_t mod;
    const bool mont = mont2_init(&mod, n);
    uint64_t xm[2];
    uint64_t ym[2];
    uint32_t kx = 0;
    uint32_t ky = 0;
#endif
    for (uint32_t idep = 0; idep < dependencies->length; idep++) {
        if (factors->length == factors->alloced) {
            break;
        }
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 1);
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
        xm[0] = 1;
        xm[1] = 0;
        ym[0] = 1;
        ym[1] = 0;
        kx    = 0;
        ky    = 0;
#endif
        //
        // Compute y   = sqrt(product(yi)) via decomposition of the yi's
        // Compute x   = product(xi)
        // Compute gcd = gcd(x - y, n)
        //
        for (uint32_t i = 0; i < nprimes; i++) {
            decomp[i] = 0;
        }
        for (uint32_t irel = 0;
             irel < dependencies->data[idep]->length;
             irel++) {

            irow = dependencies->data[idep]->data[irel];
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
            if (mont) {
                uint64_t u[2];
                mont2_residue(u, xi_array->data[irow], n, &mod, gcd);
                mont2_mul(xm, xm, u, &mod);
                kx++;
            } else
#endif
            {
                mpz_mul(x, x, xi_array->data[irow]);
                mpz_mod(x, x, n);
            }
            //
            // Get the power of each prime involved in y's decomposition on the
            // factor base
            //
            if (yi_decomp->dense != NULL) {
                const unsigned char* const row = yi_decomp->dense->data[irow];
                for (uint32_t i = 0; i < nprimes; i++) {
                    decomp[i] += row[i];
                }
            } else {
                for (uint32_t k = yi_decomp->starts[irow];
                     k < yi_decomp->starts[irow + 1]; k++) {
                    decomp[yi_decomp->sorted[k] >> 8] +=
                        yi_decomp->sorted[k] & 0xFF;
                }
            }
        }
        //
        // Compute y = sqrt(product(yi))
        //
        //
        // The exponents are small: multiply by the prime powers that fit in
        // a word, and reduce only once y outgrows n.
        //
        for (uint32_t i = 0; i < nprimes; i++) {
            const unsigned long int p = base_data[i];
            uint32_t k = decomp[i] / 2;
            while (k > 0) {
                unsigned long int pk = p;
                k--;
                while ((k > 0) && (pk <= ULONG_MAX / p)) {
                    pk *= p;
                    k--;
                }
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
                if (mont) {
                    const uint64_t b[2] = { pk, 0 };
                    mont2_mul(ym, ym, b, &mod);
                    ky++;
                    continue;
                }
#endif
                mpz_mul_ui(y, y, pk);
                if (mpz_size(y) > mpz_size(n)) {
                    mpz_mod(y, y, n);
                }
            }
        }
#if (GMP_LIMB_BITS == 64) && defined(__SIZEOF_INT128__)
        if (mont) {
            mont2_unscale(xm, kx, &mod);
            mont2_unscale(ym, ky, &mod);
            mpz_set_mont2(x, xm);
            mpz_set_mont2(y, ym);
        }
#endif
        mpz_mod(y, y, n);
        mpz_sub(x, x, y);
        mpz_gcd(gcd, x, n);

        if ( (0 != mpz_cmp(gcd, n)) && (0 != mpz_cmp_ui(gcd, 1)) ) {
            //
            // gcd is a non trivial factor of n. Add it to our list if it's
            // not already there...
            //
            if (NOT_IN_ARRAY == index_in_mpz_array(gcd, factors)) {
                append_mpz_to_array(factors, gcd);
                ecode = SOME_FACTORS_FOUND;
            }
            //
            // Early exit if n is the product of two primes
            //
            if ((factors->length == 1) && MPZ_IS_PRIME(gcd)) {
                mpz_divexact(gcd, n, gcd);
                if (MPZ_IS_PRIME(gcd)) {
                    append_mpz_to_array(factors, gcd);
                    ecode = COMPLETE_FACTORIZATION_FOUND;
                    break;
                }
            }
        }
    }
    free(decomp);
    
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(gcd);

    return ecode;
}
//-----------------------------------------------------------------------------
