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
 * \file    siqs_poly.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Structure and functions related to the polynomials used in
 *        the SIQS algorithm.
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include "macros.h"
#include "print_error.h"
#include "siqs_poly.h"

#if !defined(__TIMING__)
    #define __TIMING__ 1
#endif

#if !defined(__PREFIX__)
    #define __PREFIX__ "siqs: "
#endif

#include "timer.h"
#include "tifa_config.h"

//
// Set this to non zero to perform "safer" modular arithmetic. It involves
// extra, most of the time (always?) unneeded modular reduction.
//
#define SAFE_MODULO 0

//-----------------------------------------------------------------------------
uint32_t determine_nprimes(uint32_t bitsize);
static void choose_a(siqs_poly_t* const poly);
ecode_t  full_poly_init(siqs_poly_t* const poly);
ecode_t  fast_poly_init(siqs_poly_t* const poly);
//-----------------------------------------------------------------------------
siqs_poly_t* alloc_siqs_poly(mpz_t target_a,
                             mpz_t n,
                             uint32_array_t* const factor_base,
                             uint32_array_t* const sqrtm_pi)
{    
    siqs_poly_t* const poly = malloc(sizeof(siqs_poly_t));
        
    poly->nfullpolyinit = 0;
    poly->nprimes_in_a  = determine_nprimes(mpz_sizeinbase(target_a, 2));
    poly->idx_of_a      = malloc(poly->nprimes_in_a * sizeof(uint32_t));
    
    poly->factor_base = factor_base;
    poly->sqrtm_pi    = sqrtm_pi;
        
    poly->approximer = alloc_approximer(
                           target_a, factor_base, poly->nprimes_in_a
                       );
        
    mpz_init_set(poly->n, n);
    mpz_init(poly->a);
    mpz_init(poly->b);
    mpz_init(poly->c);
    
    poly->npolys = 1 << (poly->nprimes_in_a - 1);
    poly->polyno = 1;
    
    poly->sol1   = alloc_int32_array(poly->factor_base->length);
    poly->sol2   = alloc_int32_array(poly->factor_base->length);
    poly->Bl     = alloc_mpz_array(poly->nprimes_in_a);
    poly->Bainv2 = malloc(poly->factor_base->length * sizeof(uint32_t*));
    poly->Bainv2[0] = malloc(  poly->factor_base->length
                             * poly->nprimes_in_a * sizeof(uint32_t));
    for (uint32_t i = 1; i < poly->factor_base->length; i++) {
        poly->Bainv2[i] = poly->Bainv2[0] + i * poly->nprimes_in_a;
    }
    choose_a(poly);
    
    full_poly_init(poly);

    return poly;
}
//-----------------------------------------------------------------------------
void
free_siqs_poly(siqs_poly_t* poly)
{
    if (poly != NULL) {
        free_mpz_array(poly->Bl);
        
        free(poly->Bainv2[0]);
        free(poly->Bainv2);
                
        mpz_clear(poly->a);
        mpz_clear(poly->b);
        mpz_clear(poly->c);
        mpz_clear(poly->n);
        
        free_approximer(poly->approximer);        
        free_int32_array(poly->sol1);
        free_int32_array(poly->sol2);

        free(poly->idx_of_a);

        free(poly);
    }
}
//-----------------------------------------------------------------------------
uint32_t
determine_nprimes(uint32_t bitsize)
{   
    //
    // Determines the number of factors in the leading coefficient according
    // to the size 'bitsize' of the target in bits.
    //
    if (bitsize < 15)  { return 2; }
    if (bitsize < 32)  { return 3; }
    if (bitsize < 44)  { return 4; }
    if (bitsize < 59)  { return 5; }
    if (bitsize < 75)  { return 6; }
    if (bitsize < 88)  { return 7; }
    
    return (8 + (bitsize - 88) / 11);
}
//-----------------------------------------------------------------------------
ecode_t
update_polynomial(siqs_poly_t* const poly)
{   
    //
    // Switch to another polynomial using another 'a' coefficient if needed
    //
    ecode_t ecode;
    
    if (poly->polyno != poly->npolys) {    
        ecode = fast_poly_init(poly);
    } else {
        choose_a(poly);
        ecode = full_poly_init(poly);
    }    
    return ecode;
}
//-----------------------------------------------------------------------------
//
// Roots and Bainv2 entries stay in [0, p), so one conditional correction
// replaces a division when switching to the next polynomial.
//
static inline int32_t add_mod(int32_t r, uint32_t d, uint32_t p) {
    uint32_t s = (uint32_t) r + d;
    return (int32_t) ((s >= p) ? (s - p) : s);
}
static inline int32_t sub_mod(int32_t r, uint32_t d, uint32_t p) {
    int32_t s = r - (int32_t) d;
    return (s < 0) ? (s + (int32_t) p) : s;
}
//-----------------------------------------------------------------------------
static void
choose_a(siqs_poly_t* const poly)
{
    //
    // A prime dividing kn has the single root 0, which makes its B_l zero:
    // successive b's would then repeat, sieving the same polynomials again
    // and filling the matrix with duplicate relations (hence only trivial
    // dependencies). Skip any 'a' with such a prime. Give up after one pass
    // over all the subsets the approximer can produce.
    //
    const uint32_t* const sqrtm = poly->sqrtm_pi->data;
    const uint32_t ntries = poly->approximer->nsubsets_odd;

    for (uint32_t itry = 0; itry <= ntries; itry++) {
        random_approximation(poly->approximer, poly->a, poly->idx_of_a);

        bool divides_kn = false;
        for (uint32_t j = 0; j < poly->nprimes_in_a; j++) {
            if (sqrtm[poly->idx_of_a[j]] == 0) {
                divides_kn = true;
                break;
            }
        }
        if (!divides_kn) {
            return;
        }
    }
}
//-----------------------------------------------------------------------------
//
// The inverse of a (0 < a < p) modulo the prime p, by the extended Euclidean
// algorithm: several times faster here than mpz_invert.
//
static inline uint32_t inverse_mod(const uint32_t a, const uint32_t p) {
    uint32_t r0 = p;
    uint32_t r1 = a;
    int64_t  t0 = 0;
    int64_t  t1 = 1;
    while (r1 != 0) {
        const uint32_t q = r0 / r1;
        const uint32_t r = r0 - q * r1;
        r0 = r1;
        r1 = r;
        const int64_t t = t0 - (int64_t) q * t1;
        t0 = t1;
        t1 = t;
    }
    return (uint32_t) ((t0 < 0) ? t0 + p : t0);
}
//
// a, b and the Bl when they fit in 64 bits (n up to about 140 bits): their
// residues modulo the primes of the base then need no GMP call.
//
#define MAX_NPRIMES_IN_A_U64 8
typedef struct {
    bool     fit;
    uint64_t a;
    uint64_t b;
    uint64_t Bl[MAX_NPRIMES_IN_A_U64];
} poly_u64_t;
//
// For the i-th prime p of the base, not dividing a: the 2.Bj.inv(a) mod p of
// Bainv2[i], and the two solutions to g_{a,b}(x) = 0 mod p.
//
static ecode_t init_prime(siqs_poly_t* const poly, const uint32_t i,
                          const uint32_t curprime, const poly_u64_t* const w,
                          mpz_t inv, mpz_t mpzp) {
    const uint32_t nprimes_in_a = poly->nprimes_in_a;
    uint32_t* const bainv2 = poly->Bainv2[i];

    unsigned long int ainv;
    unsigned long int bmodp;
    uint64_t gamma;

    if (w->fit) {
        const uint32_t amodp = (uint32_t) (w->a % curprime);
        if (amodp == 0) {
            PRINTF_STDERR("\n");
            PRINT_ERROR("%Zd has no inverse mod %"PRIu32"!\n", poly->a,
                        curprime);
            PRINT_ERROR("This should not happen - Factorization aborted!\n");
            return FATAL_INTERNAL_ERROR;
        }
        ainv  = inverse_mod(amodp, curprime);
        bmodp = (unsigned long int) (w->b % curprime);
        for (uint32_t j = 0; j < nprimes_in_a; j++) {
            gamma  = 2 * (w->Bl[j] % curprime);
            gamma *= ainv;
            gamma %= curprime;
            bainv2[j] = (uint32_t) gamma;
        }
    } else {
        mpz_set_ui(mpzp, curprime);
        if (0 == mpz_invert(inv, poly->a, mpzp)) {
            PRINTF_STDERR("\n");
            PRINT_ERROR("%Zd has no inverse mod %"PRIu32"!\n", poly->a,
                        curprime);
            PRINT_ERROR("This should not happen - Factorization aborted!\n");
            return FATAL_INTERNAL_ERROR;
        }
        ainv  = mpz_get_ui(inv);
        bmodp = mpz_fdiv_ui(poly->b, curprime);
        for (uint32_t j = 0; j < nprimes_in_a; j++) {
            gamma  = 2 * mpz_fdiv_ui(poly->Bl->data[j], curprime);
            gamma *= ainv;
            gamma %= curprime;
            bainv2[j] = (uint32_t) gamma;
        }
    }
    //
    // The two solutions to the g_{a,b}(x) = 0 mod p congruence.
    //
    int64_t tmp64;

    tmp64  = (int64_t)(poly->sqrtm_pi->data[i]);
    tmp64 -= (int64_t)bmodp;
#if SAFE_MODULO
    tmp64 %= curprime;
#endif
    tmp64 *= ainv;
    tmp64 %= curprime;
    if (tmp64 < 0) {
        tmp64 += curprime;
    }
    poly->sol1->data[i] = (int32_t) tmp64;

    tmp64  = 0;
    tmp64 -= (int64_t)(poly->sqrtm_pi->data[i]);
    tmp64 -= (int64_t)bmodp;
#if SAFE_MODULO
    tmp64 %= curprime;
#endif
    tmp64 *= ainv;
    tmp64 %= curprime;
    if (tmp64 < 0) {
        tmp64 += curprime;
    }
    poly->sol2->data[i] = (int32_t) tmp64;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
ecode_t
full_poly_init(siqs_poly_t* const poly)
{
    poly->nfullpolyinit++;
    //
    // Determines the first polynomial having 'a' as its leading coefficient
    // and performs all needed precomputations. Notations are mainly that of
    // S.P. Contini's thesis "Factoring integers with the Self-initializing
    // Quadratic Sieve".
    //
    const uint32_array_t* const sqrtm_pi = poly->sqrtm_pi;
    const uint32_t*       const idx_of_a = poly->idx_of_a;
    const uint32_t*       const fb       = poly->factor_base->data;

    const uint32_t nprimes_in_a = poly->nprimes_in_a;
    const uint32_t fb_length = poly->factor_base->length;
    const uint32_t imin = poly->approximer->imin;
    const uint32_t imax = poly->approximer->imax;
    
    const mpz_ptr a = poly->a;
    mpz_ptr b       = poly->b;
    mpz_ptr c       = poly->c;

    mpz_array_t*   const Bl   = poly->Bl;
    int32_array_t* const sol1 = poly->sol1;
    int32_array_t* const sol2 = poly->sol2;

    ecode_t ecode = SUCCESS;

    mpz_t inv;
    mpz_init(inv);
    mpz_t mpzp;
    mpz_init(mpzp);

    uint64_t gamma    = 0;
    uint32_t curprime = 0;
    
    mpz_set_ui(b, 0);
    //
    // Computes the {Bl}, 1 <= l <= s such that:
    //     Bl^2 = N mod ql
    //     Bl = 0 mod qj forall j != l
    // where:
    //     qi are the primes such that a = q1 x q2 x ... x qs
    //
    for (uint32_t l = 0; l < nprimes_in_a; l++) {
  
        curprime = fb[idx_of_a[l]];
        mpz_set_ui(mpzp, curprime);
        mpz_divexact_ui(Bl->data[l], a, curprime);

        if (0 == mpz_invert(inv, Bl->data[l], mpzp)) {
            PRINTF_STDERR("\n");
            PRINT_ERROR("%Zd has no inverse mod %Zd!\n", Bl->data[l], mpzp);
            PRINT_ERROR("this should not happen - factorization aborted!\n");
            fflush(stdout);
            
            ecode = FATAL_INTERNAL_ERROR;
            goto clean_and_return;
        };
        gamma  = mpz_get_ui(inv) * sqrtm_pi->data[idx_of_a[l]];
        //
        // _NOTE_: One has to be careful when dealing with a mix of several
        //         integer types and pay attention to the possible automatic
        //         type coercion. For example a int64_t % a uint32_t will give
        //         the expected result, but not a int32_t % a uint32_t...
        //
        gamma %= curprime;
        if (gamma > (curprime >> 1)) {
            gamma = curprime - gamma;
        }
        mpz_mul_ui(Bl->data[l], Bl->data[l], (unsigned long int)gamma);
        mpz_add(b, b, Bl->data[l]);        
    }
    Bl->length = Bl->alloced;
    
    mpz_mul(c, b, b);
    mpz_sub(c, c, poly->n);
    mpz_divexact(c, c, a);

    poly->loga = mpz_sizeinbase(a, 2);
    poly->logb = mpz_sizeinbase(b, 2);
    poly->logc = mpz_sizeinbase(c, 2);
    //
    // Computes the two solutions to g_{a,b}(x) = 0 mod p for each prime
    // and stores them in the arrays sol1 and sol2.
    //
    // Also, for each prime p that does not divide the polynomial coefficient
    // 'a', computes and store the values 2 x Bj x inv(a).mod(p) mod p foreach
    // prime pj dividing 'a'. These values will be kept in the array
    // Bainv2 indexed first by the primes p, then by the indexes j.
    //
    uint32_t* avoid = (uint32_t*)&(idx_of_a[0]);

    uint32_t rimin = MAX(imin, avoid[0]);
    uint32_t rimax = MIN(imax, avoid[nprimes_in_a - 1]);
    
    avoid++;
    
    poly_u64_t w;
    w.fit = (sizeof(unsigned long int) >= 8)
            && (nprimes_in_a <= MAX_NPRIMES_IN_A_U64)
            && (mpz_sizeinbase(a, 2) <= 62) && (mpz_sizeinbase(b, 2) <= 64);
    if (w.fit) {
        w.a = mpz_get_ui(a);
        w.b = mpz_get_ui(b);
        for (uint32_t j = 0; j < nprimes_in_a; j++) {
            w.Bl[j] = mpz_get_ui(Bl->data[j]);
        }
    }
    const uint32_t* fbptr = &(fb[0]);
    
    for (uint32_t i = 0; i < rimin; i++) {
        ecode = init_prime(poly, i, *fbptr, &w, inv, mpzp);
        if (ecode != SUCCESS) {
            goto clean_and_return;
        }
        fbptr++;
    }
    fbptr++;
    
    for (uint32_t i = rimin + 1; i < rimax; i++) {
        if (i == *avoid) {
           avoid++;
           fbptr++;
           continue;
        }
        ecode = init_prime(poly, i, *fbptr, &w, inv, mpzp);
        if (ecode != SUCCESS) {
            goto clean_and_return;
        }
        fbptr++;
    }
    fbptr++;
    
    for (uint32_t i = rimax + 1; i < fb_length; i++) {
        ecode = init_prime(poly, i, *fbptr, &w, inv, mpzp);
        if (ecode != SUCCESS) {
            goto clean_and_return;
        }
        fbptr++;
    }
    sol1->length = sol1->alloced;
    sol2->length = sol2->alloced;
    
    poly->polyno = 1;
    
  clean_and_return:

    mpz_clear(inv);
    mpz_clear(mpzp);

    return ecode;
}
//------------------------------------------------------------------------------
ecode_t
fast_poly_init(siqs_poly_t* const poly)
{
    //
    // Determines the next polynomial having 'a' as its leading coefficient and
    // updates the needed data. Notations are mainly that of S.P. Contini's
    // thesis "Factoring integers with the Self-initializing Quadratic Sieve".
    //
    // The current polynomial is indexed by ipol so we are now initializing
    // the polynomial indexed by (ipol + 1)
    //
    const uint32_t ipol = poly->polyno;
    const uint32_t imin = poly->approximer->imin;
    const uint32_t imax = poly->approximer->imax;
    mpz_ptr b           = poly->b;
    mpz_ptr c           = poly->c;

    const uint32_array_t* const base = poly->factor_base;
    const mpz_array_t*    const Bl   = poly->Bl;

    const uint32_t fb_length = base->length;

    uint32_t**     const Bainv2 = poly->Bainv2;
    int32_array_t* const sol1   = poly->sol1;
    int32_array_t* const sol2   = poly->sol2;
    
    int32_t* sol1ptr = sol1->data;
    int32_t* sol2ptr = sol2->data;
    
    //
    // Compute the 'b' coefficient of the next polynomial to use
    //
    uint32_t nu = 0;
    while ( ((2 * ipol) & (1<<nu)) == 0) {
        nu++;
    }
    uint32_t io2nu = (ipol / (1<<nu)) + 1;

    if (IS_EVEN(io2nu)) {
        mpz_add(b, b, Bl->data[nu-1]);
        mpz_add(b, b, Bl->data[nu-1]);
    } else {
        mpz_sub(b, b, Bl->data[nu-1]);
        mpz_sub(b, b, Bl->data[nu-1]);
    }
    mpz_mul(c, b, b);
    mpz_sub(c, c, poly->n);
    mpz_divexact(c, c, poly->a);

    poly->logb = mpz_sizeinbase(b, 2);
    poly->logc = mpz_sizeinbase(c, 2);

    uint32_t curprime = 0;

    uint32_t* avoid_primes = poly->idx_of_a;

    uint32_t rimin = MAX(imin, avoid_primes[0]);
    uint32_t rimax = MIN(imax, avoid_primes[poly->nprimes_in_a - 1]);
    
    avoid_primes++;
    
    //
    // Compute the new sol1 and sol2 arrrays for sieving with the new polynomial
    //
    if (IS_ODD(io2nu)) {

        for (uint32_t i = 0; i < rimin; i++) {
            curprime = base->data[i];

            *sol1ptr = add_mod(*sol1ptr, Bainv2[i][nu-1], curprime);
            *sol2ptr = add_mod(*sol2ptr, Bainv2[i][nu-1], curprime);
            
            sol1ptr++;
            sol2ptr++;
        }
        sol1ptr++;
        sol2ptr++;
        
        for (uint32_t i = rimin + 1; i < rimax; i++) {
            if (i == avoid_primes[0]) {
                avoid_primes++;
                sol1ptr++;
                sol2ptr++;
                continue;
            }
            curprime = base->data[i];

            *sol1ptr = add_mod(*sol1ptr, Bainv2[i][nu-1], curprime);
            *sol2ptr = add_mod(*sol2ptr, Bainv2[i][nu-1], curprime);
            
            sol1ptr++;
            sol2ptr++;       
        }
        sol1ptr++;
        sol2ptr++;
        
        for (uint32_t i = rimax + 1; i < fb_length; i++) {
            curprime = base->data[i];

            *sol1ptr = add_mod(*sol1ptr, Bainv2[i][nu-1], curprime);
            *sol2ptr = add_mod(*sol2ptr, Bainv2[i][nu-1], curprime);
            
            sol1ptr++;
            sol2ptr++;        
        }
        
    } else {
        for (uint32_t i = 0; i < rimin; i++) {
            curprime = base->data[i];
            
            *sol1ptr = sub_mod(*sol1ptr, Bainv2[i][nu-1], curprime);
            *sol2ptr = sub_mod(*sol2ptr, Bainv2[i][nu-1], curprime);
            
            sol1ptr++;
            sol2ptr++;
        }
        sol1ptr++;
        sol2ptr++;
        
        for (uint32_t i = rimin + 1; i < rimax; i++) {
            if (i == avoid_primes[0]) {
                avoid_primes++;
                sol1ptr++;
                sol2ptr++;
                continue;                
            }
            curprime = base->data[i];
            
            *sol1ptr = sub_mod(*sol1ptr, Bainv2[i][nu-1], curprime);
            *sol2ptr = sub_mod(*sol2ptr, Bainv2[i][nu-1], curprime);
            
            sol1ptr++;
            sol2ptr++;
        }
        sol1ptr++;
        sol2ptr++;
        
        for (uint32_t i = rimax + 1; i < fb_length; i++) {
            curprime = base->data[i];
            
            *sol1ptr = sub_mod(*sol1ptr, Bainv2[i][nu-1], curprime);
            *sol2ptr = sub_mod(*sol2ptr, Bainv2[i][nu-1], curprime);
            
            sol1ptr++;
            sol2ptr++;
        }
    }
    poly->polyno++;
    
    return SUCCESS;
}
//-----------------------------------------------------------------------------
int na_used(siqs_poly_t* const poly) {
    return poly->nfullpolyinit;
}
//-----------------------------------------------------------------------------
