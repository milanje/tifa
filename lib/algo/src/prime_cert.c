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
 * \file    prime_cert.c
 *
 * \brief Primality certificates (see prime_cert.h).
 */

#include <stdlib.h>
#include <stdbool.h>

#include "tifa_config.h"
#include "prime_cert.h"
#include "funcs.h"
#include "tifa_factor.h"
#include "tifa_alloc.h"

#define NSPRP_BASES 12

static const uint32_t sprp_bases[NSPRP_BASES] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
};

//
// For a prime p, a^((p - 1) / q) = 1 for a fraction 1/q of the a: beyond
// this, p is almost surely composite, but nothing is proven.
//
#define MAX_WITNESS 1000000

//------------------------------------------------------------------------------
//                          Primes below 2^64
//------------------------------------------------------------------------------

static bool is_small_prime(const mpz_t n) {
    if ((mpz_sgn(n) <= 0) || (mpz_sizeinbase(n, 2) > 64)) {
        return false;
    }
    for (uint32_t i = 0; i < NSPRP_BASES; i++) {
        if (mpz_cmp_ui(n, sprp_bases[i]) == 0) {
            return true;
        }
        if (mpz_divisible_ui_p(n, sprp_bases[i])) {
            return false;
        }
    }
    if (mpz_cmp_ui(n, 1) == 0) {
        return false;
    }
    mpz_t nm1, d, x;
    mpz_init(nm1);
    mpz_init(d);
    mpz_init(x);
    mpz_sub_ui(nm1, n, 1);
    const mp_bitcnt_t s = mpz_scan1(nm1, 0);
    mpz_tdiv_q_2exp(d, nm1, s);

    bool prime = true;
    for (uint32_t i = 0; prime && (i < NSPRP_BASES); i++) {
        mpz_set_ui(x, sprp_bases[i]);
        mpz_powm(x, x, d, n);
        bool pass = (mpz_cmp_ui(x, 1) == 0) || (mpz_cmp(x, nm1) == 0);
        for (mp_bitcnt_t r = 1; !pass && (r < s); r++) {
            mpz_mul(x, x, x);
            mpz_mod(x, x, n);
            pass = (mpz_cmp(x, nm1) == 0);
        }
        prime = pass;
    }
    mpz_clear(nm1);
    mpz_clear(d);
    mpz_clear(x);
    return prime;
}

//------------------------------------------------------------------------------
//                              Certificates
//------------------------------------------------------------------------------

prime_cert_t* alloc_prime_cert(void) {
    prime_cert_t* const cert = tifa_malloc(sizeof(prime_cert_t));
    cert->nsteps  = 0;
    cert->alloced = 0;
    cert->steps   = NULL;
    return cert;
}
//-----------------------------------------------------------------------------
static void clear_prime_cert(prime_cert_t* const cert) {
    for (uint32_t i = 0; i < cert->nsteps; i++) {
        mpz_clear(cert->steps[i].p);
        free_mpz_array(cert->steps[i].primes);
        free_uint32_array(cert->steps[i].witnesses);
    }
    cert->nsteps = 0;
}
//-----------------------------------------------------------------------------
void free_prime_cert(prime_cert_t* cert) {
    if (cert != NULL) {
        clear_prime_cert(cert);
        free(cert->steps);
        free(cert);
    }
}
//-----------------------------------------------------------------------------
static void append_step(prime_cert_t* const cert, const mpz_t p,
                        const mpz_array_t* const primes,
                        const uint32_array_t* const witnesses) {
    if (cert->nsteps == cert->alloced) {
        cert->alloced = (cert->alloced == 0) ? 4 : 2 * cert->alloced;
        cert->steps   = tifa_realloc(cert->steps,
                                     cert->alloced * sizeof(prime_cert_step_t));
    }
    prime_cert_step_t* const step = &cert->steps[cert->nsteps++];
    mpz_init_set(step->p, p);
    step->primes    = alloc_mpz_array(primes->length);
    step->witnesses = alloc_uint32_array(witnesses->length);
    append_mpz_array(step->primes, primes);
    append_uint32_array(step->witnesses, witnesses);
}
//-----------------------------------------------------------------------------
static bool is_proven(const prime_cert_t* const cert, const uint32_t before,
                      const mpz_t q) {
    if (mpz_sizeinbase(q, 2) <= 64) {
        return is_small_prime(q);
    }
    for (uint32_t j = 0; j < before; j++) {
        if (mpz_cmp(cert->steps[j].p, q) == 0) {
            return true;
        }
    }
    return false;
}
//-----------------------------------------------------------------------------
//
// Whether c1^2 - 4 c2 is a square, with p - 1 = (c2 F + c1) F, 0 <= c1 < F.
//
static bool bls_square(const mpz_t pm1, const mpz_t f) {
    mpz_t r, c1, c2;
    mpz_init(r);
    mpz_init(c1);
    mpz_init(c2);
    mpz_divexact(r, pm1, f);
    mpz_fdiv_qr(c2, c1, r, f);
    mpz_mul(c1, c1, c1);
    mpz_submul_ui(c1, c2, 4);
    const bool square = (mpz_sgn(c1) >= 0) && mpz_perfect_square_p(c1);
    mpz_clear(r);
    mpz_clear(c1);
    mpz_clear(c2);
    return square;
}
//-----------------------------------------------------------------------------
//
// Whether (F + 1)^2 > p, or F^3 > p and the test of Brillhart, Lehmer and
// Selfridge passes.
//
static bool large_enough(const mpz_t p, const mpz_t pm1, const mpz_t f) {
    mpz_t t;
    mpz_init(t);
    mpz_add_ui(t, f, 1);
    mpz_mul(t, t, t);
    bool proven = (mpz_cmp(t, p) > 0);
    if (!proven) {
        mpz_pow_ui(t, f, 3);
        proven = (mpz_cmp(t, p) > 0) && !bls_square(pm1, f);
    }
    mpz_clear(t);
    return proven;
}
//-----------------------------------------------------------------------------
//
// 1 if a proves q for p, 0 if not, -1 if a shows p composite. x is
// a^((p - 1) / q).
//
static int try_witness(mpz_t x, const mpz_t p, const mpz_t pm1, const mpz_t q,
                       const uint32_t a) {
    mpz_t t;
    mpz_init(t);
    mpz_divexact(t, pm1, q);
    mpz_set_ui(x, a);
    mpz_powm(x, x, t, p);
    mpz_powm(t, x, q, p);

    int outcome = 1;
    if (mpz_cmp_ui(t, 1) != 0) {
        outcome = -1;
    } else if (mpz_cmp_ui(x, 1) == 0) {
        outcome = 0;
    } else {
        mpz_sub_ui(t, x, 1);
        mpz_gcd(t, t, p);
        outcome = (mpz_cmp_ui(t, 1) == 0) ? 1 : -1;
    }
    mpz_clear(t);
    return outcome;
}
//-----------------------------------------------------------------------------
static bool check_step(const prime_cert_t* const cert, const uint32_t i) {
    const prime_cert_step_t* const step = &cert->steps[i];
    const uint32_t k = step->primes->length;

    if ((mpz_cmp_ui(step->p, 3) < 0) || mpz_even_p(step->p)
            || (k == 0) || (step->witnesses->length != k)) {
        return false;
    }
    mpz_t pm1, rest, f, x;
    mpz_init(pm1);
    mpz_init(rest);
    mpz_init(f);
    mpz_init(x);
    mpz_sub_ui(pm1, step->p, 1);
    mpz_set(rest, pm1);

    bool valid = true;
    for (uint32_t j = 0; valid && (j < k); j++) {
        const mpz_ptr q = step->primes->data[j];
        valid =    is_proven(cert, i, q)
                && mpz_divisible_p(rest, q)
                && (try_witness(x, step->p, pm1, q,
                                step->witnesses->data[j]) == 1);
        if (valid) {
            mpz_remove(rest, rest, q);
        }
    }
    if (valid) {
        mpz_divexact(f, pm1, rest);
        valid = large_enough(step->p, pm1, f);
    }
    mpz_clear(pm1);
    mpz_clear(rest);
    mpz_clear(f);
    mpz_clear(x);
    return valid;
}
//-----------------------------------------------------------------------------
bool tifa_check_prime_cert(const prime_cert_t* cert, const mpz_t n) {
    if (cert->nsteps == 0) {
        return is_small_prime(n);
    }
    if (mpz_cmp(cert->steps[cert->nsteps - 1].p, n) != 0) {
        return false;
    }
    for (uint32_t i = 0; i < cert->nsteps; i++) {
        if (!check_step(cert, i)) {
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------
//                                 Proofs
//------------------------------------------------------------------------------

//
// Appends to cert the steps proving n, after those of the primes they use.
//
static prime_cert_status_t prove(prime_cert_t* const cert, const mpz_t n) {
    if (mpz_sizeinbase(n, 2) <= 64) {
        return is_small_prime(n) ? PRIME_PROVEN : NOT_PRIME;
    }
    if (is_proven(cert, cert->nsteps, n)) {
        return PRIME_PROVEN;
    }
    if (!is_prime_mpz(n)) {
        return NOT_PRIME;
    }
    mpz_t pm1, rest, f, x;
    mpz_init(pm1);
    mpz_init(rest);
    mpz_init(f);
    mpz_init(x);
    mpz_sub_ui(pm1, n, 1);

    mpz_array_t*    factors   = alloc_mpz_array(ELONGATION);
    uint32_array_t* multis    = alloc_uint32_array(ELONGATION);
    mpz_array_t*    primes    = alloc_mpz_array(ELONGATION);
    uint32_array_t* witnesses = alloc_uint32_array(ELONGATION);

    prime_cert_status_t status = PRIME_NOT_PROVEN;
    const ecode_t ecode = tifa_factor(factors, multis, pm1,
                                      FIND_COMPLETE_FACTORIZATION);
    if (   (ecode != COMPLETE_FACTORIZATION_FOUND)
        && (ecode != PARTIAL_FACTORIZATION_FOUND)) {
        goto clean_up;
    }
    const bool all_prime = (ecode == COMPLETE_FACTORIZATION_FOUND);
    //
    // The smallest primes of p - 1 first: those above 2^64 need proofs.
    //
    mpz_set(rest, pm1);
    mpz_set_ui(f, 1);
    bool enough = false;
    while (!enough) {
        mpz_ptr q = NULL;
        for (uint32_t i = 0; i < factors->length; i++) {
            const mpz_ptr c = factors->data[i];
            if (   mpz_divisible_p(rest, c)
                && ((q == NULL) || (mpz_cmp(c, q) < 0))
                && (all_prime || is_prime_mpz(c))) {
                q = c;
            }
        }
        if (q == NULL) {
            break;
        }
        mpz_remove(rest, rest, q);
        mpz_divexact(f, pm1, rest);
        append_mpz_to_array(primes, q);

        mpz_pow_ui(x, f, 3);
        enough = (mpz_cmp(x, n) > 0);
    }
    if (!enough) {
        goto clean_up;
    }
    for (uint32_t i = 0; i < primes->length; i++) {
        int outcome = 0;
        uint32_t a  = 1;
        while ((outcome == 0) && (a < MAX_WITNESS)) {
            a++;
            outcome = try_witness(x, n, pm1, primes->data[i], a);
        }
        if (outcome != 1) {
            status = (outcome == -1) ? NOT_PRIME : PRIME_NOT_PROVEN;
            goto clean_up;
        }
        append_uint32_to_array(witnesses, a);
    }
    if (!large_enough(n, pm1, f)) {
        status = NOT_PRIME;
        goto clean_up;
    }
    for (uint32_t i = 0; i < primes->length; i++) {
        if (prove(cert, primes->data[i]) != PRIME_PROVEN) {
            goto clean_up;
        }
    }
    append_step(cert, n, primes, witnesses);
    status = PRIME_PROVEN;

  clean_up:
    free_mpz_array(factors);
    free_uint32_array(multis);
    free_mpz_array(primes);
    free_uint32_array(witnesses);
    mpz_clear(pm1);
    mpz_clear(rest);
    mpz_clear(f);
    mpz_clear(x);
    return status;
}
//-----------------------------------------------------------------------------
prime_cert_status_t tifa_prove_prime(prime_cert_t* cert, const mpz_t n) {
    clear_prime_cert(cert);
    const prime_cert_status_t status = prove(cert, n);
    if (status != PRIME_PROVEN) {
        clear_prime_cert(cert);
    }
    return status;
}
