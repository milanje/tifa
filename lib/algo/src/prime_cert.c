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

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <pthread.h>

#include "tifa_config.h"
#include "prime_cert.h"
#include "ecm_schedule.h"
#include "funcs.h"
#include "mont_ecm.h"
#include "tdiv.h"
#include "tifa_factor.h"
#include "tifa_alloc.h"
#include "impl.h"

#define NSPRP_BASES 12

static const uint32_t sprp_bases[NSPRP_BASES] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
};

//
// For a prime p, a^((p - 1) / q) = 1 for a fraction 1/q of the a, and
// U_((p + 1) / r) = 0 for about 1/r of the P: beyond these, p is almost
// surely composite, but nothing is proven.
//
#define MAX_WITNESS 1000000
#define MAX_D       1000000
//
// Trial division of p - 1 and p + 1 by the primes up to 7919, then ECM on
// what is left above 100 bits (below, tifa_factor is as fast), for a
// twentieth of the time of SIQS: what it misses, the other side or SIQS on
// the smallest composite left often get for less (0.05 took 0.7 to 0.85 of
// the time of 0.2 from 192 to 256 bits).
//
#define CERT_NPRIMES_TDIV 1000
#define CERT_EASY_BITS    100
#define CERT_ECM_BUDGET   0.05

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
        prime_cert_step_t* const step = &cert->steps[i];
        mpz_clear(step->p);
        free_mpz_array(step->pm1_primes);
        free_uint32_array(step->pm1_witnesses);
        free_mpz_array(step->pp1_primes);
        free_uint32_array(step->pp1_witnesses);
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
static mpz_array_t* copy_mpz_array(const mpz_array_t* const a) {
    mpz_array_t* const copy = alloc_mpz_array(a->length);
    append_mpz_array(copy, a);
    return copy;
}
//-----------------------------------------------------------------------------
static uint32_array_t* copy_uint32_array(const uint32_array_t* const a) {
    uint32_array_t* const copy = alloc_uint32_array(a->length);
    append_uint32_array(copy, a);
    return copy;
}
//-----------------------------------------------------------------------------
static void append_step(prime_cert_t* const cert, const mpz_t p,
                        const mpz_array_t* const pm1_primes,
                        const uint32_array_t* const pm1_witnesses,
                        const int32_t d,
                        const mpz_array_t* const pp1_primes,
                        const uint32_array_t* const pp1_witnesses) {
    if (cert->nsteps == cert->alloced) {
        cert->alloced = (cert->alloced == 0) ? 4 : 2 * cert->alloced;
        cert->steps   = tifa_realloc(cert->steps,
                                     cert->alloced * sizeof(prime_cert_step_t));
    }
    prime_cert_step_t* const step = &cert->steps[cert->nsteps++];
    mpz_init_set(step->p, p);
    step->pm1_primes    = copy_mpz_array(pm1_primes);
    step->pm1_witnesses = copy_uint32_array(pm1_witnesses);
    step->pp1_d         = d;
    step->pp1_primes    = copy_mpz_array(pp1_primes);
    step->pp1_witnesses = copy_uint32_array(pp1_witnesses);
}
//-----------------------------------------------------------------------------
static bool has_step(const prime_cert_t* const cert, const mpz_t q) {
    for (uint32_t j = 0; j < cert->nsteps; j++) {
        if (mpz_cmp(cert->steps[j].p, q) == 0) {
            return true;
        }
    }
    return false;
}

//------------------------------------------------------------------------------
//                       The conditions of a step
//------------------------------------------------------------------------------

//
// 1 if a proves q for p - 1, 0 if not, -1 if a shows p composite.
//
static int try_witness(const mpz_t p, const mpz_t pm1, const mpz_t q,
                       const uint32_t a) {
    mpz_t x, t;
    mpz_init(x);
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
    mpz_clear(x);
    mpz_clear(t);
    return outcome;
}
//-----------------------------------------------------------------------------
static void half_mod(mpz_t x, const mpz_t n) {
    if (mpz_odd_p(x)) {
        mpz_add(x, x, n);
    }
    mpz_tdiv_q_2exp(x, x, 1);
}
//-----------------------------------------------------------------------------
//
// U_k mod n for the Lucas sequence of parameters P and Q = (P^2 - D) / 4,
// n odd, k >= 1.
//
static void lucas_u(mpz_t u, const mpz_t k, const uint32_t p, const mpz_t q,
                    const int32_t d, const mpz_t n) {
    mpz_t v, qk, t;
    mpz_init_set_ui(v, p);
    mpz_init(qk);
    mpz_init(t);
    mpz_mod(qk, q, n);
    mpz_set_ui(u, 1);

    for (long b = (long) mpz_sizeinbase(k, 2) - 2; b >= 0; b--) {
        mpz_mul(u, u, v);
        mpz_mod(u, u, n);
        mpz_mul(v, v, v);
        mpz_submul_ui(v, qk, 2);
        mpz_mod(v, v, n);
        mpz_mul(qk, qk, qk);
        mpz_mod(qk, qk, n);
        if (mpz_tstbit(k, (mp_bitcnt_t) b)) {
            mpz_mul_ui(t, u, p);
            mpz_add(t, t, v);
            mpz_mul_si(u, u, d);
            mpz_addmul_ui(u, v, p);
            mpz_mod(v, u, n);
            mpz_mod(u, t, n);
            half_mod(u, n);
            half_mod(v, n);
            mpz_mul(qk, qk, q);
            mpz_mod(qk, qk, n);
        }
    }
    mpz_clear(v);
    mpz_clear(qk);
    mpz_clear(t);
}
//-----------------------------------------------------------------------------
//
// 1 if P proves r for p + 1 with D, 0 if not, -1 if P shows p composite
// (when (D / p) = -1).
//
static int try_lucas(const mpz_t p, const mpz_t pp1, const mpz_t r,
                     const int32_t d, const uint32_t pw) {
    mpz_t q, u, t;
    mpz_init(q);
    mpz_init(u);
    mpz_init(t);
    mpz_set_ui(q, pw);
    mpz_mul(q, q, q);
    mpz_set_si(t, d);
    mpz_sub(q, q, t);

    int outcome = 0;
    if (mpz_divisible_2exp_p(q, 2)) {
        mpz_divexact_ui(q, q, 4);
        mpz_gcd(t, q, p);
        if (mpz_sgn(q) == 0) {
            outcome = 0;
        } else if (mpz_cmp_ui(t, 1) != 0) {
            outcome = -1;
        } else {
            lucas_u(u, pp1, pw, q, d, p);
            if (mpz_sgn(u) != 0) {
                outcome = -1;
            } else {
                mpz_divexact(t, pp1, r);
                lucas_u(u, t, pw, q, d, p);
                mpz_gcd(t, u, p);
                outcome = (mpz_sgn(u) == 0)         ? 0
                        : (mpz_cmp_ui(t, 1) == 0)   ? 1
                        :                             -1;
            }
        }
    }
    mpz_clear(q);
    mpz_clear(u);
    mpz_clear(t);
    return outcome;
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
// The three conditions of prime_cert.h, given that every prime factor of p is
// 1 mod f1 and +-1 mod f2.
//
static bool large_enough(const mpz_t p, const mpz_t pm1, const mpz_t f1,
                         const mpz_t f2) {
    mpz_t t, g;
    mpz_init(t);
    mpz_init(g);
    mpz_add_ui(t, f1, 1);
    mpz_mul(t, t, t);
    bool proven = (mpz_cmp(t, p) > 0);
    if (!proven) {
        mpz_pow_ui(t, f1, 3);
        proven = (mpz_cmp(t, p) > 0) && !bls_square(pm1, f1);
    }
    if (!proven) {
        mpz_lcm(g, f1, f2);
        mpz_mul(t, g, g);
        if (mpz_cmp(t, p) > 0) {
            mpz_mod(t, p, g);
            proven =    (mpz_cmp_ui(t, 1) == 0) || (mpz_cmp(t, p) >= 0)
                     || !mpz_divisible_p(p, t);
        }
    }
    mpz_clear(t);
    mpz_clear(g);
    return proven;
}
//-----------------------------------------------------------------------------
//
// The part f of x made of the primes, which must divide what is left of x
// once those before them are taken out.
//
static bool part_of(mpz_t f, const mpz_t x, const mpz_array_t* const primes) {
    mpz_t rest;
    mpz_init_set(rest, x);
    bool valid = true;
    for (uint32_t j = 0; valid && (j < primes->length); j++) {
        valid =    (mpz_cmp_ui(primes->data[j], 1) > 0)
                && mpz_divisible_p(rest, primes->data[j]);
        if (valid) {
            mpz_remove(rest, rest, primes->data[j]);
        }
    }
    mpz_divexact(f, x, rest);
    mpz_clear(rest);
    return valid;
}
//-----------------------------------------------------------------------------
//
// The steps' primes, sorted: a search for each q of a certificate of k steps
// would otherwise take time in k^2.
//
typedef struct {
    mpz_srcptr p;
    uint32_t   step;
} index_entry_t;

static int compare_entries(const void* const a, const void* const b) {
    return mpz_cmp(((const index_entry_t*) a)->p, ((const index_entry_t*) b)->p);
}
//-----------------------------------------------------------------------------
static bool is_proven(const index_entry_t* const index, const uint32_t nsteps,
                      const uint32_t before, const mpz_t q) {
    if (mpz_sizeinbase(q, 2) <= 64) {
        return is_small_prime(q);
    }
    const index_entry_t key = { q, 0 };
    const index_entry_t* const found =
        bsearch(&key, index, nsteps, sizeof(index_entry_t), compare_entries);
    return (found != NULL) && (found->step < before);
}
//-----------------------------------------------------------------------------
static unsigned long abs_d(const int32_t d) {
    return (d < 0) ? 0UL - (unsigned long) d : (unsigned long) d;
}
//-----------------------------------------------------------------------------
static bool check_step(const prime_cert_t* const cert,
                       const index_entry_t* const index, const uint32_t i) {
    const prime_cert_step_t* const step = &cert->steps[i];
    const mpz_array_t* const qs = step->pm1_primes;
    const mpz_array_t* const rs = step->pp1_primes;

    if (   (mpz_cmp_ui(step->p, 3) < 0) || mpz_even_p(step->p)
        || (qs->length + rs->length == 0)
        || (step->pm1_witnesses->length != qs->length)
        || (step->pp1_witnesses->length != rs->length)
        || ((rs->length != 0)
            && ((step->pp1_d == 0)
                || (mpz_gcd_ui(NULL, step->p, abs_d(step->pp1_d)) != 1)))) {
        return false;
    }
    mpz_t pm1, pp1, f1, f2;
    mpz_init(pm1);
    mpz_init(pp1);
    mpz_init(f1);
    mpz_init(f2);
    mpz_sub_ui(pm1, step->p, 1);
    mpz_add_ui(pp1, step->p, 1);

    bool valid = part_of(f1, pm1, qs) && part_of(f2, pp1, rs);
    for (uint32_t j = 0; valid && (j < qs->length); j++) {
        valid =    is_proven(index, cert->nsteps, i, qs->data[j])
                && (try_witness(step->p, pm1, qs->data[j],
                                step->pm1_witnesses->data[j]) == 1);
    }
    for (uint32_t j = 0; valid && (j < rs->length); j++) {
        valid =    is_proven(index, cert->nsteps, i, rs->data[j])
                && (try_lucas(step->p, pp1, rs->data[j], step->pp1_d,
                              step->pp1_witnesses->data[j]) == 1);
    }
    valid = valid && large_enough(step->p, pm1, f1, f2);

    mpz_clear(pm1);
    mpz_clear(pp1);
    mpz_clear(f1);
    mpz_clear(f2);
    return valid;
}
//-----------------------------------------------------------------------------
bool tifa_check_prime_cert(const prime_cert_t* cert, const mpz_t n) {
    if (cert->nsteps == 0) {
        return is_small_prime(n);
    }
    const uint32_t last = cert->nsteps - 1;
    if (mpz_cmp(cert->steps[last].p, n) != 0) {
        return false;
    }
    index_entry_t* const index =
        tifa_malloc(cert->nsteps * sizeof(index_entry_t));
    for (uint32_t i = 0; i < cert->nsteps; i++) {
        index[i].p    = cert->steps[i].p;
        index[i].step = i;
    }
    qsort(index, cert->nsteps, sizeof(index_entry_t), compare_entries);
    //
    // One step per prime, n's the largest: no step costs more than n's, so
    // the time is bounded by the length of the certificate.
    //
    bool valid = (index[last].step == last);
    for (uint32_t i = 1; valid && (i < cert->nsteps); i++) {
        valid = (mpz_cmp(index[i - 1].p, index[i].p) != 0);
    }
    for (uint32_t i = 0; valid && (i < cert->nsteps); i++) {
        valid = check_step(cert, index, i);
    }
    free(index);
    return valid;
}

//------------------------------------------------------------------------------
//                     The factors of p - 1 and p + 1
//------------------------------------------------------------------------------

#if TIFA_HAVE_MONT_ECM
static ecm_schedule_t* ecm_schedule      = NULL;
static pthread_once_t  ecm_schedule_once = PTHREAD_ONCE_INIT;

static void build_ecm_schedule(void) {
    ecm_schedule = alloc_ecm_schedule(60);
}
#endif

//
// What is known of p - 1 or p + 1: its (probable) primes found, and the
// composites not split yet.
//
struct side_struct {
    mpz_t x;
    mpz_array_t* primes;
    mpz_array_t* hard;
};
typedef struct side_struct side_t;

static void init_side(side_t* const s, const mpz_t x) {
    mpz_init_set(s->x, x);
    s->primes = alloc_mpz_array(ELONGATION);
    s->hard   = alloc_mpz_array(2);
}
//-----------------------------------------------------------------------------
static void clear_side(side_t* const s) {
    mpz_clear(s->x);
    free_mpz_array(s->primes);
    free_mpz_array(s->hard);
}
//-----------------------------------------------------------------------------
static void add_prime(side_t* const s, const mpz_t q) {
    for (uint32_t i = 0; i < s->primes->length; i++) {
        if (mpz_cmp(s->primes->data[i], q) == 0) {
            return;
        }
    }
    append_mpz_to_array(s->primes, q);
}
//-----------------------------------------------------------------------------
//
// Adds the factors of c that tifa_factor finds (all of them, unless it
// fails); the composites left are dropped.
//
static void factor_fully(side_t* const s, const mpz_t c,
                         const tifa_factor_hint_t hint) {
    mpz_array_t*    f = alloc_mpz_array(ELONGATION);
    uint32_array_t* m = alloc_uint32_array(ELONGATION);
    const ecode_t ecode = tifa_factor_impl(f, m, c,
                                           FIND_COMPLETE_FACTORIZATION, hint);
    if (ecode == NO_FACTOR_FOUND) {
        add_prime(s, c);
    }
    for (uint32_t i = 0; i < f->length; i++) {
        if ((ecode == COMPLETE_FACTORIZATION_FOUND) || is_prime_mpz(f->data[i])) {
            add_prime(s, f->data[i]);
        }
    }
    free_mpz_array(f);
    free_uint32_array(m);
}
//-----------------------------------------------------------------------------
//
// Trial division and ECM, but no SIQS on composites above CERT_EASY_BITS.
//
static void search_side(side_t* const s) {
    mpz_array_t*    f    = alloc_mpz_array(ELONGATION);
    uint32_array_t* m    = alloc_uint32_array(ELONGATION);
    mpz_array_t*    todo = alloc_mpz_array(4);

    const ecode_t ecode = tdiv(f, m, s->x, CERT_NPRIMES_TDIV);
    if (ecode == NO_FACTOR_FOUND) {
        append_mpz_to_array(todo, s->x);
    }
    for (uint32_t i = 0; i < f->length; i++) {
        if ((ecode == SOME_FACTORS_FOUND) && (i == f->length - 1)) {
            append_mpz_to_array(todo, f->data[i]);
        } else {
            add_prime(s, f->data[i]);
        }
    }
    mpz_t c, g;
    mpz_init(c);
    mpz_init(g);
    while (todo->length != 0) {
        todo->length--;
        mpz_set(c, todo->data[todo->length]);

        bool split = false;
        if (is_prime_mpz(c)) {
            add_prime(s, c);
            continue;
        }
        if (   (mpz_sizeinbase(c, 2) <= CERT_EASY_BITS)
            || mpz_perfect_power_p(c)) {
            factor_fully(s, c, TIFA_FACTOR_RANDOM);
            continue;
        }
#if TIFA_HAVE_MONT_ECM
        pthread_once(&ecm_schedule_once, build_ecm_schedule);
        split =    mont_ecm_handles(c)
                && ecm_schedule_prepass(g, c, ecm_schedule, CERT_ECM_BUDGET,
                                        INFINITY);
#endif
        if (split) {
            append_mpz_to_array(todo, g);
            mpz_divexact(g, c, g);
            append_mpz_to_array(todo, g);
        } else {
            append_mpz_to_array(s->hard, c);
        }
    }
    mpz_clear(c);
    mpz_clear(g);
    free_mpz_array(f);
    free_uint32_array(m);
    free_mpz_array(todo);
}
//-----------------------------------------------------------------------------
//
// Factors the smallest composite left on either side: by SIQS if the ECM
// pre-pass searched it, otherwise as tifa_factor does. Returns false if there
// was none.
//
static bool factor_smallest_hard(side_t* const sides, const uint32_t nsides) {
    side_t*  side = NULL;
    uint32_t at   = 0;
    for (uint32_t i = 0; i < nsides; i++) {
        for (uint32_t j = 0; j < sides[i].hard->length; j++) {
            const mpz_ptr c = sides[i].hard->data[j];
            if (   (side == NULL)
                || (mpz_cmp(c, side->hard->data[at]) < 0)) {
                side = &sides[i];
                at   = j;
            }
        }
    }
    if (side == NULL) {
        return false;
    }
    mpz_t c;
    mpz_init_set(c, side->hard->data[at]);
    mpz_array_t* const hard = side->hard;
    mpz_swap(hard->data[at], hard->data[hard->length - 1]);
    hard->length--;
#if TIFA_HAVE_MONT_ECM
    const bool searched = mont_ecm_handles(c);
#else
    const bool searched = false;
#endif
    factor_fully(side, c, searched ? TIFA_FACTOR_BALANCED : TIFA_FACTOR_RANDOM);
    mpz_clear(c);
    return true;
}

//------------------------------------------------------------------------------
//                                 Proofs
//------------------------------------------------------------------------------

//
// The smallest prime of the side not taken yet (NULL if none).
//
static mpz_ptr smallest_left(const side_t* const s,
                             const mpz_array_t* const taken) {
    mpz_ptr q = NULL;
    for (uint32_t i = 0; i < s->primes->length; i++) {
        const mpz_ptr c = s->primes->data[i];
        bool is_taken = false;
        for (uint32_t j = 0; !is_taken && (j < taken->length); j++) {
            is_taken = (mpz_cmp(taken->data[j], c) == 0);
        }
        if (!is_taken && ((q == NULL) || (mpz_cmp(c, q) < 0))) {
            q = c;
        }
    }
    return q;
}
//-----------------------------------------------------------------------------
//
// Takes the smallest primes found, of p - 1 only (nsides = 1) or of both
// sides, until p - 1 alone gives F1^3 > p or both give lcm(F1, F2)^2 > p.
//
static bool choose_primes(mpz_array_t* const qs, mpz_array_t* const rs,
                          const mpz_t p, const side_t* const sides,
                          const uint32_t nsides) {
    mpz_t f1, f2, g, t;
    mpz_init(f1);
    mpz_init(f2);
    mpz_init(g);
    mpz_init(t);
    bool enough = false;

    for (uint32_t pass = 1; !enough && (pass <= nsides); pass++) {
        reset_mpz_array(qs);
        reset_mpz_array(rs);
        mpz_set_ui(f1, 1);
        mpz_set_ui(f2, 1);
        mpz_set_ui(g, 1);
        while (!enough) {
            const mpz_ptr q = smallest_left(&sides[0], qs);
            const mpz_ptr r = (pass == 2) ? smallest_left(&sides[1], rs) : NULL;
            if ((q == NULL) && (r == NULL)) {
                break;
            }
            if ((r == NULL) || ((q != NULL) && (mpz_cmp(q, r) <= 0))) {
                append_mpz_to_array(qs, q);
                part_of(f1, sides[0].x, qs);
            } else {
                append_mpz_to_array(rs, r);
                part_of(f2, sides[1].x, rs);
            }
            mpz_lcm(g, f1, f2);
            mpz_pow_ui(t, f1, 3);
            enough = (mpz_cmp(t, p) > 0);
            if (!enough && (pass == 2)) {
                mpz_mul(t, g, g);
                enough = (mpz_cmp(t, p) > 0);
            }
        }
    }
    mpz_clear(f1);
    mpz_clear(f2);
    mpz_clear(g);
    mpz_clear(t);
    return enough;
}
//-----------------------------------------------------------------------------
//
// D for the Lucas sequences: the first of 5, -7, 9, -11, ... with
// (D / p) = -1, as in the Baillie-PSW test, or 0 if there is none below
// MAX_D (*composite is set if one shares a factor with p).
//
static int32_t find_d(const mpz_t p, bool* const composite) {
    for (int32_t a = 5; a < MAX_D; a += 2) {
        const int32_t d = ((a / 2) % 2 == 0) ? a : -a;
        const int k = mpz_si_kronecker(d, p);
        if (k == -1) {
            return d;
        }
        if (k == 0) {
            *composite = true;
            return 0;
        }
    }
    return 0;
}
//-----------------------------------------------------------------------------
static prime_cert_status_t find_witnesses(uint32_array_t* const aw,
                                          int32_t* const d,
                                          uint32_array_t* const pw,
                                          const mpz_t p, const mpz_t pm1,
                                          const mpz_t pp1,
                                          const mpz_array_t* const qs,
                                          const mpz_array_t* const rs) {
    for (uint32_t i = 0; i < qs->length; i++) {
        int outcome = 0;
        uint32_t a  = 1;
        while ((outcome == 0) && (a < MAX_WITNESS)) {
            a++;
            outcome = try_witness(p, pm1, qs->data[i], a);
        }
        if (outcome != 1) {
            return (outcome == -1) ? NOT_PRIME : PRIME_NOT_PROVEN;
        }
        append_uint32_to_array(aw, a);
    }
    *d = 0;
    if (rs->length == 0) {
        return PRIME_PROVEN;
    }
    bool composite = false;
    *d = find_d(p, &composite);
    if (*d == 0) {
        return composite ? NOT_PRIME : PRIME_NOT_PROVEN;
    }
    for (uint32_t i = 0; i < rs->length; i++) {
        int outcome = 0;
        uint32_t pp = 1;
        while ((outcome == 0) && (pp < MAX_WITNESS)) {
            outcome = try_lucas(p, pp1, rs->data[i], *d, pp);
            pp += (outcome == 0) ? 2 : 0;
        }
        if (outcome != 1) {
            return (outcome == -1) ? NOT_PRIME : PRIME_NOT_PROVEN;
        }
        append_uint32_to_array(pw, pp);
    }
    return PRIME_PROVEN;
}
//-----------------------------------------------------------------------------
//
// Appends to cert the steps proving n, after those of the primes they use.
//
static prime_cert_status_t prove(prime_cert_t* const cert, const mpz_t n) {
    if (mpz_sizeinbase(n, 2) <= 64) {
        return is_small_prime(n) ? PRIME_PROVEN : NOT_PRIME;
    }
    if (has_step(cert, n)) {
        return PRIME_PROVEN;
    }
    if (!is_prime_mpz(n)) {
        return NOT_PRIME;
    }
    mpz_t pm1, pp1, f1, f2;
    mpz_init(pm1);
    mpz_init(pp1);
    mpz_init(f1);
    mpz_init(f2);
    mpz_sub_ui(pm1, n, 1);
    mpz_add_ui(pp1, n, 1);

    side_t sides[2];
    init_side(&sides[0], pm1);
    init_side(&sides[1], pp1);
    mpz_array_t*    qs = alloc_mpz_array(ELONGATION);
    mpz_array_t*    rs = alloc_mpz_array(ELONGATION);
    uint32_array_t* aw = alloc_uint32_array(ELONGATION);
    uint32_array_t* pw = alloc_uint32_array(ELONGATION);

    search_side(&sides[0]);
    uint32_t nsides = 1;
    bool enough = choose_primes(qs, rs, n, sides, 1);
    if (!enough) {
        search_side(&sides[1]);
        nsides = 2;
        enough = choose_primes(qs, rs, n, sides, 2);
    }
    while (!enough && factor_smallest_hard(sides, nsides)) {
        enough = choose_primes(qs, rs, n, sides, nsides);
    }
    prime_cert_status_t status = PRIME_NOT_PROVEN;
    int32_t d = 0;
    if (enough) {
        status = find_witnesses(aw, &d, pw, n, pm1, pp1, qs, rs);
    }
    if (status == PRIME_PROVEN) {
        part_of(f1, pm1, qs);
        part_of(f2, pp1, rs);
        if (!large_enough(n, pm1, f1, f2)) {
            status = NOT_PRIME;
        }
    }
    for (uint32_t i = 0; (status == PRIME_PROVEN) && (i < qs->length); i++) {
        if (prove(cert, qs->data[i]) != PRIME_PROVEN) {
            status = PRIME_NOT_PROVEN;
        }
    }
    for (uint32_t i = 0; (status == PRIME_PROVEN) && (i < rs->length); i++) {
        if (prove(cert, rs->data[i]) != PRIME_PROVEN) {
            status = PRIME_NOT_PROVEN;
        }
    }
    if (status == PRIME_PROVEN) {
        append_step(cert, n, qs, aw, d, rs, pw);
    }
    clear_side(&sides[0]);
    clear_side(&sides[1]);
    free_mpz_array(qs);
    free_mpz_array(rs);
    free_uint32_array(aw);
    free_uint32_array(pw);
    mpz_clear(pm1);
    mpz_clear(pp1);
    mpz_clear(f1);
    mpz_clear(f2);
    return status;
}
//-----------------------------------------------------------------------------
prime_cert_status_t tifa_prove_prime_impl(prime_cert_t* cert, const mpz_t n) {
    clear_prime_cert(cert);
    const prime_cert_status_t status = prove(cert, n);
    if (status != PRIME_PROVEN) {
        clear_prime_cert(cert);
    }
    return status;
}
//-----------------------------------------------------------------------------
typedef struct {
    char*  str;
    size_t size;
    size_t length;
} text_t;
//-----------------------------------------------------------------------------
static void append_text(text_t* const text, const char* const format, ...) {
    const bool room = (text->length < text->size);
    va_list ap;
    va_start(ap, format);
    const int n = gmp_vsnprintf(room ? text->str + text->length : NULL,
                                room ? text->size - text->length : 0,
                                format, ap);
    va_end(ap);
    text->length += (size_t) n;
}
//-----------------------------------------------------------------------------
//
// Arrays of different lengths, which tifa_check_prime_cert refuses, give
// pairs of one number, which tifa_prime_cert_from_str refuses.
//
static void append_pairs(text_t* const text, const mpz_array_t* const primes,
                         const uint32_array_t* const witnesses) {
    append_text(text, "[");
    const uint32_t length = (primes->length > witnesses->length)
                          ? primes->length : witnesses->length;
    for (uint32_t i = 0; i < length; i++) {
        append_text(text, (i == 0) ? "[" : ", [");
        if (i < primes->length) {
            append_text(text, "%Zd", primes->data[i]);
        }
        if (i < witnesses->length) {
            append_text(text, "%s%" PRIu32, (i < primes->length) ? ", " : "",
                        witnesses->data[i]);
        }
        append_text(text, "]");
    }
    append_text(text, "]");
}
//-----------------------------------------------------------------------------
size_t tifa_prime_cert_to_str(char* str, size_t size,
                              const prime_cert_t* cert) {
    text_t text = { str, size, 0 };
    append_text(&text, "[");
    for (uint32_t i = 0; i < cert->nsteps; i++) {
        const prime_cert_step_t* const step = &cert->steps[i];
        append_text(&text, "%s[%Zd, ", (i == 0) ? "" : ",\n ", step->p);
        append_pairs(&text, step->pm1_primes, step->pm1_witnesses);
        if (   (step->pp1_d != 0) || (step->pp1_primes->length != 0)
            || (step->pp1_witnesses->length != 0)) {
            append_text(&text, ", %" PRId32 ", ", step->pp1_d);
            append_pairs(&text, step->pp1_primes, step->pp1_witnesses);
        }
        append_text(&text, "]");
    }
    append_text(&text, "]");
    return text.length;
}
//-----------------------------------------------------------------------------
typedef struct {
    const char*     s;
    mpz_t           p;
    mpz_t           x;
    mpz_array_t*    primes[2];
    uint32_array_t* witnesses[2];
} reader_t;
//-----------------------------------------------------------------------------
static const char* skip_spaces(const char* s) {
    while ((*s == ' ') || (*s == '\t') || (*s == '\n') || (*s == '\r')) {
        s++;
    }
    return s;
}
//-----------------------------------------------------------------------------
//
// Each read leaves r->s on the next character that is not a space; accept
// does not move it when c is not there.
//
static bool accept(reader_t* const r, const char c) {
    if (*r->s != c) {
        return false;
    }
    r->s = skip_spaces(r->s + 1);
    return true;
}
//-----------------------------------------------------------------------------
//
// mpz_set_str rather than digit by digit, which takes time in length^2.
//
static bool read_mpz(reader_t* const r, mpz_t x) {
    size_t length = 0;
    while (isdigit((unsigned char) r->s[length])) {
        length++;
    }
    if (length == 0) {
        return false;
    }
    char  buffer[128];
    char* const digits = (length < sizeof(buffer)) ? buffer
                                                   : tifa_malloc(length + 1);
    memcpy(digits, r->s, length);
    digits[length] = '\0';
    mpz_set_str(x, digits, 10);
    if (digits != buffer) {
        free(digits);
    }
    r->s = skip_spaces(r->s + length);
    return true;
}
//-----------------------------------------------------------------------------
static bool read_int(reader_t* const r, const int64_t min, const int64_t max,
                     int64_t* const v) {
    const bool negative = (*r->s == '-');
    r->s += negative;
    if (!isdigit((unsigned char) *r->s)) {
        return false;
    }
    int64_t x = 0;
    for (; isdigit((unsigned char) *r->s); r->s++) {
        x = 10 * x + (*r->s - '0');
        if (x > max - min) {
            return false;
        }
    }
    *v = negative ? -x : x;
    r->s = skip_spaces(r->s);
    return (*v >= min) && (*v <= max);
}
//-----------------------------------------------------------------------------
static bool read_pairs(reader_t* const r, const int side) {
    reset_mpz_array(r->primes[side]);
    reset_uint32_array(r->witnesses[side]);
    if (!accept(r, '[')) {
        return false;
    }
    if (accept(r, ']')) {
        return true;
    }
    do {
        int64_t w;
        if (   !accept(r, '[') || !read_mpz(r, r->x) || !accept(r, ',')
            || !read_int(r, 0, UINT32_MAX, &w) || !accept(r, ']')) {
            return false;
        }
        append_mpz_to_array(r->primes[side], r->x);
        append_uint32_to_array(r->witnesses[side], (uint32_t) w);
    } while (accept(r, ','));
    return accept(r, ']');
}
//-----------------------------------------------------------------------------
//
// 0 means no D.
//
static bool read_step(reader_t* const r, prime_cert_t* const cert) {
    int64_t d = 0;
    if (   !accept(r, '[') || !read_mpz(r, r->p) || !accept(r, ',')
        || !read_pairs(r, 0)) {
        return false;
    }
    if (accept(r, ',')) {
        if (   !read_int(r, -INT32_MAX, INT32_MAX, &d) || (d == 0)
            || !accept(r, ',') || !read_pairs(r, 1)) {
            return false;
        }
    } else {
        reset_mpz_array(r->primes[1]);
        reset_uint32_array(r->witnesses[1]);
    }
    if (!accept(r, ']')) {
        return false;
    }
    append_step(cert, r->p, r->primes[0], r->witnesses[0], (int32_t) d,
                r->primes[1], r->witnesses[1]);
    return true;
}
//-----------------------------------------------------------------------------
bool tifa_prime_cert_from_str(prime_cert_t* cert, const char* str) {
    clear_prime_cert(cert);

    reader_t r;
    r.s = skip_spaces(str);
    mpz_init(r.p);
    mpz_init(r.x);
    for (int i = 0; i < 2; i++) {
        r.primes[i]    = alloc_mpz_array(8);
        r.witnesses[i] = alloc_uint32_array(8);
    }
    bool ok = accept(&r, '[');
    if (ok && !accept(&r, ']')) {
        do {
            ok = read_step(&r, cert);
        } while (ok && accept(&r, ','));
        ok = ok && accept(&r, ']');
    }
    ok = ok && (*r.s == '\0');

    if (!ok) {
        clear_prime_cert(cert);
    }
    mpz_clear(r.p);
    mpz_clear(r.x);
    for (int i = 0; i < 2; i++) {
        free_mpz_array(r.primes[i]);
        free_uint32_array(r.witnesses[i]);
    }
    return ok;
}
