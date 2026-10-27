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
 * \file    cofactor.c
 *
 * \brief Cofactorization: does an integer split into primes below a bound?
 */

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "tifa_config.h"
#include "cofactor.h"
#include "factor_u64.h"
#include "first_primes.h"
#include "funcs.h"
#include "ecm_schedule.h"
#include "siqs.h"
#include "tifa_factor.h"
#include "tifa_alloc.h"

//------------------------------------------------------------------------------
//                      CONSTANT DATA, SYMBOLS AND MACROS
//------------------------------------------------------------------------------
//
// Primes divided out first when the caller knows of no factor base bound.
//
#define COFACTOR_TDIV_BOUND 1000

struct struct_cofactor_strategy_t {
    uint32_t lpb;
    mpz_t    fb_bound;
    uint32_t max_large_primes;
    uint32_t max_siqs_bits;
    //
    // No prime factor below bound (fb_bound, or COFACTOR_TDIV_BOUND when the
    // primes below it are divided out first), its square and log2.
    //
    bool     tdiv;
    mpz_t    bound;
    mpz_t    bound2;
    double   bound_bits;
#if TIFA_HAVE_MONT_ECM
    ecm_schedule_t* schedule;
#endif
};

//
// The state of one tifa_cofactor call: the primes found so far, and the
// composites left to split (each with its multiplicity in n).
//
typedef struct {
    const cofactor_strategy_t* strategy;
    mpz_array_t*    factors;
    uint32_array_t* multis;
    mpz_t*          pending;
    uint32_t*       pmultis;
    uint32_t        npending;
    uint32_t        nlarge;
    mpz_t           tmp;
} cofactor_run_t;

//------------------------------------------------------------------------------
//                          STRATEGY ALLOCATION
//------------------------------------------------------------------------------
cofactor_strategy_t* alloc_cofactor_strategy(const cofactor_params_t* params) {
    if ((params->lpb == 0) || (params->lpb > 63)) {
        return NULL;
    }
    cofactor_strategy_t* const s = tifa_malloc(sizeof(cofactor_strategy_t));

    s->lpb              = params->lpb;
    s->max_large_primes = params->max_large_primes;
    s->max_siqs_bits    = params->max_siqs_bits;
    s->tdiv             = (params->fb_bound < COFACTOR_TDIV_BOUND);

    mpz_init(s->fb_bound);
    mpz_import(s->fb_bound, 1, 1, sizeof(uint64_t), 0, 0, &params->fb_bound);
    mpz_init_set(s->bound, s->fb_bound);
    if (s->tdiv) {
        mpz_set_ui(s->bound, COFACTOR_TDIV_BOUND);
    }
    mpz_init(s->bound2);
    mpz_mul(s->bound2, s->bound, s->bound);
    s->bound_bits = log2(mpz_get_d(s->bound));

#if TIFA_HAVE_MONT_ECM
    s->schedule = alloc_ecm_schedule(s->lpb);
#endif
    return s;
}
//------------------------------------------------------------------------------
void free_cofactor_strategy(cofactor_strategy_t* strategy) {
    if (strategy == NULL) {
        return;
    }
#if TIFA_HAVE_MONT_ECM
    free_ecm_schedule(strategy->schedule);
#endif
    mpz_clear(strategy->fb_bound);
    mpz_clear(strategy->bound);
    mpz_clear(strategy->bound2);
    free(strategy);
}

//------------------------------------------------------------------------------
//                       BOOKKEEPING OF A RUN
//------------------------------------------------------------------------------
//
// Fewest primes below 2^lpb whose product can be c (c >= 2).
//
static uint32_t min_nprimes(const cofactor_run_t* const run, const mpz_t c,
                            bool composite) {
    const uint32_t lpb = run->strategy->lpb;
    const uint32_t m   = ((uint32_t) mpz_sizeinbase(c, 2) + lpb - 1) / lpb;
    return (composite && (m < 2)) ? 2 : m;
}
//------------------------------------------------------------------------------
//
// Most prime factors c can have, all of them being at least the bound.
//
static uint32_t max_nprimes(const cofactor_run_t* const run, const mpz_t c) {
    long int e;
    const double d = mpz_get_d_2exp(&e, c);
    return (uint32_t) floor((e + log2(d)) / run->strategy->bound_bits + 1e-9);
}
//------------------------------------------------------------------------------
//
// Whether the large primes found, plus the fewest the pending composites
// can hold, plus extra, exceed the limit.
//
static bool too_many_large_primes(const cofactor_run_t* const run,
                                  uint32_t extra) {
    const uint32_t limit = run->strategy->max_large_primes;
    if (limit == 0) {
        return false;
    }
    uint32_t total = run->nlarge + extra;
    for (uint32_t i = 0; i < run->npending; i++) {
        total += run->pmultis[i] * min_nprimes(run, run->pending[i], true);
    }
    return total > limit;
}
//------------------------------------------------------------------------------
static cofactor_status_t add_prime(cofactor_run_t* const run, const mpz_t p,
                                   uint32_t e) {
    const cofactor_strategy_t* const s = run->strategy;

    if (mpz_sizeinbase(p, 2) > s->lpb) {
        return COFACTOR_NOT_SMOOTH;
    }
    const uint32_t i = index_in_mpz_array(p, run->factors);
    if (i == NOT_IN_ARRAY) {
        append_mpz_to_array(run->factors, p);
        append_uint32_to_array(run->multis, e);
    } else {
        run->multis->data[i] += e;
    }
    if (mpz_cmp(p, s->fb_bound) >= 0) {
        run->nlarge += e;
        if (too_many_large_primes(run, 0)) {
            return COFACTOR_NOT_SMOOTH;
        }
    }
    return COFACTOR_SMOOTH;
}
//------------------------------------------------------------------------------
//
// Takes in a part c of n (with no prime factor below the bound) with
// multiplicity e: records it if prime, queues it if composite, and decides
// whatever can be decided without splitting it.
//
static cofactor_status_t add_part(cofactor_run_t* const run, const mpz_t c,
                                  uint32_t e) {
    const cofactor_strategy_t* const s = run->strategy;

    if (mpz_cmp_ui(c, 1) == 0) {
        return COFACTOR_SMOOTH;
    }
    const bool below_lpb = (mpz_sizeinbase(c, 2) <= s->lpb);
    if (mpz_cmp(c, s->bound2) < 0) {
        return below_lpb ? add_prime(run, c, e) : COFACTOR_NOT_SMOOTH;
    }
    //
    // Test primality only if the answer depends on it: a prime too large or
    // a composite that cannot split as allowed is not smooth either way.
    //
    const uint32_t m = min_nprimes(run, c, true);
    const bool may_split = (m <= max_nprimes(run, c))
                           && !too_many_large_primes(run, e * m);
    if (!below_lpb && !may_split) {
        return COFACTOR_NOT_SMOOTH;
    }
    if (is_prime_mpz(c)) {
        return below_lpb ? add_prime(run, c, e) : COFACTOR_NOT_SMOOTH;
    }
    if (!may_split) {
        return COFACTOR_NOT_SMOOTH;
    }
    mpz_set(run->pending[run->npending], c);
    run->pmultis[run->npending] = e;
    run->npending++;

    return COFACTOR_SMOOTH;
}
//------------------------------------------------------------------------------
#if TIFA_HAVE_FACTOR_U64 || TIFA_HAVE_MONT_ECM
static cofactor_status_t add_split(cofactor_run_t* const run, const mpz_t c,
                                   const mpz_t g, uint32_t e) {
    mpz_divexact(run->tmp, c, g);
    //
    // The larger part first: it is the likelier to be a prime too large.
    //
    const bool g_first = (mpz_cmp(g, run->tmp) > 0);
    mpz_t other;
    mpz_init_set(other, g_first ? run->tmp : g);

    cofactor_status_t status = add_part(run, g_first ? g : run->tmp, e);
    if (status == COFACTOR_SMOOTH) {
        status = add_part(run, other, e);
    }
    mpz_clear(other);
    return status;
}
#endif

//------------------------------------------------------------------------------
//                              SPLITTING
//------------------------------------------------------------------------------
//
// Factors c completely with SIQS (with tifa_factor below 64 bits, or should
// SIQS fail).
//
static cofactor_status_t split_by_siqs(cofactor_run_t* const run,
                                       const mpz_t c, uint32_t e) {
    mpz_array_t*    f = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);

    ecode_t ecode = NO_FACTOR_FOUND;
    if (mpz_sizeinbase(c, 2) > 64) {
        siqs_params_t params;
        set_siqs_params_to_default(c, &params);
        ecode = siqs(f, m, c, &params, FIND_COMPLETE_FACTORIZATION);
    }
    if (ecode != COMPLETE_FACTORIZATION_FOUND) {
        reset_mpz_array(f);
        m->length = 0;
        ecode = tifa_factor(f, m, c, FIND_COMPLETE_FACTORIZATION);
    }
    cofactor_status_t status = COFACTOR_GAVE_UP;
    if (ecode == COMPLETE_FACTORIZATION_FOUND) {
        status = COFACTOR_SMOOTH;
        for (uint32_t i = 0; (i < f->length) && (status == COFACTOR_SMOOTH);
             i++) {
            status = add_prime(run, f->data[i], e * m->data[i]);
        }
    }
    free_mpz_array(f);
    free_uint32_array(m);
    return status;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static cofactor_status_t split(cofactor_run_t* const run, const mpz_t c,
                               uint32_t e) {
    const cofactor_strategy_t* const s = run->strategy;
    const uint32_t bits = (uint32_t) mpz_sizeinbase(c, 2);

#if TIFA_HAVE_FACTOR_U64
    if (bits <= 64) {
        const uint64_t g = split_u64((uint64_t) mpz_get_ui(c));
        if (g != 0) {
            mpz_t gz;
            mpz_init_set_ui(gz, (unsigned long int) g);
            const cofactor_status_t status = add_split(run, c, gz, e);
            mpz_clear(gz);
            return status;
        }
    }
#endif
    if (mpz_perfect_power_p(c)) {
        //
        // The largest exponent first, so that the root is not a power.
        //
        mpz_t r;
        mpz_init(r);
        for (uint32_t x = (uint32_t) (bits / s->bound_bits) + 1; x >= 2; x--) {
            if (mpz_root(r, c, x)) {
                cofactor_status_t status = add_part(run, r, e * x);
                mpz_clear(r);
                return status;
            }
        }
        mpz_clear(r);
    }
#if TIFA_HAVE_MONT_ECM
    const bool no_siqs = (s->max_siqs_bits != 0) && (bits > s->max_siqs_bits);

    if (mont_ecm_handles(c)) {
        mpz_t g;
        mpz_init(g);
        //
        // If c splits as allowed, its smallest prime factor has at most
        // target bits. ECM may take as long as SIQS would: past that, SIQS
        // is cheaper, or c is given up.
        //
        const double target = fmin((double) s->lpb,
                                   (double) bits
                                   / min_nprimes(run, c, true));
        cofactor_status_t status = COFACTOR_GAVE_UP;
        if (ecm_schedule_split(g, c, s->schedule, s->bound_bits, target, 1)) {
            status = add_split(run, c, g, e);
        } else if (!no_siqs) {
            status = split_by_siqs(run, c, e);
        }
        mpz_clear(g);
        return status;
    }
    if (no_siqs) {
        return COFACTOR_GAVE_UP;
    }
#endif
    return split_by_siqs(run, c, e);
}

//------------------------------------------------------------------------------
//                              PUBLIC ENTRY
//------------------------------------------------------------------------------
cofactor_status_t tifa_cofactor(mpz_array_t* factors, uint32_array_t* multis,
                                const mpz_t n,
                                const cofactor_strategy_t* strategy) {
    reset_mpz_array(factors);
    multis->length = 0;

    if (mpz_sgn(n) <= 0) {
        return COFACTOR_GAVE_UP;
    }
    //
    // Each pending composite has two prime factors of at least
    // COFACTOR_TDIV_BOUND (about 2^10): at most b / 19 of them for b bits.
    //
    const uint32_t max_pending = (uint32_t) mpz_sizeinbase(n, 2) / 19 + 1;
    mpz_t    pending_buf[16];
    uint32_t pmultis_buf[16];

    cofactor_run_t run;
    run.strategy = strategy;
    run.factors  = factors;
    run.multis   = multis;
    run.pending  = (max_pending <= 16) ? pending_buf
                                       : tifa_malloc(max_pending
                                                     * sizeof(mpz_t));
    run.pmultis  = (max_pending <= 16) ? pmultis_buf
                                       : tifa_malloc(max_pending
                                                     * sizeof(uint32_t));
    run.npending = 0;
    run.nlarge   = 0;
    for (uint32_t i = 0; i < max_pending; i++) {
        mpz_init(run.pending[i]);
    }
    mpz_init(run.tmp);

    mpz_t c;
    mpz_init_set(c, n);

    cofactor_status_t status = COFACTOR_SMOOTH;

    if (strategy->tdiv) {
        for (uint32_t i = 0; first_primes[i] < COFACTOR_TDIV_BOUND; i++) {
            uint32_t e = 0;
            while (mpz_divisible_ui_p(c, first_primes[i])) {
                mpz_divexact_ui(c, c, first_primes[i]);
                e++;
            }
            if (e != 0) {
                mpz_set_ui(run.tmp, first_primes[i]);
                status = add_prime(&run, run.tmp, e);
                if (status != COFACTOR_SMOOTH) {
                    goto done;
                }
            }
        }
    }
    status = add_part(&run, c, 1);

    while ((status == COFACTOR_SMOOTH) && (run.npending != 0)) {
        run.npending--;
        mpz_swap(c, run.pending[run.npending]);
        status = split(&run, c, run.pmultis[run.npending]);
    }

done:
    mpz_clear(c);
    mpz_clear(run.tmp);
    for (uint32_t i = 0; i < max_pending; i++) {
        mpz_clear(run.pending[i]);
    }
    if (max_pending > 16) {
        free(run.pending);
        free(run.pmultis);
    }

    return status;
}
