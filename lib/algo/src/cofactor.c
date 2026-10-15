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
#include "mont_ecm.h"
#include "siqs.h"
#include "tifa_factor.h"

//------------------------------------------------------------------------------
//                      CONSTANT DATA, SYMBOLS AND MACROS
//------------------------------------------------------------------------------
//
// Primes divided out first when the caller knows of no factor base bound.
//
#define COFACTOR_TDIV_BOUND 1000

#if TIFA_HAVE_MONT_ECM
//
// ECM levels (as in ecm.c): the B1 for factors of fbits bits (so for the
// range (fbits - 5, fbits]), the mean number of curves that takes, and the
// time of one curve on two limbs (us, Apple M series), which only matters
// relative to siqs_cost.
//
typedef struct {
    uint32_t fbits;
    uint32_t b1;
    uint32_t ncurves;
    double   cost;
} level_t;

static const level_t levels[] = {
    { 20,   100,  2,   17.0 },
    { 25,   200,  3,   31.0 },
    { 30,   400,  4,   58.0 },
    { 35,   800,  7,  112.0 },
    { 40,  1600,  8,  220.0 },
    { 45,  3200, 10,  410.0 },
    { 50,  6400, 15,  800.0 },
    { 55, 12800, 17, 1600.0 },
    { 60, 25600, 25, 3200.0 }
};
#define NLEVELS (sizeof(levels) / sizeof(levels[0]))
//
// Curve time on n limbs relative to two limbs.
//
static const double limb_factor[MONT_ECM_MAX_LIMBS + 1] = {
    0.0, 0.6, 1.0, 2.0, 3.3
};
//
// Rounds of curves at the last level when ECM is the last resort (see
// max_siqs_bits): finds a factor of the size sought 86% of the time.
//
#define ROUNDS_AT_TOP 2
#define PM1_B1 315
#define PM1_B2 (50 * PM1_B1)
//
// Time of tifa_factor on a balanced semiprime (us, same machine).
//
static const struct { double bits; double cost; } siqs_costs[] = {
    {  64,     200 }, {  72,    320 }, {  80,    480 }, {  90,    740 },
    { 100,    1400 }, { 110,   2600 }, { 120,   4500 }, { 128,   7300 },
    { 140,   15000 }, { 160,  60000 }, { 200, 1000000 }
};
#define NSIQS_COSTS (sizeof(siqs_costs) / sizeof(siqs_costs[0]))
#endif

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
    mont_ecm_plan_t* pm1_plan;
    uint32_t         nlevels;
    mont_ecm_plan_t* plans[NLEVELS];
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
    uint32_t        sigma;
    mpz_t           tmp;
} cofactor_run_t;

//------------------------------------------------------------------------------
//                          STRATEGY ALLOCATION
//------------------------------------------------------------------------------
cofactor_strategy_t* alloc_cofactor_strategy(const cofactor_params_t* params) {
    if ((params->lpb == 0) || (params->lpb > 63)) {
        return NULL;
    }
    cofactor_strategy_t* const s = malloc(sizeof(cofactor_strategy_t));

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
    s->pm1_plan = alloc_mont_ecm_plan(PM1_B1, PM1_B2);
    //
    // Levels up to the one whose range holds factors of lpb bits: larger
    // factors are not sought.
    //
    s->nlevels = 0;
    while (   (s->nlevels < NLEVELS)
           && (levels[s->nlevels].fbits - 5 < s->lpb)) {
        s->plans[s->nlevels] = alloc_mont_ecm_plan(
                                   levels[s->nlevels].b1,
                                   50 * levels[s->nlevels].b1
                               );
        s->nlevels++;
    }
#endif
    return s;
}
//------------------------------------------------------------------------------
void free_cofactor_strategy(cofactor_strategy_t* strategy) {
    if (strategy == NULL) {
        return;
    }
#if TIFA_HAVE_MONT_ECM
    free_mont_ecm_plan(strategy->pm1_plan);
    for (uint32_t i = 0; i < strategy->nlevels; i++) {
        free_mont_ecm_plan(strategy->plans[i]);
    }
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
// can hold, exceed the limit.
//
static bool too_many_large_primes(const cofactor_run_t* const run) {
    const uint32_t limit = run->strategy->max_large_primes;
    if (limit == 0) {
        return false;
    }
    uint32_t total = run->nlarge;
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
        if (too_many_large_primes(run)) {
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
    const bool prime     = (mpz_cmp(c, s->bound2) < 0) || is_prime_mpz(c);

    if (prime) {
        return below_lpb ? add_prime(run, c, e) : COFACTOR_NOT_SMOOTH;
    }
    if (min_nprimes(run, c, true) > max_nprimes(run, c)) {
        return COFACTOR_NOT_SMOOTH;
    }
    mpz_set(run->pending[run->npending], c);
    run->pmultis[run->npending] = e;
    run->npending++;

    return too_many_large_primes(run) ? COFACTOR_NOT_SMOOTH : COFACTOR_SMOOTH;
}
//------------------------------------------------------------------------------
#if TIFA_HAVE_MONT_ECM
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
#if TIFA_HAVE_MONT_ECM
static double siqs_cost(double bits) {
    if (bits <= siqs_costs[0].bits) {
        return siqs_costs[0].cost;
    }
    for (uint32_t i = 1; i < NSIQS_COSTS; i++) {
        if (bits <= siqs_costs[i].bits) {
            const double t = (bits - siqs_costs[i - 1].bits)
                           / (siqs_costs[i].bits - siqs_costs[i - 1].bits);
            return exp(  (1 - t) * log(siqs_costs[i - 1].cost)
                       + t * log(siqs_costs[i].cost));
        }
    }
    return siqs_costs[NSIQS_COSTS - 1].cost;
}
//------------------------------------------------------------------------------
//
// p-1, then ECM curves level by level up to the size the smallest prime
// factor of c has at most if c splits as required. With a SIQS run to fall
// back on, only while the curves cost less than that run would; without,
// the last level is run ROUNDS_AT_TOP times. Returns true if g is a proper
// factor of c.
//
static bool split_by_ecm(cofactor_run_t* const run, mpz_t g, const mpz_t c,
                         bool last_resort) {
    const cofactor_strategy_t* const s = run->strategy;

    if (mont_pm1(g, c, 2, s->pm1_plan) == 1) {
        return true;
    }
    const double bits   = (double) mpz_sizeinbase(c, 2);
    const double budget = siqs_cost(bits);
    const double lfac   = limb_factor[mpz_size(c)];
    const double target = fmin((double) s->lpb,
                               bits / min_nprimes(run, c, true));
    double spent = 0;

    for (uint32_t i = 0; i < s->nlevels; i++) {
        const double hi = (double) levels[i].fbits;
        if (hi <= s->bound_bits) {
            continue;
        }
        if (hi - 5 >= target) {
            break;
        }
        spent += levels[i].ncurves * levels[i].cost * lfac;
        if (!last_resort && (spent > budget)) {
            break;
        }
        const uint32_t ncurves = (last_resort && (hi >= target))
                               ? ROUNDS_AT_TOP * levels[i].ncurves
                               : levels[i].ncurves;
        for (uint32_t j = 0; j < ncurves; j++) {
            if (mont_ecm_curve(g, c, run->sigma++, s->plans[i])) {
                return true;
            }
        }
    }
    return false;
}
#endif
//------------------------------------------------------------------------------
static cofactor_status_t split(cofactor_run_t* const run, const mpz_t c,
                               uint32_t e) {
    const cofactor_strategy_t* const s = run->strategy;
    const uint32_t bits = (uint32_t) mpz_sizeinbase(c, 2);

#if TIFA_HAVE_FACTOR_U64
    if (bits <= 64) {
        uint64_t primes[FACTOR_U64_MAX_PRIMES];
        uint32_t pmultis[FACTOR_U64_MAX_PRIMES];
        const uint32_t np = factor_u64((uint64_t) mpz_get_ui(c),
                                       primes, pmultis);
        if (np != 0) {
            cofactor_status_t status = COFACTOR_SMOOTH;
            mpz_t p;
            mpz_init(p);
            for (uint32_t i = 0; (i < np) && (status == COFACTOR_SMOOTH);
                 i++) {
                mpz_set_ui(p, (unsigned long int) primes[i]);
                status = add_prime(run, p, e * pmultis[i]);
            }
            mpz_clear(p);
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
        cofactor_status_t status = COFACTOR_GAVE_UP;
        if (split_by_ecm(run, g, c, no_siqs)) {
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
                                       : malloc(max_pending * sizeof(mpz_t));
    run.pmultis  = (max_pending <= 16) ? pmultis_buf
                                       : malloc(max_pending * sizeof(uint32_t));
    run.npending = 0;
    run.nlarge   = 0;
    run.sigma    = 6;
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
