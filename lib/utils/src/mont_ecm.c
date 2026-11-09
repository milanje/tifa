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
 * \file    mont_ecm.c
 *
 * \brief ECM curves in fixed-size Montgomery arithmetic, for moduli of one
 *        to four 64-bit limbs (five to eight: mont_ecm_wide.c).
 */

#include <stdlib.h>
#include <string.h>

#include "mont_ecm.h"

#if TIFA_HAVE_MONT_ECM

#include "first_primes.h"
#include "z2z8_curves.h"
#include "tifa_alloc.h"
#include "stopping.h"

#define MONT_MOD_LIMBS MONT_ECM_NARROW_LIMBS

#define MONT_ECM_Z12_MIN_B1 1600
#include "mont_ecm_internal.h"

//------------------------------------------------------------------------------
//                          Kernels for 1 to 4 limbs
//------------------------------------------------------------------------------

#define NL 1
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL
#define NL 2
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL
#define NL 3
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL
#define NL 4
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL

//------------------------------------------------------------------------------
//                                  Plans
//------------------------------------------------------------------------------

typedef struct {
    mont_ecm_plan_t* plan;
    uint16_t*        pairs;
    uint32_t         npairs;
    uint32_t         alloced;
    uint8_t*         needed;
    uint32_t         cur_m;
} stage2_builder_t;

static void flush_giant_step(stage2_builder_t* const sb) {
    mont_ecm_plan_t* const plan = sb->plan;
    if (sb->alloced - sb->npairs < plan->nbaby) {
        sb->alloced = 2 * sb->alloced + plan->nbaby;
        sb->pairs   = tifa_realloc(sb->pairs, sb->alloced * sizeof(uint16_t));
    }
    for (uint32_t i = 0; i < plan->nbaby; i++) {
        if (sb->needed[i]) {
            sb->pairs[sb->npairs++] = (uint16_t) i;
            sb->needed[i] = 0;
        }
    }
    plan->offsets[sb->cur_m - plan->m_first + 1] = sb->npairs;
}

static void add_stage2_prime(stage2_builder_t* const sb, const uint32_t q) {
    mont_ecm_plan_t* const plan = sb->plan;

    // In 64 bits: q may be within d of 2^32.
    const uint32_t m = (uint32_t) (((uint64_t) q + plan->d / 2) / plan->d);
    while (sb->cur_m < m) {
        flush_giant_step(sb);
        sb->cur_m++;
    }
    const uint64_t c = (uint64_t) m * plan->d;
    const uint32_t j = (uint32_t) ((q > c) ? (q - c) : (c - q));
    sb->needed[plan->jslot[j]] = 1;
}

static bool build_stage2(mont_ecm_plan_t* const plan, const bool stoppable) {
    const uint32_t d    = plan->d;
    const uint32_t half = d / 2;

    plan->jslot = tifa_malloc((half + 1) * sizeof(uint16_t));
    plan->nbaby = 0;
    for (uint32_t j = 0; j <= half; j++) {
        bool coprime = (j & 1) != 0;
        for (uint32_t p = 3; coprime && (p <= 11); p += 2) {
            if ((d % p == 0) && (j % p == 0)) {
                coprime = false;
            }
        }
        plan->jslot[j] = coprime ? (uint16_t) plan->nbaby++ : UINT16_MAX;
    }
    //
    // The primes q in (b1, b2]: q = m.d + j with |j| < d/2 (q is never
    // m.d + d/2, a multiple of d/2).
    //
    plan->m_first = (plan->b1 + 1 + half) / d;
    if (plan->m_first == 0) {
        plan->m_first = 1;
    }
    const uint32_t m_last = (uint32_t) (((uint64_t) plan->b2 + half) / d);
    plan->nm      = m_last - plan->m_first + 1;
    plan->offsets = tifa_calloc(plan->nm + 1, sizeof(uint32_t));

    stage2_builder_t sb;
    sb.plan    = plan;
    sb.pairs   = NULL;
    sb.npairs  = 0;
    sb.alloced = 0;
    sb.needed  = tifa_calloc(plan->nbaby, sizeof(uint8_t));
    sb.cur_m   = plan->m_first;

    prime_iter_t it;
    prime_iter_init(&it, plan->b1 + 1, plan->b2);
    uint32_t nprimes = 0;
    for (uint32_t q; (q = prime_iter_next(&it)) != 0; ) {
        if (   stoppable && ((++nprimes & 0xffff) == 0)
            && tifa_should_stop()) {
            prime_iter_clear(&it);
            free(sb.pairs);
            free(sb.needed);
            return false;
        }
        add_stage2_prime(&sb, q);
    }
    prime_iter_clear(&it);

    while (sb.cur_m <= m_last) {
        flush_giant_step(&sb);
        sb.cur_m++;
    }
    plan->pairs = tifa_realloc(sb.pairs, (sb.npairs + 1) * sizeof(uint16_t));
    free(sb.needed);
    return true;
}

static mont_ecm_plan_t* build_plan(uint32_t b1, uint32_t b2,
                                   const bool stoppable) {
    mont_ecm_plan_t* const plan = tifa_calloc(1, sizeof(mont_ecm_plan_t));

    //
    // d = 2310 needs fewer giant steps, but 577 baby steps and stage 1 up
    // to d/2: worth it only for wide stage 2 ranges.
    //
    uint32_t d = 0;
    if (b2 > b1) {
        d = ((b1 >= 2310 / 2) && (b2 - b1 >= 50 * 2310)) ? 2310 : 210;
        if (b1 < d / 2) {
            b1 = d / 2;
        }
        if (b2 <= b1) {
            d = 0;
        }
    }
    plan->b1 = b1;
    plan->b2 = b2;
    plan->d  = d;
    //
    // Beyond the bounds of the levels for 85-bit factors: up to seconds.
    //
    plan->long_curves = (b1 > LARGEST_PRIME) || (b2 > 50 * LARGEST_PRIME);

    mpz_t e;
    mpz_init(e);
    if (!prime_powers_product(e, b1, stoppable)) {
        mpz_clear(e);
        free_mont_ecm_plan(plan);
        return NULL;
    }
    plan->ebits = mpz_sizeinbase(e, 2);
    plan->e     = tifa_malloc((plan->ebits / 64 + 1) * sizeof(uint64_t));
    mpz_export(plan->e, &plan->ne, -1, sizeof(uint64_t), 0, 0, e);
    mpz_clear(e);

    if ((d != 0) && !build_stage2(plan, stoppable)) {
        free_mont_ecm_plan(plan);
        return NULL;
    }
    return plan;
}

mont_ecm_plan_t* alloc_mont_ecm_plan(const uint32_t b1, const uint32_t b2) {
    return build_plan(b1, b2, false);
}

mont_ecm_plan_t* alloc_mont_ecm_plan_or_stop(const uint32_t b1,
                                             const uint32_t b2) {
    return build_plan(b1, b2, true);
}

void free_mont_ecm_plan(mont_ecm_plan_t* const plan) {
    if (plan != NULL) {
        free(plan->e);
        free(plan->jslot);
        free(plan->offsets);
        free(plan->pairs);
        free(plan);
    }
}

//------------------------------------------------------------------------------
//                                  Curves
//------------------------------------------------------------------------------

bool mont_ecm_handles(const mpz_t n) {
    const size_t nl = mpz_size(n);
    return (nl >= 1) && (nl <= MONT_ECM_MAX_LIMBS) && mpz_odd_p(n)
           && (mpz_sgn(n) > 0) && (mpz_cmp_ui(n, 1) > 0);
}

//
// Sets x0 and a24 for Suyama's curve of parameter sigma: u = sigma^2 - 5,
// v = 4.sigma, starting point (u^3 : v^3) = (u^3/v^3 : 1) on the curve with
// (A + 2)/4 = (v - u)^3 (3u + v) / (16 u^3 v), both from the inverse of
// 16 u^3 v^4. Returns false, with factor = the gcd of that value with n, if
// it is not invertible.
//
static bool suyama_curve(mpz_t x0, mpz_t a24, mpz_t factor, const mpz_t n,
                         const uint32_t sigma) {
    mpz_t u, v, u3, v3, t;
    mpz_inits(u, v, u3, v3, t, NULL);

    mpz_set_ui(u, sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5);
    mpz_mod(u, u, n);
    mpz_set_ui(v, sigma);
    mpz_mul_ui(v, v, 4);
    mpz_mod(v, v, n);
    mpz_powm_ui(u3, u, 3, n);
    mpz_powm_ui(v3, v, 3, n);
    mpz_mul(x0, u3, v);
    mpz_mul_ui(x0, x0, 16);
    mpz_mod(x0, x0, n);
    mpz_mul(t, x0, v3);
    mpz_mod(t, t, n);

    const bool ok = (mpz_invert(factor, t, n) != 0);
    if (!ok) {
        mpz_gcd(factor, t, n);
    } else {
        mpz_set(t, factor);
        mpz_mul(x0, x0, u3);
        mpz_mul(x0, x0, t);
        mpz_mod(x0, x0, n);

        mpz_sub(a24, v, u);
        mpz_powm_ui(a24, a24, 3, n);
        mpz_mul(a24, a24, v3);
        mpz_mul(a24, a24, t);
        mpz_mul_ui(u, u, 3);
        mpz_add(u, u, v);
        mpz_mul(a24, a24, u);
        mpz_mod(a24, a24, n);
    }
    mpz_clears(u, v, u3, v3, t, NULL);
    return ok;
}

//
// One curve, from its starting point (x0z : 1) and a24 (residues mod n).
//
static bool run_curve(mpz_t factor, const mpz_t n, const mpz_t x0z,
                      const mpz_t a24z, const mont_ecm_edwards_t* const ed,
                      const mont_ecm_plan_t* const plan) {
    if (mpz_size(n) > MONT_ECM_NARROW_LIMBS) {
        return mont_ecm_wide_run(factor, n, x0z, a24z, ed, plan);
    }
    mod_t m;
    mod_init(&m, n);

    uint64_t x0[MONT_ECM_NARROW_LIMBS];
    uint64_t a24[MONT_ECM_NARROW_LIMBS];
    to_limbs(x0, x0z, m.nl);
    to_limbs(a24, a24z, m.nl);

    int found = 0;
    switch (m.nl) {
        case 1:
            RUN_CURVE_OR_EDWARDS(1);
            break;
        case 2:
            RUN_CURVE_OR_EDWARDS(2);
            break;
        case 3:
            RUN_CURVE_OR_EDWARDS(3);
            break;
        default:
            RUN_CURVE_OR_EDWARDS(4);
            break;
    }
    return found > 0;
}

//
// Sets up a curve with setup (then runs it), or reports the gcd found if
// setup could not invert.
//
static bool setup_and_run(mpz_t factor, const mpz_t n, const uint32_t k,
                          const mont_ecm_plan_t* const plan,
                          bool (*setup)(mpz_t, mpz_t, mpz_t, const mpz_t,
                                        uint32_t)) {
    mpz_t x0, a24;
    mpz_inits(x0, a24, NULL);
    bool found;
    if (setup(x0, a24, factor, n, k)) {
        //
        // Not in run_curve, which then no longer inlines stage 2 (1% slower
        // on three limbs).
        //
        found = plan->long_curves
                    ? mont_ecm_wide_run(factor, n, x0, a24, NULL, plan)
                    : run_curve(factor, n, x0, a24, NULL, plan);
    } else {
        found = (mpz_cmp_ui(factor, 1) > 0) && (mpz_cmp(factor, n) < 0);
    }
    mpz_clears(x0, a24, NULL);
    return found;
}

bool mont_ecm_curve(mpz_t factor, const mpz_t n, const uint32_t sigma,
                    const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n)) {
        return false;
    }
    return setup_and_run(factor, n, sigma, plan, suyama_curve);
}

//
// Curve i of a family set up by setup, its stage 1 on its Edwards form.
//
static bool run_edwards(mpz_t factor, const mpz_t n, const uint32_t i,
                        const mont_ecm_plan_t* const plan,
                        bool (*setup)(mont_ecm_edwards_t*, mpz_t, mpz_t,
                                      mpz_t, const mpz_t, uint32_t)) {
    mont_ecm_edwards_t ed;
    mpz_t x0, a24;
    mpz_inits(ed.x, ed.y, ed.z, ed.t, ed.d, x0, a24, NULL);
    bool found;
    if (setup(&ed, x0, a24, factor, n, i)) {
        found = plan->long_curves
                    ? mont_ecm_wide_run(factor, n, x0, a24, &ed, plan)
                    : run_curve(factor, n, x0, a24, &ed, plan);
    } else {
        found = (mpz_cmp_ui(factor, 1) > 0) && (mpz_cmp(factor, n) < 0);
    }
    mpz_clears(ed.x, ed.y, ed.z, ed.t, ed.d, x0, a24, NULL);
    return found;
}

bool mont_ecm_nth_curve(mpz_t factor, const mpz_t n, const uint32_t k,
                        const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n)) {
        return false;
    }
    if (k < NZ2Z8_CURVES) {
        return run_edwards(factor, n, k, plan, mont_ecm_edwards_setup);
    }
    //
    // Montgomery's Z/12 curves cost k.P on an auxiliary curve to set up:
    // only on curves long enough for it not to matter.
    //
    if (plan->b1 >= MONT_ECM_Z12_MIN_B1) {
        return run_edwards(factor, n, 2 + k - NZ2Z8_CURVES, plan,
                           mont_ecm_z12_setup);
    }
    return setup_and_run(factor, n, 6 + k - NZ2Z8_CURVES, plan, suyama_curve);
}

//------------------------------------------------------------------------------
//                              p-1 and p+1
//------------------------------------------------------------------------------


int mont_pm1(mpz_t factor, const mpz_t n, const uint32_t base,
             const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n) || (mpz_cmp_ui(n, base) <= 0)) {
        return 0;
    }
    if (mpz_size(n) > MONT_ECM_NARROW_LIMBS) {
        return mont_wide_pm1(factor, n, base, plan);
    }
    mod_t m;
    mod_init(&m, n);

    int found = 0;
    switch (m.nl) {
        case 1:
            RUN_PM1(1);
            break;
        case 2:
            RUN_PM1(2);
            break;
        case 3:
            RUN_PM1(3);
            break;
        default:
            RUN_PM1(4);
            break;
    }
    return found;
}


int mont_pp1(mpz_t factor, const mpz_t n, const uint32_t seed,
             const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n) || (mpz_cmp_ui(n, seed) <= 0)) {
        return 0;
    }
    if (mpz_size(n) > MONT_ECM_NARROW_LIMBS) {
        return mont_wide_pp1(factor, n, seed, plan);
    }
    mod_t m;
    mod_init(&m, n);

    int found = 0;
    switch (m.nl) {
        case 1:
            RUN_PP1(1);
            break;
        case 2:
            RUN_PP1(2);
            break;
        case 3:
            RUN_PP1(3);
            break;
        default:
            RUN_PP1(4);
            break;
    }
    return found;
}

#endif
