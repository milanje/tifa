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
 *        to four 64-bit limbs.
 */

#include <stdlib.h>
#include <string.h>

#include "mont_ecm.h"

#if TIFA_HAVE_MONT_ECM

#include "first_primes.h"
#include "z2z8_curves.h"
#include "tifa_alloc.h"

typedef unsigned __int128 u128;

struct struct_mont_ecm_plan_t {
    //
    // Stage 1: the product E of the prime powers up to b1 (least
    // significant limb first).
    //
    uint32_t  b1;
    uint64_t* e;
    size_t    ne;
    size_t    ebits;
    //
    // Stage 2 (d = 0 if none): the primes in (b1, b2] as k.d +/- j for the
    // odd j < d/2 coprime to d. jslot[j] is the index of j among them
    // (UINT16_MAX if j is not coprime to d). For the k-th giant step,
    // pairs[offsets[k] .. offsets[k + 1] - 1] lists the indices of the j
    // needed with (m_first + k).d.
    //
    uint32_t  b2;
    uint32_t  d;
    uint32_t  nbaby;
    uint16_t* jslot;
    uint32_t  m_first;
    uint32_t  nm;
    uint32_t* offsets;
    uint16_t* pairs;
};

typedef struct {
    uint64_t n[MONT_ECM_MAX_LIMBS];
    uint64_t ninv;                       // -n^-1 mod 2^64
    uint64_t one[MONT_ECM_MAX_LIMBS];    // 2^(64.nl) mod n
    uint64_t r2[MONT_ECM_MAX_LIMBS];     // 2^(128.nl) mod n
    uint64_t r3[MONT_ECM_MAX_LIMBS];     // 2^(192.nl) mod n
    int      nl;
} mod_t;

//------------------------------------------------------------------------------
//                          Kernels for 1 to 4 limbs
//------------------------------------------------------------------------------

#define FN__(name, nl) name##_##nl
#define FN_(name, nl)  FN__(name, nl)
#define FN(name)       FN_(name, NL)

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
    uint64_t* data;
    uint32_t  length;
    uint32_t  alloced;
} u64_vec_t;

static void push_u64(u64_vec_t* const v, const uint64_t x) {
    if (v->length == v->alloced) {
        v->alloced = (v->alloced == 0) ? 64 : 2 * v->alloced;
        v->data    = tifa_realloc(v->data, v->alloced * sizeof(uint64_t));
    }
    v->data[v->length++] = x;
}

//
// Calls visit(q, data) for every prime q in [lo, hi], in increasing order:
// from the prime table, then from a segmented sieve beyond it.
//
static void for_each_prime(const uint32_t lo, const uint32_t hi,
                           void (*visit)(uint32_t, void*), void* const data) {
    uint32_t ip = 0;
    while ((ip < NFIRST_PRIMES) && (first_primes[ip] < lo)) {
        ip++;
    }
    while ((ip < NFIRST_PRIMES) && (first_primes[ip] <= hi)) {
        visit(first_primes[ip], data);
        ip++;
    }
    if (hi <= LARGEST_PRIME) {
        return;
    }
    const uint64_t first = (lo > (uint64_t) LARGEST_PRIME + 1)
                           ? lo : (uint64_t) LARGEST_PRIME + 1;
    const uint64_t SEG   = 1 << 16;
    uint8_t* const composite = tifa_malloc(SEG);

    for (uint64_t s = first; s <= hi; s += SEG) {
        const uint64_t e = (s + SEG - 1 < hi) ? (s + SEG - 1) : hi;
        memset(composite, 0, SEG);
        for (uint32_t i = 0; i < NFIRST_PRIMES; i++) {
            const uint64_t p = first_primes[i];
            if (p * p > e) {
                break;
            }
            uint64_t x = ((s + p - 1) / p) * p;
            if (x < p * p) {
                x = p * p;
            }
            for (; x <= e; x += p) {
                composite[x - s] = 1;
            }
        }
        for (uint64_t x = s; x <= e; x++) {
            if (!composite[x - s]) {
                visit((uint32_t) x, data);
            }
        }
    }
    free(composite);
}

typedef struct {
    mont_ecm_plan_t* plan;
    u64_vec_t        pairs;
    uint8_t*         needed;
    uint32_t         cur_m;
} stage2_builder_t;

static void flush_giant_step(stage2_builder_t* const sb) {
    mont_ecm_plan_t* const plan = sb->plan;
    for (uint32_t i = 0; i < plan->nbaby; i++) {
        if (sb->needed[i]) {
            push_u64(&sb->pairs, i);
            sb->needed[i] = 0;
        }
    }
    plan->offsets[sb->cur_m - plan->m_first + 1] = sb->pairs.length;
}

static void add_stage2_prime(const uint32_t q, void* const data) {
    stage2_builder_t* const sb   = (stage2_builder_t*) data;
    mont_ecm_plan_t*  const plan = sb->plan;

    const uint32_t m = (q + plan->d / 2) / plan->d;
    while (sb->cur_m < m) {
        flush_giant_step(sb);
        sb->cur_m++;
    }
    const uint32_t c = m * plan->d;
    const uint32_t j = (q > c) ? (q - c) : (c - q);
    sb->needed[plan->jslot[j]] = 1;
}

static void build_stage2(mont_ecm_plan_t* const plan) {
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
    sb.plan   = plan;
    sb.pairs.data    = NULL;
    sb.pairs.length  = 0;
    sb.pairs.alloced = 0;
    sb.needed = tifa_calloc(plan->nbaby, sizeof(uint8_t));
    sb.cur_m  = plan->m_first;

    for_each_prime(plan->b1 + 1, plan->b2, add_stage2_prime, &sb);

    while (sb.cur_m <= m_last) {
        flush_giant_step(&sb);
        sb.cur_m++;
    }
    plan->pairs = tifa_malloc((sb.pairs.length + 1) * sizeof(uint16_t));
    for (uint32_t i = 0; i < sb.pairs.length; i++) {
        plan->pairs[i] = (uint16_t) sb.pairs.data[i];
    }
    free(sb.pairs.data);
    free(sb.needed);
}

//
// e <- scalars[lo] * ... * scalars[hi - 1], as a balanced product tree.
//
static void product_of_scalars(mpz_t e, const uint64_t* const scalars,
                               const uint32_t lo, const uint32_t hi) {
    if (hi - lo == 1) {
        mpz_import(e, 1, -1, sizeof(uint64_t), 0, 0, &scalars[lo]);
        return;
    }
    mpz_t right;
    mpz_init(right);
    product_of_scalars(e, scalars, lo, lo + (hi - lo) / 2);
    product_of_scalars(right, scalars, lo + (hi - lo) / 2, hi);
    mpz_mul(e, e, right);
    mpz_clear(right);
}

mont_ecm_plan_t* alloc_mont_ecm_plan(uint32_t b1, uint32_t b2) {
    mont_ecm_plan_t* const plan = tifa_calloc(1, sizeof(mont_ecm_plan_t));

    if (b1 > LARGEST_PRIME) {
        b1 = LARGEST_PRIME;
    }
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

    u64_vec_t scalars = { NULL, 0, 0 };
    uint64_t  k = 1;
    for (uint32_t ip = 0; (ip < NFIRST_PRIMES) && (first_primes[ip] <= b1);
         ip++) {
        const uint64_t p = first_primes[ip];
        uint64_t q = p;
        while (q <= b1 / p) {
            q *= p;
        }
        if (k > UINT64_MAX / q) {
            push_u64(&scalars, k);
            k = 1;
        }
        k *= q;
    }
    push_u64(&scalars, k);

    mpz_t e;
    mpz_init(e);
    product_of_scalars(e, scalars.data, 0, scalars.length);
    free(scalars.data);
    plan->ebits = mpz_sizeinbase(e, 2);
    plan->e     = tifa_malloc((plan->ebits / 64 + 1) * sizeof(uint64_t));
    mpz_export(plan->e, &plan->ne, -1, sizeof(uint64_t), 0, 0, e);
    mpz_clear(e);

    if (d != 0) {
        build_stage2(plan);
    }
    return plan;
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

static void to_limbs(uint64_t* const r, const mpz_t z, const int nl) {
    for (int i = 0; i < nl; i++) {
        r[i] = (uint64_t) mpz_getlimbn(z, i);
    }
}

static void mod_init(mod_t* const m, const mpz_t n) {
    m->nl = (int) mpz_size(n);
    to_limbs(m->n, n, m->nl);

    uint64_t x = (3 * m->n[0]) ^ 2;
    for (int i = 0; i < 4; i++) {
        x *= 2 - m->n[0] * x;
    }
    m->ninv = 0 - x;

    mpz_t r;
    mpz_init(r);
    mpz_setbit(r, 64 * m->nl);
    mpz_mod(r, r, n);
    to_limbs(m->one, r, m->nl);
    mpz_t one;
    mpz_init_set(one, r);
    mpz_mul(r, r, r);
    mpz_mod(r, r, n);
    to_limbs(m->r2, r, m->nl);
    mpz_mul(r, r, one);
    mpz_mod(r, r, n);
    to_limbs(m->r3, r, m->nl);
    mpz_clear(one);
    mpz_clear(r);
}

//
// Sets factor to gcd(a, n) where a is given by its nl limbs. Returns 1 if
// that is a proper factor, -1 if it is n, 0 if it is 1.
//
static int limbs_gcd(mpz_t factor, const uint64_t* const a, const int nl,
                     const mpz_t n) {
    mpz_import(factor, (size_t) nl, -1, sizeof(uint64_t), 0, 0, a);
    mpz_gcd(factor, factor, n);
    if (mpz_cmp_ui(factor, 1) == 0) {
        return 0;
    }
    return (mpz_cmp(factor, n) == 0) ? -1 : 1;
}

#define RUN_CURVE(NLIMBS)                                                     \
    do {                                                                      \
        FN_(pt_t, NLIMBS) p;                                                  \
        uint64_t x0m[NLIMBS];                                                 \
        uint64_t a24m[NLIMBS];                                                \
        FN_(mul, NLIMBS)(x0m, x0, m.r2, &m);                                  \
        FN_(mul, NLIMBS)(a24m, a24, m.r2, &m);                                \
        FN_(stage1, NLIMBS)(&p, x0m, a24m, plan, &m);                         \
        found = limbs_gcd(factor, p.z, NLIMBS, n);                            \
        if ((found == 0) && (plan->d != 0)) {                                 \
            uint64_t acc[NLIMBS];                                             \
            if (FN_(stage2, NLIMBS)(acc, &p, a24m, plan, &m, factor, n)) {    \
                found = limbs_gcd(factor, acc, NLIMBS, n);                    \
            } else {                                                          \
                found = (mpz_cmp(factor, n) < 0) ? 1 : -1;                    \
            }                                                                 \
        }                                                                     \
    } while (0)

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
// As suyama_curve, for row i of z2z8_curves: a24 = 1/(1 - d) and
// x0 = (1 + y0)/(1 - y0), from the inverse of (dd - dn)(yd - yn).
//
static bool z2z8_curve(mpz_t x0, mpz_t a24, mpz_t factor, const mpz_t n,
                       const uint32_t i) {
    mpz_t dn, dd, yn, yd, t;
    mpz_inits(dn, dd, yn, yd, t, NULL);
    mpz_import(dn, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][0]);
    mpz_import(dd, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][1]);
    mpz_import(yn, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][2]);
    mpz_import(yd, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][3]);

    mpz_sub(a24, dd, dn);
    mpz_sub(x0, yd, yn);
    mpz_mul(t, a24, x0);
    mpz_mod(t, t, n);

    const bool ok = (mpz_invert(factor, t, n) != 0);
    if (!ok) {
        mpz_gcd(factor, t, n);
    } else {
        mpz_mul(t, factor, x0);
        mpz_mul(t, t, dd);
        mpz_mod(t, t, n);
        mpz_mul(dd, factor, a24);
        mpz_add(yd, yd, yn);
        mpz_mul(dd, dd, yd);
        mpz_mod(x0, dd, n);
        mpz_set(a24, t);
    }
    mpz_clears(dn, dd, yn, yd, t, NULL);
    return ok;
}

//
// One curve, from its starting point (x0z : 1) and a24 (residues mod n).
//
static bool run_curve(mpz_t factor, const mpz_t n, const mpz_t x0z,
                      const mpz_t a24z, const mont_ecm_plan_t* const plan) {
    mod_t m;
    mod_init(&m, n);

    uint64_t x0[MONT_ECM_MAX_LIMBS];
    uint64_t a24[MONT_ECM_MAX_LIMBS];
    to_limbs(x0, x0z, m.nl);
    to_limbs(a24, a24z, m.nl);

    int found = 0;
    switch (m.nl) {
        case 1:
            RUN_CURVE(1);
            break;
        case 2:
            RUN_CURVE(2);
            break;
        case 3:
            RUN_CURVE(3);
            break;
        default:
            RUN_CURVE(4);
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
        found = run_curve(factor, n, x0, a24, plan);
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

bool mont_ecm_nth_curve(mpz_t factor, const mpz_t n, const uint32_t k,
                        const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n)) {
        return false;
    }
    if (k < NZ2Z8_CURVES) {
        return setup_and_run(factor, n, k, plan, z2z8_curve);
    }
    return setup_and_run(factor, n, 6 + k - NZ2Z8_CURVES, plan, suyama_curve);
}

//------------------------------------------------------------------------------
//                              p-1 and p+1
//------------------------------------------------------------------------------

//
// Montgomery form of the small integer v (v < n).
//
#define TO_MONT_UI(NLIMBS, r, v)                                              \
    do {                                                                      \
        uint64_t v_[NLIMBS] = { 0 };                                          \
        v_[0] = (v);                                                          \
        FN_(mul, NLIMBS)(r, v_, m.r2, &m);                                    \
    } while (0)

//
// Montgomery form of the integer z in [0, n).
//
#define TO_MONT_MPZ(NLIMBS, r, z)                                             \
    do {                                                                      \
        uint64_t z_[NLIMBS];                                                  \
        to_limbs(z_, z, NLIMBS);                                              \
        FN_(mul, NLIMBS)(r, z_, m.r2, &m);                                    \
    } while (0)

#define FROM_MONT_MPZ(NLIMBS, z, r)                                           \
    do {                                                                      \
        uint64_t one_[NLIMBS] = { 0 };                                        \
        uint64_t r_[NLIMBS];                                                  \
        one_[0] = 1;                                                          \
        FN_(mul, NLIMBS)(r_, r, one_, &m);                                    \
        mpz_import(z, NLIMBS, -1, sizeof(uint64_t), 0, 0, r_);                \
    } while (0)

#define RUN_PM1(NLIMBS)                                                       \
    do {                                                                      \
        uint64_t b[NLIMBS], x[NLIMBS], t[NLIMBS];                             \
        TO_MONT_UI(NLIMBS, b, base);                                          \
        FN_(pm1_stage1, NLIMBS)(x, b, base == 2, plan, &m);                   \
        FN_(sub, NLIMBS)(t, x, m.one, &m);                                    \
        found = limbs_gcd(factor, t, NLIMBS, n);                              \
        if ((found != 0) || (plan->d == 0)) {                                 \
            break;                                                            \
        }                                                                     \
        mpz_t a, ainv;                                                        \
        mpz_init(a);                                                          \
        mpz_init(ainv);                                                       \
        FROM_MONT_MPZ(NLIMBS, a, x);                                          \
        if (mpz_invert(ainv, a, n) == 0) {                                    \
            mpz_gcd(factor, a, n);                                            \
            found = (mpz_cmp(factor, n) < 0) ? 1 : -1;                        \
        } else {                                                              \
            uint64_t q[NLIMBS];                                               \
            mpz_add(a, a, ainv);                                              \
            mpz_mod(a, a, n);                                                 \
            TO_MONT_MPZ(NLIMBS, q, a);                                        \
            FN_(lucas_stage2, NLIMBS)(t, q, plan, &m);                        \
            found = limbs_gcd(factor, t, NLIMBS, n);                          \
        }                                                                     \
        mpz_clear(a);                                                         \
        mpz_clear(ainv);                                                      \
    } while (0)

int mont_pm1(mpz_t factor, const mpz_t n, const uint32_t base,
             const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n) || (mpz_cmp_ui(n, base) <= 0)) {
        return 0;
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

#define RUN_PP1(NLIMBS)                                                       \
    do {                                                                      \
        uint64_t p[NLIMBS], v[NLIMBS], two[NLIMBS], t[NLIMBS];               \
        TO_MONT_UI(NLIMBS, p, seed);                                          \
        FN_(pp1_stage1, NLIMBS)(v, p, plan, &m);                              \
        FN_(add, NLIMBS)(two, m.one, m.one, &m);                              \
        FN_(sub, NLIMBS)(t, v, two, &m);                                      \
        found = limbs_gcd(factor, t, NLIMBS, n);                              \
        if ((found == 0) && (plan->d != 0)) {                                 \
            FN_(lucas_stage2, NLIMBS)(t, v, plan, &m);                        \
            found = limbs_gcd(factor, t, NLIMBS, n);                          \
        }                                                                     \
    } while (0)

int mont_pp1(mpz_t factor, const mpz_t n, const uint32_t seed,
             const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n) || (mpz_cmp_ui(n, seed) <= 0)) {
        return 0;
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
