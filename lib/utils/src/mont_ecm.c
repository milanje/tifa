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

typedef unsigned __int128 u128;

struct struct_mont_ecm_plan_t {
    //
    // Stage 1: the prime powers up to b1, packed into 64-bit scalars.
    //
    uint32_t  b1;
    uint32_t  nscalars;
    uint64_t* scalars;
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
#undef NL
#define NL 2
#include "mont_ecm_impl.h"
#undef NL
#define NL 3
#include "mont_ecm_impl.h"
#undef NL
#define NL 4
#include "mont_ecm_impl.h"
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
        v->data    = realloc(v->data, v->alloced * sizeof(uint64_t));
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
    uint8_t* const composite = malloc(SEG);

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

    plan->jslot = malloc((half + 1) * sizeof(uint16_t));
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
    const uint32_t m_last = (plan->b2 + half) / d;
    plan->nm      = m_last - plan->m_first + 1;
    plan->offsets = calloc(plan->nm + 1, sizeof(uint32_t));

    stage2_builder_t sb;
    sb.plan   = plan;
    sb.pairs.data    = NULL;
    sb.pairs.length  = 0;
    sb.pairs.alloced = 0;
    sb.needed = calloc(plan->nbaby, sizeof(uint8_t));
    sb.cur_m  = plan->m_first;

    for_each_prime(plan->b1 + 1, plan->b2, add_stage2_prime, &sb);

    while (sb.cur_m <= m_last) {
        flush_giant_step(&sb);
        sb.cur_m++;
    }
    plan->pairs = malloc((sb.pairs.length + 1) * sizeof(uint16_t));
    for (uint32_t i = 0; i < sb.pairs.length; i++) {
        plan->pairs[i] = (uint16_t) sb.pairs.data[i];
    }
    free(sb.pairs.data);
    free(sb.needed);
}

mont_ecm_plan_t* alloc_mont_ecm_plan(uint32_t b1, uint32_t b2) {
    mont_ecm_plan_t* const plan = calloc(1, sizeof(mont_ecm_plan_t));

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
    plan->scalars  = scalars.data;
    plan->nscalars = scalars.length;

    if (d != 0) {
        build_stage2(plan);
    }
    return plan;
}

void free_mont_ecm_plan(mont_ecm_plan_t* const plan) {
    if (plan != NULL) {
        free(plan->scalars);
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
    mpz_mul(r, r, r);
    mpz_mod(r, r, n);
    to_limbs(m->r2, r, m->nl);
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
        uint64_t a24m[NLIMBS];                                                \
        FN_(mul, NLIMBS)(p.x, x0, m.r2, &m);                                  \
        FN_(mul, NLIMBS)(p.z, z0, m.r2, &m);                                  \
        FN_(mul, NLIMBS)(a24m, a24, m.r2, &m);                                \
        FN_(stage1, NLIMBS)(&p, a24m, plan, &m);                              \
        found = limbs_gcd(factor, p.z, NLIMBS, n);                            \
        if ((found == 0) && (plan->d != 0)) {                                 \
            uint64_t acc[NLIMBS];                                             \
            FN_(stage2, NLIMBS)(acc, &p, a24m, plan, &m);                     \
            found = limbs_gcd(factor, acc, NLIMBS, n);                        \
        }                                                                     \
    } while (0)

bool mont_ecm_curve(mpz_t factor, const mpz_t n, const uint32_t sigma,
                    const mont_ecm_plan_t* const plan) {
    if (!mont_ecm_handles(n)) {
        return false;
    }
    //
    // Suyama's parametrization: u = sigma^2 - 5, v = 4.sigma, starting
    // point (u^3 : v^3) on the curve with (A + 2)/4 =
    // (v - u)^3 (3u + v) / (16 u^3 v).
    //
    mpz_t u, v, x0z, z0z, a24z, t;
    mpz_inits(u, v, x0z, z0z, a24z, t, NULL);

    mpz_set_ui(u, sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5);
    mpz_mod(u, u, n);
    mpz_set_ui(v, sigma);
    mpz_mul_ui(v, v, 4);
    mpz_mod(v, v, n);

    mpz_powm_ui(x0z, u, 3, n);
    mpz_powm_ui(z0z, v, 3, n);

    mpz_mul(t, x0z, v);
    mpz_mul_ui(t, t, 16);
    mpz_mod(t, t, n);

    int found = 0;
    if (mpz_invert(t, t, n) == 0) {
        mpz_mul(t, x0z, v);
        mpz_mul_ui(t, t, 16);
        mpz_gcd(factor, t, n);
        found = ((mpz_cmp_ui(factor, 1) > 0) && (mpz_cmp(factor, n) < 0));
        mpz_clears(u, v, x0z, z0z, a24z, t, NULL);
        return found;
    }
    mpz_sub(a24z, v, u);
    mpz_powm_ui(a24z, a24z, 3, n);
    mpz_mul(a24z, a24z, t);
    mpz_mul_ui(t, u, 3);
    mpz_add(t, t, v);
    mpz_mul(a24z, a24z, t);
    mpz_mod(a24z, a24z, n);

    mod_t m;
    mod_init(&m, n);

    uint64_t x0[MONT_ECM_MAX_LIMBS];
    uint64_t z0[MONT_ECM_MAX_LIMBS];
    uint64_t a24[MONT_ECM_MAX_LIMBS];
    to_limbs(x0, x0z, m.nl);
    to_limbs(z0, z0z, m.nl);
    to_limbs(a24, a24z, m.nl);
    mpz_clears(u, v, x0z, z0z, a24z, t, NULL);

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

#endif
