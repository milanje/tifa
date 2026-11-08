// Copyright (C) 2026 Jerome Milan.
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

//
// What mont_ecm.c (one to four limbs) and mont_ecm_wide.c (five to eight,
// and the long curves on fewer) share: the plan, the modulus, and the runs of a curve, of p-1 and of p+1
// on NLIMBS limbs. The includer defines MONT_MOD_LIMBS, its largest modulus
// in limbs.
//

#if !defined(_TIFA_MONT_ECM_INTERNAL_H_)
#define _TIFA_MONT_ECM_INTERNAL_H_

#include "mont_ecm.h"
#include "primes.h"
#include "tifa_alloc.h"

#define MONT_ECM_NARROW_LIMBS 4

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
    //
    // Curves long enough to check the time limit and the interrupt (on all
    // limb counts, see mont_ecm_wide.c).
    //
    bool      long_curves;
};

typedef struct {
    uint64_t n[MONT_MOD_LIMBS];
    uint64_t ninv;                       // -n^-1 mod 2^64
    uint64_t one[MONT_MOD_LIMBS];    // 2^(64.nl) mod n
    uint64_t r2[MONT_MOD_LIMBS];     // 2^(128.nl) mod n
    uint64_t r3[MONT_MOD_LIMBS];     // 2^(192.nl) mod n
    int      nl;
} mod_t;

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

//
// Stage 2 normalizes its giant steps this many at a time: one inversion
// per block.
//
#define GIANT_BLOCK 256

#define FN__(name, nl) name##_##nl
#define FN_(name, nl)  FN__(name, nl)
#define FN(name)       FN_(name, NL)

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
// When a stage catches every prime factor of n at once, it is redone with a
// gcd after each prime (stage 1) or each product (stage 2).
//
#define RUN_CURVE(NLIMBS)                                                     \
    do {                                                                      \
        FN_(pt_t, NLIMBS) p;                                                  \
        uint64_t x0m[NLIMBS];                                                 \
        uint64_t a24m[NLIMBS];                                                \
        FN_(mul, NLIMBS)(x0m, x0, m.r2, &m);                                  \
        FN_(mul, NLIMBS)(a24m, a24, m.r2, &m);                                \
        FN_(stage1, NLIMBS)(&p, x0m, a24m, plan, &m);                         \
        found = limbs_gcd(factor, p.z, NLIMBS, n);                            \
        if (found == -1) {                                                    \
            found = FN_(stage1_backtrack, NLIMBS)(factor, x0m, a24m, plan,    \
                                                  &m, n);                     \
        }                                                                     \
        if ((found == 0) && (plan->d != 0)) {                                 \
            uint64_t acc[NLIMBS];                                             \
            if (FN_(stage2, NLIMBS)(acc, &p, a24m, plan, &m, factor, n)) {    \
                found = limbs_gcd(factor, acc, NLIMBS, n);                    \
                if ((found == -1)                                             \
                    && FN_(stage2_check, NLIMBS)(acc, &p, a24m, plan, &m,      \
                                                 factor, n)) {                \
                    found = limbs_gcd(factor, acc, NLIMBS, n);                \
                }                                                             \
            } else {                                                          \
                found = (mpz_cmp(factor, n) < 0) ? 1 : -1;                    \
            }                                                                 \
        }                                                                     \
    } while (0)

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
            if (found == -1) {                                                \
                FN_(lucas_stage2_check, NLIMBS)(t, q, plan, &m, factor, n);  \
                found = limbs_gcd(factor, t, NLIMBS, n);                      \
            }                                                                 \
        }                                                                     \
        mpz_clear(a);                                                         \
        mpz_clear(ainv);                                                      \
    } while (0)

#define RUN_PP1(NLIMBS)                                                       \
    do {                                                                      \
        uint64_t p[NLIMBS], v[NLIMBS], two[NLIMBS], t[NLIMBS];               \
        TO_MONT_UI(NLIMBS, p, seed);                                          \
        FN_(pp1_stage1, NLIMBS)(v, p, plan, &m);                              \
        FN_(add, NLIMBS)(two, m.one, m.one, &m);                              \
        FN_(sub, NLIMBS)(t, v, two, &m);                                      \
        found = limbs_gcd(factor, t, NLIMBS, n);                              \
        if (found == -1) {                                                    \
            found = FN_(pp1_stage1_backtrack, NLIMBS)(factor, p, plan, &m, n);\
        }                                                                     \
        if ((found == 0) && (plan->d != 0)) {                                 \
            FN_(lucas_stage2, NLIMBS)(t, v, plan, &m);                        \
            found = limbs_gcd(factor, t, NLIMBS, n);                          \
            if (found == -1) {                                                \
                FN_(lucas_stage2_check, NLIMBS)(t, v, plan, &m, factor, n);  \
                found = limbs_gcd(factor, t, NLIMBS, n);                      \
            }                                                                 \
        }                                                                     \
    } while (0)

bool mont_ecm_wide_run(mpz_t factor, const mpz_t n, const mpz_t x0,
                       const mpz_t a24, const mont_ecm_plan_t* plan);

int mont_wide_pm1(mpz_t factor, const mpz_t n, uint32_t base,
                  const mont_ecm_plan_t* plan);

int mont_wide_pp1(mpz_t factor, const mpz_t n, uint32_t seed,
                  const mont_ecm_plan_t* plan);

#endif
