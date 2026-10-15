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

//
// ECM kernels for moduli of exactly NL 64-bit limbs. Included once per limb
// count by mont_ecm.c, with NL and FN(name) (name suffixed by NL) defined:
// the loops have constant bounds, so the compiler unrolls them.
//

//------------------------------------------------------------------------------
//                     Montgomery arithmetic (CIOS), NL limbs
//------------------------------------------------------------------------------

static inline void FN(mul)(uint64_t* const r, const uint64_t* const a,
                           const uint64_t* const b, const mod_t* const m) {
    uint64_t t[NL + 2] = { 0 };

    for (int i = 0; i < NL; i++) {
        u128 c = 0;
        for (int j = 0; j < NL; j++) {
            c += (u128) a[j] * b[i] + t[j];
            t[j] = (uint64_t) c;
            c >>= 64;
        }
        c += t[NL];
        t[NL]     = (uint64_t) c;
        t[NL + 1] = (uint64_t) (c >> 64);

        const uint64_t q = t[0] * m->ninv;
        c = ((u128) q * m->n[0] + t[0]) >> 64;
        for (int j = 1; j < NL; j++) {
            c += (u128) q * m->n[j] + t[j];
            t[j - 1] = (uint64_t) c;
            c >>= 64;
        }
        c += t[NL];
        t[NL - 1] = (uint64_t) c;
        t[NL]     = t[NL + 1] + (uint64_t) (c >> 64);
    }
    //
    // t < 2n: subtract n once if needed.
    //
    uint64_t s[NL];
    u128 borrow = 0;
    for (int j = 0; j < NL; j++) {
        const u128 d = (u128) t[j] - m->n[j] - borrow;
        s[j]   = (uint64_t) d;
        borrow = (d >> 64) & 1;
    }
    const int keep_t = (t[NL] == 0) && (borrow != 0);
    for (int j = 0; j < NL; j++) {
        r[j] = keep_t ? t[j] : s[j];
    }
}

static inline void FN(add)(uint64_t* const r, const uint64_t* const a,
                           const uint64_t* const b, const mod_t* const m) {
    uint64_t t[NL];
    u128 c = 0;
    for (int j = 0; j < NL; j++) {
        c += (u128) a[j] + b[j];
        t[j] = (uint64_t) c;
        c >>= 64;
    }
    uint64_t s[NL];
    u128 borrow = 0;
    for (int j = 0; j < NL; j++) {
        const u128 d = (u128) t[j] - m->n[j] - borrow;
        s[j]   = (uint64_t) d;
        borrow = (d >> 64) & 1;
    }
    const int keep_t = (c == 0) && (borrow != 0);
    for (int j = 0; j < NL; j++) {
        r[j] = keep_t ? t[j] : s[j];
    }
}

static inline void FN(sub)(uint64_t* const r, const uint64_t* const a,
                           const uint64_t* const b, const mod_t* const m) {
    uint64_t t[NL];
    u128 borrow = 0;
    for (int j = 0; j < NL; j++) {
        const u128 d = (u128) a[j] - b[j] - borrow;
        t[j]   = (uint64_t) d;
        borrow = (d >> 64) & 1;
    }
    if (borrow) {
        u128 c = 0;
        for (int j = 0; j < NL; j++) {
            c += (u128) t[j] + m->n[j];
            t[j] = (uint64_t) c;
            c >>= 64;
        }
    }
    for (int j = 0; j < NL; j++) {
        r[j] = t[j];
    }
}

//------------------------------------------------------------------------------
//                  Montgomery curves, x-only (X : Z) coordinates
//------------------------------------------------------------------------------

typedef struct {
    uint64_t x[NL];
    uint64_t z[NL];
} FN(pt_t);

static inline void FN(xdbl)(FN(pt_t)* const r, const FN(pt_t)* const p,
                            const uint64_t* const a24,
                            const mod_t* const m) {
    uint64_t s[NL], d[NL], s2[NL], d2[NL], t[NL], u[NL];
    FN(add)(s, p->x, p->z, m);
    FN(sub)(d, p->x, p->z, m);
    FN(mul)(s2, s, s, m);
    FN(mul)(d2, d, d, m);
    FN(sub)(t, s2, d2, m);
    FN(mul)(r->x, s2, d2, m);
    FN(mul)(u, a24, t, m);
    FN(add)(u, u, d2, m);
    FN(mul)(r->z, t, u, m);
}

//
// r = p + q, given diff = p - q. r may alias p or q.
//
static inline void FN(xadd)(FN(pt_t)* const r, const FN(pt_t)* const p,
                            const FN(pt_t)* const q,
                            const FN(pt_t)* const diff,
                            const mod_t* const m) {
    uint64_t a[NL], b[NL], u[NL], v[NL], s[NL], d[NL];
    FN(sub)(a, p->x, p->z, m);
    FN(add)(b, q->x, q->z, m);
    FN(mul)(u, a, b, m);
    FN(add)(a, p->x, p->z, m);
    FN(sub)(b, q->x, q->z, m);
    FN(mul)(v, a, b, m);
    FN(add)(s, u, v, m);
    FN(sub)(d, u, v, m);
    FN(mul)(s, s, s, m);
    FN(mul)(d, d, d, m);
    FN(mul)(r->x, diff->z, s, m);
    FN(mul)(r->z, diff->x, d, m);
}

//
// r = k.p (k >= 1). r may alias p.
//
static void FN(ladder)(FN(pt_t)* const r, const FN(pt_t)* const p,
                       const uint64_t k, const uint64_t* const a24,
                       const mod_t* const m) {
    const FN(pt_t) base = *p;
    FN(pt_t) r0 = base;
    FN(pt_t) r1;
    FN(xdbl)(&r1, &base, a24, m);

    for (int b = 62 - __builtin_clzll(k); b >= 0; b--) {
        if ((k >> b) & 1) {
            FN(xadd)(&r0, &r1, &r0, &base, m);
            FN(xdbl)(&r1, &r1, a24, m);
        } else {
            FN(xadd)(&r1, &r1, &r0, &base, m);
            FN(xdbl)(&r0, &r0, a24, m);
        }
    }
    *r = r0;
}

//------------------------------------------------------------------------------
//                              Stages 1 and 2
//------------------------------------------------------------------------------

//
// p <- k.p for the product k of the prime powers up to the stage 1 bound.
//
static void FN(stage1)(FN(pt_t)* const p, const uint64_t* const a24,
                       const mont_ecm_plan_t* const plan,
                       const mod_t* const m) {
    for (uint32_t i = 0; i < plan->nscalars; i++) {
        FN(ladder)(p, p, plan->scalars[i], a24, m);
    }
}

//
// Product over the primes q = k.d +/- j of the plan of
// x(k.d.P).z(j.P) - x(j.P).z(k.d.P), which vanishes mod a prime factor p of
// n when q.P is the point at infinity mod p. One product covers both
// k.d - j and k.d + j.
//
static void FN(stage2)(uint64_t* const acc, const FN(pt_t)* const p,
                       const uint64_t* const a24,
                       const mont_ecm_plan_t* const plan,
                       const mod_t* const m) {
    //
    // Baby steps: j.P for the odd j < d/2 coprime to d, by the chain
    // (j + 2).P = j.P + 2.P (difference (j - 2).P).
    //
    const uint32_t half = plan->d / 2;
    FN(pt_t)* const baby = malloc(plan->nbaby * sizeof(FN(pt_t)));

    FN(pt_t) p2, prev, cur, next;
    FN(xdbl)(&p2, p, a24, m);
    if (plan->jslot[1] != UINT16_MAX) {
        baby[plan->jslot[1]] = *p;
    }
    prev = *p;
    FN(xadd)(&cur, &p2, p, p, m);
    for (uint32_t j = 3; ; j += 2) {
        if (plan->jslot[j] != UINT16_MAX) {
            baby[plan->jslot[j]] = cur;
        }
        if (j == half) {
            break;
        }
        FN(xadd)(&next, &cur, &p2, &prev, m);
        prev = cur;
        cur  = next;
    }
    //
    // cur = (d/2).P, so d.P = 2.cur. Giant steps G_k = k.d.P from
    // k = plan->m_first on, by G_{k+1} = G_k + d.P (difference G_{k-1}).
    //
    FN(pt_t) dp, gprev, gcur;
    FN(xdbl)(&dp, &cur, a24, m);

    if (plan->m_first == 1) {
        gprev = dp;
        gcur  = dp;
    } else {
        FN(ladder)(&gprev, &dp, plan->m_first - 1, a24, m);
        FN(ladder)(&gcur, &dp, plan->m_first, a24, m);
    }
    uint64_t t[NL], u[NL];
    for (int j = 0; j < NL; j++) {
        acc[j] = m->one[j];
    }
    for (uint32_t k = 0; k < plan->nm; k++) {
        if (k > 0) {
            if (plan->m_first + k == 2) {
                FN(xdbl)(&gcur, &dp, a24, m);
            } else {
                FN(xadd)(&next, &gcur, &dp, &gprev, m);
                gprev = gcur;
                gcur  = next;
            }
        }
        for (uint32_t i = plan->offsets[k]; i < plan->offsets[k + 1]; i++) {
            const FN(pt_t)* const b = &baby[plan->pairs[i]];
            FN(mul)(t, gcur.x, b->z, m);
            FN(mul)(u, b->x, gcur.z, m);
            FN(sub)(t, t, u, m);
            FN(mul)(acc, acc, t, m);
        }
    }
    free(baby);
}
