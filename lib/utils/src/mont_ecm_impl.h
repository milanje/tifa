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

//
// Always inlined, for the ladder step: there, on 3 and 4 limbs, the compiler
// would call mul.
//
static inline __attribute__((always_inline)) void
FN(mul_inline)(uint64_t* const r, const uint64_t* const a,
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

static inline void FN(mul)(uint64_t* const r, const uint64_t* const a,
                           const uint64_t* const b, const mod_t* const m) {
    FN(mul_inline)(r, a, b, m);
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
// (r, s) <- (2.r, r + s), given s - r = (xd : 1). The products of the
// doubling and of the addition are interleaved: as two calls, each a chain
// of dependent products, they did not overlap.
//
static inline __attribute__((always_inline)) void
FN(ladder_step)(FN(pt_t)* const r, FN(pt_t)* const s,
                const uint64_t* const xd, const uint64_t* const a24,
                const mod_t* const m) {
    uint64_t a[NL], b[NL], c[NL], d[NL], da[NL], cb[NL], aa[NL], bb[NL];
    uint64_t e[NL];
    FN(add)(a, r->x, r->z, m);
    FN(sub)(b, r->x, r->z, m);
    FN(add)(c, s->x, s->z, m);
    FN(sub)(d, s->x, s->z, m);
    FN(mul_inline)(da, d, a, m);
    FN(mul_inline)(cb, c, b, m);
    FN(mul_inline)(aa, a, a, m);
    FN(mul_inline)(bb, b, b, m);
    FN(add)(c, da, cb, m);
    FN(sub)(d, da, cb, m);
    FN(sub)(e, aa, bb, m);
    FN(mul_inline)(s->x, c, c, m);
    FN(mul_inline)(d, d, d, m);
    FN(mul_inline)(r->x, aa, bb, m);
    FN(mul_inline)(a, a24, e, m);
    FN(add)(a, a, bb, m);
    FN(mul_inline)(s->z, xd, d, m);
    FN(mul_inline)(r->z, e, a, m);
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
// p <- E.(x0 : 1), E the product of the prime powers up to the stage 1
// bound: a ladder over the bits of E, whose difference stays (x0 : 1).
//
static void FN(stage1)(FN(pt_t)* const p, const uint64_t* const x0,
                       const uint64_t* const a24,
                       const mont_ecm_plan_t* const plan,
                       const mod_t* const m) {
    FN(pt_t) r0, r1;
    for (int j = 0; j < NL; j++) {
        r0.x[j] = x0[j];
        r0.z[j] = m->one[j];
    }
    FN(xdbl)(&r1, &r0, a24, m);
    for (long b = (long) plan->ebits - 2; b >= 0; b--) {
        if ((plan->e[b / 64] >> (b % 64)) & 1) {
            FN(ladder_step)(&r1, &r0, x0, a24, m);
        } else {
            FN(ladder_step)(&r0, &r1, x0, a24, m);
        }
    }
    *p = r0;
}

//
// Stage 1 again from (x0 : 1), one prime at a time, with a gcd after each:
// for when the whole of it caught every prime factor of n at once. Returns
// as limbs_gcd (-1 if a single prime still catches them all).
//
static int FN(stage1_backtrack)(mpz_t factor, const uint64_t* const x0,
                                const uint64_t* const a24,
                                const mont_ecm_plan_t* const plan,
                                const mod_t* const m, const mpz_t n) {
    FN(pt_t) p;
    for (int j = 0; j < NL; j++) {
        p.x[j] = x0[j];
        p.z[j] = m->one[j];
    }
    for (uint32_t i = 0; (i < NFIRST_PRIMES) && (first_primes[i] <= plan->b1);
         i++) {
        const uint64_t q = first_primes[i];
        for (uint64_t qe = q; qe <= plan->b1; qe *= q) {
            FN(ladder)(&p, &p, q, a24, m);
            const int found = limbs_gcd(factor, p.z, NL, n);
            if (found != 0) {
                return found;
            }
        }
    }
    return -1;
}

//
// xs[i] = x / z of pts[i] for i < count, by a batch inversion. Returns
// false, with g the gcd of the product of the z with n, if they are not
// invertible.
//
static bool FN(normalize)(uint64_t (*const xs)[NL], const FN(pt_t)* const pts,
                          const uint32_t count, const mod_t* const m,
                          mpz_t g, const mpz_t n) {
    for (int j = 0; j < NL; j++) {
        xs[0][j] = pts[0].z[j];
    }
    for (uint32_t i = 1; i < count; i++) {
        FN(mul)(xs[i], xs[i - 1], pts[i].z, m);
    }
    mpz_import(g, NL, -1, sizeof(uint64_t), 0, 0, xs[count - 1]);
    if (!mpz_invert(g, g, n)) {
        mpz_import(g, NL, -1, sizeof(uint64_t), 0, 0, xs[count - 1]);
        mpz_gcd(g, g, n);
        return false;
    }
    uint64_t inv[NL];
    for (int j = 0; j < NL; j++) {
        inv[j] = (uint64_t) mpz_getlimbn(g, j);
    }
    FN(mul)(inv, inv, m->r3, m);
    for (uint32_t i = count - 1; i > 0; i--) {
        uint64_t iz[NL];
        FN(mul)(iz, inv, xs[i - 1], m);
        FN(mul)(inv, inv, pts[i].z, m);
        FN(mul)(xs[i], pts[i].x, iz, m);
    }
    FN(mul)(xs[0], pts[0].x, inv, m);
    return true;
}

//
// Product over the primes q = k.d +/- j of the plan of x(k.d.P) - x(j.P),
// the baby steps j.P and the giant steps k.d.P normalized to z = 1, which
// vanishes mod a prime factor p of n when q.P is the point at infinity
// mod p. One product covers both k.d - j and k.d + j. Returns false, with
// g the gcd of a product of the z with n, if they are not invertible. With
// check, it stops at the first product sharing a factor with n: for when
// the whole of it caught every prime factor at once.
//
static inline __attribute__((always_inline)) bool
FN(stage2_body)(uint64_t* const acc, const FN(pt_t)* const p,
                const uint64_t* const a24, const mont_ecm_plan_t* const plan,
                const mod_t* const m, mpz_t g, const mpz_t n,
                const bool check) {
    //
    // Baby steps: j.P for the odd j < d/2 coprime to d, by the chain
    // (j + 2).P = j.P + 2.P (difference (j - 2).P), in pts[0, nb). The giant
    // steps follow them, GIANT_BLOCK at a time, the first block normalized
    // with the baby steps.
    //
    const uint32_t half = plan->d / 2;
    const uint32_t nb   = plan->nbaby;
    const uint32_t blk  = (plan->nm < GIANT_BLOCK) ? plan->nm : GIANT_BLOCK;
    FN(pt_t)* const pts = tifa_malloc((nb + blk) * sizeof(FN(pt_t)));
    uint64_t (*const xs)[NL] = tifa_malloc((nb + blk) * sizeof(*xs));

    FN(pt_t) p2, prev, cur, next;
    FN(xdbl)(&p2, p, a24, m);
    if (plan->jslot[1] != UINT16_MAX) {
        pts[plan->jslot[1]] = *p;
    }
    prev = *p;
    FN(xadd)(&cur, &p2, p, p, m);
    for (uint32_t j = 3; ; j += 2) {
        if (plan->jslot[j] != UINT16_MAX) {
            pts[plan->jslot[j]] = cur;
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
    //
    // Four products in turn (one with check): a single one would be a chain
    // of dependent multiplications.
    //
    const uint32_t mask = check ? 0 : 3;
    uint64_t t[NL], a[4][NL];
    for (int c = 0; c < 4; c++) {
        for (int j = 0; j < NL; j++) {
            a[c][j] = m->one[j];
        }
    }
    uint32_t c = 0;
    bool ok = true;
    for (uint32_t k0 = 0; ok && (k0 < plan->nm); k0 += blk) {
        const uint32_t cnt = (plan->nm - k0 < blk) ? plan->nm - k0 : blk;
        for (uint32_t i = 0; i < cnt; i++) {
            const uint32_t k = k0 + i;
            if (k > 0) {
                if (plan->m_first + k == 2) {
                    FN(xdbl)(&gcur, &dp, a24, m);
                } else {
                    FN(xadd)(&next, &gcur, &dp, &gprev, m);
                    gprev = gcur;
                    gcur  = next;
                }
            }
            pts[nb + i] = gcur;
        }
        ok = (k0 == 0) ? FN(normalize)(xs, pts, nb + cnt, m, g, n)
                       : FN(normalize)(xs + nb, pts + nb, cnt, m, g, n);
        for (uint32_t i = 0; ok && (i < cnt); i++) {
            const uint32_t k = k0 + i;
            for (uint32_t q = plan->offsets[k]; q < plan->offsets[k + 1];
                 q++) {
                FN(sub)(t, xs[nb + i], xs[plan->pairs[q]], m);
                FN(mul)(a[c], a[c], t, m);
                if (check && (limbs_gcd(g, a[c], NL, n) != 0)) {
                    for (int j = 0; j < NL; j++) {
                        acc[j] = a[c][j];
                    }
                    free(xs);
                    free(pts);
                    return true;
                }
                c = (c + 1) & mask;
            }
        }
    }
    FN(mul)(a[0], a[0], a[1], m);
    FN(mul)(a[2], a[2], a[3], m);
    FN(mul)(acc, a[0], a[2], m);
    free(xs);
    free(pts);
    return ok;
}

static bool FN(stage2)(uint64_t* const acc, const FN(pt_t)* const p,
                       const uint64_t* const a24,
                       const mont_ecm_plan_t* const plan,
                       const mod_t* const m, mpz_t g, const mpz_t n) {
    return FN(stage2_body)(acc, p, a24, plan, m, g, n, false);
}

static bool FN(stage2_check)(uint64_t* const acc, const FN(pt_t)* const p,
                             const uint64_t* const a24,
                             const mont_ecm_plan_t* const plan,
                             const mod_t* const m, mpz_t g, const mpz_t n) {
    return FN(stage2_body)(acc, p, a24, plan, m, g, n, true);
}
