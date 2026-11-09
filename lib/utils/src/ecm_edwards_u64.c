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

/**
 * \file    ecm_edwards_u64.c
 *
 * \brief Stage 1 of ECM below 2^64 on Edwards curves x^2 + y^2 = 1 + d.x^2.y^2.
 *
 * The curves of z2z8_curves are Edwards curves: there, a doubling costs 7
 * products rather than the 10 of a step of the Montgomery ladder, and an
 * addition, needed once every few doublings with a signed window over E, 9.
 * Formulas of Hisil, Wong, Carter and Dawson (2008), projective (x : y : z)
 * between doublings and extended (x : y : z : t), t = xy/z, before additions.
 * With d a square (full 2-torsion), they have exceptional points: meeting one
 * modulo a prime factor p of n derails the computation modulo p alone,
 * mostly into z = 0, which reveals p all the same (8 numbers in 10^5, primes
 * of 25 to 32 bits).
 */

#include "factor_u64.h"

#if TIFA_HAVE_FACTOR_U64

#include "factor_u64_internal.h"

typedef struct {
    uint64_t x;
    uint64_t y;
    uint64_t z;
} proj_t;

typedef struct {
    uint64_t x;
    uint64_t y;
    uint64_t z;
    uint64_t t;
} ext_t;

//
// A multiple of the starting point, ready to be added: s = x + y, dt = d.t.
//
typedef struct {
    uint64_t x;
    uint64_t y;
    uint64_t s;
    uint64_t z;
    uint64_t dt;
} table_t;

#define MAX_TABLE 8

static inline proj_t dbl(const proj_t p, const mont_t* const m) {
    const uint64_t xy = mont_add(p.x, p.y, m);
    const uint64_t b  = mont_mul(xy, xy, m);
    const uint64_t c  = mont_mul(p.x, p.x, m);
    const uint64_t d  = mont_mul(p.y, p.y, m);
    const uint64_t h  = mont_mul(p.z, p.z, m);
    const uint64_t e  = mont_add(c, d, m);
    const uint64_t j  = mont_sub(e, mont_add(h, h, m), m);
    proj_t r;
    r.x = mont_mul(mont_sub(b, e, m), j, m);
    r.y = mont_mul(e, mont_sub(c, d, m), m);
    r.z = mont_mul(e, j, m);
    return r;
}

static inline ext_t dbl_ext(const proj_t p, const mont_t* const m) {
    const uint64_t xy = mont_add(p.x, p.y, m);
    const uint64_t a  = mont_mul(p.x, p.x, m);
    const uint64_t b  = mont_mul(p.y, p.y, m);
    const uint64_t zz = mont_mul(p.z, p.z, m);
    const uint64_t s  = mont_mul(xy, xy, m);
    const uint64_t g  = mont_add(a, b, m);
    const uint64_t e  = mont_sub(s, g, m);
    const uint64_t f  = mont_sub(g, mont_add(zz, zz, m), m);
    const uint64_t h  = mont_sub(a, b, m);
    ext_t r;
    r.x = mont_mul(e, f, m);
    r.y = mont_mul(g, h, m);
    r.t = mont_mul(e, h, m);
    r.z = mont_mul(f, g, m);
    return r;
}

//
// p + q: A = x1.x2, B = y1.y2, C = t1.d.t2, D = z1.z2,
// E = (x1 + y1)(x2 + y2) - A - B, F = D - C, G = D + C, H = B - A, and
// (x3 : y3 : z3 : t3) = (EF : GH : FG : EH).
//
static inline proj_t add(const ext_t p, const table_t* const q,
                         const mont_t* const m) {
    const uint64_t a = mont_mul(p.x, q->x, m);
    const uint64_t b = mont_mul(p.y, q->y, m);
    const uint64_t s = mont_mul(mont_add(p.x, p.y, m), q->s, m);
    const uint64_t c = mont_mul(p.t, q->dt, m);
    const uint64_t d = mont_mul(p.z, q->z, m);
    const uint64_t e = mont_sub(s, mont_add(a, b, m), m);
    const uint64_t f = mont_sub(d, c, m);
    const uint64_t g = mont_add(d, c, m);
    const uint64_t h = mont_sub(b, a, m);
    proj_t r;
    r.x = mont_mul(e, f, m);
    r.y = mont_mul(g, h, m);
    r.z = mont_mul(f, g, m);
    return r;
}

static inline ext_t add_ext(const ext_t p, const table_t* const q,
                            const mont_t* const m) {
    const uint64_t a = mont_mul(p.x, q->x, m);
    const uint64_t b = mont_mul(p.y, q->y, m);
    const uint64_t s = mont_mul(mont_add(p.x, p.y, m), q->s, m);
    const uint64_t c = mont_mul(p.t, q->dt, m);
    const uint64_t d = mont_mul(p.z, q->z, m);
    const uint64_t e = mont_sub(s, mont_add(a, b, m), m);
    const uint64_t f = mont_sub(d, c, m);
    const uint64_t g = mont_add(d, c, m);
    const uint64_t h = mont_sub(b, a, m);
    ext_t r;
    r.x = mont_mul(e, f, m);
    r.y = mont_mul(g, h, m);
    r.t = mont_mul(e, h, m);
    r.z = mont_mul(f, g, m);
    return r;
}

static inline void to_table(table_t* const r, const ext_t p, const uint64_t d,
                            const mont_t* const m) {
    r->x  = p.x;
    r->y  = p.y;
    r->s  = mont_add(p.x, p.y, m);
    r->z  = p.z;
    r->dt = mont_mul(p.t, d, m);
}

//
// -(x : y : z : t) = (-x : y : z : -t).
//
static inline void negate(table_t* const r, const table_t* const p,
                          const mont_t* const m) {
    r->x  = mont_sub(0, p->x, m);
    r->y  = p->y;
    r->s  = mont_sub(p->y, p->x, m);
    r->z  = p->z;
    r->dt = mont_sub(0, p->dt, m);
}

void edwards_stage1_u64(xz_t* const out, const edwards_start_t* const start,
                        const int8_t* const digits, const uint32_t ndigits,
                        const uint32_t w, const mont_t* const m) {
    //
    // table[c][j] = (2j + 1).P_c and table[c][nt + j] its opposite.
    //
    const uint32_t nt = 1u << (w - 2);
    table_t table[ECM_LANES][2 * MAX_TABLE];

    for (int c = 0; c < ECM_LANES; c++) {
        const ext_t p = { start[c].x, start[c].y, start[c].z, start[c].t };
        to_table(&table[c][0], p, start[c].d, m);
        if (nt > 1) {
            const proj_t pp = { p.x, p.y, p.z };
            table_t two;
            to_table(&two, dbl_ext(pp, m), start[c].d, m);
            ext_t q = p;
            for (uint32_t j = 1; j < nt; j++) {
                q = add_ext(q, &two, m);
                to_table(&table[c][j], q, start[c].d, m);
            }
        }
        for (uint32_t j = 0; j < nt; j++) {
            negate(&table[c][nt + j], &table[c][j], m);
        }
    }
    //
    // The leading digit is positive.
    //
    proj_t r[ECM_LANES];
    for (int c = 0; c < ECM_LANES; c++) {
        const table_t* const t = &table[c][(digits[0] - 1) / 2];
        r[c].x = t->x;
        r[c].y = t->y;
        r[c].z = t->z;
    }
    for (uint32_t i = 1; i < ndigits; i++) {
        const int dg = digits[i];
        if (dg == 0) {
            for (int c = 0; c < ECM_LANES; c++) {
                r[c] = dbl(r[c], m);
            }
        } else {
            const uint32_t k = (dg > 0) ? (uint32_t) (dg - 1) / 2
                                        : nt + (uint32_t) (-dg - 1) / 2;
            for (int c = 0; c < ECM_LANES; c++) {
                r[c] = add(dbl_ext(r[c], m), &table[c][k], m);
            }
        }
    }
    //
    // To the Montgomery curve of the same group, u = (1 + y)/(1 - y).
    //
    for (int c = 0; c < ECM_LANES; c++) {
        out[c].x = mont_add(r[c].z, r[c].y, m);
        out[c].z = mont_sub(r[c].z, r[c].y, m);
    }
}

#endif
