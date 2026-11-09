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
// Shared by factor_u64.c and ecm_edwards_u64.c: residues mod a 64-bit n in
// Montgomery form, and points of ECM's curves.
//

#if !defined(_TIFA_FACTOR_U64_INTERNAL_H_)
#define _TIFA_FACTOR_U64_INTERNAL_H_

#include <stdint.h>

typedef unsigned __int128 u128;

//------------------------------------------------------------------------------
//                         Montgomery arithmetic mod n
//------------------------------------------------------------------------------

//
// Residues are kept in Montgomery form (a.2^64 mod n) and in [0, n). n odd.
//
typedef struct {
    uint64_t n;
    uint64_t ninv;   // n^-1 mod 2^64
    uint64_t one;    // 2^64 mod n
    uint64_t r2;     // 2^128 mod n
} mont_t;

static inline uint64_t mont_mul(const uint64_t a, const uint64_t b,
                                const mont_t* const m) {
    const u128     t  = (u128) a * b;
    const uint64_t lo = (uint64_t) t;
    const uint64_t hi = (uint64_t) (t >> 64);
    const uint64_t q  = lo * m->ninv;
    const uint64_t qn = (uint64_t) (((u128) q * m->n) >> 64);
    return (hi < qn) ? (hi - qn + m->n) : (hi - qn);
}

static inline uint64_t mont_add(const uint64_t a, const uint64_t b,
                                const mont_t* const m) {
    const uint64_t s = a + b;
    return ((s < a) || (s >= m->n)) ? (s - m->n) : s;
}

static inline uint64_t mont_sub(const uint64_t a, const uint64_t b,
                                const mont_t* const m) {
    return (a < b) ? (a - b + m->n) : (a - b);
}

static inline uint64_t to_mont(const uint64_t a, const mont_t* const m) {
    return mont_mul(a % m->n, m->r2, m);
}

typedef struct {
    uint64_t x;
    uint64_t z;
} xz_t;

//
// Curves run ECM_LANES at a time: they share E and the stage 2 schedule, so
// they proceed in lockstep, and their independent multiplications overlap
// (1.5x to 2x faster than one at a time, even counting the curves run past
// the one that succeeds).
//
#define ECM_LANES 3

typedef struct {
    uint64_t x;
    uint64_t y;
    uint64_t z;
    uint64_t t;
    uint64_t d;
} edwards_start_t;

//
// out[c] = E.P_c, as Montgomery (X : Z), for the ECM_LANES curves
// x^2 + y^2 = 1 + d.x^2.y^2 and points P_c = (x : y : z : t) of start, E
// given by its signed digits of width w, most significant first.
//
void edwards_stage1_u64(xz_t* out, const edwards_start_t* start,
                        const int8_t* digits, uint32_t ndigits, uint32_t w,
                        const mont_t* m);

#endif
