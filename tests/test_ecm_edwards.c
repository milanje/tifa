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

// The prime iterator (primes.h) against mpz_nextprime, across the end of the
// table of first_primes.h, the segments of the sieve and 2^32 - 1, and the
// stage-1 exponent against a direct product.


//
// edwards_stage1_u64 against k.P computed by affine double-and-add, on
// the curves of z2z8_curves, modulo a prime: u = (1 + y)/(1 - y) must agree.
//

#include <stdlib.h>
#include <gmp.h>

#include "tifa_test.h"
#include "factor_u64.h"
#include "z2z8_curves.h"

#if TIFA_HAVE_FACTOR_U64

#include "factor_u64_internal.h"

static mpz_t P, D;

static void frac(mpz_t r, uint64_t num, uint64_t den) {
    mpz_t t;
    mpz_init_set_ui(t, den);
    mpz_invert(t, t, P);
    mpz_mul_ui(r, t, num);
    mpz_mod(r, r, P);
    mpz_clear(t);
}

//
// (x1, y1) += (x2, y2) on x^2 + y^2 = 1 + d.x^2.y^2. False if a denominator
// vanishes (exceptional points; not met with these curves and primes).
//
static bool affine_add(mpz_t x1, mpz_t y1, const mpz_t x2, const mpz_t y2) {
    mpz_t t, num, den, x3;
    mpz_inits(t, num, den, x3, NULL);
    mpz_mul(t, x1, x2);
    mpz_mul(t, t, y1);
    mpz_mul(t, t, y2);
    mpz_mul(t, t, D);
    mpz_mod(t, t, P);
    bool ok = true;
    mpz_add_ui(den, t, 1);
    mpz_mul(num, x1, y2);
    mpz_addmul(num, y1, x2);
    ok = ok && mpz_invert(den, den, P);
    mpz_mul(x3, num, den);
    mpz_mod(x3, x3, P);
    mpz_ui_sub(den, 1, t);
    mpz_mod(den, den, P);
    mpz_mul(num, y1, y2);
    mpz_submul(num, x1, x2);
    ok = ok && mpz_invert(den, den, P);
    mpz_mul(y1, num, den);
    mpz_mod(y1, y1, P);
    mpz_set(x1, x3);
    mpz_clears(t, num, den, x3, NULL);
    return ok;
}

//
// (x, y) = k.(x0, y0), k >= 1: the addition law also doubles.
//
static bool affine_mul(mpz_t x, mpz_t y, const mpz_t x0, const mpz_t y0,
                       uint64_t k) {
    mpz_set(x, x0);
    mpz_set(y, y0);
    bool ok = true;
    for (int b = 62 - __builtin_clzll(k); b >= 0; b--) {
        ok = ok && affine_add(x, y, x, y);
        if ((k >> b) & 1) {
            ok = ok && affine_add(x, y, x0, y0);
        }
    }
    return ok;
}

static uint64_t mont_of(const mpz_t a, const mont_t* m) {
    return to_mont(mpz_get_ui(a), m);
}

//
// Width-w signed digits of k, most significant first.
//
static uint32_t recode(int8_t* digits, uint64_t k, uint32_t w) {
    int8_t tmp[80];
    uint32_t n = 0;
    while (k != 0) {
        int d = 0;
        if (k & 1) {
            d = (int) (k & ((1u << w) - 1));
            if (d >= (1 << (w - 1))) {
                d -= 1 << w;
            }
            k -= (uint64_t) (int64_t) d;
        }
        tmp[n++] = (int8_t) d;
        k >>= 1;
    }
    for (uint32_t i = 0; i < n; i++) {
        digits[i] = tmp[n - 1 - i];
    }
    return n;
}

int main(void) {
    TIFA_TEST_BEGIN("ecm_edwards");

    const uint64_t p = 4611686018427387847ULL;   // 2^62 - 57, prime
    mpz_init_set_ui(P, p);
    mpz_init(D);
    mont_t m;
    m.n = p;
    uint64_t x = (3 * p) ^ 2;
    for (int i = 0; i < 4; i++) {
        x *= 2 - p * x;
    }
    m.ninv = x;
    m.one  = (uint64_t) ((((u128) 1) << 64) % p);
    m.r2   = (uint64_t) (((u128) m.one * m.one) % p);

    srand(20261109);
    uint32_t nbad = 0, nchecks = 0;
    for (uint32_t first = 0; first + ECM_LANES <= NZ2Z8_CURVES;
         first += ECM_LANES) {
        mpz_t x0[ECM_LANES], y0[ECM_LANES], xk[ECM_LANES], yk[ECM_LANES];
        mpz_t d[ECM_LANES];
        edwards_start_t start[ECM_LANES];
        for (int c = 0; c < ECM_LANES; c++) {
            const uint64_t* const row = z2z8_curves[first + c];
            mpz_inits(x0[c], y0[c], xk[c], yk[c], d[c], NULL);
            frac(d[c], row[0], row[1]);
            frac(y0[c], row[2], row[3]);
            frac(x0[c], row[4], row[5]);
            start[c].x = mont_of(x0[c], &m);
            start[c].y = mont_of(y0[c], &m);
            start[c].z = m.one;
            start[c].t = mont_mul(start[c].x, start[c].y, &m);
            start[c].d = mont_of(d[c], &m);
        }
        for (uint32_t w = 2; w <= 4; w++) {
            for (uint32_t trial = 0; trial < 40; trial++) {
                uint64_t r = (uint64_t) rand() << 20;
                r ^= (uint64_t) rand();
                const uint64_t k = (trial < 20) ? trial + 1
                                                : 1 + r % (1ULL << 40);
                int8_t digits[80];
                const uint32_t nd = recode(digits, k, w);
                xz_t out[ECM_LANES];
                edwards_stage1_u64(out, start, digits, nd, w, &m);
                for (int c = 0; c < ECM_LANES; c++) {
                    mpz_set(D, d[c]);
                    const bool ok = affine_mul(xk[c], yk[c], x0[c], y0[c], k);
                    // u = (1 + y)/(1 - y): X.(1 - y) = Z.(1 + y).
                    mpz_t one_p, one_m;
                    mpz_inits(one_p, one_m, NULL);
                    mpz_add_ui(one_p, yk[c], 1);
                    mpz_ui_sub(one_m, 1, yk[c]);
                    mpz_mod(one_m, one_m, P);
                    mpz_mod(one_p, one_p, P);
                    const uint64_t lhs = mont_mul(out[c].x, mont_of(one_m, &m),
                                                  &m);
                    const uint64_t rhs = mont_mul(out[c].z, mont_of(one_p, &m),
                                                  &m);
                    nbad += !ok || (lhs != rhs);
                    nchecks++;
                    mpz_clears(one_p, one_m, NULL);
                }
            }
        }
        for (int c = 0; c < ECM_LANES; c++) {
            mpz_clears(x0[c], y0[c], xk[c], yk[c], d[c], NULL);
        }
    }
    CHECK(nchecks == 9 * 3 * 40 * 3, "every curve, width and multiplier");
    CHECK_EQ_UI(nbad, 0, "Edwards stage 1 matches affine k.P");
    mpz_clears(P, D, NULL);

    TIFA_TEST_END();
}

#else

int main(void) {
    TIFA_TEST_BEGIN("ecm_edwards");
    TIFA_TEST_END();
}

#endif
