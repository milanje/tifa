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
 * \file    mont_ecm_edwards.c
 *
 * \brief Stage 1 of the curves of z2z8_curves on their Edwards form, on 1
 *        to 4 limbs (see mont_ecm_impl.h); mont_ecm.c runs the rest.
 */

#include "mont_ecm.h"

#if TIFA_HAVE_MONT_ECM

#include "z2z8_curves.h"

#define MONT_MOD_LIMBS MONT_ECM_NARROW_LIMBS
#include "mont_ecm_internal.h"

#define MONT_ECM_EDWARDS

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

//
// Curve i: the Edwards starting point (x0 : y0 : 1) = (xn.yd : yn.xd : xd.yd),
// t = xn.yn, d = dn/dd, and for the Montgomery form, a24 = 1/(1 - d) and
// u0 = (1 + y0)/(1 - y0), from one inversion.
//
bool mont_ecm_edwards_setup(mont_ecm_edwards_t* const ed, mpz_t u0,
                            mpz_t a24, mpz_t factor, const mpz_t n,
                            const uint32_t i) {
    mpz_t dn, dd, yn, yd, xn, xd, a, y, t;
    mpz_inits(dn, dd, yn, yd, xn, xd, a, y, t, NULL);
    mpz_t* const row[6] = { &dn, &dd, &yn, &yd, &xn, &xd };
    for (int k = 0; k < 6; k++) {
        mpz_import(*row[k], 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][k]);
    }
    mpz_mul(ed->x, xn, yd);
    mpz_mod(ed->x, ed->x, n);
    mpz_mul(ed->y, yn, xd);
    mpz_mod(ed->y, ed->y, n);
    mpz_mul(ed->z, xd, yd);
    mpz_mod(ed->z, ed->z, n);
    mpz_mul(ed->t, xn, yn);
    mpz_mod(ed->t, ed->t, n);
    //
    // With a = dd - dn, y = yd - yn and inv = 1/(a.y.dd): a24 = dd.y.dd.inv,
    // d = dn.a.y.inv, u0 = (yd + yn).a.dd.inv.
    //
    mpz_sub(a, dd, dn);
    mpz_sub(y, yd, yn);
    mpz_mul(t, a, y);
    mpz_mul(t, t, dd);
    mpz_mod(t, t, n);
    const bool ok = (mpz_invert(factor, t, n) != 0);
    if (!ok) {
        mpz_gcd(factor, t, n);
    } else {
        mpz_mul(a24, factor, y);
        mpz_mod(a24, a24, n);
        mpz_mul(a24, a24, dd);
        mpz_mod(a24, a24, n);
        mpz_mul(a24, a24, dd);
        mpz_mod(a24, a24, n);
        mpz_mul(ed->d, factor, a);
        mpz_mod(ed->d, ed->d, n);
        mpz_mul(ed->d, ed->d, y);
        mpz_mod(ed->d, ed->d, n);
        mpz_mul(ed->d, ed->d, dn);
        mpz_mod(ed->d, ed->d, n);
        mpz_mul(u0, factor, a);
        mpz_mod(u0, u0, n);
        mpz_mul(u0, u0, dd);
        mpz_mod(u0, u0, n);
        mpz_add(t, yd, yn);
        mpz_mul(u0, u0, t);
        mpz_mod(u0, u0, n);
    }
    mpz_clears(dn, dd, yn, yd, xn, xd, a, y, t, NULL);
    return ok;
}


//
// k.P on v^2 = u^3 - 12.u (P = (-2, 4)), Jacobian (X : Y : Z), u = X/Z^2,
// v = Y/Z^3: doublings, and additions of P. k >= 1.
//
static void z12_point(mpz_t x, mpz_t y, mpz_t z, const mpz_t n,
                      const uint32_t k) {
    mpz_t a, b, c, d, e, f;
    mpz_inits(a, b, c, d, e, f, NULL);
    mpz_set_si(x, -2);
    mpz_mod(x, x, n);
    mpz_set_ui(y, 4);
    mpz_set_ui(z, 1);
    for (int i = 30 - __builtin_clz(k); i >= 0; i--) {
        //
        // Doubling: M = 3X^2 - 12Z^4, S = 4XY^2, X' = M^2 - 2S,
        // Y' = M(S - X') - 8Y^4, Z' = 2YZ.
        //
        mpz_mul(a, x, x);
        mpz_mul(b, z, z);
        mpz_mul(b, b, b);
        mpz_mul_ui(a, a, 3);
        mpz_submul_ui(a, b, 12);
        mpz_mod(a, a, n);
        mpz_mul(c, y, y);
        mpz_mod(c, c, n);
        mpz_mul(d, x, c);
        mpz_mul_2exp(d, d, 2);
        mpz_mod(d, d, n);
        mpz_mul(z, z, y);
        mpz_mul_2exp(z, z, 1);
        mpz_mod(z, z, n);
        mpz_mul(x, a, a);
        mpz_submul_ui(x, d, 2);
        mpz_mod(x, x, n);
        mpz_sub(e, d, x);
        mpz_mul(y, a, e);
        mpz_mul(c, c, c);
        mpz_submul_ui(y, c, 8);
        mpz_mod(y, y, n);
        if ((k >> i) & 1) {
            //
            // Adding P = (-2, 4): U2 = -2Z^2, S2 = 4Z^3, H = U2 - X,
            // R = S2 - Y, X' = R^2 - H^3 - 2X.H^2,
            // Y' = R(X.H^2 - X') - Y.H^3, Z' = Z.H.
            //
            mpz_mul(a, z, z);
            mpz_mod(a, a, n);
            mpz_mul(b, a, z);
            mpz_mul_2exp(b, b, 2);
            mpz_mod(b, b, n);
            mpz_mul_2exp(a, a, 1);
            mpz_neg(a, a);
            mpz_sub(e, a, x);
            mpz_mod(e, e, n);
            mpz_sub(f, b, y);
            mpz_mod(f, f, n);
            mpz_mul(c, e, e);
            mpz_mod(c, c, n);
            mpz_mul(d, c, e);
            mpz_mod(d, d, n);
            mpz_mul(c, c, x);
            mpz_mod(c, c, n);
            mpz_mul(x, f, f);
            mpz_sub(x, x, d);
            mpz_submul_ui(x, c, 2);
            mpz_mod(x, x, n);
            mpz_sub(c, c, x);
            mpz_mul(c, c, f);
            mpz_submul(c, y, d);
            mpz_mod(y, c, n);
            mpz_mul(z, z, e);
            mpz_mod(z, z, n);
        }
    }
    mpz_clears(a, b, c, d, e, f, NULL);
}

//
// Curve k (k >= 2) of Montgomery's family with torsion Z/12, (u, v) = k.P on
// v^2 = u^3 - 12.u (Montgomery's thesis, 6.1), in Edwards form: with
// d1 = v^2 + 12u^2, d2 = v^2 - 4u^2, the Montgomery curve has
// A = (4(32vu^3)^2 - 2d1.d2^3)/(d1.d2^3), a24 = (32vu^3)^2/(d1.d2^3) and the
// starting point (x0 : y0 : z0) = (2(u(u^2 + 12))^2 : d1.u(u^2 + 12) : 2d1.d2);
// (A + 2)/B = s^2, s = 4uv/(v^2 + 4u^2), so that x^2 + y^2 = 1 + d.x^2.y^2,
// d = (A - 2)/(A + 2), holds the point (s.x0/y0, (x0 - z0)/(x0 + z0)).
//
bool mont_ecm_z12_setup(mont_ecm_edwards_t* const ed, mpz_t u0, mpz_t a24,
                        mpz_t factor, const mpz_t n, const uint32_t k) {
    mpz_t u, v, w, uu, vv, d1, d2, x0, y0, z0, an, bn, ad, t;
    mpz_inits(u, v, w, uu, vv, d1, d2, x0, y0, z0, an, bn, ad, t, NULL);
    bool ok = false;

    z12_point(u, v, w, n, k);
    if (!mpz_invert(t, w, n)) {
        mpz_gcd(factor, w, n);
        goto done;
    }
    mpz_mul(w, t, t);
    mpz_mul(u, u, w);
    mpz_mod(u, u, n);
    mpz_mul(w, w, t);
    mpz_mul(v, v, w);
    mpz_mod(v, v, n);

    mpz_mul(uu, u, u);
    mpz_mod(uu, uu, n);
    mpz_mul(vv, v, v);
    mpz_mod(vv, vv, n);
    mpz_set(d1, vv);
    mpz_addmul_ui(d1, uu, 12);
    mpz_mod(d1, d1, n);
    mpz_set(d2, vv);
    mpz_submul_ui(d2, uu, 4);
    mpz_mod(d2, d2, n);
    mpz_add_ui(t, uu, 12);
    mpz_mul(t, t, u);
    mpz_mod(t, t, n);                       // u(u^2 + 12)
    mpz_mul(x0, t, t);
    mpz_mul_2exp(x0, x0, 1);
    mpz_mod(x0, x0, n);
    mpz_mul(y0, d1, t);
    mpz_mod(y0, y0, n);
    mpz_mul(z0, d1, d2);
    mpz_mul_2exp(z0, z0, 1);
    mpz_mod(z0, z0, n);

    mpz_mul(bn, uu, u);
    mpz_mul(bn, bn, v);
    mpz_mul_2exp(bn, bn, 5);
    mpz_mod(bn, bn, n);
    mpz_mul(bn, bn, bn);
    mpz_mod(bn, bn, n);                     // (32vu^3)^2
    mpz_mul(ad, d2, d2);
    mpz_mod(ad, ad, n);
    mpz_mul(ad, ad, d2);
    mpz_mul(ad, ad, d1);
    mpz_mod(ad, ad, n);                     // d1.d2^3
    //
    // a24 = bn/ad, d = (A - 2)/(A + 2) = (bn - ad)/bn, u0 = x0/z0: one
    // inversion of ad.bn.z0.
    //
    mpz_mul(t, ad, bn);
    mpz_mod(t, t, n);
    mpz_mul(t, t, z0);
    mpz_mod(t, t, n);
    if (!mpz_invert(factor, t, n)) {
        mpz_gcd(factor, t, n);
        goto done;
    }
    mpz_mul(t, factor, bn);
    mpz_mod(t, t, n);                       // 1/(ad.z0)
    mpz_mul(a24, t, z0);
    mpz_mul(a24, a24, bn);
    mpz_mod(a24, a24, n);
    mpz_mul(u0, t, ad);
    mpz_mul(u0, u0, x0);
    mpz_mod(u0, u0, n);
    mpz_mul(t, factor, ad);
    mpz_mul(t, t, z0);
    mpz_mod(t, t, n);                       // 1/bn
    mpz_sub(ed->d, bn, ad);
    mpz_mul(ed->d, ed->d, t);
    mpz_mod(ed->d, ed->d, n);
    //
    // (X : Y : Z) = (4uv.x0(x0 + z0) : (x0 - z0)y0(v^2 + 4u^2) :
    // y0(x0 + z0)(v^2 + 4u^2)), then extended (XZ : YZ : Z^2 : XY).
    //
    mpz_add(an, x0, z0);
    mpz_set(w, vv);
    mpz_addmul_ui(w, uu, 4);
    mpz_mod(w, w, n);
    mpz_mul(t, u, v);
    mpz_mul_2exp(t, t, 2);
    mpz_mul(t, t, x0);
    mpz_mod(t, t, n);
    mpz_mul(t, t, an);
    mpz_mod(t, t, n);                       // X
    mpz_sub(uu, x0, z0);
    mpz_mul(uu, uu, y0);
    mpz_mod(uu, uu, n);
    mpz_mul(uu, uu, w);
    mpz_mod(uu, uu, n);                     // Y
    mpz_mul(vv, y0, an);
    mpz_mod(vv, vv, n);
    mpz_mul(vv, vv, w);
    mpz_mod(vv, vv, n);                     // Z
    mpz_mul(ed->x, t, vv);
    mpz_mod(ed->x, ed->x, n);
    mpz_mul(ed->y, uu, vv);
    mpz_mod(ed->y, ed->y, n);
    mpz_mul(ed->z, vv, vv);
    mpz_mod(ed->z, ed->z, n);
    mpz_mul(ed->t, t, uu);
    mpz_mod(ed->t, ed->t, n);
    ok = true;
done:
    mpz_clears(u, v, w, uu, vv, d1, d2, x0, y0, z0, an, bn, ad, t, NULL);
    return ok;
}

#endif
