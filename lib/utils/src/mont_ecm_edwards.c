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
 * \brief ECM curves of z2z8_curves on 1 to 4 limbs, stage 1 on their
 *        Edwards form (see mont_ecm_impl.h), stage 2 as in mont_ecm.c.
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
static bool z2z8_edwards(mpz_t ex, mpz_t ey, mpz_t ez, mpz_t et, mpz_t d,
                         mpz_t u0, mpz_t a24, mpz_t factor, const mpz_t n,
                         const uint32_t i) {
    mpz_t dn, dd, yn, yd, xn, xd, a, y, t;
    mpz_inits(dn, dd, yn, yd, xn, xd, a, y, t, NULL);
    mpz_t* const row[6] = { &dn, &dd, &yn, &yd, &xn, &xd };
    for (int k = 0; k < 6; k++) {
        mpz_import(*row[k], 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][k]);
    }
    mpz_mul(ex, xn, yd);
    mpz_mod(ex, ex, n);
    mpz_mul(ey, yn, xd);
    mpz_mod(ey, ey, n);
    mpz_mul(ez, xd, yd);
    mpz_mod(ez, ez, n);
    mpz_mul(et, xn, yn);
    mpz_mod(et, et, n);
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
        mpz_mul(d, factor, a);
        mpz_mod(d, d, n);
        mpz_mul(d, d, y);
        mpz_mod(d, d, n);
        mpz_mul(d, d, dn);
        mpz_mod(d, d, n);
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

#define RUN_EDWARDS(NLIMBS)                                                   \
    do {                                                                      \
        FN_(pt_t, NLIMBS) p;                                                  \
        FN_(ext_t, NLIMBS) st;                                                \
        uint64_t dm[NLIMBS];                                                  \
        uint64_t x0m[NLIMBS];                                                 \
        uint64_t a24m[NLIMBS];                                                \
        TO_MONT_MPZ(NLIMBS, st.x, ex);                                        \
        TO_MONT_MPZ(NLIMBS, st.y, ey);                                        \
        TO_MONT_MPZ(NLIMBS, st.z, ez);                                        \
        TO_MONT_MPZ(NLIMBS, st.t, et);                                        \
        TO_MONT_MPZ(NLIMBS, dm, d);                                           \
        TO_MONT_MPZ(NLIMBS, x0m, u0);                                         \
        TO_MONT_MPZ(NLIMBS, a24m, a24);                                       \
        FN_(stage1_edwards, NLIMBS)(&p, &st, dm, plan, &m);                   \
        AFTER_STAGE1(NLIMBS);                                                 \
    } while (0)

bool mont_ecm_edwards_curve(mpz_t factor, const mpz_t n, const uint32_t i,
                            const mont_ecm_plan_t* const plan) {
    mpz_t ex, ey, ez, et, d, u0, a24;
    mpz_inits(ex, ey, ez, et, d, u0, a24, NULL);
    int found = 0;
    if (z2z8_edwards(ex, ey, ez, et, d, u0, a24, factor, n, i)) {
        mod_t m;
        mod_init(&m, n);
        switch (m.nl) {
            case 1:
                RUN_EDWARDS(1);
                break;
            case 2:
                RUN_EDWARDS(2);
                break;
            case 3:
                RUN_EDWARDS(3);
                break;
            default:
                RUN_EDWARDS(4);
                break;
        }
    } else {
        found = (mpz_cmp_ui(factor, 1) > 0) && (mpz_cmp(factor, n) < 0);
    }
    mpz_clears(ex, ey, ez, et, d, u0, a24, NULL);
    return found > 0;
}

#endif
