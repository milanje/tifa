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

#endif
