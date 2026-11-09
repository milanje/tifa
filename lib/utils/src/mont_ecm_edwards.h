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
// The curves whose stage 1 runs on their Edwards form, set up as mpz
// (mont_ecm_edwards.c); also used by the tests.
//

#if !defined(_TIFA_MONT_ECM_EDWARDS_H_)
#define _TIFA_MONT_ECM_EDWARDS_H_

#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

//
// Curve i of z2z8_curves in Edwards form (mont_ecm_edwards.c): the starting
// point (x : y : z : t) and d, and for its Montgomery form u0 and a24. False,
// with the gcd in factor, if a value was not invertible mod n.
//
typedef struct {
    mpz_t x;
    mpz_t y;
    mpz_t z;
    mpz_t t;
    mpz_t d;
} mont_ecm_edwards_t;

bool mont_ecm_edwards_setup(mont_ecm_edwards_t* ed, mpz_t u0, mpz_t a24,
                            mpz_t factor, const mpz_t n, uint32_t i);

//
// As mont_ecm_edwards_setup, for curve k >= 2 of Montgomery's family with
// torsion Z/12 (mont_ecm_edwards.c).
//
bool mont_ecm_z12_setup(mont_ecm_edwards_t* ed, mpz_t u0, mpz_t a24,
                        mpz_t factor, const mpz_t n, uint32_t k);

#endif
