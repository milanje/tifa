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
 * \file    mont_ecm_wide.c
 *
 * \brief ECM curves, p-1 and p+1 in fixed-size Montgomery arithmetic, for
 *        moduli of five to eight 64-bit limbs, and the long curves on fewer
 *        (see mont_ecm.c).
 */

#include "mont_ecm.h"

#if TIFA_HAVE_MONT_ECM

#define MONT_MOD_LIMBS MONT_ECM_MAX_LIMBS
#include "mont_ecm_internal.h"
#include "stopping.h"
//
// A curve takes up to seconds here: it checks the time limit and the
// interrupt now and then (not on fewer limbs, whose code stays the same,
// but for the long curves).
//
#define MONT_ECM_STOP_CHECKS

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
#define NL 5
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL
#define NL 6
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL
#define NL 7
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL
#define NL 8
#include "mont_ecm_impl.h"
#include "mont_pm1_impl.h"
#undef NL

bool mont_ecm_wide_run(mpz_t factor, const mpz_t n, const mpz_t x0z,
                       const mpz_t a24z, const mont_ecm_plan_t* const plan) {
    mod_t m;
    mod_init(&m, n);

    uint64_t x0[MONT_MOD_LIMBS];
    uint64_t a24[MONT_MOD_LIMBS];
    to_limbs(x0, x0z, m.nl);
    to_limbs(a24, a24z, m.nl);

    int found = 0;
    switch (m.nl) {
        case 1:
            RUN_CURVE(1);
            break;
        case 2:
            RUN_CURVE(2);
            break;
        case 3:
            RUN_CURVE(3);
            break;
        case 4:
            RUN_CURVE(4);
            break;
        case 5:
            RUN_CURVE(5);
            break;
        case 6:
            RUN_CURVE(6);
            break;
        case 7:
            RUN_CURVE(7);
            break;
        default:
            RUN_CURVE(8);
            break;
    }
    return found > 0;
}

int mont_wide_pm1(mpz_t factor, const mpz_t n, const uint32_t base,
                  const mont_ecm_plan_t* const plan) {
    mod_t m;
    mod_init(&m, n);

    int found = 0;
    switch (m.nl) {
        case 5:
            RUN_PM1(5);
            break;
        case 6:
            RUN_PM1(6);
            break;
        case 7:
            RUN_PM1(7);
            break;
        default:
            RUN_PM1(8);
            break;
    }
    return found;
}

int mont_wide_pp1(mpz_t factor, const mpz_t n, const uint32_t seed,
                  const mont_ecm_plan_t* const plan) {
    mod_t m;
    mod_init(&m, n);

    int found = 0;
    switch (m.nl) {
        case 5:
            RUN_PP1(5);
            break;
        case 6:
            RUN_PP1(6);
            break;
        case 7:
            RUN_PP1(7);
            break;
        default:
            RUN_PP1(8);
            break;
    }
    return found;
}

#endif
