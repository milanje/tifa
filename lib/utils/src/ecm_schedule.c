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

/**
 * \file    ecm_schedule.c
 *
 * \brief A bounded search for a factor of up to some size: p-1, then ECM.
 */

#include <stdlib.h>
#include <math.h>

#include "ecm_schedule.h"

#if TIFA_HAVE_MONT_ECM

//
// ECM levels (as in ecm.c): the B1 for factors of fbits bits (so for the
// range (fbits - 5, fbits]), the mean number of curves that takes, and the
// time of one curve on two limbs (us, Apple M series), which only matters
// relative to siqs_cost.
//
static const struct {
    uint32_t fbits;
    uint32_t b1;
    uint32_t ncurves;
    double   cost;
} levels[] = {
    { 15,    50,  2,    9.0 },
    { 20,   100,  2,   17.0 },
    { 25,   200,  3,   31.0 },
    { 30,   400,  4,   58.0 },
    { 35,   800,  7,  112.0 },
    { 40,  1600,  8,  220.0 },
    { 45,  3200, 10,  410.0 },
    { 50,  6400, 15,  800.0 },
    { 55, 12800, 17, 1600.0 },
    { 60, 25600, 25, 3200.0 }
};
#define NLEVELS (sizeof(levels) / sizeof(levels[0]))
//
// Curve time on n limbs relative to two limbs.
//
static const double limb_factor[MONT_ECM_MAX_LIMBS + 1] = {
    0.0, 0.6, 1.0, 2.0, 3.3
};
//
// Time of tifa_factor on a balanced semiprime (us, same machine).
//
static const struct {
    double bits;
    double cost;
} siqs_costs[] = {
    {  64,     200 }, {  72,    320 }, {  80,    480 }, {  90,    740 },
    { 100,    1400 }, { 110,   2600 }, { 120,   4500 }, { 128,   7300 },
    { 140,   15000 }, { 160,  60000 }, { 200, 1000000 }
};
#define NSIQS_COSTS (sizeof(siqs_costs) / sizeof(siqs_costs[0]))

#define PM1_B1 315

struct struct_ecm_schedule_t {
    mont_ecm_plan_t* pm1_plan;
    uint32_t         nlevels;
    mont_ecm_plan_t* plans[NLEVELS];
};

//------------------------------------------------------------------------------
ecm_schedule_t* alloc_ecm_schedule(const uint32_t max_bits) {
    ecm_schedule_t* const s = malloc(sizeof(ecm_schedule_t));

    s->pm1_plan = alloc_mont_ecm_plan(PM1_B1, 50 * PM1_B1);
    s->nlevels  = 0;
    while (   (s->nlevels < NLEVELS)
           && (levels[s->nlevels].fbits - 5 < max_bits)) {
        s->plans[s->nlevels] = alloc_mont_ecm_plan(levels[s->nlevels].b1,
                                                   50 * levels[s->nlevels].b1);
        s->nlevels++;
    }
    return s;
}
//------------------------------------------------------------------------------
void free_ecm_schedule(ecm_schedule_t* const schedule) {
    if (schedule == NULL) {
        return;
    }
    free_mont_ecm_plan(schedule->pm1_plan);
    for (uint32_t i = 0; i < schedule->nlevels; i++) {
        free_mont_ecm_plan(schedule->plans[i]);
    }
    free(schedule);
}
//------------------------------------------------------------------------------
static double siqs_cost(const double bits) {
    if (bits <= siqs_costs[0].bits) {
        return siqs_costs[0].cost;
    }
    for (uint32_t i = 1; i < NSIQS_COSTS; i++) {
        if (bits <= siqs_costs[i].bits) {
            const double t = (bits - siqs_costs[i - 1].bits)
                           / (siqs_costs[i].bits - siqs_costs[i - 1].bits);
            return exp(  (1 - t) * log(siqs_costs[i - 1].cost)
                       + t * log(siqs_costs[i].cost));
        }
    }
    return siqs_costs[NSIQS_COSTS - 1].cost;
}
//------------------------------------------------------------------------------
bool ecm_schedule_split(mpz_t g, const mpz_t n,
                        const ecm_schedule_t* const schedule,
                        const double min_bits, const double target_bits,
                        const double budget) {
    if (mont_pm1(g, n, 2, schedule->pm1_plan) == 1) {
        return true;
    }
    const double limit = budget * siqs_cost((double) mpz_sizeinbase(n, 2));
    const double lfac  = limb_factor[mpz_size(n)];
    double   spent = 0;
    uint32_t sigma = 6;

    for (uint32_t i = 0; i < schedule->nlevels; i++) {
        const double hi = (double) levels[i].fbits;
        if (hi <= min_bits) {
            continue;
        }
        if (hi - 5 >= target_bits) {
            break;
        }
        spent += levels[i].ncurves * levels[i].cost * lfac;
        if (spent > limit) {
            break;
        }
        for (uint32_t j = 0; j < levels[i].ncurves; j++) {
            if (mont_ecm_curve(g, n, sigma++, schedule->plans[i])) {
                return true;
            }
        }
    }
    return false;
}

#endif
