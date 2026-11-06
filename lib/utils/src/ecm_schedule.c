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
#include <pthread.h>

#include "ecm_schedule.h"
#include "tifa_alloc.h"
#include "stopping.h"

#if TIFA_HAVE_MONT_ECM

//
// ECM levels, tuned for Suyama's curves, though the curves run are those of
// mont_ecm_nth_curve (the Z/2 x Z/8 table first): the B1 for factors of
// fbits bits (so for the range (fbits - 5, fbits]), the mean number of
// curves that takes, and the time of one curve on two limbs (us, Apple M
// series), which only matters relative to siqs_cost. From 65 bits, B1
// doubles from level to level and the curve counts, from 20 to 80 factors
// per level, are rounded up.
//
static const struct {
    uint32_t fbits;
    uint32_t b1;
    uint32_t ncurves;
    double   cost;
} levels[] = {
    { 15,     50,  2,    11.7 },
    { 20,    100,  2,    13.5 },
    { 25,    200,  3,    22.9 },
    { 30,    400,  4,    41.3 },
    { 35,    800,  7,    77.4 },
    { 40,   1600,  8,   152.0 },
    { 45,   3200, 10,   280.0 },
    { 50,   6400, 15,   535.0 },
    { 55,  12800, 17,  1030.0 },
    { 60,  25600, 25,  2050.0 },
    { 65,  51200, 32,  3930.0 },
    { 70, 102400, 38,  7850.0 },
    { 75, 204800, 45, 15700.0 },
    { 80, 409600, 48, 30900.0 },
    { 85, 819200, 50, 63000.0 }
};
#define NLEVELS (sizeof(levels) / sizeof(levels[0]))
//
// Curve time on n limbs relative to two limbs.
//
static const double limb_factor[MONT_ECM_MAX_LIMBS + 1] = {
    0.0, 0.37, 1.0, 1.82, 2.85, 3.95, 6.19, 8.29, 10.9
};
//
// Time of SIQS on a balanced semiprime (us, same machine; from 264 bits,
// version 0.4.7's).
//
static const struct {
    double bits;
    double cost;
} siqs_costs[] = {
    {  64,      35 }, {  72,      54 }, {  80,      89 }, {  90,     177 },
    { 100,     378 }, { 110,     797 }, { 120,    1490 }, { 128,    2880 },
    { 140,    7240 }, { 160,   30400 }, { 180,  123000 }, { 200,  628000 },
    { 220, 3300000 }, { 264, 31400000 }, { 272, 61000000 },
    { 280, 98000000 }, { 288, 185000000 }
};
#define NSIQS_COSTS (sizeof(siqs_costs) / sizeof(siqs_costs[0]))

#define PM1_B1 315
//
// The plans of the levels up to EAGER_FBITS bits are built with the
// schedule, those above (7 to 70 ms each, about one curve) when a search
// first reaches them, under the lock.
//
#define EAGER_FBITS 60

struct struct_ecm_schedule_t {
    mont_ecm_plan_t* pm1_plan;
    uint32_t         nlevels;
    uint32_t         neager;
    pthread_mutex_t  lock;
    mont_ecm_plan_t* plans[NLEVELS];
};

//------------------------------------------------------------------------------
ecm_schedule_t* alloc_ecm_schedule(const uint32_t max_bits) {
    ecm_schedule_t* const s = tifa_malloc(sizeof(ecm_schedule_t));

    s->pm1_plan = alloc_mont_ecm_plan(PM1_B1, 50 * PM1_B1);
    s->nlevels  = 0;
    s->neager   = 0;
    pthread_mutex_init(&s->lock, NULL);
    while (   (s->nlevels < NLEVELS)
           && (levels[s->nlevels].fbits - 5 < max_bits)) {
        s->plans[s->nlevels] = NULL;
        if (levels[s->nlevels].fbits <= EAGER_FBITS) {
            s->plans[s->nlevels] = alloc_mont_ecm_plan(
                                       levels[s->nlevels].b1,
                                       50 * levels[s->nlevels].b1
                                   );
            s->neager++;
        }
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
    pthread_mutex_destroy(&schedule->lock);
    free(schedule);
}
//------------------------------------------------------------------------------
static const mont_ecm_plan_t* level_plan(ecm_schedule_t* const schedule,
                                         const uint32_t i) {
    if (i < schedule->neager) {
        return schedule->plans[i];
    }
    pthread_mutex_lock(&schedule->lock);
    if (schedule->plans[i] == NULL) {
        schedule->plans[i] = alloc_mont_ecm_plan(levels[i].b1,
                                                 50 * levels[i].b1);
    }
    const mont_ecm_plan_t* const plan = schedule->plans[i];
    pthread_mutex_unlock(&schedule->lock);
    return plan;
}
//------------------------------------------------------------------------------
static double siqs_cost(const double bits) {
    if (bits <= siqs_costs[0].bits) {
        return siqs_costs[0].cost;
    }
    uint32_t i = 1;
    while ((i < NSIQS_COSTS - 1) && (bits > siqs_costs[i].bits)) {
        i++;
    }
    //
    // Log-linear between two sizes, and past the last one.
    //
    const double t = (bits - siqs_costs[i - 1].bits)
                   / (siqs_costs[i].bits - siqs_costs[i - 1].bits);
    return exp((1 - t) * log(siqs_costs[i - 1].cost)
               + t * log(siqs_costs[i].cost));
}
//------------------------------------------------------------------------------
static bool run_levels(mpz_t g, const mpz_t n,
                       ecm_schedule_t* const schedule,
                       const double min_bits, const double target_bits,
                       const double budget, const bool whole_first_level) {
    const double limit = budget * siqs_cost((double) mpz_sizeinbase(n, 2));
    const double lfac  = limb_factor[mpz_size(n)];
    double   spent = 0;
    uint32_t k     = 0;

    for (uint32_t i = 0; i < schedule->nlevels; i++) {
        const double hi = (double) levels[i].fbits;
        if (hi <= min_bits) {
            continue;
        }
        if (hi - 5 >= target_bits) {
            break;
        }
        const bool first = whole_first_level && (spent == 0)
                           && (levels[i].cost * lfac <= limit);
        spent += levels[i].ncurves * levels[i].cost * lfac;
        if ((spent > limit) && !first) {
            break;
        }
        const mont_ecm_plan_t* const plan = level_plan(schedule, i);
        for (uint32_t j = 0; j < levels[i].ncurves; j++) {
            if (mont_ecm_nth_curve(g, n, k++, plan)) {
                return true;
            }
            if (tifa_should_stop()) {
                return false;
            }
        }
    }
    return false;
}
//------------------------------------------------------------------------------
bool ecm_schedule_split(mpz_t g, const mpz_t n,
                        ecm_schedule_t* const schedule,
                        const double min_bits, const double target_bits,
                        const double budget) {
    if (mont_pm1(g, n, 2, schedule->pm1_plan) == 1) {
        return true;
    }
    return run_levels(g, n, schedule, min_bits, target_bits, budget, false);
}
//------------------------------------------------------------------------------
bool ecm_schedule_prepass(mpz_t g, const mpz_t n,
                          ecm_schedule_t* const schedule,
                          const double budget) {
    if (mont_pm1(g, n, 2, schedule->pm1_plan) == 1) {
        return true;
    }
    return run_levels(g, n, schedule, 0, (double) mpz_sizeinbase(n, 2) / 2,
                      budget, true);
}

#endif
