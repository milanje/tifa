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
#include "macros.h"
#include "stopping.h"

#if TIFA_HAVE_MONT_ECM

//
// ECM levels, tuned for Suyama's curves, though the curves run are those of
// mont_ecm_nth_curve (the Z/2 x Z/8 table first): the B1 for factors of
// fbits bits (so for the range (fbits - 5, fbits]), the mean number of
// curves that takes, and the time of one curve on two limbs (us, Apple M
// series), which only matters relative to siqs_cost. From 65 bits, B1
// doubles from level to level and the curve counts, from 20 to 80 factors
// per level, are rounded up; from 90, they are Dickman's estimates scaled to
// those counted on 219 factors of 85 to 100 bits.
//
static const struct {
    uint32_t fbits;
    uint32_t b1;
    uint32_t ncurves;
    double   cost;
} levels[] = {
    {  15,      50,  2,     11.7 },
    {  20,     100,  2,     13.5 },
    {  25,     200,  3,     22.9 },
    {  30,     400,  4,     41.3 },
    {  35,     800,  7,     77.4 },
    {  40,    1600,  8,    152.0 },
    {  45,    3200, 10,    280.0 },
    {  50,    6400, 15,    535.0 },
    {  55,   12800, 17,   1030.0 },
    {  60,   25600, 25,   2050.0 },
    {  65,   51200, 32,   3930.0 },
    {  70,  102400, 38,   7850.0 },
    {  75,  204800, 45,  15700.0 },
    {  80,  409600, 48,  30900.0 },
    {  85,  819200, 50,  63000.0 },
    {  90, 1638400, 65, 123000.0 },
    {  95, 3276800, 75, 246000.0 },
    { 100, 6553600, 85, 489000.0 }
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
// schedule, those above (7 to 380 ms each, about one curve) when a search
// first reaches them, under the lock. Those above KEEP_FBITS bits (the
// default ecm_max_bits), 57 MB in all, are freed when the last search using
// them ends.
//
#define EAGER_FBITS 60
#define KEEP_FBITS  85

struct struct_ecm_schedule_t {
    mont_ecm_plan_t* pm1_plan;
    uint32_t         nlevels;
    uint32_t         neager;
    pthread_mutex_t  lock;
    mont_ecm_plan_t* plans[NLEVELS];
    uint32_t         users[NLEVELS];
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
        s->users[s->nlevels] = 0;
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
static const mont_ecm_plan_t* acquire_plan(ecm_schedule_t* const schedule,
                                           const uint32_t i) {
    if (i < schedule->neager) {
        return schedule->plans[i];
    }
    pthread_mutex_lock(&schedule->lock);
    if ((schedule->plans[i] == NULL) && (levels[i].fbits > KEEP_FBITS)) {
        //
        // Outside the lock: up to 0.4 s, which the other searches would
        // wait, each time the last search using the plan has freed it.
        //
        pthread_mutex_unlock(&schedule->lock);
        mont_ecm_plan_t* const built = alloc_mont_ecm_plan(levels[i].b1,
                                                           50 * levels[i].b1);
        pthread_mutex_lock(&schedule->lock);
        if (schedule->plans[i] == NULL) {
            schedule->plans[i] = built;
        } else {
            free_mont_ecm_plan(built);
        }
    }
    if (schedule->plans[i] == NULL) {
        schedule->plans[i] = alloc_mont_ecm_plan(levels[i].b1,
                                                 50 * levels[i].b1);
    }
    if (levels[i].fbits > KEEP_FBITS) {
        schedule->users[i]++;
    }
    const mont_ecm_plan_t* const plan = schedule->plans[i];
    pthread_mutex_unlock(&schedule->lock);
    return plan;
}

static void release_plan(ecm_schedule_t* const schedule, const uint32_t i) {
    if (levels[i].fbits <= KEEP_FBITS) {
        return;
    }
    pthread_mutex_lock(&schedule->lock);
    if (--schedule->users[i] == 0) {
        free_mont_ecm_plan(schedule->plans[i]);
        schedule->plans[i] = NULL;
    }
    pthread_mutex_unlock(&schedule->lock);
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
typedef struct {
    double   limit;
    double   lfac;
    double   min_bits;
    double   target_bits;
    double   spent;
    uint32_t next;
    bool     whole_first_level;
} level_iter_t;

static void init_level_iter(level_iter_t* const it, const mpz_t n,
                            const double min_bits, const double target_bits,
                            const double budget, const bool whole_first_level) {
    it->limit       = budget * siqs_cost((double) mpz_sizeinbase(n, 2));
    it->lfac        = limb_factor[mpz_size(n)];
    it->min_bits    = min_bits;
    it->target_bits = target_bits;
    it->spent       = 0;
    it->next        = 0;
    it->whole_first_level = whole_first_level;
}

static bool next_level(level_iter_t* const it,
                       const ecm_schedule_t* const schedule,
                       uint32_t* const level) {
    for (uint32_t i = it->next; i < schedule->nlevels; i++) {
        const double hi = (double) levels[i].fbits;
        if (hi <= it->min_bits) {
            continue;
        }
        if (hi - 5 >= it->target_bits) {
            break;
        }
        const bool first = it->whole_first_level && (it->spent == 0)
                           && (levels[i].cost * it->lfac <= it->limit);
        it->spent += levels[i].ncurves * levels[i].cost * it->lfac;
        if ((it->spent > it->limit) && !first) {
            break;
        }
        it->next = i + 1;
        *level   = i;
        return true;
    }
    it->next = schedule->nlevels;
    return false;
}
//------------------------------------------------------------------------------
static bool run_levels(mpz_t g, const mpz_t n,
                       ecm_schedule_t* const schedule,
                       level_iter_t* const it) {
    uint32_t k = 0;
    uint32_t i;
    while (next_level(it, schedule, &i)) {
        const mont_ecm_plan_t* const plan = acquire_plan(schedule, i);
        for (uint32_t j = 0; j < levels[i].ncurves; j++) {
            if (mont_ecm_nth_curve(g, n, k++, plan)) {
                release_plan(schedule, i);
                return true;
            }
            if (tifa_should_stop()) {
                release_plan(schedule, i);
                return false;
            }
        }
        release_plan(schedule, i);
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
    level_iter_t it;
    init_level_iter(&it, n, min_bits, target_bits, budget, false);
    return run_levels(g, n, schedule, &it);
}
//------------------------------------------------------------------------------
static bool run_levels_threaded(mpz_t g, const mpz_t n,
                                ecm_schedule_t* schedule, level_iter_t* it);

bool ecm_schedule_prepass(mpz_t g, const mpz_t n,
                          ecm_schedule_t* const schedule,
                          const double budget, const double max_bits) {
    if (mont_pm1(g, n, 2, schedule->pm1_plan) == 1) {
        return true;
    }
    const double half = (double) mpz_sizeinbase(n, 2) / 2;
    level_iter_t it;
    init_level_iter(&it, n, 0, MIN(half, max_bits), budget, true);
    if (tifa_threads() > 1) {
        return run_levels_threaded(g, n, schedule, &it);
    }
    return run_levels(g, n, schedule, &it);
}
//------------------------------------------------------------------------------
//
// The time limit and the interrupt check are the calling thread's.
//
// At least PREPASS_US_PER_THREAD of estimated curve time per thread: starting
// one takes tens of microseconds.
//
#define PREPASS_US_PER_THREAD 100
#define PREPASS_MAX_THREADS   64

typedef struct {
    pthread_mutex_t lock;
    mpz_srcptr      n;
    ecm_schedule_t* schedule;
    uint32_t        level[NLEVELS];
    uint32_t        first[NLEVELS + 1];
    //
    // The plans this search acquired, by position in level.
    //
    const mont_ecm_plan_t* plans[NLEVELS];
    uint32_t        nlevels;
    uint32_t        next;
    bool            done;
    bool            found;
    mpz_t           g;
} prepass_t;

typedef struct {
    prepass_t* p;
    pthread_t  thread;
    bool       caller;
    bool       running;
} prepass_worker_t;

static const mont_ecm_plan_t* prepass_plan(prepass_t* const p,
                                           const uint32_t r) {
    const mont_ecm_plan_t* plan = __atomic_load_n(&p->plans[r],
                                                  __ATOMIC_ACQUIRE);
    if (plan == NULL) {
        pthread_mutex_lock(&p->lock);
        if (p->plans[r] == NULL) {
            __atomic_store_n(&p->plans[r],
                             acquire_plan(p->schedule, p->level[r]),
                             __ATOMIC_RELEASE);
        }
        plan = p->plans[r];
        pthread_mutex_unlock(&p->lock);
    }
    return plan;
}

static bool prepass_done(void* const arg) {
    return __atomic_load_n(&((prepass_t*) arg)->done, __ATOMIC_RELAXED);
}

static void* run_prepass_curves(void* const arg) {
    const prepass_worker_t* const w = (const prepass_worker_t*) arg;
    prepass_t* const p = w->p;
    const uint32_t ncurves = p->first[p->nlevels];
    uint32_t r = 0;
    mpz_t g;
    mpz_init(g);
    if (!w->caller) {
        //
        // So that the stop checks of the wide curves end this thread's curve
        // when the calling thread stops.
        //
        tifa_set_thread_interrupt(prepass_done, p);
    }

    while (!__atomic_load_n(&p->done, __ATOMIC_RELAXED)) {
        const uint32_t k = __atomic_fetch_add(&p->next, 1, __ATOMIC_RELAXED);
        if (k >= ncurves) {
            break;
        }
        while (k >= p->first[r + 1]) {
            r++;
        }
        const mont_ecm_plan_t* const plan = prepass_plan(p, r);
        if (mont_ecm_nth_curve(g, p->n, k, plan)) {
            pthread_mutex_lock(&p->lock);
            if (!p->found) {
                p->found = true;
                mpz_set(p->g, g);
            }
            pthread_mutex_unlock(&p->lock);
            __atomic_store_n(&p->done, true, __ATOMIC_RELAXED);
        } else if (w->caller && tifa_should_stop()) {
            __atomic_store_n(&p->done, true, __ATOMIC_RELAXED);
        }
    }
    mpz_clear(g);
    return NULL;
}

static bool run_levels_threaded(mpz_t g, const mpz_t n,
                                ecm_schedule_t* const schedule,
                                level_iter_t* const it) {
    prepass_t p;
    p.nlevels  = 0;
    p.first[0] = 0;
    double cost = 0;
    uint32_t i;
    while (next_level(it, schedule, &i)) {
        p.level[p.nlevels] = i;
        p.plans[p.nlevels] = NULL;
        p.first[p.nlevels + 1] = p.first[p.nlevels] + levels[i].ncurves;
        p.nlevels++;
        cost = it->spent;
    }
    const double most = MIN(cost / PREPASS_US_PER_THREAD, PREPASS_MAX_THREADS);
    uint32_t nthreads = MIN(tifa_threads(), (uint32_t) most);
    nthreads = MAX(MIN(nthreads, p.first[p.nlevels]), 1);

    pthread_mutex_init(&p.lock, NULL);
    p.n        = n;
    p.schedule = schedule;
    p.next     = 0;
    p.done     = false;
    p.found    = false;
    mpz_init(p.g);

    prepass_worker_t* const w = tifa_malloc(nthreads * sizeof(prepass_worker_t));
    for (uint32_t t = 0; t < nthreads; t++) {
        w[t].p       = &p;
        w[t].caller  = (t == 0);
        w[t].running = false;
    }
    for (uint32_t t = 1; t < nthreads; t++) {
        w[t].running = (0 == pthread_create(&w[t].thread, NULL,
                                            run_prepass_curves, &w[t]));
    }
    run_prepass_curves(&w[0]);
    for (uint32_t t = 1; t < nthreads; t++) {
        if (w[t].running) {
            pthread_join(w[t].thread, NULL);
        }
    }
    for (uint32_t r = 0; r < p.nlevels; r++) {
        if (p.plans[r] != NULL) {
            release_plan(schedule, p.level[r]);
        }
    }
    if (p.found) {
        mpz_set(g, p.g);
    }
    const bool found = p.found;
    free(w);
    mpz_clear(p.g);
    pthread_mutex_destroy(&p.lock);
    return found;
}

#endif
