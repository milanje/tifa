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


// The time limit and the interrupt check of the options: SIQS and ECM stop
// soon after them, which hold for the call alone.

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include <pthread.h>
#include <time.h>

#include "tifa_test.h"
#include "ecm.h"
#include "mont_ecm.h"
#include "tifa_factor.h"
#include "tifa_options.h"
#include "stopping.h"

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

static void random_prime(mpz_t p, uint32_t bits, gmp_randstate_t rng) {
    mpz_urandomb(p, rng, bits);
    mpz_setbit(p, bits - 1);
    mpz_nextprime(p, p);
}

//
// Factors n under a limit of 50 ms; returns the time taken, and whether the
// factorization was complete in *complete.
//
static double factor_limited(const mpz_t n, tifa_factor_hint_t hint,
                             uint32_t threads, bool* complete) {
    mpz_array_t* f    = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    const tifa_options_t options = { .time_limit = 0.05, .hint = hint,
                                     .threads = threads };
    const double start = now();
    const ecode_t e = tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION,
                                  &options);
    const double elapsed = now() - start;
    *complete = (e == COMPLETE_FACTORIZATION_FOUND);
    free_mpz_array(f);
    free_uint32_array(m);
    return elapsed;
}

static void* not_stopped(void* arg) {
    *(bool*) arg = tifa_should_stop();
    return NULL;
}

static volatile bool stop_flag = false;

static bool read_flag(void* arg) {
    (*(unsigned*) arg)++;
    return stop_flag;
}

static void* raise_flag(void* arg) {
    (void) arg;
    const struct timespec wait = { 0, 50 * 1000 * 1000 };
    nanosleep(&wait, NULL);
    stop_flag = true;
    return NULL;
}

//
// An interrupt check that factors a number of its own, with options of its
// own, the first time it is called.
//
static bool nested_call(void* arg) {
    unsigned* const ncalls = arg;
    if ((*ncalls)++ == 0) {
        mpz_t k;
        mpz_init_set_ui(k, 15);
        mpz_array_t* f    = alloc_mpz_array(4);
        uint32_array_t* m = alloc_uint32_array(4);
        const tifa_options_t own = { .threads = 1 };
        tifa_factor(f, m, k, FIND_COMPLETE_FACTORIZATION, &own);
        free_mpz_array(f);
        free_uint32_array(m);
        mpz_clear(k);
    }
    return false;
}

int main(void) {
    TIFA_TEST_BEGIN("time limit");

    CHECK(!tifa_should_stop(), "no limit at first");

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 1);
    mpz_t p, q, n;
    mpz_inits(p, q, n, NULL);

    // A 200-bit semiprime: about a second of SIQS.
    random_prime(p, 100, rng);
    random_prime(q, 100, rng);
    mpz_mul(n, p, q);
    bool complete;
    double elapsed = factor_limited(n, TIFA_FACTOR_BALANCED, 1, &complete);
    CHECK(!complete, "SIQS stops at the limit");
    CHECK(elapsed < 1.0, "SIQS stops soon after the limit");
    CHECK(!tifa_should_stop(), "the limit ends with the call");

    // In several threads: the calling one checks the limit for all.
    random_prime(p, 100, rng);
    random_prime(q, 100, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_BALANCED, 4, &complete);
    CHECK(!complete, "SIQS in 4 threads stops at the limit");
    CHECK(elapsed < 1.0, "SIQS in 4 threads stops soon after the limit");

    // The options of a call hold in its thread alone.
    tifa_call_t call;
    const tifa_options_t passed = { .time_limit = 1e-9 };
    tifa_call_begin(&call, &passed);
    while (!tifa_should_stop()) {
    }
    bool other = true;
    pthread_t thread;
    CHECK(pthread_create(&thread, NULL, not_stopped, &other) == 0,
          "thread started");
    pthread_join(thread, NULL);
    CHECK(!other, "the limit is the call's thread's alone");
    tifa_call_end(&call);
    CHECK(!tifa_should_stop(), "the limit ends with the call");

    // A 64-bit factor of a 240-bit number: ECM searches for long.
    random_prime(p, 64, rng);
    random_prime(q, 176, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_UNBALANCED, 1, &complete);
    CHECK(!complete, "ECM stops at the limit");
    CHECK(elapsed < 1.0, "ECM stops soon after the limit");

    random_prime(p, 80, rng);
    random_prime(q, 160, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_UNBALANCED, 4, &complete);
    CHECK(!complete, "ECM in 4 threads stops at the limit");
    CHECK(elapsed < 1.0, "ECM in 4 threads stops soon after the limit");

#if TIFA_HAVE_MONT_ECM
    // On eight limbs, the other threads' curves stop within themselves, and
    // SIQS does not start (0.2 s to start four threads at 500 bits).
    random_prime(p, 250, rng);
    random_prime(q, 250, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_RANDOM, 4, &complete);
    CHECK(!complete, "tifa_factor in 4 threads stops at the limit");
    CHECK(elapsed < 0.2, "tifa_factor in 4 threads stops soon at 500 bits");
#endif

    // The longest curves stop within themselves: in GMP arithmetic above
    // 512 bits (2.5 s each at 600 bits on an Apple M5 Pro), and on eight
    // limbs (0.7 s).
    random_prime(p, 300, rng);
    random_prime(q, 300, rng);
    mpz_mul(n, p, q);
    ecm_params_t ecm_params = { 576000, 50 * 576000, 4, 1 };
    const tifa_options_t short_limit = { .time_limit = 0.05 };
    mpz_array_t* found = alloc_mpz_array(4);
    double start = now();
    ecm(found, NULL, n, &ecm_params, SINGLE_RUN, &short_limit);
    elapsed = now() - start;
    CHECK(elapsed < 0.5, "a GMP curve stops soon after the limit");
#if TIFA_HAVE_MONT_ECM
    random_prime(p, 250, rng);
    random_prime(q, 250, rng);
    mpz_mul(n, p, q);
    mont_ecm_plan_t* const plan = alloc_mont_ecm_plan(819200, 50 * 819200);
    start = now();
    tifa_call_begin(&call, &short_limit);
    mont_ecm_curve(p, n, 6, plan);
    tifa_call_end(&call);
    elapsed = now() - start;
    CHECK(elapsed < 0.3, "an eight-limb curve stops soon after the limit");
    free_mont_ecm_plan(plan);
#endif
    free_mpz_array(found);

    // Another thread asks SIQS to stop, through a flag.
    random_prime(p, 100, rng);
    random_prime(q, 100, rng);
    mpz_mul(n, p, q);
    unsigned ncalls = 0;
    const tifa_options_t flagged = { .interrupt     = read_flag,
                                     .interrupt_arg = &ncalls,
                                     .hint          = TIFA_FACTOR_BALANCED };
    CHECK(pthread_create(&thread, NULL, raise_flag, NULL) == 0,
          "thread started");
    mpz_array_t* g    = alloc_mpz_array(4);
    uint32_array_t* h = alloc_uint32_array(4);
    start = now();
    ecode_t e = tifa_factor(g, h, n, FIND_COMPLETE_FACTORIZATION, &flagged);
    elapsed = now() - start;
    pthread_join(thread, NULL);
    CHECK(e != COMPLETE_FACTORIZATION_FOUND, "SIQS stops when interrupted");
    CHECK(elapsed < 1.0, "SIQS stops soon after the interrupt");
    CHECK(ncalls > 1, "the check is called again and again");
    CHECK(!tifa_should_stop(), "the interrupt ends with the call");

    // A call made from the interrupt check, with options of its own, gives
    // the outer call back its time limit.
    random_prime(p, 100, rng);
    random_prime(q, 100, rng);
    mpz_mul(n, p, q);
    ncalls = 0;
    const tifa_options_t nesting = { .time_limit    = 0.05,
                                     .interrupt     = nested_call,
                                     .interrupt_arg = &ncalls,
                                     .hint          = TIFA_FACTOR_BALANCED };
    g->length = 0;
    h->length = 0;
    start = now();
    e = tifa_factor(g, h, n, FIND_COMPLETE_FACTORIZATION, &nesting);
    elapsed = now() - start;
    CHECK(e != COMPLETE_FACTORIZATION_FOUND, "the outer call stops");
    CHECK(elapsed < 1.0,
          "the outer call stops at its limit after a nested call");
    CHECK(ncalls > 1, "its check holds after the nested call");
    free_mpz_array(g);
    free_uint32_array(h);

    // Stopped early, a composite is never said to have no factor.
    static const char* const composites[] = {
        "78664248049168147517", "1552322315513162798881",
        "23981845589681487731813794469"
    };
    unsigned nfound = 0;
    for (unsigned r = 0; r < 100; r++) {
        for (unsigned i = 0; i < 3; i++) {
            mpz_set_str(n, composites[i], 10);
            mpz_array_t* cf    = alloc_mpz_array(4);
            uint32_array_t* cm = alloc_uint32_array(4);
            const tifa_options_t tiny = { .time_limit = 1e-6 * (20 + 6 * r) };
            nfound += (tifa_factor(cf, cm, n, FIND_COMPLETE_FACTORIZATION,
                                   &tiny) == NO_FACTOR_FOUND);
            free_mpz_array(cf);
            free_uint32_array(cm);
        }
    }
    CHECK(nfound == 0, "a composite stopped early: never NO_FACTOR_FOUND");

    mpz_set_str(n, "19409887069306060631113", 10);
    mpz_array_t* f    = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    CHECK(tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION, NULL)
          == COMPLETE_FACTORIZATION_FOUND, "factors without a limit");
    free_mpz_array(f);
    free_uint32_array(m);

    mpz_clears(p, q, n, NULL);
    gmp_randclear(rng);

    TIFA_TEST_END();
}
