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


// tifa_set_time_limit and tifa_set_interrupt_check: SIQS and ECM stop soon
// after the limit, or the interrupt, which are the calling thread's alone.

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include <pthread.h>
#include <time.h>

#include "tifa_test.h"
#include "ecm.h"
#include "mont_ecm.h"
#include "tifa_factor.h"
#include "time_limit.h"

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
                             bool* complete) {
    mpz_array_t* f    = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    const double start = now();
    tifa_set_time_limit(0.05);
    const ecode_t e = tifa_factor_with_hint(f, m, n,
                                            FIND_COMPLETE_FACTORIZATION, hint);
    const double elapsed = now() - start;
    *complete = (e == COMPLETE_FACTORIZATION_FOUND);
    free_mpz_array(f);
    free_uint32_array(m);
    return elapsed;
}

static void* limit_seen(void* arg) {
    *(bool*) arg = tifa_time_limit_reached();
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

int main(void) {
    TIFA_TEST_BEGIN("time limit");

    CHECK(!tifa_time_limit_reached(), "no limit at first");

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
    double elapsed = factor_limited(n, TIFA_FACTOR_BALANCED, &complete);
    CHECK(!complete, "SIQS stops at the limit");
    CHECK(elapsed < 1.0, "SIQS stops soon after the limit");
    CHECK(tifa_time_limit_reached(), "the limit has passed");

    // In several threads: the calling one checks the limit for all.
    tifa_set_threads(4);
    random_prime(p, 100, rng);
    random_prime(q, 100, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_BALANCED, &complete);
    CHECK(!complete, "SIQS in 4 threads stops at the limit");
    CHECK(elapsed < 1.0, "SIQS in 4 threads stops soon after the limit");
    tifa_set_threads(1);

    bool other = true;
    pthread_t thread;
    CHECK(pthread_create(&thread, NULL, limit_seen, &other) == 0,
          "thread started");
    pthread_join(thread, NULL);
    CHECK(!other, "the limit is the thread's alone");

    // A 64-bit factor of a 240-bit number: ECM searches for long.
    random_prime(p, 64, rng);
    random_prime(q, 176, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_UNBALANCED, &complete);
    CHECK(!complete, "ECM stops at the limit");
    CHECK(elapsed < 1.0, "ECM stops soon after the limit");

    tifa_set_threads(4);
    random_prime(p, 80, rng);
    random_prime(q, 160, rng);
    mpz_mul(n, p, q);
    elapsed = factor_limited(n, TIFA_FACTOR_UNBALANCED, &complete);
    CHECK(!complete, "ECM in 4 threads stops at the limit");
    CHECK(elapsed < 1.0, "ECM in 4 threads stops soon after the limit");
    tifa_set_threads(1);

    // The longest curves stop within themselves: in GMP arithmetic above
    // 512 bits (2.5 s each at 600 bits on an Apple M5 Pro), and on eight
    // limbs (0.7 s).
    random_prime(p, 300, rng);
    random_prime(q, 300, rng);
    mpz_mul(n, p, q);
    ecm_params_t ecm_params = { 576000, 50 * 576000, 4, 1 };
    mpz_array_t* found = alloc_mpz_array(4);
    double start = now();
    tifa_set_time_limit(0.05);
    ecm(found, NULL, n, &ecm_params, SINGLE_RUN);
    elapsed = now() - start;
    CHECK(elapsed < 0.5, "a GMP curve stops soon after the limit");
#if TIFA_HAVE_MONT_ECM
    random_prime(p, 250, rng);
    random_prime(q, 250, rng);
    mpz_mul(n, p, q);
    mont_ecm_plan_t* const plan = alloc_mont_ecm_plan(819200, 50 * 819200);
    start = now();
    tifa_set_time_limit(0.05);
    mont_ecm_curve(p, n, 6, plan);
    elapsed = now() - start;
    CHECK(elapsed < 0.3, "an eight-limb curve stops soon after the limit");
    free_mont_ecm_plan(plan);
#endif
    free_mpz_array(found);

    tifa_set_time_limit(0.0);
    CHECK(!tifa_time_limit_reached(), "the limit is removed");

    // Another thread asks SIQS to stop, through a flag.
    random_prime(p, 100, rng);
    random_prime(q, 100, rng);
    mpz_mul(n, p, q);
    unsigned ncalls = 0;
    tifa_set_interrupt_check(read_flag, &ncalls);
    CHECK(pthread_create(&thread, NULL, raise_flag, NULL) == 0,
          "thread started");
    mpz_array_t* g    = alloc_mpz_array(4);
    uint32_array_t* h = alloc_uint32_array(4);
    start = now();
    const ecode_t e = tifa_factor_with_hint(g, h, n,
                                            FIND_COMPLETE_FACTORIZATION,
                                            TIFA_FACTOR_BALANCED);
    elapsed = now() - start;
    pthread_join(thread, NULL);
    CHECK(e != COMPLETE_FACTORIZATION_FOUND, "SIQS stops when interrupted");
    CHECK(elapsed < 1.0, "SIQS stops soon after the interrupt");
    CHECK(tifa_interrupted(), "interrupted");
    CHECK(!tifa_time_limit_reached(), "not timed out");
    CHECK(ncalls > 1, "the check is called again and again");
    tifa_set_interrupt_check(NULL, NULL);
    CHECK(!tifa_interrupted(), "no longer interrupted");
    free_mpz_array(g);
    free_uint32_array(h);
    mpz_set_str(n, "19409887069306060631113", 10);
    mpz_array_t* f    = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    CHECK(tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION)
          == COMPLETE_FACTORIZATION_FOUND, "factors without a limit");
    free_mpz_array(f);
    free_uint32_array(m);

    mpz_clears(p, q, n, NULL);
    gmp_randclear(rng);

    TIFA_TEST_END();
}
