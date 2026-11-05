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


// tifa_set_time_limit: SIQS and ECM stop soon after the limit, which is the
// calling thread's alone.

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include <pthread.h>
#include <time.h>

#include "tifa_test.h"
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

    tifa_set_time_limit(0.0);
    CHECK(!tifa_time_limit_reached(), "the limit is removed");
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
