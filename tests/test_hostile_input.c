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


// Inputs from the security review: each crashed, hung or ran for minutes.

#include <pthread.h>
#include <stdlib.h>

#include "tifa_test.h"
#include "factor_u64.h"
#include "factoring_machine.h"
#include "fermat.h"
#include "rho.h"
#include "siqs.h"
#include "tifa_factor.h"

static void* factor_power(void* arg) {
    mpz_t n;
    mpz_init(n);
    mpz_ui_pow_ui(n, 1009, 5001);
    mpz_mul_ui(n, n, 1013);
    mpz_array_t* f    = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    const ecode_t e = tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION);
    *(bool*) arg = (e == COMPLETE_FACTORIZATION_FOUND) && (f->length == 2)
                && (mpz_cmp_ui(f->data[0], 1009) == 0) && (m->data[0] == 5001)
                && (mpz_cmp_ui(f->data[1], 1013) == 0) && (m->data[1] == 1);
    free_mpz_array(f);
    free_uint32_array(m);
    mpz_clear(n);
    return NULL;
}

static ecode_t found(const char* const* digits, uint32_t count, uint32_t nn,
                     factoring_mode_t mode, mpz_array_t* f) {
    reset_mpz_array(f);
    mpz_t x;
    mpz_init(x);
    for (uint32_t i = 0; i < count; i++) {
        mpz_set_str(x, digits[i], 10);
        append_mpz_to_array(f, x);
    }
    mpz_set_ui(x, nn);
    uint32_array_t* m = alloc_uint32_array(4);
    const ecode_t e = process_found_factors(f, m, x, mode, tifa_factor);
    free_uint32_array(m);
    mpz_clear(x);
    return e;
}

int main(void) {
    TIFA_TEST_BEGIN("hostile input");

    // 1009^5001 * 1013 in a thread with 256 KiB of stack: the coprime base
    // recursed once per power of 1009.
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, 256 * 1024);
    pthread_t thread;
    bool power_ok = false;
    CHECK(pthread_create(&thread, &attr, factor_power, &power_ok) == 0,
          "thread started");
    pthread_join(thread, NULL);
    pthread_attr_destroy(&attr);
    CHECK(power_ok, "1009^5001 * 1013 factored in a small stack");

    // Factors found by a caller: 0 was divided by, -3 brought -1 into the
    // coprime base (an endless loop), a non-divisor gave garbage.
    mpz_array_t* f = alloc_mpz_array(4);
    const char* zero[]  = { "0" };
    const char* neg[]   = { "-3", "3" };
    const char* trivial[] = { "1", "15" };
    const char* other[] = { "7" };
    CHECK(found(zero, 1, 15, FIND_SOME_COPRIME_FACTORS, f) == NO_FACTOR_FOUND,
          "factor 0 dropped");
    CHECK((found(neg, 2, 45, FIND_COMPLETE_FACTORIZATION, f)
           == COMPLETE_FACTORIZATION_FOUND) && (f->length == 2),
          "factors -3 and 3 of 45");
    CHECK(found(trivial, 2, 15, FIND_SOME_COPRIME_FACTORS, f)
          == NO_FACTOR_FOUND, "factors 1 and n dropped");
    CHECK(found(other, 1, 15, FIND_SOME_COPRIME_FACTORS, f) == FAILURE,
          "a non-divisor refused");
    free_mpz_array(f);

    // Searches past the end of a sorted array.
    uint32_array_t* sorted = alloc_uint32_array(4);
    for (uint32_t i = 1; i <= 4; i++) {
        append_uint32_to_array(sorted, 10 * i);
    }
    CHECK(index_in_sorted_uint32_array(25, sorted, 0, 100) == NOT_IN_ARRAY,
          "25 not found, max_index beyond the array");
    CHECK_EQ_UI(index_in_sorted_uint32_array(40, sorted, 0, 100), 3,
                "40 found, max_index beyond the array");
    CHECK(index_in_sorted_uint32_array(10, sorted, 3, 1) == NOT_IN_ARRAY,
          "min_index above max_index");
    free_uint32_array(sorted);

#if TIFA_HAVE_FACTOR_U64
    // split_u64 on what is not composite (32 s on 2^64 - 59).
    CHECK((split_u64(18446744073709551557ULL) == 0) && (split_u64(0) == 0)
          && (split_u64(1) == 0) && (split_u64(3) == 0),
          "no split of 0, 1, primes");
    const uint64_t s = split_u64(15);
    CHECK((s == 3) || (s == 5), "15 split");
#endif

    mpz_array_t*    fa = alloc_mpz_array(4);
    uint32_array_t* mu = alloc_uint32_array(4);
    mpz_t n;
    mpz_init(n);

    // Fermat on a prime went through 2^24 primes (45 s for 3).
    fermat_params_t fp;
    set_fermat_params_to_default(&fp);
    mpz_set_ui(n, 2305843009213693951ULL);
    CHECK(fermat(fa, mu, n, &fp, FIND_COMPLETE_FACTORIZATION) == NO_FACTOR_FOUND,
          "fermat on 2^61 - 1");

    // rho with more accumulated than allowed iterations (258 s).
    rho_params_t rp;
    mpz_set_ui(n, 15);
    set_rho_params_to_default(n, &rp);
    rp.acc_niters = UINT32_MAX;
    rp.max_niters = 1000;
    reset_mpz_array(fa);
    reset_uint32_array(mu);
    CHECK(rho(fa, mu, n, &rp, FIND_COMPLETE_FACTORIZATION)
          == COMPLETE_FACTORIZATION_FOUND, "rho on 15, acc_niters = UINT32_MAX");

    // SIQS with the smallest factor base, above 140 bits: the candidates
    // were all dropped, but the fills not counted as idle (no end).
    mpz_t p, q;
    mpz_init(p);
    mpz_init(q);
    mpz_ui_pow_ui(p, 2, 89);
    mpz_nextprime(p, p);
    mpz_ui_pow_ui(q, 2, 72);
    mpz_nextprime(q, q);
    mpz_mul(n, p, q);
    siqs_params_t sp;
    set_siqs_params_to_default(n, &sp);
    sp.nprimes_in_base = 0;
    reset_mpz_array(fa);
    reset_uint32_array(mu);
    const ecode_t e = siqs(fa, mu, n, &sp, FIND_COMPLETE_FACTORIZATION);
    CHECK((e == GIVING_UP) || (e == COMPLETE_FACTORIZATION_FOUND),
          "siqs with nprimes_in_base = 0 returns");

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n);
    free_mpz_array(fa);
    free_uint32_array(mu);

    TIFA_TEST_END();
}
