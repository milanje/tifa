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

#include <stdlib.h>

#include "tifa_test.h"
#include "prime_cert.h"

static prime_cert_status_t prove_and_check(prime_cert_t* c, const mpz_t n,
                                           int* checked) {
    const prime_cert_status_t s = tifa_prove_prime(c, n);
    *checked = tifa_check_prime_cert(c, n);
    return s;
}

static void expect_proven(prime_cert_t* c, const char* digits) {
    mpz_t n;
    mpz_init_set_str(n, digits, 10);
    int checked;
    char msg[200];
    snprintf(msg, sizeof msg, "%s proven", digits);
    CHECK(prove_and_check(c, n, &checked) == PRIME_PROVEN, msg);
    snprintf(msg, sizeof msg, "%s: certificate checked", digits);
    CHECK(checked, msg);
    mpz_clear(n);
}

static void expect_not_prime(prime_cert_t* c, const char* digits) {
    mpz_t n;
    mpz_init_set_str(n, digits, 10);
    int checked;
    char msg[200];
    snprintf(msg, sizeof msg, "%s not prime", digits);
    CHECK(prove_and_check(c, n, &checked) == NOT_PRIME, msg);
    snprintf(msg, sizeof msg, "%s: certificate refused", digits);
    CHECK(!checked && (c->nsteps == 0), msg);
    mpz_clear(n);
}

int main(void) {
    TIFA_TEST_BEGIN("prime_cert");

    prime_cert_t* c = alloc_prime_cert();
    mpz_t n, m;
    mpz_inits(n, m, NULL);
    int checked;
    char msg[160];

    // ---- Below 2^64: no steps.
    const char* small_primes[] = {
        "2", "3", "5", "37", "41", "2305843009213693951",
        "18446744073709551557"
    };
    for (size_t i = 0; i < sizeof small_primes / sizeof *small_primes; i++) {
        expect_proven(c, small_primes[i]);
        CHECK_EQ_UI(c->nsteps, 0, "no steps below 2^64");
    }
    // Strong pseudoprimes to the bases up to 7 and up to 31.
    const char* small_composites[] = {
        "0", "1", "4", "-7", "18446744073709551615", "3215031751",
        "3825123056546413051"
    };
    for (size_t i = 0; i < sizeof small_composites / sizeof *small_composites;
         i++) {
        expect_not_prime(c, small_composites[i]);
    }

    // ---- Above 2^64.
    expect_proven(c, "18446744073709551629");
    CHECK_EQ_UI(c->nsteps, 1, "2^64 + 13: one step");
    expect_proven(c, "618970019642690137449562111");
    expect_proven(c, "170141183460469231731687303715884105727");
    expect_not_prime(c, "18446744073709551617");
    expect_not_prime(c, "340282366920938463463374607431768211457");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261104UL);
    for (unsigned long bits = 65; bits <= 225; bits += 16) {
        for (int i = 0; i < 3; i++) {
            do {
                mpz_urandomb(n, state, bits);
                mpz_setbit(n, bits - 1);
            } while (!mpz_probab_prime_p(n, 30));
            snprintf(msg, sizeof msg, "random %lu-bit prime proven", bits);
            CHECK(prove_and_check(c, n, &checked) == PRIME_PROVEN, msg);
            snprintf(msg, sizeof msg, "random %lu-bit prime checked", bits);
            CHECK(checked && (c->nsteps >= 1), msg);
            CHECK(mpz_cmp(c->steps[c->nsteps - 1].p, n) == 0,
                  "last step proves n");

            mpz_nextprime(m, n);
            CHECK(!tifa_check_prime_cert(c, m), "not a certificate of another prime");
            mpz_mul_ui(m, n, 3);
            CHECK(!tifa_check_prime_cert(c, m), "not a certificate of 3n");
        }
    }

    // ---- Altered certificates are refused. A 160-bit prime with a step
    //      for a prime above 2^64 is searched for.
    bool found = false;
    while (!found) {
        do {
            mpz_urandomb(n, state, 160);
            mpz_setbit(n, 159);
        } while (!mpz_probab_prime_p(n, 30));
        found = (tifa_prove_prime(c, n) == PRIME_PROVEN) && (c->nsteps >= 2);
    }
    prime_cert_step_t* last = &c->steps[c->nsteps - 1];
    CHECK(tifa_check_prime_cert(c, n), "two-step certificate checked");

    const uint32_t a = last->witnesses->data[0];
    last->witnesses->data[0] = 1;
    CHECK(!tifa_check_prime_cert(c, n), "witness 1 refused");
    last->witnesses->data[0] = 0;
    CHECK(!tifa_check_prime_cert(c, n), "witness 0 refused");
    last->witnesses->data[0] = a;
    CHECK(tifa_check_prime_cert(c, n), "witness restored");

    mpz_add_ui(last->primes->data[0], last->primes->data[0], 2);
    CHECK(!tifa_check_prime_cert(c, n), "altered prime refused");
    mpz_sub_ui(last->primes->data[0], last->primes->data[0], 2);

    last->witnesses->length--;
    CHECK(!tifa_check_prime_cert(c, n), "missing witness refused");
    last->witnesses->length++;

    c->nsteps--;
    CHECK(!tifa_check_prime_cert(c, n), "missing last step refused");
    c->nsteps++;

    prime_cert_step_t first = c->steps[0];
    c->steps[0] = c->steps[c->nsteps - 1];
    c->steps[c->nsteps - 1] = first;
    CHECK(!tifa_check_prime_cert(c, n), "steps out of order refused");
    c->steps[c->nsteps - 1] = c->steps[0];
    c->steps[0] = first;
    CHECK(tifa_check_prime_cert(c, n), "order restored");

    // A step for 2^64 + 1 = 274177 * 67280421310721, from 2 and 3.
    mpz_set_str(n, "18446744073709551617", 10);
    mpz_set(last->p, n);
    mpz_set_ui(last->primes->data[0], 2);
    last->primes->length = 1;
    last->witnesses->data[0] = 3;
    last->witnesses->length = 1;
    CHECK(!tifa_check_prime_cert(c, n), "certificate of 2^64 + 1 refused");

    // ---- Proving again replaces the certificate.
    mpz_set_str(n, "618970019642690137449562111", 10);
    CHECK(tifa_prove_prime(c, n) == PRIME_PROVEN, "2^89 - 1 proven again");
    mpz_set_str(m, "170141183460469231731687303715884105727", 10);
    CHECK(tifa_prove_prime(c, m) == PRIME_PROVEN, "2^127 - 1 proven again");
    CHECK(!tifa_check_prime_cert(c, n) && tifa_check_prime_cert(c, m),
          "certificate replaced");

    gmp_randclear(state);
    mpz_clears(n, m, NULL);
    free_prime_cert(c);

    TIFA_TEST_END();
}
