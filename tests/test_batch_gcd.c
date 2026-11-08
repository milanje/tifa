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

// tifa_batch_gcd (batch_gcd.h) against direct gcds: random moduli of 64 to
// 1024 bits, some sharing a prime or repeated, in sets of 1 to 300.

#include <stdlib.h>

#include "tifa_test.h"
#include "batch_gcd.h"
#include "tifa_options.h"

static void random_prime(mpz_t p, gmp_randstate_t state, unsigned long bits) {
    mpz_urandomb(p, state, bits);
    mpz_setbit(p, bits - 1);
    mpz_nextprime(p, p);
}

// Whether gcds->data[first + i] is gcd(n_i, prod_{j != i} n_j) for all i,
// the products from prefix and suffix products.
static bool matches_direct(const mpz_array_t* const gcds, const uint32_t first,
                           const mpz_array_t* const moduli) {
    const uint32_t count = moduli->length;
    if (gcds->length != first + count) {
        return false;
    }
    mpz_t* const pre = malloc((count + 1) * sizeof(mpz_t));
    mpz_t* const suf = malloc((count + 1) * sizeof(mpz_t));
    mpz_init_set_ui(pre[0], 1);
    mpz_init_set_ui(suf[count], 1);
    for (uint32_t i = 0; i < count; i++) {
        mpz_init(pre[i + 1]);
        mpz_mul(pre[i + 1], pre[i], moduli->data[i]);
    }
    for (uint32_t i = count; i > 0; i--) {
        mpz_init(suf[i - 1]);
        mpz_mul(suf[i - 1], suf[i], moduli->data[i - 1]);
    }
    mpz_t g;
    mpz_init(g);
    bool ok = true;
    for (uint32_t i = 0; ok && (i < count); i++) {
        mpz_mul(g, pre[i], suf[i + 1]);
        mpz_gcd(g, g, moduli->data[i]);
        ok = (mpz_cmp(g, gcds->data[first + i]) == 0);
    }
    for (uint32_t i = 0; i <= count; i++) {
        mpz_clears(pre[i], suf[i], NULL);
    }
    free(pre);
    free(suf);
    mpz_clear(g);
    return ok;
}

int main(void) {
    TIFA_TEST_BEGIN("batch_gcd");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261106UL);

    mpz_t p, q;
    mpz_inits(p, q, NULL);
    mpz_array_t* const moduli = alloc_mpz_array(4);
    mpz_array_t* const gcds   = alloc_mpz_array(4);
    char msg[160];

    CHECK(tifa_batch_gcd(gcds, moduli, NULL) == SUCCESS && (gcds->length == 0),
          "no moduli, no gcd");

    mpz_set_ui(p, 15);
    append_mpz_to_array(moduli, p);
    CHECK(   (tifa_batch_gcd(gcds, moduli, NULL) == SUCCESS)
          && (gcds->length == 1) && (mpz_cmp_ui(gcds->data[0], 1) == 0),
          "a single modulus: 1");

    reset_mpz_array(gcds);
    static const unsigned long small[] = { 15, 21, 22, 1, 35, 22 };
    reset_mpz_array(moduli);
    for (uint32_t i = 0; i < 6; i++) {
        mpz_set_ui(p, small[i]);
        append_mpz_to_array(moduli, p);
    }
    CHECK(   (tifa_batch_gcd(gcds, moduli, NULL) == SUCCESS)
          && matches_direct(gcds, 0, moduli),
          "15, 21, 22, 1, 35, 22: shared primes, 1 and a duplicate");

    static const uint32_t counts[] = { 2, 3, 7, 64, 65, 300 };
    static const uint32_t bits[]   = { 64, 128, 512, 1024, 256, 512 };
    for (uint32_t t = 0; t < 6; t++) {
        //
        // Every fifth modulus shares a prime with the one before, and the
        // last repeats the first.
        //
        reset_mpz_array(moduli);
        for (uint32_t i = 0; i < counts[t]; i++) {
            if (i % 5 != 1) {
                random_prime(p, state, bits[t] / 2);
            }
            random_prime(q, state, bits[t] / 2);
            mpz_mul(q, p, q);
            append_mpz_to_array(moduli, q);
        }
        if (counts[t] > 2) {
            mpz_set(moduli->data[counts[t] - 1], moduli->data[0]);
        }
        reset_mpz_array(gcds);
        mpz_set_ui(p, 7);
        append_mpz_to_array(gcds, p);
        snprintf(msg, sizeof(msg), "%u moduli of %u bits match direct "
                 "gcds, appended", counts[t], bits[t]);
        CHECK(   (tifa_batch_gcd(gcds, moduli, NULL) == SUCCESS)
              && (mpz_cmp_ui(gcds->data[0], 7) == 0)
              && matches_direct(gcds, 1, moduli), msg);
    }

    reset_mpz_array(gcds);
    mpz_set_si(p, -15);
    append_mpz_to_array(moduli, p);
    CHECK(   (tifa_batch_gcd(gcds, moduli, NULL) == FAILURE)
          && (gcds->length == 0), "a negative modulus fails");
    mpz_set_ui(moduli->data[moduli->length - 1], 0);
    CHECK(   (tifa_batch_gcd(gcds, moduli, NULL) == FAILURE)
          && (gcds->length == 0), "a zero modulus fails");

    const tifa_options_t limited = { .time_limit = 1e-9 };
    moduli->length--;
    CHECK(   (tifa_batch_gcd(gcds, moduli, &limited) == GIVING_UP)
          && (gcds->length == 0), "the time limit gives up, appends nothing");
    const tifa_options_t bad = { .reserved = { 1 } };
    CHECK(   (tifa_batch_gcd(gcds, moduli, &bad) == FAILURE)
          && (gcds->length == 0), "options with a reserved field set fail");

    free_mpz_array(moduli);
    free_mpz_array(gcds);
    mpz_clears(p, q, NULL);
    gmp_randclear(state);

    TIFA_TEST_END();
}
