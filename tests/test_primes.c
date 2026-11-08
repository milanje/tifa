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

// The prime iterator (primes.h) against mpz_nextprime, across the end of the
// table of first_primes.h, the segments of the sieve and 2^32 - 1, and the
// stage-1 exponent against a direct product.

#include "tifa_test.h"
#include "first_primes.h"
#include "primes.h"

// 0 if the primes of [lo, hi] are those of mpz_nextprime, else the first
// prime that differs (or UINT32_MAX for a missing one).
static uint64_t mismatch(const uint32_t lo, const uint32_t hi) {
    mpz_t x;
    mpz_init_set_ui(x, lo);
    if (lo > 0) {
        mpz_sub_ui(x, x, 1);
    }
    prime_iter_t it;
    prime_iter_init(&it, lo, hi);
    uint64_t bad = 0;
    for (;;) {
        mpz_nextprime(x, x);
        const uint64_t want = (mpz_cmp_ui(x, hi) <= 0) ? mpz_get_ui(x) : 0;
        const uint32_t got  = prime_iter_next(&it);
        if (got != want) {
            bad = (got != 0) ? got : UINT32_MAX;
            break;
        }
        if (got == 0) {
            bad = (prime_iter_next(&it) != 0) ? UINT32_MAX : 0;
            break;
        }
    }
    prime_iter_clear(&it);
    mpz_clear(x);
    return bad;
}

int main(void) {
    TIFA_TEST_BEGIN("primes");

    const uint32_t lp  = LARGEST_PRIME;
    const uint32_t top = UINT32_MAX;
    const struct { uint32_t lo, hi; } windows[] = {
        { 0, 2000 }, { 2, 2 }, { 4, 4 }, { 10, 9 }, { 0, 0 },
        { lp - 2000, lp }, { lp, lp }, { lp + 1, lp + 1 },
        { lp - 1000, lp + 300000 }, { lp + 2, lp + 140000 },
        { 1000000000, 1000400000 }, { top - 5000, top }, { top, top },
        { top - 4, top - 4 },
    };
    char msg[120];
    for (uint32_t i = 0; i < sizeof(windows) / sizeof(windows[0]); i++) {
        snprintf(msg, sizeof(msg), "primes of [%u, %u]", windows[i].lo,
                 windows[i].hi);
        CHECK_EQ_UI(mismatch(windows[i].lo, windows[i].hi), 0, msg);
    }
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261108UL);
    uint32_t nbad = 0;
    for (uint32_t t = 0; t < 200; t++) {
        const uint32_t lo = (uint32_t) gmp_urandomm_ui(state, top);
        const uint32_t w  = (uint32_t) gmp_urandomm_ui(state, 5000);
        nbad += (mismatch(lo, (lo > top - w) ? top : lo + w) != 0);
    }
    CHECK_EQ_UI(nbad, 0, "primes of random windows");

    mpz_t e, ref, x;
    mpz_inits(e, ref, x, NULL);
    const uint32_t bounds[] = { 0, 1, 2, 3, 8, 9, 50, 125, lp, lp + 1,
                                823543, 900000 };
    for (uint32_t i = 0; i < sizeof(bounds) / sizeof(bounds[0]); i++) {
        const uint32_t b = bounds[i];
        mpz_set_ui(ref, 1);
        mpz_set_ui(x, 1);
        for (mpz_nextprime(x, x); mpz_cmp_ui(x, b) <= 0; mpz_nextprime(x, x)) {
            const unsigned long q = mpz_get_ui(x);
            unsigned long qe = q;
            while (qe <= b / q) {
                qe *= q;
            }
            mpz_mul_ui(ref, ref, qe);
        }
        prime_powers_product(e, b);
        snprintf(msg, sizeof(msg), "prime powers up to %u", b);
        CHECK(mpz_cmp(e, ref) == 0, msg);
    }
    mpz_clears(e, ref, x, NULL);
    gmp_randclear(state);

    TIFA_TEST_END();
}
