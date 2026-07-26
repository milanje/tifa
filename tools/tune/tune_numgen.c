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

/**
 * \file    tune_numgen.c
 */

#include <stdbool.h>

#include "tune_numgen.h"

//
// Draws a random prime of exactly `bits` bits into `out`. The top bit is forced
// so the pre-nextprime candidate already spans the range, and the result is
// retried on the rare occasion mpz_nextprime carries it past `bits` bits.
//
static void random_prime_of_size(mpz_t out, unsigned int bits,
                                  gmp_randstate_t rng) {
    do {
        mpz_urandomb(out, rng, bits);
        mpz_setbit(out, bits - 1);      // force exactly `bits` significant bits
        mpz_nextprime(out, out);
    } while (mpz_sizeinbase(out, 2) != bits);
}

//------------------------------------------------------------------------------
void tune_gen_composite(mpz_t n, mpz_t smallest,
                        unsigned int nfactors, unsigned int n_bits,
                        gmp_randstate_t rng) {
    if (nfactors < 1) {
        nfactors = 1;
    }
    if (nfactors > TUNE_MAX_NFACTORS) {
        nfactors = TUNE_MAX_NFACTORS;
    }

    // Split n_bits as evenly as possible across the factors; the first `extra`
    // factors take one additional bit.
    unsigned int base  = n_bits / nfactors;
    unsigned int extra = n_bits % nfactors;
    if (base < 2) {
        base = 2;
    }

    mpz_t primes[TUNE_MAX_NFACTORS];
    for (unsigned int i = 0; i < nfactors; i++) {
        mpz_init(primes[i]);
    }

    bool done = false;
    while (!done) {
        mpz_set_ui(n, 1UL);
        for (unsigned int i = 0; i < nfactors; i++) {
            unsigned int b = base + ((i < extra) ? 1u : 0u);
            random_prime_of_size(primes[i], b, rng);
            mpz_mul(n, n, primes[i]);
        }
        if (mpz_sizeinbase(n, 2) != n_bits) {
            continue;   // product missed the target bit-size; retry
        }
        // Reject repeated primes (which would give n a repeated factor).
        bool distinct = true;
        for (unsigned int i = 0; i < nfactors && distinct; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (mpz_cmp(primes[i], primes[j]) == 0) {
                    distinct = false;
                    break;
                }
            }
        }
        done = distinct;
    }

    mpz_set(smallest, primes[0]);
    for (unsigned int i = 1; i < nfactors; i++) {
        if (mpz_cmp(primes[i], smallest) < 0) {
            mpz_set(smallest, primes[i]);
        }
    }

    for (unsigned int i = 0; i < nfactors; i++) {
        mpz_clear(primes[i]);
    }
}
