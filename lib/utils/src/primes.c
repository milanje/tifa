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
 * \file    primes.c
 *
 * \brief The primes up to 2^32 - 1, and the stage-1 exponent over them.
 */

#include <stdlib.h>
#include <string.h>

#include "primes.h"
#include "first_primes.h"
#include "tifa_alloc.h"
#include "stopping.h"

//
// Odd numbers per segment of the sieve beyond the table.
//
#define SEG (1U << 16)

void prime_iter_init(prime_iter_t* const it, const uint32_t lo,
                     const uint32_t hi) {
    it->hi = hi;
    it->ip = 0;
    while ((it->ip < NFIRST_PRIMES) && (first_primes[it->ip] < lo)) {
        it->ip++;
    }
    const uint64_t first = (uint64_t) LARGEST_PRIME + 1;
    it->start     = ((lo > first) ? lo : first) | 1;
    it->seg_lo    = 0;
    it->pos       = 0;
    it->len       = 0;
    it->done      = false;
    it->composite = NULL;
}

static bool sieve_next_segment(prime_iter_t* const it) {
    const uint64_t lo = (it->composite == NULL)
                      ? it->start : it->seg_lo + 2 * (uint64_t) SEG;
    if (lo > it->hi) {
        return false;
    }
    uint32_t n = SEG;
    if ((it->hi - lo) / 2 + 1 < SEG) {
        n = (uint32_t) ((it->hi - lo) / 2 + 1);
    }
    if (it->composite == NULL) {
        it->composite = tifa_malloc(SEG);
    }
    uint8_t* const composite = it->composite;
    memset(composite, 0, n);
    const uint64_t last = lo + 2 * (uint64_t) (n - 1);
    //
    // lo is beyond the table, so no multiple of p marked is p itself.
    //
    for (uint32_t i = 1; i < NFIRST_PRIMES; i++) {
        const uint64_t p = first_primes[i];
        if (p * p > last) {
            break;
        }
        uint64_t x = ((lo + p - 1) / p) * p;
        if ((x & 1) == 0) {
            x += p;
        }
        for (uint64_t j = (x - lo) / 2; j < n; j += p) {
            composite[j] = 1;
        }
    }
    it->seg_lo = lo;
    it->pos    = 0;
    it->len    = n;
    return true;
}

uint32_t prime_iter_next(prime_iter_t* const it) {
    if (it->done) {
        return 0;
    }
    if (it->ip < NFIRST_PRIMES) {
        const uint32_t p = first_primes[it->ip];
        if (p <= it->hi) {
            it->ip++;
            return p;
        }
        it->done = true;
        return 0;
    }
    for (;;) {
        while (it->pos < it->len) {
            const uint32_t i = it->pos++;
            if (!it->composite[i]) {
                return (uint32_t) (it->seg_lo + 2 * (uint64_t) i);
            }
        }
        if (!sieve_next_segment(it)) {
            it->done = true;
            return 0;
        }
    }
}

void prime_iter_clear(prime_iter_t* const it) {
    free(it->composite);
    it->composite = NULL;
}

//
// e <- scalars[lo] * ... * scalars[hi - 1], as a balanced product tree.
//
static bool product_of_scalars(mpz_t e, const uint64_t* const scalars,
                               const uint32_t lo, const uint32_t hi,
                               const bool stoppable) {
    if (hi - lo == 1) {
        mpz_import(e, 1, -1, sizeof(uint64_t), 0, 0, &scalars[lo]);
        return true;
    }
    if (stoppable && (hi - lo >= 256) && tifa_should_stop()) {
        return false;
    }
    mpz_t right;
    mpz_init(right);
    const bool done
        =    product_of_scalars(e, scalars, lo, lo + (hi - lo) / 2, stoppable)
          && product_of_scalars(right, scalars, lo + (hi - lo) / 2, hi,
                                stoppable);
    if (done) {
        mpz_mul(e, e, right);
    }
    mpz_clear(right);
    return done;
}

typedef struct {
    uint64_t* data;
    uint32_t  length;
    uint32_t  alloced;
} u64_vec_t;

static void push_u64(u64_vec_t* const v, const uint64_t x) {
    if (v->length == v->alloced) {
        v->alloced = (v->alloced == 0) ? 64 : 2 * v->alloced;
        v->data    = tifa_realloc(v->data, v->alloced * sizeof(uint64_t));
    }
    v->data[v->length++] = x;
}

bool prime_powers_product(mpz_t e, const uint32_t b, const bool stoppable) {
    u64_vec_t scalars = { NULL, 0, 0 };
    uint64_t  k = 1;
    uint32_t  nprimes = 0;
    prime_iter_t it;
    prime_iter_init(&it, 2, b);
    for (uint64_t p; (p = prime_iter_next(&it)) != 0; ) {
        if (   stoppable && ((++nprimes & 0xffff) == 0)
            && tifa_should_stop()) {
            prime_iter_clear(&it);
            free(scalars.data);
            return false;
        }
        uint64_t q = p;
        while (q <= b / p) {
            q *= p;
        }
        if (k > UINT64_MAX / q) {
            push_u64(&scalars, k);
            k = 1;
        }
        k *= q;
    }
    prime_iter_clear(&it);
    push_u64(&scalars, k);
    const bool done = product_of_scalars(e, scalars.data, 0, scalars.length,
                                         stoppable);
    free(scalars.data);
    return done;
}
