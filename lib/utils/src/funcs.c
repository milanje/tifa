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
 * \file    funcs.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <pthread.h>

#include "first_primes.h"

#include "funcs.h"
#include "factor_u64.h"
#include "macros.h"
#include "tifa_alloc.h"

/*
 *-----------------------------------------------------------------------------
 *                      Miscellaneous functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static const uint8_t log2table[] = {
    //
    // From public domain code from the Bit Twiddling Hacks web page:
    // http://graphics.stanford.edu/~seander/bithacks.html
    //
  0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};
//-----------------------------------------------------------------------------
uint32_t most_significant_bit(uint32_t n) {
    //
    // Returns the most significant bit of an integer in (more or less)
    // constant time. Bit 0 is the least significant bit, bit 31 is the
    // most significant.
    //
    // From public domain code from the Bit Twiddling Hacks web page:
    // http://graphics.stanford.edu/~seander/bithacks.html
    //
    uint32_t tmp1;
    uint32_t tmp2;

    tmp2 = (n >> 16);

    if (tmp2 != 0) {
        tmp1 = (n >> 24);
        if (tmp1 != 0) {
            return (24 + log2table[tmp1]);
        } else {
            return (16 + log2table[tmp2 & 0xFF]);
        }
    } else {
        tmp1 = (n >> 8);
        if (tmp1 != 0) {
            return (8 + log2table[tmp1]);
        } else {
            return log2table[n];
        }
    }
}
//-----------------------------------------------------------------------------
void augment_coprime_base(mpz_t f, mpz_array_t* base) {
    //
    // Try to augment the coprime base with 'f':
    //
    //     - If 'f' is coprime with all other integers in the 'base' array:
    //          - add 'f' in the base
    //
    //     - Else if y is the first element in the base for which
    //       gcd(y, 'f') != 1:
    //          - If gcd(y, 'f') == 'f':
    //                - add 'f' in the base
    //                - remove y from the base
    //                - call augment_coprime_base for y/gcd(y, 'f')
    //          - If gcd(y, 'f') == y:
    //                - keep y in the base
    //                - call augment_coprime_base for 'f'/gcd(y, 'f')
    //          - Otherwise:
    //                - add gcd in the base
    //                - remove y from the base
    //                - call augment_coprime_base for f/gcd(y, 'f')
    //                - call augment_coprime_base for y/gcd(y, 'f')
    //
    if (is_in_mpz_array(f, base)) {
        return;
    }
    if (mpz_cmp_ui(f, 1) == 0) {
        return;
    }
    if (base->length == 0) {
        append_mpz_to_array(base, f);
        return;
    }

    uint32_t size_f = mpz_sizeinbase(f, 2);
    uint32_t len    = base->length;

    bool coprime_to_all_others = true;

    mpz_t cofactor_1;
    mpz_t cofactor_2;
    mpz_t gcd;

    mpz_init2(cofactor_1, size_f);
    mpz_init2(cofactor_2, size_f);
    mpz_init2(gcd, size_f);

    for (uint32_t i = 0; i < len; i++) {

        mpz_gcd(gcd, base->data[i], f);

        if ((mpz_cmp_ui(gcd, 1) == 0) ) {
            continue;
        }
        coprime_to_all_others = false;

        if ((mpz_cmp(gcd, f) == 0) ) {
            //
            // Keep the gcd in the base and call recursively
            // augment_coprime_base with base->data[i]'s cofactor.
            //
            mpz_divexact(cofactor_2, base->data[i], gcd);
            mpz_set(base->data[i], gcd);
            augment_coprime_base(cofactor_2, base);
            break;
        }
        if ((mpz_cmp(gcd, base->data[i]) == 0) ) {
            //
            // Keep the gcd in the base and call recursively
            // augment_coprime_base with f's cofactor.
            //
            mpz_divexact(cofactor_2, f, gcd);
            augment_coprime_base(cofactor_2, base);
            break;
        }
        //
        // Here, we have the following inequalities:
        //     1) gcd != 1
        //     2) gcd != f
        //     3) gcd != base->data[i]
        //
        // Now, keep gcd in the base array and call recursively
        // augment_coprime_base for the found cofactors
        //
        mpz_divexact(cofactor_1, f, gcd);
        mpz_divexact(cofactor_2, base->data[i], gcd);
        mpz_set(base->data[i], gcd);
        augment_coprime_base(cofactor_1, base);
        augment_coprime_base(cofactor_2, base);
        break;
    }
    if (coprime_to_all_others) {
        //
        // Add the integer f in the base if it is coprime with all other
        // integers in the base.
        //
        append_mpz_to_array(base, f);
    }
    mpz_clear(gcd);
    mpz_clear(cofactor_1);
    mpz_clear(cofactor_2);
}
//-----------------------------------------------------------------------------
void find_coprime_base(mpz_array_t* const base, const mpz_t n,
                       const mpz_array_t* const factors) {

    mpz_t cofactor;
    mpz_init2(cofactor, mpz_sizeinbase(n, 2));

    uint32_t len = factors->length;

    for (uint32_t i = 0 ; i < len; i++) {
        augment_coprime_base(factors->data[i], base);
        mpz_divexact(cofactor, n, factors->data[i]);
        augment_coprime_base(cofactor, base);
    }
    mpz_clear(cofactor);
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                      Number theoretical functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
const int8_t kronecker_data[8] = {0, 1, 0, -1, 0, -1, 0, 1};
//-----------------------------------------------------------------------------
int8_t  kronecker_ui(uint32_t a, uint32_t b) {
    //
    // Algorithm 1.4.10 from the book "A Course in Computational Algebraic
    // Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (b == 0) {
        if (a != 1) {
            return 0;
        } else {
            return 1;
        }
    }
    if (ARE_EVEN(a, b)) {
        return 0;
    }
    uint32_t v = 0;
    while (IS_EVEN(b)) {
        v++;
        b = b >> 1;
    }
    int32_t k = 1;
    if (IS_ODD(v)) {
        k = kronecker_data[a & 7];      // k = (-1)^((a^2-1)/8)
    }
    uint32_t r = 0;
    while (1) {
        if (a == 0) {
            if (b == 1) {
                return k;
            } else {
                return 0;
            }
        }
        v = 0;
        while (IS_EVEN(a)) {
            v++;
            a = a >> 1;
        }
        if (IS_ODD(v)) {
            k *= kronecker_data[b & 7]; // k = (-1)^((b^2-1)/8) * k
        }
        if (a & b & 2) {                // k = (-1)^((a-1)(b-1)/4) * k
            k = -k;
        }
        r = a;
        a = b % r;
        b = r;
    }
}
//-----------------------------------------------------------------------------
uint32_t powm(uint32_t base, uint32_t power, uint32_t modulus) {
    //
    // Left-right binary algorithm to compute (base^power) mod modulus.
    //
    // See for example algorithm 1.2.3 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (power == 0) {
        return 1;
    }
    uint32_t n = power;
    uint32_t z = base % modulus;

    uint64_t y = z;
    uint32_t f = most_significant_bit(power);

    while (f != 0) {
        f--;
        y *= y;
        y %= modulus;
        if (BIT(n, f) == 1) {
            y *= z;
            y %= modulus;
        }
    }
    return (uint32_t)y;
}
//-----------------------------------------------------------------------------
uint32_t sqrtm(uint32_t a, uint32_t p) {
    //
    // Shanks's algorithm for modular square roots.
    //
    // See for example algorithm 1.5.1 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (a == 0) {
        //
        // a does not have a modular square root!
        //
        return NO_SQRT_MOD_P;
    }
    if (kronecker_ui(a, p) != 1) {
        return NO_SQRT_MOD_P;
    }
    //
    // For p = 3 (mod 4) and p = 5 (mod 8), a single exponentiation, giving
    // the root the general algorithm below would (2 is the non-residue it
    // finds for p = 5 mod 8).
    //
    if ((p & 3) == 3) {
        return p - powm(a, (p + 1) / 4, p);
    }
    if ((p & 7) == 5) {
        const uint32_t q  = (p - 1) / 4;
        const uint64_t y  = powm(a, q / 2, p);
        const uint64_t x0 = (y * a) % p;
        if ((((y * y) % p) * a) % p == 1) {
            return p - (uint32_t) x0;
        }
        return (uint32_t) ((x0 * powm(2, q, p)) % p);
    }
    //
    // Find n such that p = q.2^n + 1 with q odd
    //
    int n = ffs(p - 1) - 1;
    int k = n;

    uint32_t q = (p - 1) >> n;
    //
    // Find u, a quadratic non-residue mod p.
    //
    // _WARNING_: Here p should be an odd prime!
    //
    // _NOTE_: In his book "A Course in Computational Algebraic Number Theory",
    //         Henri Cohen suggests to choose u at random until a non-residue
    //         is found instead of a sequential search. In practice however, it
    //         should not make any significant difference so we follow here the
    //         path of least resistance.
    //
    int u = 1;
    while (kronecker_ui(u, p) != -1) {
        u++;
    }
    uint64_t t   = 0;
    uint64_t b2m = 0;

    uint64_t z   = powm(u, q, p);
    uint64_t x   = powm(a, q / 2, p);
    uint64_t b   = (x * x) % p;

    if (IS_ODD(q)) {
        b = (b * a) % p;
        x = (x * a) % p;
    }
    //
    // Now x = powm(a, (q + 1) / 2, p) and b = powm(a, q, p)
    //
    while (1) {
        int m = 0;
        b2m   = b;
        //
        // Find the least integer m such that b^(2^m) = 1 (mod p).
        //
        while (b2m != 1) {
            m++;
            b2m = (b2m * b2m) % p;
        }
        if (m == k) {
            //
            // a does not have a modular square root!
            //
            return NO_SQRT_MOD_P;
        }
        t = powm(z, 1 << (k - m - 1), p);
        z = (t * t) % p;
        b = (b * z) % p;
        x = (x * t) % p;

        if (b == 1) {
            return (uint32_t)x;
        }
        k = m;
    }
}
//-----------------------------------------------------------------------------
const unsigned short qres_mod_315[315] = {
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0
};
const unsigned short qres_mod_256[256] = {
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0
};
const unsigned short qres_mod_221[221] = {
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 1, 0, 1, 0,
    0, 0, 1, 1, 0, 0, 0, 0, 0, 1,
    0, 1, 1, 1, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 1, 0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
    1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 1, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
    1, 0, 1, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
    0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    1
};
extern unsigned long int is_square(unsigned long int x);
//-----------------------------------------------------------------------------
//
// The following threshold values were experimentally determined on a
// PowerPC 7455 (for the 32 bit case) and on an Opteron 250 (64 bit) using
// respectively GCC 4.0.0 and GCC 4.1.2. Be aware that this threshold is
// probably very different on an 32 bit x86 processor.
//
#if TIFA_WORDSIZE == 64
    #define PRIMALITY_TDIV_THRESHOLD 2500000
#elif TIFA_WORDSIZE == 32
    #define PRIMALITY_TDIV_THRESHOLD 80000000
#else
    //
    // This is picked up out of the blue and should be determined... if we
    // really care about very (old|odd) architectures.
    //
    #define PRIMALITY_TDIV_THRESHOLD 15000
#endif
//-----------------------------------------------------------------------------
bool is_prime(uint32_t n) {
    //
    // A very basic (Dubois-Selfridge-) Miller-Rabin composition test.
    //
    // See for example algorithm 8.2.2 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993. The
    // description given in this book has been straitforwardly adapted.
    //
    // We also roughly follow the strategy used in the mpz_probab_prime_p()
    // function from the GMP library with some trial divisions up to a certain
    // threshold.
    //
    if (n == 1) { return false; }
    if (n == 2) { return true; }

    if (IS_EVEN(n)) {
        return false;
    }
    //
    // Trial division if n is smaller than PRIMALITY_TDIV_THRESHOLD. See above
    // comment in this symbol's definition.
    //
    uint32_t q;
    uint32_t r = 1;
    uint32_t d = 3;

    if (n < PRIMALITY_TDIV_THRESHOLD) {
        while (r != 0) {
            q = n / d;
            r = n - q * d;
            if (q < d) {
   	            return true;
            }
            d += 2;
        }
        return false;
    }
    r = n % 255255; // r = n mod (3 * 5 * 7 * 11 * 13 * 17)

    if (   !(r % 3) || !(r % 5) || !(r % 7) || !(r % 11) || !(r % 13)
        || !(r % 17)) {
        return false;
    }

    r = n % 392863; // r = n mod (19 * 23 * 29 * 31)

    if (!(r % 19) || !(r % 23) || !(r % 29) || !(r % 31)) {
        return false;
    }

    //
    // Strong probable prime tests to the bases 2, 7 and 61, which are
    // conclusive below 4759123141 (Jaeschke).
    //
    static const uint32_t bases[] = { 2, 7, 61 };

    q = n - 1;
    int t = 0;

    while (IS_EVEN(q)) {
        q >>= 1;
        t++;
    }
    for (uint32_t ib = 0; ib < 3; ib++) {
        int      e = 0;
        uint32_t a = bases[ib];
        uint64_t b = powm(a, q, n);

        if (b != 1) {
            while ( (b != 1) && (b != (n - 1)) && (e <= (t - 2)) ) {
                b *= b;
                b %= n;
                e += 1;
            }
            if (b != n - 1) {
                return false;
            }
        }
    }
    return true;
}
//-----------------------------------------------------------------------------
bool is_prime_mpz(const mpz_t n) {
#if TIFA_HAVE_FACTOR_U64
    if (mpz_sizeinbase(n, 2) <= 64) {
        return is_prime_u64((uint64_t) mpz_get_ui(n));
    }
#endif
    return (0 != mpz_probab_prime_p(n, NMILLER_RABIN));
}
//-----------------------------------------------------------------------------
unsigned long int gcd_ulint(unsigned long int a, unsigned long int b) {
    //
    // Standard "right-shift" binary gcd algorithm.
    //
    // See for example algorithm 1.3.5 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    unsigned long int tmp = 0;
    if (a < b) {
        tmp = a;
        a   = b;
        b   = tmp;
    }
    tmp = a % b;
    a   = b;
    b   = tmp;

    if (b == 0) {
        return a;
    }
    unsigned long int shift = 0;

    while (ARE_EVEN(a, b)) {
        shift++;
        a >>= 1;
        b >>= 1;
    }
    while (IS_EVEN(a)) {
        a >>= 1;
    }
    while (IS_EVEN(b)) {
        b >>= 1;
    }
    unsigned long int t = 0;

    do {
        if (a > b) {
            t   = a;
            t  -= b;
            t >>= 1;
            if (t == 0) {
                return (a << shift);
            }
            while (IS_EVEN(t)) {
                t >>= 1;
            }
            a = t;
        } else {
            t   = b;
            t  -= a;
            t >>= 1;
            if (t == 0) {
                return (a << shift);
            }
            while (IS_EVEN(t)) {
                t >>= 1;
            }
            b = t;
        }
    } while (1);

    return (a << shift);
}
//-----------------------------------------------------------------------------
//
// The following macros are used by the modinv_ui function, where the
// *_bit_mask variables are declared, and are thus not intended to be reused
// elsewhere.
//
#define BIT_N(x)           ( (x) & (nth_bit_mask))
#define BIT_N_MINUS_ONE(x) ( (x) & (nth_minus_one_bit_mask) )
#define HIGH_BITS(x)       ( (x) & (high_bit_mask) )
#define LOW_BITS(x)        ( (x) & (low_bit_mask) )
//-----------------------------------------------------------------------------
unsigned long int modinv_ui(unsigned long int num, unsigned long int p) {
    //
    // Reference:
    //      "New Algorithm for Classical Modular Inverse",
    //      Robert Lorencz,
    //      Lecture Notes In Computer Science, Volume 2523/2003,
    //      Revised Papers from the 4th International Workshop on
    //      Cryptographic Hardware and Embedded Systems
    //
    unsigned long int a = num % p;
    long int u = p;
    long int v = a;
    long int r = 0;
    long int s = 1;
    unsigned long int cu = 0;
    unsigned long int cv = 0;
    long int exp_cu = 1;
    long int exp_cv = 1;
    unsigned long int n  = 0;
    for (unsigned long int t = p; t != 0; t >>= 1) {
        n++;
    }
    //
    // If n is the number of bits of p, compute masks to be used by the
    // previously defined macros, to:
    //   - keep the nth bit (i.e. the most significant one);
    //   - keep the nth-1 bit (i.e. the next most significant one);
    //   - keep the 2 most significant bits;
    //   - keep the nth-2 least significant bits.
    //
    const unsigned long int high_bit_mask = 3UL << (n-1);
    const unsigned long int low_bit_mask  = (1UL << (n-1)) - 1;
    const unsigned long int nth_bit_mask  = 1UL << n;
    const unsigned long int nth_minus_one_bit_mask = 1UL << (n-1);

    //
    // u, v, r and s may be negative: they are doubled and halved, not
    // shifted. r and s are always even when halved.
    //
    while ((u != exp_cu) && (u != -exp_cu) && (v != exp_cv) && (v != -exp_cv)) {

        if (!HIGH_BITS(u) || (BIT_N(u) && BIT_N_MINUS_ONE(u) && LOW_BITS(u))) {
            if (cu >= cv) {
                u *= 2;
                r *= 2;
                cu++;
                exp_cu <<= 1;
            } else {
                u *= 2;
                s /= 2;
                cu++;
                exp_cu <<= 1;
            }
        } else {
            if (   !HIGH_BITS(v)
                || (BIT_N(v) && BIT_N_MINUS_ONE(v) && LOW_BITS(v)) ) {

               if (cv >= cu) {
                   v *= 2;
                   s *= 2;
                   cv++;
                   exp_cv <<= 1;
               } else {
                   v *= 2;
                   r /= 2;
                   cv++;
                   exp_cv <<= 1;
               }
            } else {
                if (BIT_N(v) == BIT_N(u)) {
                    if (cu <= cv) {
                        u -= v;
                        r -= s;
                    } else {
                        v -= u;
                        s -= r;
                    }
                } else {
                    if (cu <= cv) {
                        u += v;
                        r += s;
                    } else {
                        v += u;
                        s += r;
                    }
                }
            }
        }
    }
    if ((v == exp_cv) || (v == -exp_cv)) {
        r = s;
        u = v;
    }
    if (BIT_N(u) != 0) {
        if (r < 0) {
            r = -r;
        } else {
            r = p - r;
            if (r < 0) {
                r += p;
            }
        }
    } else {
        if (r < 0) {
            r += p;
        }
    }
    return r;
}
//-----------------------------------------------------------------------------
unsigned long int sqrtm_p2(uint32_t n, uint32_t p) {
    //
    // _NOTE_: Now this is getting a bit messy. Indeed, the code now mixes
    //         apparently randomly (unsigned) long ints with C99 integer types
    //         like uint32_t. This can be seen as bad pratice and indeed, it is
    //         not very pretty. The reason for this is to convey an idea of
    //         the size of the numbers involved. Maybe only using long ints and
    //         clearly documenting the expected range of the variables would
    //         have been better.
    //
    long int s = 0;

    uint32_t np = n % p;
    uint32_t sp = sqrtm(np, p);

    if (sp == NO_SQRT_MOD_P) {
        //
        // No solution
        //
        return NO_SQRT_MOD_P2;
    }

    s = sp * sp;
    s = n - s;

    uint32_t inv = modinv_ui(sp << 1, p);

    s /= (long int)p;
    s *= inv;

    if (s < 0) {
        s = p - ((-s) % p);
    } else {
        s = s % p;
    }

    s *= p;
    s += sp;

    return (unsigned long int)s;
}
//-----------------------------------------------------------------------------
#define KS_NPRIMES 300
#define KS_ANY_RESIDUE_BITSIZE 100
//
// The candidate multipliers: the odd squarefree k up to LARGEST_MULTIPLIER.
//
#define KS_MAX_MULTS ((LARGEST_MULTIPLIER + 1) / 2)
//
// What does not depend on n, computed once: for each candidate k, the
// product of the (q/p) over the primes q dividing k (0 if p divides k), and
// log(p)/(p-1), for the first KS_NPRIMES primes p.
//
static uint32_t       ks_mults[KS_MAX_MULTS];
static uint32_t       ks_nmults;
static int8_t         ks_legendre[KS_MAX_MULTS][KS_NPRIMES];
static double         ks_logp[KS_NPRIMES];
static double         ks_logp2[KS_NPRIMES];
static pthread_once_t ks_once = PTHREAD_ONCE_INIT;

static void ks_init(void) {
    //
    // (q/p) for q <= LARGEST_MULTIPLIER, from quadratic reciprocity and the
    // squares mod q.
    //
    uint32_t nq = 0;
    while (first_primes[nq] <= LARGEST_MULTIPLIER) {
        nq++;
    }
    int8_t (*legendre_q)[KS_NPRIMES] = tifa_malloc(nq * sizeof(*legendre_q));
    uint8_t isqr[LARGEST_MULTIPLIER + 1];

    for (uint32_t j = 0; j < nq; j++) {
        const uint32_t q = first_primes[j];

        if (q > 2) {
            memset(isqr, 0, q);
            for (uint32_t x = 1; x < q; x++) {
                isqr[(x * x) % q] = 1;
            }
        }
        for (uint32_t i = 1; i < KS_NPRIMES; i++) {
            const uint32_t p = first_primes[i];
            int8_t l;

            if (p == q) {
                l = 0;
            } else if (q == 2) {
                l = (((p & 7) == 1) || ((p & 7) == 7)) ? 1 : -1;
            } else {
                l = isqr[p % q] ? 1 : -1;
                if (((p & 3) == 3) && ((q & 3) == 3)) {
                    l = -l;
                }
            }
            legendre_q[j][i] = l;
        }
    }
    for (uint32_t i = 1; i < KS_NPRIMES; i++) {
        const uint32_t p = first_primes[i];
        ks_logp[i]  = log((double) p) / (p - 1);
        ks_logp2[i] = 2.0 * ks_logp[i];
    }
    ks_nmults = 0;
    for (uint32_t k = 1; k <= LARGEST_MULTIPLIER; k += 2) {
        if ((k % 9 == 0) || (k % 25 == 0) || (k % 49 == 0)) {
            continue;
        }
        int8_t* const row = ks_legendre[ks_nmults];
        for (uint32_t i = 1; i < KS_NPRIMES; i++) {
            row[i] = 1;
        }
        for (uint32_t j = 0; j < nq; j++) {
            if (k % first_primes[j] == 0) {
                for (uint32_t i = 1; i < KS_NPRIMES; i++) {
                    row[i] *= legendre_q[j][i];
                }
            }
        }
        ks_mults[ks_nmults++] = k;
    }
    free(legendre_q);
}

uint32_t ks_multiplier(const mpz_t n, const uint32_t size_base) {
    //
    // Knuth-Schroeppel function, as in msieve and YAFU: for each odd
    // squarefree k, add the expected log contribution of the small primes to
    // the sieve values of k.n (2.log(p)/(p-1) when k.n is a quadratic residue
    // mod p, log(p)/(p-1) when p divides k, and a 2-adic term depending on
    // k.n mod 8), minus the (1/2).log(k) the larger values cost.
    //
    // Below KS_ANY_RESIDUE_BITSIZE bits only k.n = 1 mod 8 is allowed: the
    // SIQS thresholds for small numbers were tuned with that restriction, and
    // lifting it there costs 5-25%, against 2-12% gained above 120 bits.
    //
    // Primes beyond the first few hundred change the score very little, so
    // the sum stops at KS_NPRIMES primes (fewer for small factor bases).
    //
    pthread_once(&ks_once, ks_init);

    const uint32_t nprimes = MIN(2 * size_base + 1, KS_NPRIMES);
    const bool     any_residue
        = (mpz_sizeinbase(n, 2) >= KS_ANY_RESIDUE_BITSIZE);

    int8_t legendre_n[KS_NPRIMES];
    for (uint32_t i = 1; i < nprimes; i++) {
        const uint32_t p = first_primes[i];
        legendre_n[i] = kronecker_ui(mpz_fdiv_ui(n, p), p);
    }
    const uint32_t nmod8 = mpz_fdiv_ui(n, 8);
    const double   log2  = log(2.0);

    //
    // The candidates are scored four at a time: one sum is bound by the
    // latency of its additions, four independent ones overlap. Each still
    // adds its terms in the same order.
    //
    uint32_t cand[KS_MAX_MULTS + 3];
    double   score[KS_MAX_MULTS + 3];
    uint32_t ncand = 0;
    for (uint32_t c = 0; c < ks_nmults; c++) {
        const uint32_t k = ks_mults[c];
        if (!any_residue && (((k * nmod8) & 7) != 1)) {
            continue;
        }
        double sc = -0.5 * log((double) k);

        switch ((k * nmod8) & 7) {
            case 1:
                sc += 2.0 * log2;
                break;
            case 5:
                sc += log2;
                break;
            case 3:
            case 7:
                sc += 0.5 * log2;
                break;
            default:
                break;
        }
        cand[ncand]  = c;
        score[ncand] = sc;
        ncand++;
    }
    const uint32_t nblocks = (ncand + 3) / 4;
    for (uint32_t c = ncand; c < 4 * nblocks; c++) {
        cand[c]  = cand[0];
        score[c] = 0.0;
    }
    for (uint32_t b = 0; b < nblocks; b++) {
        const int8_t* const r0 = ks_legendre[cand[4 * b]];
        const int8_t* const r1 = ks_legendre[cand[4 * b + 1]];
        const int8_t* const r2 = ks_legendre[cand[4 * b + 2]];
        const int8_t* const r3 = ks_legendre[cand[4 * b + 3]];
        double s0 = score[4 * b];
        double s1 = score[4 * b + 1];
        double s2 = score[4 * b + 2];
        double s3 = score[4 * b + 3];
        for (uint32_t i = 1; i < nprimes; i++) {
            const int8_t  ln = legendre_n[i];
            const double  w1 = ks_logp[i];
            const double  w2 = ks_logp2[i];
            s0 += (r0[i] == 0) ? w1 : ((r0[i] * ln == 1) ? w2 : 0.0);
            s1 += (r1[i] == 0) ? w1 : ((r1[i] * ln == 1) ? w2 : 0.0);
            s2 += (r2[i] == 0) ? w1 : ((r2[i] * ln == 1) ? w2 : 0.0);
            s3 += (r3[i] == 0) ? w1 : ((r3[i] * ln == 1) ? w2 : 0.0);
        }
        score[4 * b]     = s0;
        score[4 * b + 1] = s1;
        score[4 * b + 2] = s2;
        score[4 * b + 3] = s3;
    }
    uint32_t mult      = 1;
    double   bestscore = 0.0;
    for (uint32_t c = 0; c < ncand; c++) {
        if ((c == 0) || (score[c] > bestscore)) {
            bestscore = score[c];
            mult      = ks_mults[cand[c]];
        }
    }
    //
    // k is squarefree: kn is a square, which would make polynomial values
    // vanish, only when k divides n and n/k is a square.
    //
    if ((mult != 1) && mpz_divisible_ui_p(n, mult)) {
        mpz_t q;
        mpz_init(q);
        mpz_divexact_ui(q, n, mult);
        if (mpz_perfect_square_p(q)) {
            mult = 1;
        }
        mpz_clear(q);
    }
    return mult;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                          Hash functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
uint32_t hash_rj_32(const void* const keyptr) {

    uint32_t hash = *((uint32_t*)keyptr);
    //
    // Robert Jenkins' 32 bit mix function.
    // http://www.concentric.net/~Ttwang/tech/inthash.htm
    //
    hash += (hash << 12);
    hash ^= (hash >> 22);
    hash += (hash << 4);
    hash ^= (hash >> 9);
    hash += (hash << 10);
    hash ^= (hash >> 2);
    hash += (hash << 7);
    hash ^= (hash >> 12);
    return hash;
}
//-----------------------------------------------------------------------------
uint32_t hash_pjw(const void* const keyptr) {
    //
    // Found in the Dragon book:
    // "Compilers: Principles, Techniques and Tools", Aho, Sethi, & Ullman.
    //
    // Attributed to P. J. Weinberger.
    // Supposed to perform well on strings.
    //
    char *p;
    uint32_t hash = 0;
    uint32_t g = 0;

    for (p = (char*)keyptr; *p != '\0'; p++) {
        hash = (hash << 4) + (*p);
        g = (hash & 0xf0000000);
        if (0 != g) {
            hash ^= (g >> 24);
            hash ^= g;
        }
    }
    return hash;
}
//-----------------------------------------------------------------------------
#define GET_16_BITS(d) (*((const uint16_t *) (d)))
//-----------------------------------------------------------------------------
uint32_t hash_sfh_ph(const void* const keyptr) {
    //
    // The so-called "SuperFastHash" function By Paul Hsieh.
    // http://www.azillionmonkeys.com/qed/hash.html
    //
    // The original code has been slightly altered, but the logic is
    // unchanged. Refer to the aforementioned URL for P. Hsieh's original code.
    //
    if ((keyptr == NULL)) {
        return 0;
    }
    const char* data = (const char*)keyptr;

    uint32_t length = strlen(data)*sizeof(char);

    uint32_t hash = length;
    uint32_t tmp;

    uint32_t rem = length & 3;
    length >>= 2;

    while (length > 0) {
        //
        // Use the next 32 bits of *data
        //
        length--;
        hash  += GET_16_BITS(data);
        tmp    = (GET_16_BITS(data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2 * sizeof(uint16_t);
        hash  += hash >> 11;
    }
    switch (rem) {
        //
        // If the number of bits of *data is not a multiple of 32
        // i.e. (length % 4) != 0
        //
        case 3:
            hash += GET_16_BITS(data);
            hash ^= hash << 16;
            hash ^= data[sizeof(uint16_t)] << 18;
            hash += hash >> 11;
            break;
        case 2:
            hash += GET_16_BITS(data);
            hash ^= hash << 11;
            hash += hash >> 17;
            break;
        case 1:
            hash += *data;
            hash ^= hash << 10;
            hash += hash >> 1;
    }
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 2;
    hash += hash >> 15;
    hash ^= hash << 10;

    return hash;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                          Comparison functions
 *-----------------------------------------------------------------------------
 */

//
// The following functions return an int instead of an int32_t only
// for compatibility with pre-C99 code such as the standard library
// qsort function...
//

//-----------------------------------------------------------------------------
int mpz_cmp_func(const void* const mpza, const void* const mpzb) {

    int cmp = mpz_cmp(*((mpz_t*)mpza), *((mpz_t*)mpzb));
    if (cmp == 0) return 0;
    if (cmp > 0) {
        return 1;
    } else {
        return -1;
    }
}
//-----------------------------------------------------------------------------
int uint32_cmp_func(const void* const uinta, const void* const uintb) {

    if (*((uint32_t*)uinta) == *((uint32_t*)uintb)) return 0;
    if (*((uint32_t*)uinta) > *((uint32_t*)uintb)) {
        return 1;
    } else {
        return -1;
    }
}
//-----------------------------------------------------------------------------
int string_cmp_func(const void* const stra, const void* const strb) {
    const char* a = *(const char**)stra;
    const char* b = *(const char**)strb;
    int cmp = strcmp(a, b);
    if (cmp == 0) return 0;
    if (cmp > 0) {
        return 1;
    } else {
        return -1;
    }
}
//-----------------------------------------------------------------------------
int cmp_mult_data(const void* mda, const void* mdb) {

    const mult_data_t* const a = (mult_data_t*)mda;
    const mult_data_t* const b = (mult_data_t*)mdb;

    if (a->count > b->count) {
        return 1;
    }
    if (a->count < b->count) {
        return -1;
    }
    //
    // Here, a->count == b->count
    //
    if (a->sum_inv_pi > b->sum_inv_pi) {
        return 1;
    }
    if (a->sum_inv_pi < b->sum_inv_pi) {
        return -1;
    }
    //
    // Here, a->count      == b->count
    // and   a->sum_inv_pi == b->sum_inv_pi
    //
    // The better (i.e. greater) mult_data_t is the one with the smaller
    // multiplier.
    //
    if (a->multiplier < b->multiplier) {
        return 1;
    }
    if (a->multiplier > b->multiplier) {
        return -1;
    }
    return 0;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                        Combinatorial functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
uint32_t n_choose_k(uint8_t n, uint8_t k) {
    if (k > n) {
        return 0;
    }
    if (k > (n / 2)) {
        return n_choose_k(n, n - k);
    }

    uint64_t nck = 1;
    for (uint8_t d = 1; d <= k; d++) {
      nck *= n;
      n--;
      nck /= d;
    }
    return (uint32_t)nck;
}
//-----------------------------------------------------------------------------
void next_subset_lex(uint32_t n, uint32_t k, uint32_t* subset, bool* end) {
    //
    // This is Algorithm 2.6 from the book "Combinatorial Algorithms -
    // Generation, Enumeration, and Search" by Donald L. Kreher and Douglas
    // Stinson.
    //
    // This function is adapted from the implementation provided by the
    // authors at:
    //
    //      http://www.math.mtu.edu/~kreher/cages/Src.html
    //
    int32_t i;
    uint32_t tmp[k];

    *end = false;

    for(i = 0; i < (int32_t)k; i++) {
        tmp[i] = subset[i];
    }
    i = (int32_t)k - 1;

    while ((i >= 0) && (subset[i] == (n - k + i + 1)) ) {
        i--;
    }

    if (i == -1) {
        *end = true;
    } else {
        for (uint32_t j = i; j < k; j++) {
            tmp[j] = subset[i] + j - i + 1;
        }
        for (uint32_t j = 0; j < k ; j++) {
            subset[j] = tmp[j];
        }
    }
}
//-----------------------------------------------------------------------------
void unrank_subset_lex(uint32_t n, uint32_t k, uint32_t r, uint32_t* subset) {
    //
    // This is Algorithm 2.8 from the book "Combinatorial Algorithms -
    // Generation, Enumeration, and Search" by Donald L. Kreher and Douglas
    // Stinson.
    //
    // This function is adapted from the implementation provided by the
    // authors at:
    //
    //      http://www.math.mtu.edu/~kreher/cages/Src.html
    //
    uint32_t x = 1;
    uint32_t y;

    for (uint32_t i = 0; i < k; i++) {

        y = n_choose_k(n - x, k - i - 1);

        while (y <= r) {
            r = r - y;
            x++;
            y = n_choose_k(n - x, k - i - 1);
        }
        subset[i] = x;
        x++;
    }
}
//-----------------------------------------------------------------------------
void unrank_subset_revdoor(uint32_t n, uint32_t k,
                           uint32_t r, uint32_t* subset) {
    //
    // This is Algorithm 2.12 from the book "Combinatorial Algorithms -
    // Generation, Enumeration, and Search" by Donald L. Kreher and Douglas
    // Stinson.
    //
    // This function is adapted from the implementation provided by the
    // authors at:
    //
    //      http://www.math.mtu.edu/~kreher/cages/Src.html
    //
    uint32_t x = n;

    for (int32_t i = k - 1; i >= 0; i = i - 1) {
        uint32_t y = n_choose_k(x, i + 1);

        while (y > r) {
            x = x - 1;
            y = n_choose_k(x, i + 1);
        }
        subset[i] = x+1;
        r = n_choose_k(x + 1, i + 1) - r - 1;
    }
}
//-----------------------------------------------------------------------------

