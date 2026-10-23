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

// Fixed-limb Montgomery p-1 and p+1 (mont_ecm.h) against plain GMP
// references: identical stage 1 outcomes, and every factor the reference
// finds in stage 2 is found too. Moduli of 1 to 4 limbs, with factors p
// chosen so that p - 1 (or p + 1) is smooth half of the time.

#include <stdlib.h>

#include "tifa_test.h"
#include "first_primes.h"
#include "mont_ecm.h"

#if TIFA_HAVE_MONT_ECM

#define B1 200
#define B2 5000

static void exponent(mpz_t e, uint32_t b1) {
    mpz_set_ui(e, 1);
    for (uint32_t i = 0; first_primes[i] <= b1; i++) {
        unsigned long q = first_primes[i];
        while (q <= b1 / first_primes[i]) {
            q *= first_primes[i];
        }
        mpz_mul_ui(e, e, q);
    }
}

// 1 if gcd(a, n) is a proper factor (stored in g), -1 if it is n, 0 if 1.
static int classify(mpz_t g, const mpz_t a, const mpz_t n) {
    mpz_gcd(g, a, n);
    if (mpz_cmp_ui(g, 1) == 0) {
        return 0;
    }
    return (mpz_cmp(g, n) == 0) ? -1 : 1;
}

// V_k(p) mod n.
static void lucas_v(mpz_t v, const mpz_t p, const mpz_t k, const mpz_t n) {
    mpz_t v1, t;
    mpz_init(v1);
    mpz_init(t);
    mpz_set(v, p);
    mpz_mul(v1, p, p);
    mpz_sub_ui(v1, v1, 2);
    mpz_mod(v1, v1, n);
    for (long b = (long) mpz_sizeinbase(k, 2) - 2; b >= 0; b--) {
        mpz_mul(t, v, v1);
        mpz_sub(t, t, p);
        mpz_mod(t, t, n);
        if (mpz_tstbit(k, b)) {
            mpz_set(v, t);
            mpz_mul(v1, v1, v1);
            mpz_sub_ui(v1, v1, 2);
            mpz_mod(v1, v1, n);
        } else {
            mpz_set(v1, t);
            mpz_mul(v, v, v);
            mpz_sub_ui(v, v, 2);
            mpz_mod(v, v, n);
        }
    }
    mpz_clear(v1);
    mpz_clear(t);
}

// Reference p-1 (pm1 != 0) or p+1 (seed s): stage 1 outcome in *stage1,
// stage 1 then prime-by-prime stage 2 outcome returned.
static int ref_run(mpz_t g, int* stage1, const mpz_t n, int pm1,
                   unsigned long s) {
    mpz_t e, x, acc, t, k;
    mpz_inits(e, x, acc, t, k, NULL);
    exponent(e, B1);
    mpz_set_ui(t, s);
    if (pm1) {
        mpz_powm(x, t, e, n);
        mpz_sub_ui(t, x, 1);
    } else {
        lucas_v(x, t, e, n);
        mpz_sub_ui(t, x, 2);
    }
    int r = classify(g, t, n);
    *stage1 = r;
    if (r == 0) {
        mpz_set_ui(acc, 1);
        for (uint32_t i = 0; first_primes[i] <= B2; i++) {
            if (first_primes[i] <= B1) {
                continue;
            }
            mpz_set_ui(k, first_primes[i]);
            if (pm1) {
                mpz_powm(t, x, k, n);
                mpz_sub_ui(t, t, 1);
            } else {
                lucas_v(t, x, k, n);
                mpz_sub_ui(t, t, 2);
            }
            mpz_mul(acc, acc, t);
            mpz_mod(acc, acc, n);
        }
        r = classify(g, acc, n);
    }
    mpz_clears(e, x, acc, t, k, NULL);
    return r;
}

static void random_prime(mpz_t p, gmp_randstate_t state, unsigned long bits) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);
    } while (!mpz_probab_prime_p(p, 30));
}

// A prime p of about 24 bits with p - delta (delta = 1 or -1) a product of
// primes below B1, or below B2 for one of them when wide is set.
static void smooth_prime(mpz_t p, gmp_randstate_t state, int delta,
                         int wide) {
    do {
        mpz_set_ui(p, 2);
        if (wide) {
            uint32_t i = 0;
            while (first_primes[i] <= B1) {
                i++;
            }
            uint32_t j = i + gmp_urandomm_ui(state, 600);
            while (first_primes[j] > B2) {
                j--;
            }
            mpz_mul_ui(p, p, first_primes[j]);
        }
        while (mpz_sizeinbase(p, 2) < 24) {
            mpz_mul_ui(p, p, first_primes[1 + gmp_urandomm_ui(state, 40)]);
        }
        if (delta > 0) {
            mpz_add_ui(p, p, 1);
        } else {
            mpz_sub_ui(p, p, 1);
        }
    } while (!mpz_probab_prime_p(p, 30));
}

#endif

int main(void) {
    TIFA_TEST_BEGIN("mont_pm1");

#if TIFA_HAVE_MONT_ECM
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261015UL);

    mpz_t n, p, q, g, gref;
    mpz_inits(n, p, q, g, gref, NULL);

    mont_ecm_plan_t* const plan1 = alloc_mont_ecm_plan(B1, 0);
    mont_ecm_plan_t* const plan2 = alloc_mont_ecm_plan(B1, B2);
    char msg[160];

    for (int pm1 = 1; pm1 >= 0; pm1--) {
        const char* name = pm1 ? "p-1" : "p+1";
        for (uint32_t nl = 1; nl <= 4; nl++) {
            uint32_t mismatch1 = 0, missed2 = 0, bad = 0, nfound = 0;
            for (uint32_t t = 0; t < 30; t++) {
                switch (t % 3) {
                    case 0:
                        random_prime(p, state, 24);
                        break;
                    case 1:
                        smooth_prime(p, state, pm1 ? 1 : -1, 0);
                        break;
                    default:
                        smooth_prime(p, state, pm1 ? 1 : -1, 1);
                        break;
                }
                const unsigned long qbits = 64 * nl - mpz_sizeinbase(p, 2)
                                            - (nl == 1 ? 1 : 5);
                do {
                    random_prime(q, state, qbits);
                    mpz_mul(n, p, q);
                } while (mpz_size(n) != nl);

                const unsigned long s = pm1 ? 2 : 7;
                int ref1;
                const int ref2 = ref_run(gref, &ref1, n, pm1, s);

                const int mine1 = pm1 ? mont_pm1(g, n, s, plan1)
                                      : mont_pp1(g, n, s, plan1);
                if ((mine1 != ref1) || ((mine1 == 1) && mpz_cmp(g, gref))) {
                    mismatch1++;
                }
                const int mine2 = pm1 ? mont_pm1(g, n, s, plan2)
                                      : mont_pp1(g, n, s, plan2);
                if ((ref2 == 1) && (mine2 == 0)) {
                    missed2++;
                }
                if ((mine2 == 1) && !mpz_divisible_p(n, g)) {
                    bad++;
                }
                nfound += (mine2 == 1);
            }
            snprintf(msg, sizeof(msg), "%s, %u limb(s): stage 1 matches",
                     name, nl);
            CHECK_EQ_UI(mismatch1, 0, msg);
            snprintf(msg, sizeof(msg), "%s, %u limb(s): stage 2 finds what "
                     "the reference finds", name, nl);
            CHECK_EQ_UI(missed2, 0, msg);
            snprintf(msg, sizeof(msg), "%s, %u limb(s): factors divide n",
                     name, nl);
            CHECK_EQ_UI(bad, 0, msg);
            snprintf(msg, sizeof(msg), "%s, %u limb(s): finds factors",
                     name, nl);
            CHECK(nfound >= 5, msg);
        }
    }

    free_mont_ecm_plan(plan1);
    free_mont_ecm_plan(plan2);
    mpz_clears(n, p, q, g, gref, NULL);
    gmp_randclear(state);
#endif

    TIFA_TEST_END();
}
