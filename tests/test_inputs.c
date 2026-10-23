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

// Inputs at the edges of the public functions: n below 2, perfect squares,
// output arrays that are full or already hold factors, and parameters out of
// their useful range. Each call must return (no crash, no endless loop), and
// any factorization it reports must be right.

#include "tifa_test.h"
#include "tifa_factor.h"
#include "cfrac.h"
#include "ecm.h"
#include "fermat.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "rho.h"
#include "siqs.h"
#include "squfof.h"
#include "tdiv.h"

// p60 * p68 (128 bits), 1000003 * 1000033 (40 bits), and an odd 300-bit
// composite with a 20-bit factor: 1000003 * (p30 * p226).
static const char* N128 = "105550643491250938224431408772864108703";
static const char* N40  = "1000036000099";
static const char* N300 =
    "368467674441974501388690716582945795262606216917154768567183016512125596391233025592773064113917";

// The factors found from index 'from' on multiply back to n (when the mode
// gives multiplicities) or divide it.
static int found_ok(const mpz_t n, const mpz_array_t* f,
                    const uint32_array_t* m, uint32_t from) {
    if (f->length <= from) {
        return 0;
    }
    mpz_t prod, t;
    mpz_init_set_ui(prod, 1);
    mpz_init(t);
    int ok = 1;
    for (uint32_t i = from; i < f->length; i++) {
        if (!mpz_divisible_p(n, f->data[i]) || (mpz_cmp_ui(f->data[i], 1) <= 0)
            || (mpz_cmp(f->data[i], n) == 0)) {
            ok = 0;
        }
        if (m != NULL) {
            mpz_pow_ui(t, f->data[i], m->data[i]);
            mpz_mul(prod, prod, t);
        }
    }
    if ((m != NULL) && (mpz_cmp(prod, n) != 0)) {
        ok = 0;
    }
    mpz_clear(prod);
    mpz_clear(t);
    return ok;
}

typedef ecode_t (*algo_t)(mpz_array_t*, uint32_array_t*, const mpz_t);

#define WRAP(name, type, setter)                                            \
    static ecode_t run_##name(mpz_array_t* f, uint32_array_t* m,            \
                              const mpz_t n) {                              \
        type p;                                                             \
        setter;                                                             \
        return name(f, m, n, &p, FIND_COMPLETE_FACTORIZATION);              \
    }
WRAP(siqs,      siqs_params_t,      set_siqs_params_to_default(n, &p))
WRAP(cfrac,     cfrac_params_t,     set_cfrac_params_to_default(n, &p))
WRAP(ecm,       ecm_params_t,       set_ecm_params_to_default(n, &p))
WRAP(p_minus_1, p_minus_1_params_t, set_p_minus_1_params_to_default(n, &p))
WRAP(p_plus_1,  p_plus_1_params_t,  set_p_plus_1_params_to_default(n, &p))
WRAP(rho,       rho_params_t,       set_rho_params_to_default(n, &p))
WRAP(squfof,    squfof_params_t,    set_squfof_params_to_default(&p))
WRAP(fermat,    fermat_params_t,    set_fermat_params_to_default(&p))

static void test_small_and_square_n(void) {
    static const algo_t algos[] = {
        run_siqs, run_cfrac, run_ecm, run_p_minus_1, run_p_plus_1, run_rho,
        run_squfof, run_fermat
    };
    static const char* names[] = {
        "siqs", "cfrac", "ecm", "p_minus_1", "p_plus_1", "rho", "squfof",
        "fermat"
    };
    static const char* below2[] = { "0", "1", "-1", "-15" };
    static const char* squares[] = {
        "4", "25", "18446744073709551616", "340282366920938463463374607431768211456"
    };
    mpz_t n;
    mpz_init(n);
    mpz_array_t*    f = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);

    for (unsigned a = 0; a < sizeof(algos) / sizeof(algos[0]); a++) {
        for (unsigned i = 0; i < sizeof(below2) / sizeof(below2[0]); i++) {
            mpz_set_str(n, below2[i], 10);
            f->length = 0;
            m->length = 0;
            CHECK_EQ_UI(algos[a](f, m, n), NO_FACTOR_FOUND, names[a]);
            CHECK_EQ_UI(f->length, 0, names[a]);
        }
        for (unsigned i = 0; i < sizeof(squares) / sizeof(squares[0]); i++) {
            mpz_set_str(n, squares[i], 10);
            f->length = 0;
            m->length = 0;
            CHECK_EQ_UI(algos[a](f, m, n), FAILURE, names[a]);
        }
    }
    for (unsigned i = 0; i < sizeof(below2) / sizeof(below2[0]); i++) {
        mpz_set_str(n, below2[i], 10);
        CHECK_EQ_UI(tdiv(f, m, n, 0), NO_FACTOR_FOUND, "tdiv below 2");
        CHECK_EQ_UI(tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION),
                    NO_FACTOR_FOUND, "tifa_factor below 2");
    }
    free_mpz_array(f);
    free_uint32_array(m);
    mpz_clear(n);
}

static void test_output_arrays(void) {
    mpz_t n, t;
    mpz_init_set_str(n, N128, 10);
    mpz_init(t);

    // An array with no room left grows.
    mpz_array_t*    f = alloc_mpz_array(1);
    uint32_array_t* m = alloc_uint32_array(1);
    CHECK_EQ_UI(tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION),
                COMPLETE_FACTORIZATION_FOUND, "full arrays");
    CHECK(found_ok(n, f, m, 0), "full arrays: factorization");

    // Factors already in the arrays are kept and left out of the work.
    for (unsigned i = 0; i < 5; i++) {
        f->length = 0;
        m->length = 0;
        for (unsigned k = 0; k < 5; k++) {
            mpz_set_ui(t, (k == i) ? 1 : 1000 + k);
            append_mpz_to_array(f, t);
            append_uint32_to_array(m, 7);
        }
        CHECK_EQ_UI(tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION),
                    COMPLETE_FACTORIZATION_FOUND, "filled arrays");
        CHECK(found_ok(n, f, m, 5), "filled arrays: factorization");
        CHECK_EQ_UI(mpz_get_ui(f->data[0]), (i == 0) ? 1 : 1000,
                    "filled arrays: kept");
    }
    free_mpz_array(f);
    free_uint32_array(m);
    mpz_clear(n);
    mpz_clear(t);
}

static void test_siqs_params(void) {
    mpz_t n;
    mpz_init_set_str(n, N128, 10);
    mpz_array_t*    f = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);

    // The last case needs more than 255 sieve chunks and bucket sieving.
    struct { uint32_t npb, hw, tdiv, nrel, extra; } cases[] = {
        { 0,    0,          0,    0,          0          },
        { 1,    1,          1,    1,          1          },
        { 400,  0xFFFFFFFF, 400,  0xFFFFFFFF, 0xFFFFFFFF },
        { 2500, 17000000,   2500, 24,         2          },
    };
    for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        siqs_params_t p;
        set_siqs_params_to_default(n, &p);
        p.nprimes_in_base  = cases[i].npb;
        p.sieve_half_width = cases[i].hw;
        p.nprimes_tdiv     = cases[i].tdiv;
        p.nrelations       = cases[i].nrel;
        p.extra_primes_in_a = cases[i].extra;
        f->length = 0;
        m->length = 0;
        if (siqs(f, m, n, &p, FIND_COMPLETE_FACTORIZATION)
            == COMPLETE_FACTORIZATION_FOUND) {
            CHECK(found_ok(n, f, m, 0), "siqs extreme parameters");
        }
    }
    // Numbers too small for the sieve's usual 'a', and n = k.m^2 for a
    // multiplier k, with kn a square.
    static const unsigned long small[] = { 15015, 3, 27, 75, 1009 };
    for (unsigned i = 0; i < sizeof(small) / sizeof(small[0]); i++) {
        mpz_set_ui(n, small[i]);
        siqs_params_t p;
        set_siqs_params_to_default(n, &p);
        f->length = 0;
        m->length = 0;
        if (siqs(f, m, n, &p, FIND_COMPLETE_FACTORIZATION)
            == COMPLETE_FACTORIZATION_FOUND) {
            CHECK(found_ok(n, f, m, 0), "siqs on a small n");
        }
    }

    free_mpz_array(f);
    free_uint32_array(m);
    mpz_clear(n);
}

static void test_cfrac_params(void) {
    mpz_t n;
    mpz_init_set_str(n, N40, 10);
    mpz_array_t*    f = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);

    struct { uint32_t npb, tdiv, nrel; } cases[] = {
        { 0, 0, 0 }, { 7, 192, 32 }, { 10, 65536, 32 }, { 65536, 65536, 0xFFFFFFFF }
    };
    for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        cfrac_params_t p;
        set_cfrac_params_to_default(n, &p);
        p.nprimes_in_base = cases[i].npb;
        p.nprimes_tdiv    = cases[i].tdiv;
        p.nrelations      = cases[i].nrel;
        f->length = 0;
        m->length = 0;
        if (cfrac(f, m, n, &p, FIND_COMPLETE_FACTORIZATION)
            == COMPLETE_FACTORIZATION_FOUND) {
            CHECK(found_ok(n, f, m, 0), "cfrac extreme parameters");
        }
    }
    free_mpz_array(f);
    free_uint32_array(m);
    mpz_clear(n);
}

static void test_bounds(void) {
    mpz_t n128, n300;
    mpz_init_set_str(n128, N128, 10);
    mpz_init_set_str(n300, N300, 10);
    mpz_array_t* f = alloc_mpz_array(4);

    struct { uint32_t b1, b2; } bounds[] = {
        { 0, 0 }, { 1, 1000 }, { 10, 50 }, { 50, 100 }, { 1000, 0xFFFFFFFF }
    };
    for (unsigned i = 0; i < sizeof(bounds) / sizeof(bounds[0]); i++) {
        const mpz_ptr ns[] = { n128, n300 };
        for (unsigned k = 0; k < 2; k++) {
            ecm_params_t e;
            set_ecm_params_to_default(ns[k], &e);
            e.b1      = bounds[i].b1;
            e.b2      = bounds[i].b2;
            e.ncurves = 1;
            f->length = 0;
            if (ecm(f, NULL, ns[k], &e, SINGLE_RUN) == SOME_FACTORS_FOUND) {
                CHECK(found_ok(ns[k], f, NULL, 0), "ecm extreme bounds");
            }
            p_minus_1_params_t pm;
            set_p_minus_1_params_to_default(ns[k], &pm);
            pm.b1 = bounds[i].b1;
            pm.b2 = bounds[i].b2;
            f->length = 0;
            if (p_minus_1(f, NULL, ns[k], &pm, SINGLE_RUN)
                == SOME_FACTORS_FOUND) {
                CHECK(found_ok(ns[k], f, NULL, 0), "p-1 extreme bounds");
            }
            p_plus_1_params_t pp;
            set_p_plus_1_params_to_default(ns[k], &pp);
            pp.b1 = bounds[i].b1;
            pp.b2 = bounds[i].b2;
            f->length = 0;
            if (p_plus_1(f, NULL, ns[k], &pp, SINGLE_RUN)
                == SOME_FACTORS_FOUND) {
                CHECK(found_ok(ns[k], f, NULL, 0), "p+1 extreme bounds");
            }
        }
    }
    free_mpz_array(f);
    mpz_clear(n128);
    mpz_clear(n300);
}

static void test_giving_up(void) {
    // A prime, which the modes other than SINGLE_RUN retry on: rho must give
    // up rather than retry forever.
    mpz_t n;
    mpz_init_set_ui(n, 1000003);
    mpz_array_t* f = alloc_mpz_array(4);

    rho_params_t r;
    set_rho_params_to_default(n, &r);
    CHECK(rho(f, NULL, n, &r, FIND_SOME_FACTORS) != SOME_FACTORS_FOUND,
          "rho on a prime");
    r.acc_niters = 0;
    CHECK(rho(f, NULL, n, &r, FIND_SOME_FACTORS) != SOME_FACTORS_FOUND,
          "rho without accumulation");

    free_mpz_array(f);
    mpz_clear(n);
}

int main(void) {
    TIFA_TEST_BEGIN("test_inputs");
    test_small_and_square_n();
    test_output_arrays();
    test_siqs_params();
    test_cfrac_params();
    test_bounds();
    test_giving_up();
    TIFA_TEST_END();
}
