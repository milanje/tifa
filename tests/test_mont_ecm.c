// Fixed-limb Montgomery ECM (mont_ecm.h) against a plain GMP reference with
// the same curves: identical stage 1 outcomes for every sigma, and every
// factor the reference finds in stage 2 is found too. Moduli of 1 to 4 limbs.

#include <stdlib.h>

#include "tifa_test.h"
#include "first_primes.h"
#include "mont_ecm.h"
#include "z2z8_curves.h"

#if TIFA_HAVE_MONT_ECM

typedef struct {
    mpz_t x;
    mpz_t z;
} ref_pt_t;

static void ref_xdbl(ref_pt_t* r, const ref_pt_t* p, const mpz_t a24,
                     const mpz_t n) {
    mpz_t s, d, t;
    mpz_inits(s, d, t, NULL);
    mpz_add(s, p->x, p->z);
    mpz_mul(s, s, s);
    mpz_sub(d, p->x, p->z);
    mpz_mul(d, d, d);
    mpz_sub(t, s, d);
    mpz_mul(r->x, s, d);
    mpz_mod(r->x, r->x, n);
    mpz_mul(s, a24, t);
    mpz_add(s, s, d);
    mpz_mul(r->z, t, s);
    mpz_mod(r->z, r->z, n);
    mpz_clears(s, d, t, NULL);
}

static void ref_xadd(ref_pt_t* r, const ref_pt_t* p, const ref_pt_t* q,
                     const ref_pt_t* diff, const mpz_t n) {
    mpz_t u, v, a, b;
    mpz_inits(u, v, a, b, NULL);
    mpz_sub(a, p->x, p->z);
    mpz_add(b, q->x, q->z);
    mpz_mul(u, a, b);
    mpz_add(a, p->x, p->z);
    mpz_sub(b, q->x, q->z);
    mpz_mul(v, a, b);
    mpz_add(a, u, v);
    mpz_mul(a, a, a);
    mpz_sub(b, u, v);
    mpz_mul(b, b, b);
    mpz_mul(r->x, diff->z, a);
    mpz_mod(r->x, r->x, n);
    mpz_mul(r->z, diff->x, b);
    mpz_mod(r->z, r->z, n);
    mpz_clears(u, v, a, b, NULL);
}

static void ref_ladder(ref_pt_t* r, const ref_pt_t* p, const mpz_t k,
                       const mpz_t a24, const mpz_t n) {
    ref_pt_t r0, r1;
    mpz_init_set(r0.x, p->x);
    mpz_init_set(r0.z, p->z);
    mpz_inits(r1.x, r1.z, NULL);
    ref_xdbl(&r1, p, a24, n);
    for (long b = (long) mpz_sizeinbase(k, 2) - 2; b >= 0; b--) {
        if (mpz_tstbit(k, b)) {
            ref_xadd(&r0, &r1, &r0, p, n);
            ref_xdbl(&r1, &r1, a24, n);
        } else {
            ref_xadd(&r1, &r1, &r0, p, n);
            ref_xdbl(&r0, &r0, a24, n);
        }
    }
    mpz_set(r->x, r0.x);
    mpz_set(r->z, r0.z);
    mpz_clears(r0.x, r0.z, r1.x, r1.z, NULL);
}

// 1 if gcd(a, n) is a proper factor (stored in g), 0 otherwise.
static int proper_gcd(mpz_t g, const mpz_t a, const mpz_t n) {
    mpz_gcd(g, a, n);
    return (mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0);
}

// Reference ECM curve: same Suyama curve, stage 1 to b1, then (if b2 > b1)
// q.Q for every prime q in (b1, b2]. Returns 1 if a proper factor is found.
static int ref_curve(mpz_t g, const mpz_t n, uint32_t sigma, uint32_t b1,
                     uint32_t b2) {
    mpz_t u, v, a24, t, k;
    mpz_inits(u, v, a24, t, k, NULL);
    ref_pt_t p;
    mpz_inits(p.x, p.z, NULL);

    mpz_set_ui(u, sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5);
    mpz_set_ui(v, sigma);
    mpz_mul_ui(v, v, 4);
    mpz_powm_ui(p.x, u, 3, n);
    mpz_powm_ui(p.z, v, 3, n);
    mpz_mul(t, p.x, v);
    mpz_mul_ui(t, t, 16);
    int found = 0;
    if (!mpz_invert(t, t, n)) {
        mpz_mul(t, p.x, v);
        mpz_mul_ui(t, t, 16);
        found = proper_gcd(g, t, n);
        goto done;
    }
    mpz_sub(a24, v, u);
    mpz_powm_ui(a24, a24, 3, n);
    mpz_mul(a24, a24, t);
    mpz_mul_ui(t, u, 3);
    mpz_add(t, t, v);
    mpz_mul(a24, a24, t);
    mpz_mod(a24, a24, n);

    mpz_set_ui(k, 1);
    for (uint32_t i = 0; first_primes[i] <= b1; i++) {
        unsigned long q = first_primes[i];
        while (q <= b1 / first_primes[i]) {
            q *= first_primes[i];
        }
        mpz_mul_ui(k, k, q);
    }
    ref_ladder(&p, &p, k, a24, n);
    found = proper_gcd(g, p.z, n);
    mpz_gcd(t, p.z, n);
    if (found || (mpz_cmp(t, n) == 0) || (b2 <= b1)) {
        goto done;
    }
    ref_pt_t q;
    mpz_inits(q.x, q.z, NULL);
    mpz_set_ui(t, 1);
    for (uint32_t i = 0; first_primes[i] <= b2; i++) {
        if (first_primes[i] <= b1) {
            continue;
        }
        mpz_set_ui(k, first_primes[i]);
        ref_ladder(&q, &p, k, a24, n);
        mpz_mul(t, t, q.z);
        mpz_mod(t, t, n);
    }
    found = proper_gcd(g, t, n);
    mpz_clears(q.x, q.z, NULL);
done:
    mpz_clears(u, v, a24, t, k, p.x, p.z, NULL);
    return found;
}

// Stage 1 of the curve with starting point (x0 : 1) and a24 (residues), to
// b1: 1 if gcd(z, n) is a proper factor (in g).
static int ref_stage1(mpz_t g, const mpz_t n, const mpz_t x0, const mpz_t a24,
                      uint32_t b1) {
    ref_pt_t p;
    mpz_init_set(p.x, x0);
    mpz_init_set_ui(p.z, 1);
    mpz_t k;
    mpz_init_set_ui(k, 1);
    for (uint32_t i = 0; first_primes[i] <= b1; i++) {
        unsigned long q = first_primes[i];
        while (q <= b1 / first_primes[i]) {
            q *= first_primes[i];
        }
        mpz_mul_ui(k, k, q);
    }
    ref_ladder(&p, &p, k, a24, n);
    const int found = proper_gcd(g, p.z, n);
    mpz_clears(p.x, p.z, k, NULL);
    return found;
}

// a24 = 1/(1 - d), x0 = (1 + y0)/(1 - y0) for row i of z2z8_curves; 0 if
// not invertible mod n.
static int z2z8_ref(mpz_t x0, mpz_t a24, const mpz_t n, uint32_t i) {
    mpz_t d, y, t;
    mpz_inits(d, y, t, NULL);
    mpz_set_ui(d, 0);
    mpz_import(t, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][1]);
    mpz_import(d, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][0]);
    int ok = mpz_invert(t, t, n);
    mpz_mul(d, d, t);
    mpz_import(t, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][3]);
    mpz_import(y, 1, 1, sizeof(uint64_t), 0, 0, &z2z8_curves[i][2]);
    ok = ok && mpz_invert(t, t, n);
    mpz_mul(y, y, t);
    mpz_ui_sub(t, 1, d);
    ok = ok && mpz_invert(a24, t, n);
    mpz_ui_sub(t, 1, y);
    ok = ok && mpz_invert(t, t, n);
    mpz_add_ui(x0, y, 1);
    mpz_mul(x0, x0, t);
    mpz_mod(x0, x0, n);
    mpz_clears(d, y, t, NULL);
    return ok;
}

static void random_prime(mpz_t p, gmp_randstate_t state, unsigned long bits) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);
    } while (!mpz_probab_prime_p(p, 30));
}

#endif

int main(void) {
    TIFA_TEST_BEGIN("mont_ecm");

#if TIFA_HAVE_MONT_ECM
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261015UL);

    mpz_t n, p, q, g, gref;
    mpz_inits(n, p, q, g, gref, NULL);

    mont_ecm_plan_t* const plan1 = alloc_mont_ecm_plan(300, 0);
    mont_ecm_plan_t* const plan2 = alloc_mont_ecm_plan(200, 5000);
    char msg[160];

    for (uint32_t nl = 1; nl <= 4; nl++) {
        // n = p q with a 22-bit p and exactly nl limbs.
        random_prime(p, state, 22);
        do {
            random_prime(q, state, 64 * nl - 22 - (nl == 1 ? 1 : 5));
            mpz_mul(n, p, q);
        } while (mpz_size(n) != nl);

        uint32_t mismatches = 0;
        uint32_t bad_factor = 0;
        uint32_t found1     = 0;
        for (uint32_t sigma = 6; sigma < 6 + 120; sigma++) {
            const int mine = mont_ecm_curve(g, n, sigma, plan1);
            const int ref  = ref_curve(gref, n, sigma, 300, 0);
            if ((mine != ref) || (mine && (mpz_cmp(g, gref) != 0))) {
                mismatches++;
            }
            if (mine && !mpz_divisible_p(n, g)) {
                bad_factor++;
            }
            found1 += mine;
        }
        snprintf(msg, sizeof(msg), "%u limb(s): stage 1 matches reference",
                 nl);
        CHECK_EQ_UI(mismatches, 0, msg);
        snprintf(msg, sizeof(msg), "%u limb(s): stage 1 finds some", nl);
        CHECK(found1 > 0, msg);

        uint32_t missed = 0;
        uint32_t found2 = 0;
        for (uint32_t sigma = 6; sigma < 6 + 40; sigma++) {
            const int mine = mont_ecm_curve(g, n, sigma, plan2);
            const int ref  = ref_curve(gref, n, sigma, 200, 5000);
            if (ref && !mine) {
                missed++;
            }
            if (mine && !mpz_divisible_p(n, g)) {
                bad_factor++;
            }
            found2 += mine;
        }
        snprintf(msg, sizeof(msg), "%u limb(s): stage 2 finds what the "
                 "reference finds", nl);
        CHECK_EQ_UI(missed, 0, msg);
        snprintf(msg, sizeof(msg), "%u limb(s): factors divide n", nl);
        CHECK_EQ_UI(bad_factor, 0, msg);
        snprintf(msg, sizeof(msg), "%u limb(s): stage 2 finds more", nl);
        CHECK(found2 > found1 * 40 / 120, msg);
    }

    // The curve sequence of mont_ecm_nth_curve: the Z/2 x Z/8 curves match
    // the reference in stage 1, then come Suyama's curves.
    {
        mpz_t x0, a24;
        mpz_inits(x0, a24, NULL);
        for (uint32_t nl = 1; nl <= 4; nl++) {
            random_prime(p, state, 20);
            do {
                random_prime(q, state, 64 * nl - 20 - (nl == 1 ? 1 : 5));
                mpz_mul(n, p, q);
            } while (mpz_size(n) != nl);
            uint32_t mismatches = 0;
            uint32_t found = 0;
            for (uint32_t k = 0; k < NZ2Z8_CURVES; k++) {
                const int mine = mont_ecm_nth_curve(g, n, k, plan1);
                const int ref  = z2z8_ref(x0, a24, n, k)
                                 && ref_stage1(gref, n, x0, a24, 300);
                mismatches += (mine != ref) || (mine && mpz_cmp(g, gref));
                found += mine;
            }
            snprintf(msg, sizeof(msg), "%u limb(s): Z/2xZ/8 curves match "
                     "the reference", nl);
            CHECK_EQ_UI(mismatches, 0, msg);
            snprintf(msg, sizeof(msg), "%u limb(s): Z/2xZ/8 curves find "
                     "some", nl);
            CHECK(found > 0, msg);
            mismatches = 0;
            for (uint32_t k = NZ2Z8_CURVES; k < NZ2Z8_CURVES + 30; k++) {
                const int mine = mont_ecm_nth_curve(g, n, k, plan1);
                const int ref  = mont_ecm_curve(gref, n, 6 + k - NZ2Z8_CURVES,
                                                plan1);
                mismatches += (mine != ref) || (mine && mpz_cmp(g, gref));
            }
            snprintf(msg, sizeof(msg), "%u limb(s): then Suyama's curves",
                     nl);
            CHECK_EQ_UI(mismatches, 0, msg);
        }
        mpz_clears(x0, a24, NULL);
    }

    // A stage 2 reaching past the prime table (segmented sieve): a 36-bit
    // factor of a 3-limb number.
    mont_ecm_plan_t* const plan3 = alloc_mont_ecm_plan(2000, 900000);
    random_prime(p, state, 36);
    random_prime(q, state, 150);
    mpz_mul(n, p, q);
    int found = 0;
    for (uint32_t sigma = 6; (sigma < 6 + 300) && !found; sigma++) {
        found = mont_ecm_curve(g, n, sigma, plan3);
    }
    CHECK(found && (mpz_cmp(g, p) == 0 || mpz_cmp(g, q) == 0),
          "36-bit factor with b2 beyond the prime table");

    free_mont_ecm_plan(plan1);
    free_mont_ecm_plan(plan2);
    free_mont_ecm_plan(plan3);
    mpz_clears(n, p, q, g, gref, NULL);
    gmp_randclear(state);
#endif

    TIFA_TEST_END();
}
