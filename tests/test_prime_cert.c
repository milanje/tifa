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
#include <string.h>

#include "tifa_test.h"
#include "prime_cert.h"

static prime_cert_status_t prove_and_check(prime_cert_t* c, const mpz_t n,
                                           int* checked) {
    const prime_cert_status_t s = tifa_prove_prime(c, n, NULL);
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

// Each line of the text: "[[p, [" for the first step, " [p, [" then.
static bool text_matches(const prime_cert_t* c) {
    const size_t length = tifa_prime_cert_to_str(NULL, 0, c);
    char* const text = malloc(length + 1);
    bool ok = (tifa_prime_cert_to_str(text, length + 1, c) == length)
           && (strlen(text) == length);
    const char* line = text;
    for (uint32_t i = 0; ok && (i < c->nsteps); i++) {
        char start[200];
        gmp_snprintf(start, sizeof start, "%s[%Zd, [", (i == 0) ? "[" : " ",
                     c->steps[i].p);
        const char* const end = strchr(line, '\n');
        ok = (strncmp(line, start, strlen(start)) == 0)
          && ((end == NULL) == (i == c->nsteps - 1));
        line = (end == NULL) ? line : end + 1;
    }
    ok = ok && (strcmp(line + strlen(line) - 2, "]]") == 0);
    free(text);
    return ok;
}

// The text of c, read into d, is that of a certificate of n, and the same
// text again.
static bool reads_back(const prime_cert_t* c, prime_cert_t* d,
                       const mpz_t n) {
    const size_t length = tifa_prime_cert_to_str(NULL, 0, c);
    char* const text  = malloc(length + 1);
    char* const again = malloc(length + 1);
    tifa_prime_cert_to_str(text, length + 1, c);
    const bool ok = tifa_prime_cert_from_str(d, text)
                 && tifa_check_prime_cert(d, n)
                 && (tifa_prime_cert_to_str(again, length + 1, d) == length)
                 && (strcmp(text, again) == 0);
    free(text);
    free(again);
    return ok;
}

static char* text_of(const prime_cert_t* c) {
    const size_t length = tifa_prime_cert_to_str(NULL, 0, c);
    char* const text = malloc(length + 1);
    tifa_prime_cert_to_str(text, length + 1, c);
    return text;
}

// The steps of a, then those of b, as text: "[" a's "," b's "]".
static char* joined(const char* a, const char* b) {
    const size_t la = strlen(a);
    const size_t lb = strlen(b);
    char* const text = malloc(la + lb + 2);
    memcpy(text, a, la - 1);
    text[la - 1] = ',';
    memcpy(text + la, b + 1, lb);
    text[la + lb - 1] = '\0';
    return text;
}

static bool check_text(prime_cert_t* c, const char* text, const char* digits) {
    mpz_t n;
    mpz_init_set_str(n, digits, 10);
    const bool ok = tifa_prime_cert_from_str(c, text)
                 && tifa_check_prime_cert(c, n);
    mpz_clear(n);
    return ok;
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
        found =    (tifa_prove_prime(c, n, NULL) == PRIME_PROVEN)
                && (c->nsteps >= 2);
    }
    prime_cert_step_t* last = &c->steps[c->nsteps - 1];
    CHECK(tifa_check_prime_cert(c, n), "two-step certificate checked");

    const uint32_t a = last->pm1_witnesses->data[0];
    last->pm1_witnesses->data[0] = 1;
    CHECK(!tifa_check_prime_cert(c, n), "witness 1 refused");
    last->pm1_witnesses->data[0] = 0;
    CHECK(!tifa_check_prime_cert(c, n), "witness 0 refused");
    last->pm1_witnesses->data[0] = a;
    CHECK(tifa_check_prime_cert(c, n), "witness restored");

    mpz_add_ui(last->pm1_primes->data[0], last->pm1_primes->data[0], 2);
    CHECK(!tifa_check_prime_cert(c, n), "altered prime refused");
    mpz_sub_ui(last->pm1_primes->data[0], last->pm1_primes->data[0], 2);

    last->pm1_witnesses->length--;
    CHECK(!tifa_check_prime_cert(c, n), "missing witness refused");
    last->pm1_witnesses->length++;

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
    mpz_set_ui(last->pm1_primes->data[0], 2);
    last->pm1_primes->length = 1;
    last->pm1_witnesses->data[0] = 3;
    last->pm1_witnesses->length = 1;
    CHECK(!tifa_check_prime_cert(c, n), "certificate of 2^64 + 1 refused");

    // ---- Steps with primes of p + 1: searched among 200-bit primes.
    found = false;
    for (int tries = 0; !found && (tries < 400); tries++) {
        do {
            mpz_urandomb(n, state, 200);
            mpz_setbit(n, 199);
        } while (!mpz_probab_prime_p(n, 30));
        CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN,
              "200-bit prime proven");
        CHECK(tifa_check_prime_cert(c, n), "200-bit prime checked");
        last  = &c->steps[c->nsteps - 1];
        found = (last->pp1_primes->length != 0);
    }
    CHECK(found, "a step with primes of p + 1");
    if (found) {
        const uint32_t pw = last->pp1_witnesses->data[0];
        last->pp1_witnesses->data[0] = pw + 1;
        CHECK(!tifa_check_prime_cert(c, n), "even P refused");
        last->pp1_witnesses->data[0] = pw;

        const int32_t d = last->pp1_d;
        last->pp1_d = 0;
        CHECK(!tifa_check_prime_cert(c, n), "D = 0 refused");
        last->pp1_d = d;
        CHECK(tifa_check_prime_cert(c, n), "D restored");
        CHECK(mpz_si_kronecker(d, n) == -1, "D is a non-residue");

        mpz_add_ui(last->pp1_primes->data[0], last->pp1_primes->data[0], 2);
        CHECK(!tifa_check_prime_cert(c, n), "altered prime of p + 1 refused");
        mpz_sub_ui(last->pp1_primes->data[0], last->pp1_primes->data[0], 2);
    }

    // A step from p + 1 alone: 2^89 - 1, with p + 1 = 2^89. The checker
    // accepts some D and P, and those D are non-residues mod p.
    mpz_set_str(n, "618970019642690137449562111", 10);
    CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN, "2^89 - 1 proven");
    CHECK_EQ_UI(c->nsteps, 1, "2^89 - 1: one step");
    last = &c->steps[0];
    last->pm1_primes->length    = 0;
    last->pm1_witnesses->length = 0;
    reset_mpz_array(last->pp1_primes);
    reset_uint32_array(last->pp1_witnesses);
    mpz_set_ui(m, 2);
    append_mpz_to_array(last->pp1_primes, m);
    append_uint32_to_array(last->pp1_witnesses, 1);
    int accepted = 0;
    int wrong    = 0;
    for (int32_t a = 5; a < 60; a += 2) {
        for (uint32_t pp = 1; pp < 40; pp += 2) {
            last->pp1_d = ((a / 2) % 2 == 0) ? a : -a;
            last->pp1_witnesses->data[0] = pp;
            if (tifa_check_prime_cert(c, n)) {
                accepted++;
                wrong += (mpz_si_kronecker(last->pp1_d, n) != -1);
            }
        }
    }
    CHECK(accepted > 0, "2^89 - 1 from p + 1 alone");
    CHECK_EQ_UI(wrong, 0, "only non-residues accepted");
    mpz_set_str(n, "618970019642690137449562113", 10);
    mpz_set(last->p, n);
    CHECK(!tifa_check_prime_cert(c, n), "2^89 + 1 refused");

    // ---- Text.
    char text[200];
    mpz_set_ui(n, 1000003);
    CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN, "1000003 proven");
    CHECK_EQ_UI(tifa_prime_cert_to_str(text, sizeof text, c), 2,
                "text length, no step");
    CHECK(strcmp(text, "[]") == 0, "text, no step");

    mpz_set_str(n, "618970019642690137449562111", 10);
    CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN, "2^89 - 1 proven");
    last = &c->steps[0];
    reset_mpz_array(last->pm1_primes);
    reset_uint32_array(last->pm1_witnesses);
    reset_mpz_array(last->pp1_primes);
    reset_uint32_array(last->pp1_witnesses);
    mpz_set_ui(m, 2);
    append_mpz_to_array(last->pm1_primes, m);
    append_uint32_to_array(last->pm1_witnesses, 3);
    mpz_set_ui(m, 3);
    append_mpz_to_array(last->pm1_primes, m);
    append_uint32_to_array(last->pm1_witnesses, 5);
    mpz_set_ui(m, 2);
    append_mpz_to_array(last->pp1_primes, m);
    append_uint32_to_array(last->pp1_witnesses, 1);
    last->pp1_d = -7;
    const char* const expected =
        "[[618970019642690137449562111, [[2, 3], [3, 5]], -7, [[2, 1]]]]";
    CHECK_EQ_UI(tifa_prime_cert_to_str(text, sizeof text, c),
                strlen(expected), "text length, one step");
    CHECK(strcmp(text, expected) == 0, "text, one step");
    CHECK_EQ_UI(tifa_prime_cert_to_str(NULL, 0, c), strlen(expected),
                "text length, no room");
    memset(text, 'x', sizeof text);
    CHECK_EQ_UI(tifa_prime_cert_to_str(text, 10, c), strlen(expected),
                "text length, truncated");
    CHECK((strlen(text) == 9) && (strncmp(text, expected, 9) == 0)
          && (text[10] == 'x'), "text truncated");

    prime_cert_t* const d = alloc_prime_cert();
    uint32_t most = 0;
    bool with_pp1 = false;
    bool matched  = true;
    bool read     = true;
    for (int i = 0; i < 40; i++) {
        do {
            mpz_urandomb(n, state, 200);
            mpz_setbit(n, 199);
        } while (!mpz_probab_prime_p(n, 30));
        CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN,
              "200-bit prime proven");
        matched = matched && text_matches(c);
        read    = read && reads_back(c, d, n);
        most = (c->nsteps > most) ? c->nsteps : most;
        for (uint32_t j = 0; j < c->nsteps; j++) {
            with_pp1 = with_pp1 || (c->steps[j].pp1_d != 0);
        }
    }
    CHECK(matched, "text of 200-bit primes, a line per step");
    CHECK(read, "text of 200-bit primes read back and checked");
    CHECK((most >= 2) && with_pp1, "texts with several steps, and p + 1");

    // ---- Reading text.
    mpz_set_str(n, "5704689200685129054721", 10);
    CHECK(tifa_prime_cert_from_str(d, "[[5704689200685129054721, [[2, 7], "
                                      "[3, 3], [5, 3], [12497, 3]]]]")
          && tifa_check_prime_cert(d, n), "certificate of 2^128 + 1's factor");
    CHECK(tifa_prime_cert_from_str(d, " [ ]\n") && (d->nsteps == 0),
          "no step");
    CHECK(tifa_prime_cert_from_str(d, " [ [ 618970019642690137449562111 ,\t"
                                      "[ ] , -7 ,\n[ [ 2 , 1 ] ] ] ]\r\n"),
          "spaces, tabs and newlines");
    CHECK_EQ_UI(tifa_prime_cert_to_str(text, sizeof text, d), 49,
                "read: length");
    CHECK(strcmp(text, "[[618970019642690137449562111, [], -7, [[2, 1]]]]")
          == 0, "read: text");
    CHECK(tifa_prime_cert_from_str(d, "[[5, [[2, 4294967295]], 2147483647, "
                                      "[]], [7, [], -2147483647, [[2, 0]]]]")
          && (d->nsteps == 2)
          && (d->steps[0].pm1_witnesses->data[0] == UINT32_MAX)
          && (d->steps[0].pp1_d == INT32_MAX)
          && (d->steps[1].pp1_d == -INT32_MAX)
          && (d->steps[1].pp1_witnesses->data[0] == 0), "bounds");

    static const char* const malformed[] = {
        "", "x", "[", "]", "[,]", "[[]]", "[5]", "[[5]]", "[[5, []",
        "[[5, [[2]]]]", "[[5, [[2, 3, 4]]]]", "[[5, [[2, 3]], 7]]",
        "[[5, [[2, 3]], 0, [[2, 1]]]]", "[[5, [[2, 3]], -7, [[2, 1]], 1]]",
        "[[5, [[2, 3]]]] x", "[[5, [[2, 3]]]][]", "[[5, [[2, 3.0]]]]",
        "[[5, [[2, -3]]]]", "[[-5, [[2, 3]]]]", "[[+5, [[2, 3]]]]",
        "[[5, [[-2, 3]]]]", "[[5, [[2, 4294967296]]]]",
        "[[5, [[2, 99999999999999999999999]]]]", "[[5, [], 2147483648, []]]",
        "[[5, [], -2147483648, []]]", "[[5, [[2, 3]]],]", "[[5, [[2, 3],]]]",
        "[[5, [[2, 3]]] [7, []]]", "[[5, [[2, 3]], -, []]]",
    };
    const uint32_t nmalformed = sizeof malformed / sizeof malformed[0];
    uint32_t refused = 0;
    for (uint32_t i = 0; i < nmalformed; i++) {
        tifa_prime_cert_from_str(d, "[[5, [[2, 3]]]]");
        if (!tifa_prime_cert_from_str(d, malformed[i]) && (d->nsteps == 0)) {
            refused++;
        } else {
            printf("    accepted: \"%s\"\n", malformed[i]);
        }
    }
    CHECK_EQ_UI(refused, nmalformed, "malformed texts refused, no steps left");

    // ---- Hostile certificates: refused, without crash or long checks.
    static const char* const zeros[] = {
        "[[5, [[0, 2]]]]", "[[5, [], 3, [[0, 1]]]]",
        "[[5, [[2, 2]], 3, [[0, 1]]]]", "[[5, [[1, 2]]]]",
        "[[5, [[2, 2], [1, 2]]]]", "[[7, [[2, 3]], 3, [[1, 1]]]]",
    };
    bool none = true;
    for (uint32_t i = 0; i < sizeof zeros / sizeof zeros[0]; i++) {
        none = none && !check_text(d, zeros[i], (zeros[i][3] == '7') ? "7" : "5");
    }
    CHECK(none, "primes 0 and 1 refused");

    mpz_set_str(n, "170141183460469231731687303715884105727", 10);
    CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN, "2^127 - 1 proven");
    char* const m127 = text_of(c);
    mpz_set_str(n, "618970019642690137449562111", 10);
    CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN, "2^89 - 1 proven");
    char* const m89 = text_of(c);
    char* const big_first   = joined(m127, m89);
    char* const small_first = joined(m89, m127);
    CHECK(!check_text(d, big_first, "618970019642690137449562111"),
          "step for a prime above n refused");
    CHECK(check_text(d, small_first,
                     "170141183460469231731687303715884105727"),
          "unused step for a prime below n accepted");
    char* const twice = joined(m89, m89);
    CHECK(tifa_prime_cert_from_str(d, twice) && (d->nsteps == 2)
          && !check_text(d, twice, "618970019642690137449562111"),
          "two steps for the same prime refused");

    // 20000 steps, all for 2^89 - 1, before that of 2^127 - 1.
    const size_t l89  = strlen(m89) - 2;
    const size_t l127 = strlen(m127) - 2;
    char* const many = malloc(19999 * (l89 + 1) + l127 + 3);
    char* end = many;
    *end++ = '[';
    for (int i = 0; i < 19999; i++) {
        memcpy(end, m89 + 1, l89);
        end += l89;
        *end++ = ',';
    }
    memcpy(end, m127 + 1, l127);
    end += l127;
    *end++ = ']';
    *end   = '\0';
    CHECK(tifa_prime_cert_from_str(d, many) && (d->nsteps == 20000)
          && !check_text(d, many, "170141183460469231731687303715884105727"),
          "20000 steps for one prime refused");
    free(many);
    free(m127);
    free(m89);
    free(big_first);
    free(small_first);
    free(twice);

    // Arrays of different lengths, or primes of p + 1 with D = 0: refused,
    // and their text too.
    CHECK(tifa_prime_cert_from_str(d, "[[618970019642690137449562111, "
                                      "[], -7, [[2, 1]]]]"), "a step to alter");
    append_uint32_to_array(d->steps[0].pp1_witnesses, 3);
    char* const longer = text_of(d);
    CHECK(strcmp(longer, "[[618970019642690137449562111, [], -7, "
                         "[[2, 1], [3]]]]") == 0, "text, a witness too many");
    d->steps[0].pp1_witnesses->length = 0;
    char* const shorter = text_of(d);
    CHECK(strcmp(shorter, "[[618970019642690137449562111, [], -7, [[2]]]]")
          == 0, "text, a witness too few");
    d->steps[0].pp1_witnesses->length = 1;
    d->steps[0].pp1_d = 0;
    char* const no_d = text_of(d);
    CHECK(strcmp(no_d, "[[618970019642690137449562111, [], 0, [[2, 1]]]]")
          == 0, "text, primes of p + 1 and D = 0");
    CHECK(!tifa_prime_cert_from_str(c, longer)
          && !tifa_prime_cert_from_str(c, shorter)
          && !tifa_prime_cert_from_str(c, no_d), "their texts refused");
    free(longer);
    free(shorter);
    free(no_d);
    free_prime_cert(d);

    // ---- Proving again replaces the certificate.
    mpz_set_str(n, "618970019642690137449562111", 10);
    CHECK(tifa_prove_prime(c, n, NULL) == PRIME_PROVEN,
          "2^89 - 1 proven again");
    mpz_set_str(m, "170141183460469231731687303715884105727", 10);
    CHECK(tifa_prove_prime(c, m, NULL) == PRIME_PROVEN,
          "2^127 - 1 proven again");
    CHECK(!tifa_check_prime_cert(c, n) && tifa_check_prime_cert(c, m),
          "certificate replaced");

    gmp_randclear(state);
    mpz_clears(n, m, NULL);
    free_prime_cert(c);

    TIFA_TEST_END();
}
