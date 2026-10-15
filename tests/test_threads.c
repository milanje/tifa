// Concurrent factorizations, one per thread, as in a cofactorization farm:
// tifa_factor, tifa_cofactor, and the algorithms called directly, with
// parameter structs and cofactorization strategy shared by all the threads. Phase 1 starts on an empty tuning cache, so the
// threads race to load it; in phase 2 another thread keeps resetting it.
// Built with -fsanitize=thread, this catches the data races; otherwise it
// checks that every result is right.

#include <pthread.h>
#include <sched.h>
#include <stdlib.h>
#include <string.h>

#include "tifa_test.h"
#include "ecm.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "siqs.h"
#include "cofactor.h"
#include "tifa_factor.h"
#include "tuning.h"

#if !defined(TIFA_TEST_FIXTURE_DIR)
    #error "TIFA_TEST_FIXTURE_DIR must be defined by the build"
#endif

#define NTHREADS 8

typedef struct {
    ecm_params_t       ecm;
    p_minus_1_params_t pm1;
    p_plus_1_params_t  pp1;
    siqs_params_t      siqs;
    cofactor_strategy_t* cofactor;
} shared_params_t;

typedef struct {
    unsigned int           index;
    const shared_params_t* shared;
    uint32_t               pm1_b1;
    uint32_t               nrun;
    uint32_t               nfailed;
} worker_t;

static pthread_mutex_t toggle_lock = PTHREAD_MUTEX_INITIALIZER;
static int             workers_done = 0;

static void random_prime(mpz_t p, gmp_randstate_t state, unsigned long bits) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);
    } while (!mpz_probab_prime_p(p, 30));
}

typedef ecode_t (*run_t)(mpz_array_t*, uint32_array_t*, const mpz_t,
                         const shared_params_t*);

static ecode_t run_tifa_factor(mpz_array_t* f, uint32_array_t* m,
                               const mpz_t n, const shared_params_t* s) {
    (void) s;
    return tifa_factor(f, m, n, FIND_COMPLETE_FACTORIZATION);
}

static ecode_t run_ecm(mpz_array_t* f, uint32_array_t* m, const mpz_t n,
                       const shared_params_t* s) {
    return ecm(f, m, n, &s->ecm, FIND_COMPLETE_FACTORIZATION);
}

static ecode_t run_pm1(mpz_array_t* f, uint32_array_t* m, const mpz_t n,
                       const shared_params_t* s) {
    return p_minus_1(f, m, n, &s->pm1, FIND_COMPLETE_FACTORIZATION);
}

static ecode_t run_pp1(mpz_array_t* f, uint32_array_t* m, const mpz_t n,
                       const shared_params_t* s) {
    return p_plus_1(f, m, n, &s->pp1, FIND_COMPLETE_FACTORIZATION);
}

static ecode_t run_siqs(mpz_array_t* f, uint32_array_t* m, const mpz_t n,
                        const shared_params_t* s) {
    return siqs(f, m, n, &s->siqs, FIND_COMPLETE_FACTORIZATION);
}

static void factor_and_check(worker_t* w, run_t run, const mpz_t n) {
    mpz_array_t*    f = alloc_mpz_array(8);
    uint32_array_t* m = alloc_uint32_array(8);
    const ecode_t rc = run(f, m, n, w->shared);
    if (   (rc != COMPLETE_FACTORIZATION_FOUND)
        || !tifa_is_complete_factorization(n, f, m)) {
        w->nfailed++;
        gmp_fprintf(stderr, "  thread %u failed on %Zd\n", w->index, n);
    }
    w->nrun++;
    free_mpz_array(f);
    free_uint32_array(m);
}

static void* worker(void* arg) {
    worker_t* const w = arg;
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20261015UL + w->index);
    mpz_t n, p, q;
    mpz_inits(n, p, q, NULL);

    random_prime(p, state, 32);
    random_prime(q, state, 32);
    mpz_mul(n, p, q);
    p_minus_1_params_t pm1;
    set_p_minus_1_params_to_default(n, &pm1);
    w->pm1_b1 = pm1.b1;

    static const unsigned long sizes[] = { 40, 64, 80, 100, 120, 140 };
    for (uint32_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
        random_prime(p, state, sizes[i] / 2);
        random_prime(q, state, sizes[i] - sizes[i] / 2);
        mpz_mul(n, p, q);
        factor_and_check(w, run_tifa_factor, n);
    }
    // p^2.q.r: repeated factors.
    random_prime(p, state, 30);
    random_prime(q, state, 40);
    mpz_mul(n, p, p);
    mpz_mul(n, n, q);
    random_prime(q, state, 50);
    mpz_mul(n, n, q);
    factor_and_check(w, run_tifa_factor, n);

    // A 24-bit factor, the size the shared ECM parameters are set for.
    random_prime(p, state, 24);
    random_prime(q, state, 100);
    mpz_mul(n, p, q);
    factor_and_check(w, run_ecm, n);

    // p - 1 = 2.5003 and p + 1 = 4.1013: caught in stage 2.
    random_prime(q, state, 100);
    mpz_mul_ui(n, q, 10007);
    factor_and_check(w, run_pm1, n);
    mpz_mul_ui(n, q, 4051);
    factor_and_check(w, run_pp1, n);

    random_prime(p, state, 50);
    random_prime(q, state, 50);
    mpz_mul(n, p, q);
    factor_and_check(w, run_siqs, n);

    // Three large primes, above 64 bits.
    mpz_array_t*    f = alloc_mpz_array(4);
    uint32_array_t* m = alloc_uint32_array(4);
    for (uint32_t i = 0; i < 3; i++) {
        random_prime(p, state, 28);
        random_prime(q, state, 30);
        mpz_mul(n, p, q);
        random_prime(p, state, 31);
        mpz_mul(n, n, p);
        if (tifa_cofactor(f, m, n, w->shared->cofactor) != COFACTOR_SMOOTH) {
            w->nfailed++;
            gmp_fprintf(stderr, "  thread %u: %Zd not smooth\n", w->index, n);
        }
        w->nrun++;
    }
    free_mpz_array(f);
    free_uint32_array(m);

    mpz_clears(n, p, q, NULL);
    gmp_randclear(state);
    return NULL;
}

static void* toggler(void* arg) {
    (void) arg;
    for (uint32_t i = 0; ; i++) {
        pthread_mutex_lock(&toggle_lock);
        const int done = workers_done;
        pthread_mutex_unlock(&toggle_lock);
        if (done) {
            break;
        }
        tifa_tuning_set_dir((i & 1) ? NULL : TIFA_TEST_FIXTURE_DIR);
        sched_yield();
    }
    return NULL;
}

static void run_workers(const shared_params_t* shared, worker_t* workers,
                        bool with_toggler) {
    pthread_t threads[NTHREADS];
    pthread_t toggle_thread;
    workers_done = 0;
    if (with_toggler) {
        pthread_create(&toggle_thread, NULL, toggler, NULL);
    }
    for (unsigned int i = 0; i < NTHREADS; i++) {
        workers[i] = (worker_t) { .index = i, .shared = shared };
        pthread_create(&threads[i], NULL, worker, &workers[i]);
    }
    for (unsigned int i = 0; i < NTHREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    if (with_toggler) {
        pthread_mutex_lock(&toggle_lock);
        workers_done = 1;
        pthread_mutex_unlock(&toggle_lock);
        pthread_join(toggle_thread, NULL);
    }
}

int main(void) {
    TIFA_TEST_BEGIN("threads");

    shared_params_t shared;
    memset(&shared, 0, sizeof(shared));
    mpz_t n;
    mpz_init(n);
    mpz_ui_pow_ui(n, 2, 63);

    set_ecm_params_for_factor_size(24, &shared.ecm);
    set_p_minus_1_params_to_default(n, &shared.pm1);
    shared.pm1.b1 = 1000;
    shared.pm1.b2 = 10000;
    set_p_plus_1_params_to_default(n, &shared.pp1);
    shared.pp1.b1 = 1000;
    shared.pp1.b2 = 10000;
    mpz_ui_pow_ui(n, 2, 99);
    set_siqs_params_to_default(n, &shared.siqs);
    const cofactor_params_t cofactor = { .lpb = 31, .fb_bound = 1UL << 22,
                                         .max_large_primes = 3 };
    shared.cofactor = alloc_cofactor_strategy(&cofactor);
    shared_params_t before;
    memcpy(&before, &shared, sizeof(shared));
    mpz_clear(n);

    worker_t workers[NTHREADS];
    char msg[160];

    tifa_tuning_set_dir(TIFA_TEST_FIXTURE_DIR);
    run_workers(&shared, workers, false);
    for (unsigned int i = 0; i < NTHREADS; i++) {
        snprintf(msg, sizeof(msg), "phase 1, thread %u: all factored", i);
        CHECK_EQ_UI(workers[i].nfailed, 0, msg);
        snprintf(msg, sizeof(msg), "phase 1, thread %u: p-1 b1 from the "
                 "fixture tuning file", i);
        CHECK_EQ_UI(workers[i].pm1_b1, 12000, msg);
    }

    run_workers(&shared, workers, true);
    for (unsigned int i = 0; i < NTHREADS; i++) {
        snprintf(msg, sizeof(msg), "phase 2, thread %u: all factored", i);
        CHECK_EQ_UI(workers[i].nfailed, 0, msg);
    }
    tifa_tuning_set_dir(NULL);

    CHECK(memcmp(&shared, &before, sizeof(shared)) == 0,
          "shared parameters left untouched");
    free_cofactor_strategy(shared.cofactor);

    TIFA_TEST_END();
}
