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
 * \file    time_limit.c
 */

// clock_gettime, which strict C99 hides.
#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>

#include "stopping.h"

static pthread_key_t  state_key;
static pthread_once_t state_once = PTHREAD_ONCE_INIT;
static bool           state_key_ok = false;

static void make_state_key(void) {
    const bool ok = (0 == pthread_key_create(&state_key, free));
    __atomic_store_n(&state_key_ok, ok, __ATOMIC_RELEASE);
}
//
// Read without pthread_once by the threads that never set anything: the
// acquire makes the key they then use the one created.
//
static bool key_ok(void) {
    return __atomic_load_n(&state_key_ok, __ATOMIC_ACQUIRE);
}

//
// The calling thread's state, created if asked for (NULL if it can not be).
//
static tifa_stop_state_t* state(const bool create) {
    if (create) {
        pthread_once(&state_once, make_state_key);
    }
    if (!key_ok()) {
        return NULL;
    }
    tifa_stop_state_t* s = pthread_getspecific(state_key);
    if ((s == NULL) && create) {
        s = calloc(1, sizeof(tifa_stop_state_t));
        if ((s != NULL) && (pthread_setspecific(state_key, s) != 0)) {
            free(s);
            s = NULL;
        }
    }
    return s;
}

static bool passed(const tifa_stop_state_t* const s) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return ((int64_t) now.tv_sec > s->deadline_sec)
           || (   ((int64_t) now.tv_sec == s->deadline_sec)
               && ((int64_t) now.tv_nsec >= s->deadline_nsec));
}

static bool valid(const tifa_options_t* const options) {
    for (uint32_t i = 0; i < sizeof(options->reserved) / sizeof(uint64_t);
         i++) {
        if (options->reserved[i] != 0) {
            return false;
        }
    }
    return true;
}

static bool is_default(const tifa_stop_state_t* const s) {
    return    !s->limited && (s->check == NULL) && (s->threads <= 1)
           && (s->siqs_max_bits == 0) && (s->ecm_max_bits == 0);
}

static const tifa_options_t defaults;

//------------------------------------------------------------------------------
bool tifa_call_begin(tifa_call_t* const call,
                     const tifa_options_t* options) {
    call->installed = false;
    if (options == NULL) {
        //
        // The defaults, which the thread has unless within a call with
        // options (from its interrupt check): nothing to install then.
        //
        const tifa_stop_state_t* const s = key_ok()
                                         ? pthread_getspecific(state_key)
                                         : NULL;
        if ((s == NULL) || is_default(s)) {
            return true;
        }
        options = &defaults;
    }
    if (!valid(options)) {
        return false;
    }
    tifa_stop_state_t* const s = state(true);
    if (s == NULL) {
        //
        // Options that can not hold are refused rather than dropped.
        //
        return    !(isfinite(options->time_limit) && (options->time_limit > 0))
               && (options->interrupt == NULL) && (options->threads <= 1)
               && (options->siqs_max_bits == 0) && (options->ecm_max_bits == 0);
    }
    call->saved     = *s;
    call->installed = true;

    double seconds = options->time_limit;
    s->limited = isfinite(seconds) && (seconds > 0.0);
    if (s->limited) {
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        // Past about 30 years, as good as no limit, and no overflow.
        if (seconds > 1e9) {
            seconds = 1e9;
        }
        const double whole = floor(seconds);
        s->deadline_sec  = (int64_t) now.tv_sec + (int64_t) whole;
        s->deadline_nsec = (int64_t) now.tv_nsec
                         + (int64_t) ((seconds - whole) * 1e9);
        if (s->deadline_nsec >= 1000000000) {
            s->deadline_sec++;
            s->deadline_nsec -= 1000000000;
        }
    }
    s->check         = options->interrupt;
    s->arg           = options->interrupt_arg;
    s->interrupted   = false;
    s->threads       = options->threads;
    s->siqs_max_bits = options->siqs_max_bits;
    s->ecm_max_bits  = options->ecm_max_bits;
    return true;
}
//------------------------------------------------------------------------------
void tifa_call_end(tifa_call_t* const call) {
    if (!call->installed) {
        return;
    }
    tifa_stop_state_t* const s = state(false);
    if (s != NULL) {
        *s = call->saved;
    }
}
//------------------------------------------------------------------------------
bool tifa_should_stop(void) {
    if (!key_ok()) {
        return false;
    }
    tifa_stop_state_t* const s = pthread_getspecific(state_key);
    if (s == NULL) {
        return false;
    }
    if (s->interrupted || (s->limited && passed(s))) {
        return true;
    }
    if (s->check != NULL) {
        //
        // Not called again by the calls the check makes.
        //
        const tifa_interrupt_check_t check = s->check;
        s->check = NULL;
        const bool stop = check(s->arg);
        s->check = check;
        if (stop) {
            s->interrupted = true;
        }
    }
    return s->interrupted;
}
//------------------------------------------------------------------------------
static const tifa_stop_state_t* current(void) {
    if (!key_ok()) {
        return NULL;
    }
    return pthread_getspecific(state_key);
}
//------------------------------------------------------------------------------
uint32_t tifa_threads(void) {
    const tifa_stop_state_t* const s = current();
    return ((s == NULL) || (s->threads < 1)) ? 1 : s->threads;
}
//------------------------------------------------------------------------------
uint32_t tifa_siqs_max_bits(void) {
    const tifa_stop_state_t* const s = current();
    return ((s == NULL) || (s->siqs_max_bits == 0)) ? TIFA_SIQS_MAX_BITS
                                                    : s->siqs_max_bits;
}
//------------------------------------------------------------------------------
uint32_t tifa_ecm_max_bits(void) {
    const tifa_stop_state_t* const s = current();
    const uint32_t bits = ((s == NULL) || (s->ecm_max_bits == 0))
                        ? TIFA_ECM_MAX_BITS : s->ecm_max_bits;
    return (bits > TIFA_ECM_TOP_BITS) ? TIFA_ECM_TOP_BITS : bits;
}
//------------------------------------------------------------------------------
void tifa_set_thread_interrupt(tifa_interrupt_check_t check, void* arg) {
    tifa_stop_state_t* const s = state(true);
    if (s != NULL) {
        s->check       = check;
        s->arg         = arg;
        s->interrupted = false;
    }
}
