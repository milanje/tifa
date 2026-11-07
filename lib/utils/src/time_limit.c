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

#include "time_limit.h"
#include "stopping.h"

typedef struct {
    bool                   limited;
    struct timespec        deadline;
    tifa_interrupt_check_t check;
    void*                  arg;
    bool                   interrupted;
    uint32_t               threads;
} stop_state_t;

static pthread_key_t  state_key;
static pthread_once_t state_once = PTHREAD_ONCE_INIT;
static bool           state_key_ok = false;

static void make_state_key(void) {
    state_key_ok = (0 == pthread_key_create(&state_key, free));
}

//
// The calling thread's state, created if asked for (NULL if it can not be).
//
static stop_state_t* state(const bool create) {
    pthread_once(&state_once, make_state_key);
    if (!state_key_ok) {
        return NULL;
    }
    stop_state_t* s = pthread_getspecific(state_key);
    if ((s == NULL) && create) {
        s = calloc(1, sizeof(stop_state_t));
        if ((s != NULL) && (pthread_setspecific(state_key, s) != 0)) {
            free(s);
            s = NULL;
        }
    }
    return s;
}

static bool passed(const struct timespec* const deadline) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec > deadline->tv_sec)
           || (   (now.tv_sec == deadline->tv_sec)
               && (now.tv_nsec >= deadline->tv_nsec));
}

//------------------------------------------------------------------------------
void tifa_set_time_limit(double seconds) {
    const bool limited = isfinite(seconds) && (seconds > 0.0);
    stop_state_t* const s = state(limited);
    if (s == NULL) {
        return;
    }
    s->limited = limited;
    if (!limited) {
        return;
    }
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    // Past about 30 years, as good as no limit, and no overflow of time_t.
    if (seconds > 1e9) {
        seconds = 1e9;
    }
    const double whole = floor(seconds);
    s->deadline.tv_sec  = now.tv_sec + (time_t) whole;
    s->deadline.tv_nsec = now.tv_nsec + (long) ((seconds - whole) * 1e9);
    if (s->deadline.tv_nsec >= 1000000000L) {
        s->deadline.tv_sec++;
        s->deadline.tv_nsec -= 1000000000L;
    }
}
//------------------------------------------------------------------------------
bool tifa_time_limit_reached(void) {
    if (!state_key_ok) {
        return false;
    }
    const stop_state_t* const s = pthread_getspecific(state_key);
    return (s != NULL) && s->limited && passed(&s->deadline);
}
//------------------------------------------------------------------------------
void tifa_set_interrupt_check(tifa_interrupt_check_t check, void* arg) {
    stop_state_t* const s = state(check != NULL);
    if (s == NULL) {
        return;
    }
    s->check       = check;
    s->arg         = arg;
    s->interrupted = false;
}
//------------------------------------------------------------------------------
bool tifa_interrupted(void) {
    if (!state_key_ok) {
        return false;
    }
    const stop_state_t* const s = pthread_getspecific(state_key);
    return (s != NULL) && s->interrupted;
}
//------------------------------------------------------------------------------
bool tifa_should_stop(void) {
    if (!state_key_ok) {
        return false;
    }
    stop_state_t* const s = pthread_getspecific(state_key);
    if (s == NULL) {
        return false;
    }
    if (s->interrupted || (s->limited && passed(&s->deadline))) {
        return true;
    }
    if ((s->check != NULL) && s->check(s->arg)) {
        s->interrupted = true;
    }
    return s->interrupted;
}
//------------------------------------------------------------------------------
void tifa_set_threads(uint32_t threads) {
    stop_state_t* const s = state(threads > 1);
    if (s != NULL) {
        s->threads = threads;
    }
}
//------------------------------------------------------------------------------
uint32_t tifa_threads(void) {
    if (!state_key_ok) {
        return 1;
    }
    const stop_state_t* const s = pthread_getspecific(state_key);
    return ((s == NULL) || (s->threads < 1)) ? 1 : s->threads;
}
