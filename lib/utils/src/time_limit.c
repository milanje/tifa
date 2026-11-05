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

typedef struct {
    bool            active;
    struct timespec deadline;
} time_limit_t;

static pthread_key_t  limit_key;
static pthread_once_t limit_once = PTHREAD_ONCE_INIT;
static bool           limit_key_ok = false;

static void make_limit_key(void) {
    limit_key_ok = (0 == pthread_key_create(&limit_key, free));
}

//------------------------------------------------------------------------------
void tifa_set_time_limit(double seconds) {
    pthread_once(&limit_once, make_limit_key);
    if (!limit_key_ok) {
        return;
    }
    time_limit_t* limit = pthread_getspecific(limit_key);
    if (!isfinite(seconds) || (seconds <= 0.0)) {
        if (limit != NULL) {
            limit->active = false;
        }
        return;
    }
    if (limit == NULL) {
        limit = malloc(sizeof(time_limit_t));
        if ((limit == NULL) || (pthread_setspecific(limit_key, limit) != 0)) {
            free(limit);
            return;
        }
    }
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    // Past about 30 years, as good as no limit, and no overflow of time_t.
    if (seconds > 1e9) {
        seconds = 1e9;
    }
    const double whole = floor(seconds);
    limit->deadline.tv_sec  = now.tv_sec + (time_t) whole;
    limit->deadline.tv_nsec = now.tv_nsec + (long) ((seconds - whole) * 1e9);
    if (limit->deadline.tv_nsec >= 1000000000L) {
        limit->deadline.tv_sec++;
        limit->deadline.tv_nsec -= 1000000000L;
    }
    limit->active = true;
}
//------------------------------------------------------------------------------
bool tifa_time_limit_reached(void) {
    if (!limit_key_ok) {
        return false;
    }
    const time_limit_t* const limit = pthread_getspecific(limit_key);
    if ((limit == NULL) || !limit->active) {
        return false;
    }
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec > limit->deadline.tv_sec)
           || (   (now.tv_sec == limit->deadline.tv_sec)
               && (now.tv_nsec >= limit->deadline.tv_nsec));
}
