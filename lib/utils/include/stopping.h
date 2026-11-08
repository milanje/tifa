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
 * \file    stopping.h
 * \brief   Whether a thread's factorizations should stop (internal).
 */

#if !defined(_TIFA_STOPPING_H_)
#define _TIFA_STOPPING_H_

#include <stdbool.h>
#include <stdint.h>

#include "tifa_options.h"

//
// What the options of the calling thread's current call became: zero is the
// defaults (no call, or a call without options).
//
typedef struct {
    bool                   limited;
    int64_t                deadline_sec;
    int64_t                deadline_nsec;
    tifa_interrupt_check_t check;
    void*                  arg;
    bool                   interrupted;
    uint32_t               threads;
    uint32_t               siqs_max_bits;
    uint32_t               ecm_max_bits;
} tifa_stop_state_t;

typedef struct {
    tifa_stop_state_t saved;
    bool              installed;
} tifa_call_t;

//
// Around the work of a public function: its options hold until
// tifa_call_end, which gives back those of the call it was made from. NULL
// keeps those, and costs nothing. Returns false, installing nothing, if the
// options are invalid (reserved fields not zero, unknown hint).
//
bool tifa_call_begin(tifa_call_t* call, const tifa_options_t* options);
void tifa_call_end(tifa_call_t* call);

   /**
    * \brief Whether the calling thread's factorizations should give up.
    *
    * \return true if the current call's time limit has passed or its
    *         interrupt check has returned true (see tifa_options.h).
    */
bool tifa_should_stop(void);
//
// The current call's options, with their defaults.
//
uint32_t tifa_threads(void);
uint32_t tifa_siqs_max_bits(void);
uint32_t tifa_ecm_max_bits(void);
//
// For the threads a call starts: an interrupt check of their own, which
// holds until they end.
//
void tifa_set_thread_interrupt(tifa_interrupt_check_t check, void* arg);

#endif
