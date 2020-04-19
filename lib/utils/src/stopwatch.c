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
 * \file    stopwatch.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief A very basic stopwatch-like timer.
 *
 * This file implements a very basic stopwatch-like timer (with microsecond
 * precision) based on the \c rusage structure and using the \c getrusage
 * function.
 */

#include <stdlib.h>
#include "stopwatch.h"

//----------------------------------------------------------------------------
void init_stopwatch(stopwatch_t* const watch) {    
    watch->elapsed_usec = 0;
    watch->is_running   = false;
}
//----------------------------------------------------------------------------
void start_stopwatch(stopwatch_t* const watch) {    
    if (! watch->is_running) {
        getrusage(RUSAGE_SELF, watch->rsg);
        
        watch->started_usec  = (uint64_t) watch->rsg->ru_utime.tv_sec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_sec;
        watch->started_usec *= 1000000;
        watch->started_usec += (uint64_t) watch->rsg->ru_utime.tv_usec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_usec;
        
        watch->is_running = true;
    }
}
//----------------------------------------------------------------------------
void stop_stopwatch(stopwatch_t* const watch) {
    if (watch->is_running) {
        getrusage(RUSAGE_SELF, watch->rsg);
        
        watch->elapsed_usec += (   (uint64_t) watch->rsg->ru_utime.tv_sec
                                 + (uint64_t) watch->rsg->ru_stime.tv_sec
                               ) * 1000000;                                         
        watch->elapsed_usec += (uint64_t) watch->rsg->ru_utime.tv_usec;
        watch->elapsed_usec += (uint64_t) watch->rsg->ru_stime.tv_usec;
        watch->elapsed_usec -= watch->started_usec;
        
        watch->is_running  = false;
    }
}
//----------------------------------------------------------------------------
void reset_stopwatch(stopwatch_t* const watch) {
    watch->elapsed_usec = 0;
    if (watch->is_running) {
        getrusage(RUSAGE_SELF, watch->rsg);
        
        watch->started_usec  = (uint64_t) watch->rsg->ru_utime.tv_sec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_sec;
        watch->started_usec *= 1000000;        
        watch->started_usec += (uint64_t) watch->rsg->ru_utime.tv_usec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_usec;
    }
}
//----------------------------------------------------------------------------
double get_stopwatch_elapsed(stopwatch_t* const watch) {
    return  (watch->elapsed_usec / 1000000.0);
}
//----------------------------------------------------------------------------
