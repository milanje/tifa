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
 * \file    timer.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief This file defines some macros used to perform timing measurements.
 *
 * \warning The \c __TIMING__ symbol should be defined in the file including
 * this header. It is under the responsability of the including file to check
 * for its definition or to define it, if needed.
 *
 * If \c __TIMING__ is set to 0, then the macros defined in this file do
 * nothing.
 */

#if !defined(_TIFA_TIMER_H_)
   /**
    * \def _TIFA_TIMER_H_
    * Standard include guard.
    */
#define _TIFA_TIMER_H_

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__TIMING__)
#error "__TIMING__ must be defined before including this header."
#endif

#if __TIMING__
    #include "stopwatch.h"
    
       /**
        * \def TIMING_FORMAT
        * Format used to print timing measurements.
        */
    #define TIMING_FORMAT "%8.3f"

       /**
        * \def INIT_NAMED_TIMER(NAME)
        * Initialize a timer named NAME.
        *
        * \warning The timer name should not be enclosed in quotes, e.g.
        * <tt>INIT_NAMED_TIMER(my_timer)</tt> is correct, but
        * <tt>INIT_NAMED_TIMER("my_timer")</tt> is wrong and will result in a
        * compilation error.
        *
        * This warning holds for all of the <tt>*_NAMED_TIMER</tt> macros.
        */
    #define INIT_NAMED_TIMER(NAME)                  \
        stopwatch_t __TIFA_STOPWATCH_ ##NAME;       \
        init_stopwatch(&__TIFA_STOPWATCH_ ##NAME);

       /**
        * \def RESET_NAMED_TIMER(NAME)
        * Reset the timer named NAME to zero.
        */
    #define RESET_NAMED_TIMER(NAME)                     \
        do {                                            \
            reset_stopwatch(&__TIFA_STOPWATCH_ ##NAME); \
        } while (0)

       /**
        * \def START_NAMED_TIMER(NAME)
        * Start the timer named NAME.
        *
        * \note Consecutive "calls" to START_NAMED_TIMER are without effect.
        */
    #define START_NAMED_TIMER(NAME)                     \
        do {                                            \
            start_stopwatch(&__TIFA_STOPWATCH_ ##NAME); \
        } while (0)

       /**
        * \def STOP_NAMED_TIMER(NAME)
        * Stop the timer named NAME. Successive "calls" to START_NAMED_TIMER
        * and STOP_NAMED_TIMER are cumulative. In other words, the timer's state
        * holds the time elapsed during all previous time intervals defined by
        * a "call" to START_NAMED_TIMER followed by a "call" to
        * STOP_NAMED_TIMER (provided that the timer was not reset via
        * RESET_NAMED_TIMER).
        *
        * \note Consecutive "calls" to STOP_NAMED_TIMER are without side
        * effects.
        */
    #define STOP_NAMED_TIMER(NAME)                      \
        do {                                            \
            stop_stopwatch(&__TIFA_STOPWATCH_ ##NAME);  \
        } while (0)

       /**
        * \def GET_NAMED_TIMING(NAME)
        * Return the timing of the timer named NAME as seconds.
        *
        * \warning The returned result is only meaningful if the timer is
        * not running (i.e. it has been stopped via STOP_NAMED_TIMER).
        */
    #define GET_NAMED_TIMING(NAME)                          \
        get_stopwatch_elapsed(&__TIFA_STOPWATCH_ ##NAME)
        
       /**
        * \def INIT_TIMER
        * Initialize an unamed timer.
        */
    #define INIT_TIMER      INIT_NAMED_TIMER()
    
       /**
        * \def RESET_TIMER
        * Reset an unamed timer.
        */
    #define RESET_TIMER     RESET_NAMED_TIMER()
    
       /**
        * \def START_TIMER
        * Start an unamed timer.
        */
    #define START_TIMER     START_NAMED_TIMER()
    
       /**
        * \def STOP_TIMER
        * Stop an unamed timer.
        */
    #define STOP_TIMER      STOP_NAMED_TIMER()
    
       /**
        * \def GET_TIMING
        * Get timing from an unamed timer.
        */
    #define GET_TIMING      GET_NAMED_TIMING()
#else
    #define INIT_NAMED_TIMER(NAME)    /* intentionally left empty */
    #define RESET_NAMED_TIMER(NAME)   /* intentionally left empty */
    #define START_NAMED_TIMER(NAME)   /* intentionally left empty */
    #define STOP_NAMED_TIMER(NAME)    /* intentionally left empty */
    #define GET_NAMED_TIMING(NAME)    (-1.0)

    #define INIT_TIMER    /* intentionally left empty */
    #define RESET_TIMER   /* intentionally left empty */
    #define START_TIMER   /* intentionally left empty */
    #define STOP_TIMER    /* intentionally left empty */
    #define GET_TIMING    (-1.0)
#endif

#ifdef __cplusplus
}
#endif

#endif
