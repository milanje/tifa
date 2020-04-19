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
 * \file    print_error.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Error printing macro
 *
 * This file defines a macro used to output critical error messages on stderr
 * if the global symbol \c TIFA_PRINT_ERROR is set to non-zero.
 *
 */

#if !defined(_TIFA_PRINT_ERROR_H_)
   /**
    * \def _TIFA_PRINT_ERROR_H_
    * Standard include guard.
    */
#define _TIFA_PRINT_ERROR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "tifa_config.h"

#if !defined(__PREFIX__)
    #define __PREFIX__ ""
#endif

#if !defined(__VERBOSE__)
    #define __VERBOSE__ 0
#endif

//
// No error messages are printed if TIFA_PRINT_ERROR is set to 0.
//
#if TIFA_PRINT_ERROR
    #if __VERBOSE__
            /**
             * \def _TIFA_ERR_PRFX_
             * String prefixed to all error messages.
             */
        #define _TIFA_ERR_PRFX_ "ERROR: "
        
            /**
             * \def PRINTF_STDERR(...)
             * Print on the standard error. Messages are prefixed by
             * the name of running algorithm.
             */
        #define PRINTF_STDERR(...) do {                 \
            fprintf(stderr, __PREFIX__ __VA_ARGS__);    \
            fflush(stdout);                             \
        } while (0)
            
            /**
             * \def PRINT_ERROR(...)
             * Print on the standard error. Messages are prefixed by
             * the name of running algorithm and the function where
             * the error occured.
             */
        #define PRINT_ERROR(...) do {                                      \
            fprintf(stderr,__PREFIX__ "%s%s: ",_TIFA_ERR_PRFX_,__func__);  \
            gmp_fprintf(stderr, __VA_ARGS__);                              \
        } while (0)
    #else
        #define _TIFA_ERR_PRFX_ "TIFA_ERROR: "
        
            /**
             * \def PRINTF_STDERR(...)
             * Print on the standard error. Messages are prefixed by
             * \c _TIFA_ERR_PRFX_.
             */
        #define PRINTF_STDERR(...) fprintf(stderr, __VA_ARGS__); fflush(stdout)
        
            /**
             * \def PRINTF_ERROR(...)
             * Print on the standard error. Messages are prefixed by
             * \c _TIFA_ERR_PRFX_ and the name of the function where
             * the error occured.
             */
        #define PRINT_ERROR(...) do {                              \
            fprintf(stderr, "%s%s: ", _TIFA_ERR_PRFX_, __func__);  \
            gmp_fprintf(stderr, __VA_ARGS__);                      \
        } while (0)
    #endif
#else
    #define PRINTF_STDERR(...) /* intentionally left empty */
    #define PRINT_ERROR(...)   /* intentionally left empty */
#endif

#ifdef __cplusplus
}
#endif

#endif
