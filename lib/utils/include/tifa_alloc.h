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
 * \file    tifa_alloc.h
 * \author  Jerome Milan
 *
 * \brief Memory allocation that stops the program when memory runs out.
 *
 * As GMP's own allocation does: the library cannot go on without the memory
 * it asks for, and a NULL pointer used further on is worse than a stop.
 */

#if !defined(_TIFA_ALLOC_H_)
   /**
    * \def _TIFA_ALLOC_H_
    * Standard include guard.
    */
#define _TIFA_ALLOC_H_

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \brief Prints an error message and aborts.
    */
static inline void tifa_out_of_memory(void) {
    fprintf(stderr, "TIFA: out of memory\n");
    abort();
}

   /**
    * \brief \c malloc, aborting if the memory cannot be allocated.
    */
static inline void* tifa_malloc(size_t size) {
    void* const ptr = malloc(size);
    if ((ptr == NULL) && (size != 0)) {
        tifa_out_of_memory();
    }
    return ptr;
}

   /**
    * \brief \c calloc, aborting if the memory cannot be allocated.
    */
static inline void* tifa_calloc(size_t nmemb, size_t size) {
    void* const ptr = calloc(nmemb, size);
    if ((ptr == NULL) && (nmemb != 0) && (size != 0)) {
        tifa_out_of_memory();
    }
    return ptr;
}

   /**
    * \brief \c realloc, aborting if the memory cannot be allocated.
    */
static inline void* tifa_realloc(void* ptr, size_t size) {
    void* const newptr = realloc(ptr, size);
    if ((newptr == NULL) && (size != 0)) {
        tifa_out_of_memory();
    }
    return newptr;
}

#ifdef __cplusplus
}
#endif

#endif
