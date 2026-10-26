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
 * \file    tool_utils.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#if HAVE_CONFIG_H
#    include <tifa_config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"

//------------------------------------------------------------------------------
bool is_a_number(const char* const str_n, uint32_t length) {

    if (0 == isdigit((unsigned char) str_n[0])) {
        return false;
    }
    for (uint32_t i = 1; i < length; i++) {
        if (str_n[i] == '\0') {
            return true;
        }
        if (0 == isdigit((unsigned char) str_n[i])) {
            return false;
        }
    }
    // Longer than length: the rest was never checked.
    return false;
}
//------------------------------------------------------------------------------
char* read_number(char* const buffer, uint32_t length) {
    if (fgets(buffer, (int) length, stdin) == NULL) {
        fprintf(stderr, "\nERROR: no number to read!\n");
        exit(-1);
    }
    char* const newline = strchr(buffer, '\n');
    if (newline != NULL) {
        *newline = '\0';
        if ((newline != buffer) && (newline[-1] == '\r')) {
            newline[-1] = '\0';
        }
    } else if (!feof(stdin)) {
        fprintf(stderr, "\nERROR: the number has more than %" PRIu32
                        " digits!\n", length - 2);
        exit(-1);
    }
    return buffer;
}
//------------------------------------------------------------------------------

