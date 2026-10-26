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
 * \file    tool_utils.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Miscellaneous helpful functions.
 *
 * Miscellaneous functions used by the "tool programs".
 */

#if !defined(_TIFA_TOOL_UTILS_H_)
   /**
    * \def _TIFA_TOOL_UTILS_H_
    * Standard include guard.
    */
#define _TIFA_TOOL_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <stdbool.h>

   /**
    * \brief Does a string \c str_n read as a number?
    *
    * Returns \c true if the string \c str_n represents a number in the
    * decimal base, of fewer than \c length digits, \c false otherwise.
    *
    * \param[in] str_n The string to check.
    * \param[in] length The bound on the length of the string.
    */
bool is_a_number(const char* const str_n, uint32_t length);

   /**
    * \brief Reads a line of the standard input, without its newline.
    *
    * Exits the program with an error message if nothing can be read or if
    * the line does not fit in \c length - 1 characters.
    *
    * \param[out] buffer The line read.
    * \param[in]  length The size of \c buffer.
    * \return \c buffer.
    */
char* read_number(char* const buffer, uint32_t length);

#ifdef __cplusplus
}
#endif

#endif
