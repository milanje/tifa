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
 * \file    x_array_list.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include "stdlib.h"
#include "stdio.h"
#include "x_array_list.h"

/*
 *-----------------------------------------------------------------------------
 *              uint32_array_list_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
uint32_array_list_t* alloc_uint32_array_list(uint32_t alloced) {

    uint32_array_list_t* list = malloc(sizeof(uint32_array_list_t));

    list->alloced = alloced;
    list->length  = 0;
    list->data    = malloc(list->alloced*sizeof(uint32_array_t*));

    return list;
}
//-----------------------------------------------------------------------------
void free_uint32_array_list(uint32_array_list_t* const list) {
    for (uint32_t i = 0U; i < list->length; i++) {
        if (list->data[i] != NULL) {
            free_uint32_array(list->data[i]);
        }
    }
    free(list->data);
    free(list);
}
//-----------------------------------------------------------------------------
void print_uint32_array_list(const uint32_array_list_t* const list) {
    //
    // Mostly for debugging purposes...
    //
    if (list->length == 0U) {
        printf("<empty uint32_array_list>\n");
    }
    for (uint32_t i = 0U; i != list->length; i++) {
        if (list->data[i] != NULL) {
            print_uint32_array(list->data[i]);
        } else {
            printf("<null entry>\n");
        }
        printf("------------------------------\n");
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              mpz_array_list_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
mpz_array_list_t* alloc_mpz_array_list(uint32_t alloced) {

    mpz_array_list_t* list = malloc(sizeof(mpz_array_list_t));

    list->alloced = alloced;
    list->length  = 0;
    list->data    = malloc(list->alloced*sizeof(mpz_array_t*));

    return list;
}
//-----------------------------------------------------------------------------
void free_mpz_array_list(mpz_array_list_t* const list) {
    for (uint32_t i = 0U; i != list->length; i++) {
        if (list->data[i] != NULL) {
            free_mpz_array(list->data[i]);
        }
    }
    free(list->data);
    free(list);
}
//-----------------------------------------------------------------------------
void print_mpz_array_list(const mpz_array_list_t* const list) {
    //
    // Mostly for debugging purposes...
    //
    if (list->length == 0U) {
        printf("<empty mpz_array_list>\n");
    }
    for (uint32_t i = 0U; i < list->length; i++) {
        if (list->data[i] != NULL) {
            print_mpz_array(list->data[i]);
        } else {
            printf("\t<null entry>\n");
        }
        printf("------------------------------\n");
    }
}
//-----------------------------------------------------------------------------
