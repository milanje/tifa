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
 * \file    batch_gcd_program.c
 * \author  Jerome Milan
 *
 * The integers of files (RSA moduli, say) that share a prime with another,
 * from tifa_batch_gcd.
 */

#include <tifa_config.h>

#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "array.h"
#include "batch_gcd.h"
#include "exit_codes.h"
#include "tifa_alloc.h"
#include "tool_utils.h"

//------------------------------------------------------------------------------
static void print_usage(const char* const name) {
    fprintf(stderr, "Usage: %s [--hex] [--all] [<file>...]\n\n", name);
    fprintf(stderr,
        "Reads positive integers, separated by spaces or newlines, from the\n"
        "files (from the standard input if there are none, or for '-'),\n"
        "'#' starting a comment to the end of its line. Prints, for each\n"
        "integer sharing a prime with another, its position among them\n"
        "(from 1) and that gcd: the gcd of the integer with the product of\n"
        "the others. It is the integer itself if all its primes are shared,\n"
        "as with a duplicate; otherwise, a factor of it.\n\n"
        "--hex: the integers are hexadecimal (with or without 0x), and so\n"
        "       are the gcds printed.\n"
        "--all: prints the gcd of every integer instead, one per line, in\n"
        "       order (1 for those sharing no prime).\n\n"
        "Exit status: %d some integers share primes, %d none do, "
        "%d invalid input.\n\n",
        TIFA_EXIT_COMPLETE, TIFA_EXIT_NO_FACTOR, TIFA_EXIT_USAGE);
}
//------------------------------------------------------------------------------
//
// Appends the token to ints, if a positive integer in the base.
//
static bool append_token(mpz_array_t* const ints, char* token,
                         const int base) {
    if ((base == 16) && (token[0] == '0')
            && ((token[1] == 'x') || (token[1] == 'X'))) {
        token += 2;
    }
    if (token[0] == '\0') {
        return false;
    }
    for (const char* c = token; *c != '\0'; c++) {
        if ((base == 16) ? !isxdigit((unsigned char) *c)
                         : !isdigit((unsigned char) *c)) {
            return false;
        }
    }
    //
    // Doubled when full: append_mpz_to_array grows by ELONGATION only.
    //
    if (ints->length == ints->alloced) {
        resize_mpz_array(ints, 2 * ints->alloced);
    }
    mpz_set_str(ints->data[ints->length], token, base);
    if (mpz_sgn(ints->data[ints->length]) <= 0) {
        return false;
    }
    ints->length++;
    return true;
}
//------------------------------------------------------------------------------
static bool read_integers(mpz_array_t* const ints, FILE* const file,
                          const char* const name, const int base) {
    size_t alloced = 1024;
    size_t length  = 0;
    char*  token   = tifa_malloc(alloced);
    unsigned long line = 1;
    bool ok = true;

    int c = getc(file);
    while (ok && (c != EOF)) {
        if (c == '#') {
            while ((c != EOF) && (c != '\n')) {
                c = getc(file);
            }
            continue;
        }
        if (isspace(c)) {
            line += (c == '\n');
            c = getc(file);
            continue;
        }
        length = 0;
        while ((c != EOF) && !isspace(c)) {
            if (length + 1 == alloced) {
                alloced *= 2;
                token = tifa_realloc(token, alloced);
            }
            token[length++] = (char) c;
            c = getc(file);
        }
        token[length] = '\0';
        if (!append_token(ints, token, base)) {
            fprintf(stderr, "%s:%lu: not a positive integer: %.40s%s\n",
                    name, line, token, (length > 40) ? "..." : "");
            ok = false;
        }
    }
    if (ok && ferror(file)) {
        fprintf(stderr, "%s: read error\n", name);
        ok = false;
    }
    free(token);
    return ok;
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
    int  base  = 10;
    bool all   = false;
    int  nargs = 0;

    for (int i = 1; i < argc; i++) {
        if (0 == strcmp(argv[i], "--hex")) {
            base = 16;
        } else if (0 == strcmp(argv[i], "--all")) {
            all = true;
        } else if (   (0 == strcmp(argv[i], "--help"))
                   || (0 == strcmp(argv[i], "-h"))) {
            print_usage(argv[0]);
            return TIFA_EXIT_COMPLETE;
        } else if ((argv[i][0] == '-') && (argv[i][1] != '\0')) {
            fprintf(stderr, "Unknown option: %s\n\n", argv[i]);
            print_usage(argv[0]);
            return TIFA_EXIT_USAGE;
        } else {
            argv[++nargs] = argv[i];
        }
    }
    mpz_array_t* const ints = alloc_mpz_array(1024);
    bool ok = true;

    if (nargs == 0) {
        ok = read_integers(ints, stdin, "<stdin>", base);
    }
    for (int i = 1; ok && (i <= nargs); i++) {
        if (0 == strcmp(argv[i], "-")) {
            ok = read_integers(ints, stdin, "<stdin>", base);
            continue;
        }
        FILE* const file = fopen(argv[i], "r");
        if (file == NULL) {
            perror(argv[i]);
            ok = false;
            break;
        }
        ok = read_integers(ints, file, argv[i], base);
        fclose(file);
    }
    if (!ok) {
        free_mpz_array(ints);
        return TIFA_EXIT_USAGE;
    }
    mpz_array_t* const gcds = alloc_mpz_array(ints->length);
    tifa_batch_gcd(gcds, ints);

    uint32_t nshared = 0;
    for (uint32_t i = 0; i < gcds->length; i++) {
        const bool shared = (mpz_cmp_ui(gcds->data[i], 1) != 0);
        nshared += shared;
        if (all || shared) {
            if (!all) {
                printf("%" PRIu32 " ", i + 1);
            }
            mpz_out_str(stdout, base, gcds->data[i]);
            putchar('\n');
        }
    }
    fprintf(stderr, "%" PRIu32 " integers, %" PRIu32
                    " sharing a prime with another\n",
            ints->length, nshared);

    free_mpz_array(gcds);
    free_mpz_array(ints);

    return (nshared != 0) ? TIFA_EXIT_COMPLETE : TIFA_EXIT_NO_FACTOR;
}
