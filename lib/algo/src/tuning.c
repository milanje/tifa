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
 * \file    tuning.c
 * \author  Jerome Milan
 * \date    2026
 * \version 1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdbool.h>
#include <pthread.h>

#include <gmp.h>

#include "tifa_config.h"
#include "tuning.h"

#if defined(TIFA_USE_BAKED_TUNING) && TIFA_USE_BAKED_TUNING
#include "tuning_baked.h"
#endif

#define __PREFIX__ "tuning: "

#include "print_error.h"

//
// Hard limits kept deliberately small: a tuning file has one row per bit-size
// step, so a few dozen rows is typical and 256 is comfortably above any real
// need. A line long enough to overflow this buffer cannot be a valid data row.
//
#define TUNING_MAX_ROWS  256
#define TUNING_LINE_MAX  1024

//------------------------------------------------------------------------------
//                             FILE PARSING
//------------------------------------------------------------------------------

//
// Skips leading blanks and returns a pointer to the first non-blank character.
//
static const char* skip_blanks(const char* p) {
    while ((*p == ' ') || (*p == '\t') || (*p == '\r') || (*p == '\n')) {
        p++;
    }
    return p;
}

//
// Returns true if the line (already stripped of leading blanks) is empty or a
// comment.
//
static bool is_ignorable(const char* p) {
    return (*p == '\0') || (*p == '#');
}

//
// Counts the whitespace-separated tokens in a string.
//
static uint32_t count_tokens(const char* p) {
    uint32_t n = 0;
    p = skip_blanks(p);
    while (*p != '\0') {
        n++;
        while ((*p != '\0') && (*p != ' ') && (*p != '\t')
               && (*p != '\r') && (*p != '\n')) {
            p++;
        }
        p = skip_blanks(p);
    }
    return n;
}

//
// Validates a "columns:" directive line. The tokens after the "columns:"
// keyword must be "min_bitsize" followed by exactly expect_ncols names.
// Returns true if the directive is consistent with expect_ncols.
//
static bool validate_columns_line(const char* p, uint32_t expect_ncols) {
    // Skip the "columns" keyword and an optional ':'.
    p += strlen("columns");
    p = skip_blanks(p);
    if (*p == ':') {
        p++;
    }
    return (count_tokens(p) == (expect_ncols + 1));
}

//
// Parses a data row into out_row. Requires exactly (1 + expect_ncols) unsigned
// decimal integers: min_bitsize followed by the parameter values. Returns true
// on success.
//
static bool parse_row(const char* p, uint32_t expect_ncols,
                      tuning_row_t* const out_row) {
    uint32_t nvals = 1 + expect_ncols;
    uint32_t values[1 + TIFA_TUNING_MAX_COLS];

    for (uint32_t i = 0; i < nvals; i++) {
        char* end = NULL;
        unsigned long v = strtoul(p, &end, 10);
        if (end == p) {
            // Not enough integers on the line.
            return false;
        }
        values[i] = (uint32_t) v;
        p = end;
    }
    // Nothing but blanks may follow the expected integers.
    p = skip_blanks(p);
    if (*p != '\0') {
        return false;
    }

    out_row->min_bitsize = values[0];
    for (uint32_t i = 0; i < expect_ncols; i++) {
        out_row->cols[i] = values[1 + i];
    }
    for (uint32_t i = expect_ncols; i < TIFA_TUNING_MAX_COLS; i++) {
        out_row->cols[i] = 0;
    }
    return true;
}

//
// qsort comparator ordering rows by ascending min_bitsize.
//
static int cmp_rows(const void* a, const void* b) {
    const tuning_row_t* ra = (const tuning_row_t*) a;
    const tuning_row_t* rb = (const tuning_row_t*) b;
    if (ra->min_bitsize < rb->min_bitsize) {
        return -1;
    }
    if (ra->min_bitsize > rb->min_bitsize) {
        return 1;
    }
    return 0;
}

//------------------------------------------------------------------------------
int tifa_tuning_load_file(const char* path, uint32_t expect_ncols,
                          tuning_table_t** out_table) {
    *out_table = NULL;

    if ((expect_ncols == 0) || (expect_ncols > TIFA_TUNING_MAX_COLS)) {
        return 1;
    }

    FILE* file = fopen(path, "r");
    if (file == NULL) {
        // A missing file is the common, expected case: no error message.
        return 1;
    }

    tuning_row_t* rows = malloc(TUNING_MAX_ROWS * sizeof(tuning_row_t));
    if (rows == NULL) {
        fclose(file);
        return 1;
    }

    uint32_t nrows   = 0;
    bool     bad     = false;
    char     line[TUNING_LINE_MAX];

    while (fgets(line, (int) sizeof(line), file) != NULL) {
        const char* p = skip_blanks(line);

        if (is_ignorable(p)) {
            continue;
        }
        if (strncmp(p, "columns", strlen("columns")) == 0) {
            if (!validate_columns_line(p, expect_ncols)) {
                PRINT_ERROR("%s: inconsistent 'columns' directive\n", path);
                bad = true;
                break;
            }
            continue;
        }
        if (nrows == TUNING_MAX_ROWS) {
            PRINT_ERROR("%s: too many rows (max %d)\n", path, TUNING_MAX_ROWS);
            bad = true;
            break;
        }
        if (!parse_row(p, expect_ncols, &rows[nrows])) {
            PRINT_ERROR("%s: malformed data row\n", path);
            bad = true;
            break;
        }
        nrows++;
    }

    fclose(file);

    if (bad || (nrows == 0)) {
        free(rows);
        return 1;
    }

    tuning_table_t* table = malloc(sizeof(tuning_table_t));
    if (table == NULL) {
        free(rows);
        return 1;
    }

    qsort(rows, nrows, sizeof(tuning_row_t), cmp_rows);

    table->ncols = expect_ncols;
    table->nrows = nrows;
    table->rows  = rows;
    table->algo[0] = '\0';

    *out_table = table;
    return 0;
}

//------------------------------------------------------------------------------
void tifa_tuning_free(tuning_table_t* table) {
    if (table != NULL) {
        free(table->rows);
        free(table);
    }
}

//------------------------------------------------------------------------------
const tuning_row_t* tifa_tuning_lookup(const tuning_table_t* table,
                                       uint32_t bitsize) {
    if ((table == NULL) || (table->nrows == 0)) {
        return NULL;
    }
    // Rows are sorted ascending: pick the last row whose min_bitsize <= bitsize,
    // clamping to the first row when bitsize is below the whole table.
    const tuning_row_t* chosen = &table->rows[0];
    for (uint32_t i = 0; i < table->nrows; i++) {
        if (table->rows[i].min_bitsize <= bitsize) {
            chosen = &table->rows[i];
        } else {
            break;
        }
    }
    return chosen;
}

//------------------------------------------------------------------------------
bool tifa_tuning_interpolate(const tuning_table_t* table, uint32_t bitsize,
                             uint32_t interp_cols, tuning_row_t* out) {
    const tuning_row_t* lo = tifa_tuning_lookup(table, bitsize);
    if (lo == NULL) {
        return false;
    }
    *out = *lo;
    const tuning_row_t* hi = lo + 1;
    if ((bitsize <= lo->min_bitsize) || (hi == table->rows + table->nrows)) {
        return true;
    }
    const double r = (double) (bitsize - lo->min_bitsize)
                   / (double) (hi->min_bitsize - lo->min_bitsize);
    for (uint32_t c = 0; c < table->ncols; c++) {
        if (interp_cols & (1U << c)) {
            out->cols[c] = (uint32_t) (lo->cols[c]
                         + r * ((double) hi->cols[c] - lo->cols[c]) + 0.5);
        }
    }
    return true;
}

//------------------------------------------------------------------------------
//                     DIRECTORY-BASED LAZY RESOLVER
//------------------------------------------------------------------------------

//
// Small per-algorithm cache. Each entry records the parsed table (NULL on a
// miss) and whether a load has already been attempted, so each file is opened
// at most once. Sized above the number of tunable algorithms.
//
#define TUNING_CACHE_SIZE 16

typedef struct {
    char            algo[TIFA_TUNING_NAME_MAX];
    tuning_table_t* table;
    bool            in_use;
    bool            tried;
} tuning_cache_entry_t;

static tuning_cache_entry_t tuning_cache[TUNING_CACHE_SIZE];
static char*                tuning_dir_override = NULL;
//
// Guards tuning_cache and tuning_dir_override.
//
static pthread_mutex_t      tuning_lock = PTHREAD_MUTEX_INITIALIZER;

//
// Drops all cached tables (called when the search directory changes).
//
static void tuning_cache_reset(void) {
    for (uint32_t i = 0; i < TUNING_CACHE_SIZE; i++) {
        if (tuning_cache[i].in_use) {
            tifa_tuning_free(tuning_cache[i].table);
        }
        tuning_cache[i].table  = NULL;
        tuning_cache[i].in_use = false;
        tuning_cache[i].tried  = false;
        tuning_cache[i].algo[0] = '\0';
    }
}

//
// Returns the cache entry for algo, allocating a slot on first use. Returns NULL
// only if the cache is full (cannot happen for the built-in algorithm set).
//
static tuning_cache_entry_t* tuning_cache_entry(const char* algo) {
    tuning_cache_entry_t* free_slot = NULL;
    for (uint32_t i = 0; i < TUNING_CACHE_SIZE; i++) {
        if (tuning_cache[i].in_use) {
            if (strcmp(tuning_cache[i].algo, algo) == 0) {
                return &tuning_cache[i];
            }
        } else if (free_slot == NULL) {
            free_slot = &tuning_cache[i];
        }
    }
    if (free_slot != NULL) {
        snprintf(free_slot->algo, sizeof(free_slot->algo), "%s", algo);
        free_slot->in_use = true;
        free_slot->tried  = false;
        free_slot->table  = NULL;
    }
    return free_slot;
}

//
// Resolves the tuning directory: the programmatic override wins, then the
// TIFA_TUNING_DIR environment variable, then the compiled-in default (if any).
// Returns NULL when no directory is configured.
//
static const char* resolve_tuning_dir(void) {
    if (tuning_dir_override != NULL) {
        return tuning_dir_override;
    }
    const char* env = getenv("TIFA_TUNING_DIR");
    if ((env != NULL) && (env[0] != '\0')) {
        return env;
    }
#if defined(TIFA_TUNING_DIR)
    if (TIFA_TUNING_DIR[0] != '\0') {
        return TIFA_TUNING_DIR;
    }
#endif
    return NULL;
}

//------------------------------------------------------------------------------
void tifa_tuning_set_dir(const char* dir) {
    pthread_mutex_lock(&tuning_lock);
    free(tuning_dir_override);
    tuning_dir_override = NULL;
    if (dir != NULL) {
        size_t len = strlen(dir);
        tuning_dir_override = malloc(len + 1);
        if (tuning_dir_override != NULL) {
            memcpy(tuning_dir_override, dir, len + 1);
        }
    }
    // A directory change invalidates every cached table.
    tuning_cache_reset();
    pthread_mutex_unlock(&tuning_lock);
}

//------------------------------------------------------------------------------
//
// The table for algo: its runtime file, else its baked table, else NULL. The
// caller holds the lock.
//
static const tuning_table_t* table_for(const char* algo,
                                       uint32_t expect_ncols) {
    tuning_cache_entry_t* entry = tuning_cache_entry(algo);
    if (entry == NULL) {
        return NULL;
    }

    if (!entry->tried) {
        entry->tried = true;

        const char* dir = resolve_tuning_dir();
        if (dir != NULL) {
            char path[TUNING_LINE_MAX];
            snprintf(path, sizeof(path), "%s/%s.tuning", dir, algo);
            // Failure leaves entry->table == NULL, i.e. a cached miss.
            (void) tifa_tuning_load_file(path, expect_ncols, &entry->table);
        }
    }

    // Highest priority: a runtime tuning file.
    if (entry->table != NULL) {
        return entry->table;
    }

    // Middle priority (only when compiled in): baked tables. The caller's
    // own hardcoded defaults remain the final fallback (we return NULL).
#if defined(TIFA_USE_BAKED_TUNING) && TIFA_USE_BAKED_TUNING
    return tifa_tuning_baked_table(algo, expect_ncols);
#else
    return NULL;
#endif
}

//------------------------------------------------------------------------------
const tuning_row_t* tifa_tuning_row_for(const char* algo, uint32_t expect_ncols,
                                        uint32_t bitsize) {
    pthread_mutex_lock(&tuning_lock);
    const tuning_row_t* row =
        tifa_tuning_lookup(table_for(algo, expect_ncols), bitsize);
    pthread_mutex_unlock(&tuning_lock);
    return row;
}

//------------------------------------------------------------------------------
bool tifa_tuning_get_row(const char* algo, uint32_t expect_ncols,
                         uint32_t bitsize, uint32_t interp_cols,
                         tuning_row_t* const out) {
    pthread_mutex_lock(&tuning_lock);
    const bool found = tifa_tuning_interpolate(
                           table_for(algo, expect_ncols), bitsize,
                           interp_cols, out
                       );
    pthread_mutex_unlock(&tuning_lock);
    return found;
}
