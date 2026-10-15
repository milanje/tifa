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
 * \file    tuning.h
 * \author  Jerome Milan
 * \date    2026
 * \version 1.0
 *
 * \brief Runtime loader for machine-tuned algorithm parameters.
 *
 * The factoring algorithms depend on numeric parameters whose best values
 * depend on the bit-size of the number being factored. This module lets those
 * values be read at runtime from per-algorithm text files (produced by the
 * \c tifa_tune tool) instead of the compiled-in defaults.
 *
 * Each <tt>set_&lt;algo&gt;_params_to_default</tt> function first fills the
 * parameters with its usual defaults, then consults \c tifa_tuning_row_for. If
 * a tuning file provides a row for the number's bit-size, the tunable fields
 * are overridden; otherwise the defaults are left untouched. When no tuning
 * file is present the library behaves exactly as if this module did not exist.
 *
 * \section format Tuning file format
 *
 * A tuning file is a line-based text file. Lines beginning with \c '#' and
 * blank lines are ignored (the tuner emits a metadata block as comments). An
 * optional directive line <tt>columns: min_bitsize &lt;name...&gt;</tt> may
 * describe the columns. Data rows are whitespace-separated unsigned decimal
 * integers: the first field is \c min_bitsize, followed by exactly \c ncols
 * parameter values. Rows may appear in any order; they are sorted internally.
 *
 * Lookup for a bit-size \c s returns the row with the largest
 * <tt>min_bitsize &lt;= s</tt>; if \c s is below the first row, the first row
 * is used.
 */

#if !defined(_TIFA_TUNING_H_)
#define _TIFA_TUNING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <stdbool.h>

   /**
    * \def TIFA_TUNING_MAX_COLS
    * Maximum number of tunable parameter columns per row (excluding
    * \c min_bitsize).
    */
#define TIFA_TUNING_MAX_COLS 12

   /**
    * \def TIFA_TUNING_NAME_MAX
    * Maximum length (including the terminating NUL) of an algorithm name.
    */
#define TIFA_TUNING_NAME_MAX 32

//
// Per-algorithm column layout. These macros are the single source of truth for
// each tuning file's schema; they are shared by the loader (validation), the
// set_<algo>_params_to_default fillers (column -> field mapping), and the
// tifa_tune tool (grid and file output). Adding a new algorithm means adding
// one macro pair here, one filler tail, and one tuner descriptor: the parser
// itself never changes.
//
#define TIFA_TUNING_NCOLS_P_MINUS_1 3
#define TIFA_TUNING_COLS_P_MINUS_1  "b1 b2 base"

#define TIFA_TUNING_NCOLS_P_PLUS_1  3
#define TIFA_TUNING_COLS_P_PLUS_1   "b1 b2 seed"

#define TIFA_TUNING_NCOLS_RHO       4
#define TIFA_TUNING_COLS_RHO        "a b acc_niters max_niters"

#define TIFA_TUNING_NCOLS_ECM       3
#define TIFA_TUNING_COLS_ECM        "b1 b2 ncurves"

#define TIFA_TUNING_NCOLS_CFRAC     5
#define TIFA_TUNING_COLS_CFRAC      "nprimes_in_base nprimes_tdiv nrelations use_large_primes nsteps_early_abort"

#define TIFA_TUNING_NCOLS_SIQS      6
#define TIFA_TUNING_COLS_SIQS       "nprimes_in_base sieve_half_width threshold nprimes_tdiv nrelations use_large_primes"

   /**
    * \struct struct_tuning_row_t tuning.h lib/algo/include/tuning.h
    * \brief  A single row of a tuning table.
    *
    * A row associates a minimum bit-size with a set of tunable parameter
    * values. All values are stored as \c uint32_t (booleans and enumerations
    * fit in a \c uint32_t) so the on-disk format is a uniform grid of unsigned
    * integers.
    */
struct struct_tuning_row_t {
        /**
         * The smallest bit-size of \c n for which this row applies.
         */
    uint32_t min_bitsize;
        /**
         * The tunable parameter values, in the column order defined by the
         * algorithm's \c TIFA_TUNING_COLS_* macro. Only the first \c ncols
         * entries (see \c tuning_table_t) are meaningful.
         */
    uint32_t cols[TIFA_TUNING_MAX_COLS];
};

   /**
    * \typedef tuning_row_t
    * \brief Equivalent to <tt>struct struct_tuning_row_t</tt>.
    */
typedef struct struct_tuning_row_t tuning_row_t;

   /**
    * \struct struct_tuning_table_t tuning.h lib/algo/include/tuning.h
    * \brief  A parsed tuning table for one algorithm.
    */
struct struct_tuning_table_t {
        /**
         * The algorithm name this table was loaded for.
         */
    char algo[TIFA_TUNING_NAME_MAX];
        /**
         * Number of meaningful parameter columns per row (excluding
         * \c min_bitsize).
         */
    uint32_t ncols;
        /**
         * Number of rows in the \c rows array.
         */
    uint32_t nrows;
        /**
         * The rows, heap-allocated and sorted ascending by \c min_bitsize.
         */
    tuning_row_t* rows;
};

   /**
    * \typedef tuning_table_t
    * \brief Equivalent to <tt>struct struct_tuning_table_t</tt>.
    */
typedef struct struct_tuning_table_t tuning_table_t;

   /**
    * \brief Loads a tuning table from a file.
    *
    * Parses the tuning file at \c path, requiring each data row to hold exactly
    * <tt>expect_ncols</tt> parameter values (in addition to \c min_bitsize). On
    * success a newly allocated \c tuning_table_t is stored in <tt>*out_table</tt>
    * and the caller is responsible for freeing it with \c tifa_tuning_free.
    *
    * \param[in]  path         Path to the tuning file.
    * \param[in]  expect_ncols The expected number of parameter columns.
    * \param[out] out_table    Receives the allocated table on success.
    * \return 0 on success, a non-zero value on failure (file not found, parse
    *         error, column-count mismatch, or no data rows). On failure
    *         <tt>*out_table</tt> is set to NULL.
    */
int tifa_tuning_load_file(const char* path, uint32_t expect_ncols,
                          tuning_table_t** out_table);

   /**
    * \brief Frees a tuning table.
    *
    * \param[in] table The table to free (may be NULL).
    */
void tifa_tuning_free(tuning_table_t* table);

   /**
    * \brief Looks up the row applicable to a given bit-size.
    *
    * Returns the row with the largest <tt>min_bitsize &lt;= bitsize</tt>. If
    * \c bitsize is below the first row, the first row is returned.
    *
    * \param[in] table   The tuning table (must be non-NULL and non-empty).
    * \param[in] bitsize The bit-size to look up.
    * \return The applicable row, or NULL if \c table is NULL or empty.
    */
const tuning_row_t* tifa_tuning_lookup(const tuning_table_t* table,
                                       uint32_t bitsize);

   /**
    * \brief Sets the directory searched for tuning files.
    *
    * Overrides, for the rest of the process, the directory in which
    * \c tifa_tuning_row_for looks for <tt>&lt;algo&gt;.tuning</tt> files. This
    * takes precedence over the \c TIFA_TUNING_DIR environment variable and the
    * compiled-in default directory. Passing NULL clears the override.
    *
    * \note Calling this function invalidates the internal cache so a subsequent
    * lookup reloads from the new directory.
    *
    * \param[in] dir The directory to search, or NULL to clear the override.
    */
void tifa_tuning_set_dir(const char* dir);

   /**
    * \brief Returns the tuning row for an algorithm and bit-size, if any.
    *
    * Resolves the tuning directory (the value set by \c tifa_tuning_set_dir,
    * else the \c TIFA_TUNING_DIR environment variable, else the compiled-in
    * default), lazily loads <tt>&lt;dir&gt;/&lt;algo&gt;.tuning</tt> the first
    * time it is needed, and returns the row applicable to \c bitsize.
    *
    * Loaded tables are cached per algorithm for the lifetime of the process, so
    * each file is parsed at most once. Any failure (no directory, missing file,
    * parse error, or column-count mismatch) yields NULL and is cached as a miss.
    *
    * \warning The returned row lives in the cache: it is freed by the next
    * call to \c tifa_tuning_set_dir. Use \c tifa_tuning_get_row when that
    * call may happen in another thread.
    *
    * \param[in] algo         The algorithm name (matches the file base name).
    * \param[in] expect_ncols The expected number of parameter columns.
    * \param[in] bitsize      The bit-size of the number to factor.
    * \return The applicable tuning row, or NULL if no tuning data applies.
    */
const tuning_row_t* tifa_tuning_row_for(const char* algo, uint32_t expect_ncols,
                                        uint32_t bitsize);

   /**
    * \brief Copies the tuning row for an algorithm and bit-size, if any.
    *
    * Same as \c tifa_tuning_row_for, but copies the row into \c out, so the
    * result stays valid whatever other threads do. The tuning functions may
    * be called from several threads at once.
    *
    * \param[in]  algo         The algorithm name (matches the file base name).
    * \param[in]  expect_ncols The expected number of parameter columns.
    * \param[in]  bitsize      The bit-size of the number to factor.
    * \param[out] out          The applicable tuning row, if any.
    * \return true if a tuning row applies (and was copied into \c out).
    */
bool tifa_tuning_get_row(const char* algo, uint32_t expect_ncols,
                         uint32_t bitsize, tuning_row_t* out);

#ifdef __cplusplus
}
#endif

#endif
