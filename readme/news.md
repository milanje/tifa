# TIFA news

## 27 October 2026 &mdash; six primes in `a`

Six primes in `a` were never clearly best, and even slower than both five and
seven from 145 to 155 bits. The approximer, which picks the primes of each
`a`, matched the last part of its target with a product of three primes, and
for six it always precomputed and sorted all C(49, 3) = 18424 triples of its
pool: 0.5 ms per run, 2.4% at 150 bits, against 1176 pairs for five. With
pairs, six primes run 1.1% to 3.8% faster from 145 to 160 bits:

| SIQS with six primes in `a`, bits | 145 | 150 | 155 | 160 |
|---|---:|---:|---:|---:|
| time, pairs over triples | 0.99 | 0.97 | 0.96 | 0.97 |

The defaults are unchanged: six primes were then best at 150 bits only, by
1.3% over five. Six primes come from `extra_primes_in_a` only, in tuning files
or parameters. Pairs from seven primes up gave slower `a`
(2 to 4% at 170 and 190 bits), so these keep triples.

## 27 October 2026 &mdash; SIQS collects fewer extra relations

SIQS collected 24 relations beyond the size of its factor base, for 24
dependencies: with k of them, none splits n = pq with probability 2^-k. It now
collects 8 (`SIQS_DFLT_NRELATIONS`), and when none of the dependencies splits
n, 8 more, twice at most: 24 in the worst case, as before, which happens once
in 65536 runs. The relations already collected are kept. At 70 bits, 9 of
2000 runs needed more relations, and all 2000 were complete. Against the
previous version (100 numbers per size up to 120 bits, 30 up to 150 and 8
above):

| SIQS, bits | 66 | 70 | 75 | 80 | 85 | 90 | 95 | 100 | 105 | 110 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.87 | 0.87 | 0.85 | 0.89 | 0.90 | 0.93 | 0.92 | 0.92 | 0.96 | 0.96 |

| SIQS, bits | 115 | 120 | 130 | 140 | 150 | 160 | 170 | 180 | 200 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.97 | 0.97 | 0.97 | 1.00 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 |

(At 180 bits, 1.02 in a first run, 1.00 in a second.)

| `tifa_factor`, bits | 66 | 80 | 100 | 128 | 160 |
|---|---:|---:|---:|---:|---:|
| time | 0.88 | 0.91 | 0.94 | 0.96 | 1.00 |

Starting from 6 was 1.5% faster at 66 and 70 bits only, from 4 not clearly
faster than from 6. The other defaults were measured again after the
overhead cuts below, at 70, 90 and 110 bits: a factor base 10% to 30% larger
or 10% smaller, a sieve 20% narrower or 25% wider, and a threshold 1 or 2
higher or lower were all as fast or slower, and so was the large prime
variation turned on at 95 and 100 bits or off at 105 and 110. The relations' values are no longer
multiplied by their `a` after the trial division: nothing read them since
the change below.

## 27 October 2026 &mdash; SIQS overheads

Outside the sieve, a small SIQS run spent much of its time in four places,
now cheaper:

- The Gaussian elimination keeps the rows in lists by their rightmost column
  left, so that each pivot meets only the rows that have it, instead of
  scanning every row for every column: 2.1 to 2.6 times faster on the
  matrices of 70 to 110-bit runs.
- The trial division of the relations works on g(x)/a, one limb shorter at
  70 bits, and adds the primes of `a`, known, instead of finding them again.
  It tests 64 primes at a time without branches, then divides out the
  divisors only: a branch per prime mispredicted on the few that divide.
- Each new `a` needs its inverse modulo every prime of the base: the
  Euclidean algorithms run four at a time, interleaved, twice as fast per
  inverse.
- The size of the L1 data cache, which sets the size of the sieve chunks, is
  looked up once per process: `sysctlbyname` took 0.7% of a 66-bit run.

Against the previous version (time after over time before, balanced
semiprimes, one thread of an Apple M5 Pro; 100 numbers per size up to 120
bits, 30 up to 150 and 8 above):

| SIQS, bits | 66 | 70 | 75 | 80 | 85 | 90 | 95 | 100 | 105 | 110 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.82 | 0.79 | 0.83 | 0.82 | 0.85 | 0.82 | 0.86 | 0.89 | 0.93 | 0.93 |

| SIQS, bits | 115 | 120 | 130 | 140 | 150 | 160 | 170 | 180 | 190 | 200 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.94 | 0.92 | 0.97 | 0.93 | 0.96 | 0.96 | 0.94 | 0.98 | 0.95 | 0.96 |

Each change on its own, against the one before: the elimination 0.92 to 0.93
from 66 to 90 bits; the trial division of g(x)/a 0.98 from 70 to 90 bits,
then without branches 0.93 at 66 and 70 bits and 0.95 to 0.98 up to 110; the
inverses 0.96 to 0.975 from 66 to 110 bits. The library's interface is
unchanged.

## 27 October 2026 &mdash; wider CI

CI ran the tests with GCC and Clang only. It now also runs them under ASan,
UBSan and LeakSanitizer, without 128-bit integers (the GMP code paths the
Montgomery and machine-word code bypass), with `TIFA_USE_CALLOC_MEMSET=OFF`,
and on macOS arm64. Run first on Linux, the new configurations found that
three programs leaked the number they factored, and that some builds did not
compile with `-Werror`: without 128-bit integers, a function and a parameter
used only with them; in Debug with GCC, a format truncation in the programs'
banner.

## 27 October 2026 &mdash; version 0.3.0

Version 0.3.0 gathers the work below since 0.2.0: faster and retuned SIQS,
fixes from a review of latent defects (prime powers, hangs, wrong roots,
memory errors), ECM, p+1 and p-1 splitting small composites, and renamed
programs with Unix exit statuses. Against 0.2.0, on one thread of an Apple
M5 Pro (time of 0.3.0 over that of 0.2.0, balanced semiprimes):

| SIQS, bits | 66 | 70 | 75 | 80 | 85 | 90 | 95 | 100 | 105 | 110 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.89 | 0.88 | 0.87 | 0.80 | 0.86 | 0.91 | 0.89 | 0.89 | 0.92 | 0.93 |

| SIQS, bits | 115 | 120 | 130 | 140 | 150 | 160 | 170 | 180 | 190 | 200 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.95 | 0.96 | 0.97 | 0.97 | 0.94 | 0.95 | 0.94 | 0.98 | 0.90 | 0.99 |

| `tifa_factor`, bits | 40 | 50 | 64 | 66 | 80 | 100 | 128 | 160 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 1.00 | 1.00 | 1.00 | 0.90 | 0.85 | 0.90 | 0.97 | 0.96 |

The ABI changes (`cfrac_params_t` lost `nsteps_early_abort`, `siqs_params_t`
gained `extra_primes_in_a`, `sqrtm_p2` takes an `unsigned long`, among
others): the shared library is now `libtifa.so.0.3` (`libtifa.0.3.dylib` on
macOS), and `find_package(Tifa 0.3)` accepts 0.3.x only. The programs are
renamed `tifa-<name>` and exit with 0, 1, 2 or 3 (complete, partial, none,
invalid input), and tuning files are read more strictly.

## 27 October 2026 &mdash; program names

The programs are now named `tifa-<name>`, with dashes as is usual for Unix
commands: `tifa-factor`, `tifa-tune`, `tifa-cfrac`, `tifa-ecm`, `tifa-fermat`,
`tifa-siqs`, `tifa-squfof`, `tifa-rho`, `tifa-tdiv`, `tifa-p-minus-1` and
`tifa-p-plus-1`. Plain names such as `ecm` clashed with other tools (GMP-ECM's
`ecm`), and the others were too generic for a shared `bin` directory. The
CMake targets keep their names (`tifa_factor`, `tifa_tune`, `ecm`...).

## 27 October 2026 &mdash; exit statuses of the programs

The programs exited with TIFA's internal codes: 5 for a complete
factorization, 4 for a partial one, 6 when no factor was found, 0 for a
prime and 255 for invalid input, so that `tifa_factor N && ...` took success
for failure. They now exit with 0 for a complete factorization (or a prime,
or 1), 1 for a partial one, 2 when no factor was found and 3 for invalid
input, as their usage messages say. Scripts that tested the old values need
updating; the library's `ecode_t` values are unchanged.

## 27 October 2026 &mdash; ECM, p+1 and p-1 on small composites

When a stage caught every prime factor of n at once (a gcd equal to n), ECM
and p+1 gave up: on 18% of the odd composites near 10^5, 11% near 10^6 and
4% near 10^8 for ECM, and on 25% of those near 10^5 for p+1. Like p-1, they
now redo such a stage one prime at a time (stage 1) or one product at a time
(stage 2), with a gcd after each, on both the Montgomery and the GMP paths.

p-1 already did so, but still gave up on 0.5% of them near 10^5 (2% below
5000): the orders of 2 modulo both primes ended in the same prime. Its
backtracking also goes through the primes in decreasing order, and when no
order splits n (equal orders), it tries the bases 3, 5, 7 and 11.

All three now split every composite tried, from 10^3 to 10^10; other runs
take the same time.

## 26 October 2026 &mdash; fixes and documentation

A review for latent defects, with a test for most fixes.

Wrong results, hangs and missed factors:

- `tifa_factor` takes apart any perfect power, not only squares: p^3 and
  most higher powers of a prime above 64 bits gave up (no congruence of
  squares splits a prime power), and so did `tifa_cofactor` on them.
- `cfrac` gives up instead of looping forever when k.n = m^2 + 1, or when
  the continued fraction has a short period without smooth residues.
- `squfof` no longer hangs on even n = k.t^2 (a racer divided by 0), and
  no longer returns n itself as a factor of a few small primes; `ecm`
  no longer returns n and 1 when n divides its curve's discriminant.
- `sqrtm_p2` was wrong for every p > 65535: a single run of `fermat` now
  factors 20 of 20 balanced semiprimes of 70 bits instead of 1.
- `cfrac`'s trial division filters dropped every negative smooth residue,
  and the early-abort one also its base's largest prime and composite
  cofactors; the SIQS approximer could put 2 in `a`, which repeats
  polynomials; p+1's GMP stage 2 skipped the primes in (b1, 104].
- Without `__int128`, `tifa_factor` reported a prime below 8192 as its
  own factor.

Memory and undefined behaviour: an ECM leak at each larger bound, a SQUFOF
clear of an uninitialised integer, an out-of-bounds read in the early-abort
filter, a sieve threshold above 127 (now capped), and latent bugs in helpers
the library rarely or never calls (bit clearing, byte sorting, hash tables,
lists, product trees, binary array sizes, `MPN_ADD`, `MPN_TDIV_QR`,
`MPZ_SWAP`, a ceiling log2).

Programs, tuning and build:

- The programs print the right exponents for perfect powers, fit their
  default parameters to what trial division leaves, and check their
  arguments (`cfrac`'s method numbers crashed it); `rho`'s usage lists its
  `<method>` argument.
- Tuning files with negative, out-of-range or duplicate rows, or overlong
  lines, are rejected; `tifa_tuning_set_dir("")` selects the built-in
  defaults only, which `tifa_tune` now centres its search on; `tifa_tune`
  checks its options and output directory before tuning.
- `-DTIFA_USE_CALLOC_MEMSET=OFF` compiles; `tifa.pc` is relocatable (and
  works on macOS); the exported CMake targets name `m` instead of a path;
  source packages build from a checkout under `/build...`; `rho.h`,
  `p_minus_1.h` and `p_plus_1.h` compile on their own.
- The comments and documentation were checked against the code throughout.

## 26 October 2026 &mdash; number of primes in `a` retuned

`determine_nprimes`, which picks the number of primes in `a` from the size of
its target, now follows the best number measured for each size: two primes
below 15 bits, three below 26, four up to 35, five up to 59, seven up to 69
and eight up to 98 (six were never clearly best). The defaults no longer add
primes for some sizes of n, which the thresholds now cover. Against the
previous defaults:

| Bits | 77 | 83 | 85 | 86 | 88 | 95 | 96 | 98 | 145 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.97 | 0.97 | 0.92 | 0.95 | 0.97 | 0.97 | 0.97 | 0.94 | 0.97 |

The other sizes from 66 to 200 bits are within 1.5% of before.

## 26 October 2026 &mdash; SIQS defaults retuned

`tifa_tune` retuned SIQS from 65 to 200 bits, the number of primes in `a`
included. The defaults take its values only where they measured faster on
other numbers: narrower sieves from 75 to 89 bits, at 110 and 190 bits, and
one more prime in `a` from 75 to 84, 100 to 119 and 150 to 209 bits. Against
the previous defaults:

| Bits | 76 | 78 | 80 | 82 | 84 | 100 | 103 | 110 | 150 | 160 | 170 | 185 | 190 | 195 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| time | 0.97 | 0.91 | 0.89 | 0.91 | 0.95 | 0.95 | 0.96 | 0.99 | 0.97 | 0.97 | 0.95 | 0.86 | 0.91 | 0.90 |

The other sizes from 66 to 200 bits are within 1.5% of before.

## 23 October 2026 &mdash; faster SIQS setup

Outside the sieve, a SIQS run spent a fifth of its time at 70 bits on its
setup and per-run overheads. The factor base takes its residues from the
multiplier choice (bitmaps of quadratic residues instead of Jacobi symbols,
no GMP division or `log()` per prime), `sqrtm` rejects non-residues first and
needs one or two exponentiations for three primes in four, the relations'
trial division stops once the residue is below p^2, and the polynomial
coefficient approximer and the batch smoothness test lost most of their libm
calls, sorting and GMP allocations. The results are unchanged.

| `tifa_factor`, balanced semiprimes | 72 | 80 | 90 | 100 | 110 | 128 | 160 bits |
|---|---:|---:|---:|---:|---:|---:|---:|
| time, against the cleanups below | 0.90 | 0.93 | 0.94 | 0.95 | 0.96 | 0.97 | 0.98 |

`tifa_cofactor` in its exact mode is 2 to 3.5% faster on the sets of the
October tables.

Relations and candidates then keep an index of their polynomial's `a`
instead of a copy (a quarter fewer GMP allocations per run, 2.7% faster at
66-70 bits), and the number of primes in `a` can be raised by one or two
(`extra_primes_in_a`): `tifa_tune` tunes it, in a seventh column of
`siqs.tuning` (files with six columns are refused). Tuned from 85 to 105
bits, it runs 2.4 to 6.2% faster than the defaults, which do not use it.

## 23 October 2026 &mdash; cleanups

The version is now 0.3.0 (`libtifa.so.0.3`), not released yet: builds of
this code no longer pass for 0.2.0, whose ABI differs.

- Every allocation is checked: out of memory, the library prints an error
  and aborts, as GMP does, instead of writing through a NULL pointer
  (`tifa_alloc.h`).
- Arrays allocated or resized to 0 elements no longer leak, shrinking an
  `mpz_array_t` no longer exposes unused entries, the sorted-array searches
  no longer read one element past the end, and the early-abort filter frees
  its arrays' data.
- Dead code is gone, which changes the ABI (hence 0.3.0 next): the
  `nsteps_early_abort` field of `cfrac_params_t` and its tuning column
  (`cfrac.tuning` files now have 4 columns), the unused `bern_21*`
  functions, and the unused SIQS variant of `smooth_filter_t`.
- CI builds with `-Werror` under GCC too; source packages leave out every
  build directory and local files; the tests have license headers.

## 22 October 2026 &mdash; version 0.2.0

Version 0.2.0 gathers all the October 2026 work below: the rework for
speed, correctness and thread safety, `tifa_cofactor`, the review of
invalid inputs and parameters, the `tifa_factor` program and the updated
documentation. It requires GMP 6.0 or later.

The ABI changes (`tifa_rand` and `tifa_srand` are gone, among others), so
the shared library's soname now carries the minor version while the major
one is 0: `libtifa.so.0.2` (`libtifa.0.2.dylib` on macOS). Programs built
against 0.1.5 must be rebuilt. Likewise, `find_package(Tifa 0.2)` accepts
0.2.x only. `tifa.pc` now lists GMP as a public dependency: programs using
TIFA call GMP, and `pkg-config --libs tifa` left `-lgmp` out.

## 22 October 2026 &mdash; invalid inputs and parameters

A review of the library and its programs for memory errors and invalid
inputs (AddressSanitizer, UndefinedBehaviorSanitizer and ThreadSanitizer
builds, edge-case sweeps of every public function, hostile tuning files, a
fuzzed tuning-file parser) found and fixed:

- `tifa_factor` and `tifa_cofactor` read past the table of primes on hard
  composites from 260 bits (the SIQS factor base outgrew it).
- The algorithms called directly crashed or looped forever on n below 2,
  on perfect squares and on some tiny n (SIQS on 3), and wrote out of
  bounds or hung with parameters out of their useful range (factor base
  of 0 or beyond the table, sieve half-width of 0 or above 256 chunks,
  stage 2 bounds near 2^32). They now bound their parameters, return
  `NO_FACTOR_FOUND` below 2 and `FAILURE` on squares; the bounds are in
  the headers. A tuning file can no longer crash them either.
- Factors already in the caller's arrays were taken as divisors of n, and
  full arrays got no factor at all from SIQS and CFRAC.
- rho retried forever on numbers it could not split; Fermat read
  uninitialized memory in its second pass; the GMP ECM stage 2 was wrong
  below b1 = 105.
- The programs crashed on an empty standard input, cut long numbers
  silently, and looped forever on 0.

`test_inputs` covers these cases. Factorizations and speed are unchanged.

## 22 October 2026 &mdash; command-line tools and documentation

A `tifa_factor` program factors a number with the methods the `tifa_factor()`
function chooses, and the `fermat` and `tdiv` programs now take the number
alone, as the others do. The Perl scripts of `tools/scripts` (a wrapper for the
programs, benchmark and plotting scripts, and their modules) are removed: they
had not worked as installed for a long time, and `tifa_tune` measures the
algorithms' parameters instead.

The documentation (the readme files and the headers) was brought in line with
the code. The unused `SIQS_DFLT_*` and `CFRAC_DFLT_*` macros for the sizes of
the factor base, of the trial division and of the sieve are removed: these
defaults depend on the size of the number.

## October 2026 &mdash; faster, correct, thread-safe

A series of 56 commits (`aac2494..6e6deaf`) reworked most of the library for
speed, correctness and use inside a sieve (one factorization per thread, as in
the cofactorization step of the number field sieve). The tables below compare
the code before (`aac2494`) and after (`6e6deaf`) on the same inputs.

All times are means per call on one thread of an Apple M5 Pro (Apple clang 21,
GMP 6.3, release builds, `-O3`). "Failed" counts the outputs that are not a
complete factorization into primes with multiplicities.

### `tifa_factor`

Balanced semiprimes:

| Bits | Before | After | Speed-up | Failed before |
|-----:|-------:|------:|---------:|--------------:|
|   32 |   10.7 µs |    1.80 µs |   5.9x | 0 / 300 |
|   40 |   34.6 µs |    3.84 µs |   9.0x | 0 / 300 |
|   48 |   75.9 µs |    5.56 µs |  13.7x | 0 / 300 |
|   56 |  249 µs   |    9.22 µs |  27x   | 0 / 300 |
|   64 |  856 µs   |   20.9 µs  |  41x   | 0 / 300 |
|   72 |  17.1 ms  |  166 µs    | 103x   | 133 / 300 |
|   80 |  23.0 ms  |  289 µs    |  80x   | 196 / 300 |
|   90 |  24.8 ms  |  511 µs    |  48x   | 30 / 40 |
|  100 |  40.4 ms  |  1.05 ms   |  38x   | 38 / 40 |
|  110 |  37.8 ms  |  1.94 ms   |  19x   | 36 / 40 |
|  128 |  44.7 ms  |  5.53 ms   |  8.1x  | 12 / 12 |
|  140 |  68.5 ms  |  12.1 ms   |  5.7x  | 12 / 12 |
|  160 |   139 ms  |  62.2 ms   |  2.2x  | 4 / 4 |
|  180 |   313 ms  |   173 ms   |  1.8x  | 2 / 2 |
|  200 |  1.39 s   |   963 ms   |  1.4x  | 2 / 2 |

Other inputs (the "after" code had no failure on any of them):

| Input | Bits | Before | After | Speed-up | Failed before |
|-------|-----:|-------:|------:|---------:|--------------:|
| Composites without factors below 1000 |  64 |   847 µs |  9.95 µs |  85x   | 0 / 200 |
|                                        |  80 |  6.06 ms |  99.3 µs |  61x   | 30 / 200 |
|                                        | 100 |  12.2 ms |   339 µs |  36x   | 10 / 40 |
|                                        | 128 |  11.9 ms |  1.27 ms |  9.4x  | 3 / 12 |
|                                        | 160 |  45.3 ms |  11.9 ms |  3.8x  | 1 / 4 |
| Random composites                      |  64 |  17.2 ms |   3.3 µs | 5221x  | 11 / 40 |
|                                        | 100 |  9.41 ms |   174 µs |  54x   | 9 / 40 |
|                                        | 128 |  17.7 ms |  1.01 ms |  17x   | 16 / 40 |
| Primes                                 |  64 |   383 ms |   1.7 µs | about 230000x | 0 / 3 |
|                                        | 100 |  40.3 ms |   8.0 µs | 5036x  | 0 / 3 |
|                                        | 128 |  46.0 ms |  11.3 µs | 4056x  | 0 / 3 |
|                                        | 160 |   128 ms |  17.7 µs | 7260x  | 0 / 3 |

### The algorithms

| Algorithm | Input | Before | After | Speed-up |
|-----------|-------|-------:|------:|---------:|
| SIQS | 80-bit semiprimes  |   492 µs |   262 µs | 1.9x |
|      | 100-bit semiprimes |  1.39 ms |   902 µs | 1.5x |
|      | 128-bit semiprimes |  7.83 ms |  5.27 ms | 1.5x |
|      | 160-bit semiprimes |  83.6 ms |  57.9 ms | 1.4x |
|      | 200-bit semiprimes |  1.32 s  |   950 ms | 1.4x |
| ECM, default bounds | 25-bit factor, 90-bit n  |   875 µs |  54.4 µs |  16x |
|                     | 30-bit factor, 100-bit n |  2.37 ms |   168 µs |  14x |
|                     | 35-bit factor, 120-bit n |  9.96 ms |   499 µs |  20x |
|                     | 40-bit factor, 128-bit n |  21.6 ms |  1.47 ms |  15x |
|                     | 40-bit factor, 190-bit n |  3.37 s  |  1.82 ms | 1851x |
| p-1, B1 = 10000 | 128 bits | 2.45 ms |   417 µs | 5.9x |
|                 | 200 bits | 4.92 ms |   794 µs | 6.2x |
| p+1, B1 = 10000 | 128 bits | 2.24 ms |   464 µs | 4.8x |
|                 | 200 bits | 3.80 ms |  1.03 ms | 3.7x |

The SIQS "before" outputs all failed the check: `siqs()` did not fill in the
multiplicities.

### Cofactorization

Sets of integers as a sieve leaves them: windows of consecutive integers
without prime factors below the factor base bound, of lpb + 1 bits up to the
cofactor bound, each with its true answer. A cofactor is smooth when its prime
factors are all below 2^lpb, with at most the given number of large primes.
Before, there was no cofactorization function: the test is a size and
primality filter, then `tifa_factor`. After, `tifa_cofactor` answers exactly
(`max_siqs_bits = 0`), or gives up on the composites above 2 lpb bits that ECM
does not split within the time SIQS would take (`max_siqs_bits = 2 lpb`).

| Set | Before | After, same filter | `tifa_cofactor`, exact | `tifa_cofactor`, giving up |
|-----|-------:|-------------------:|-----------------------:|---------------------------:|
| lpb 32, 2 large primes, base bound 2^23, 33-64 bits (3000, 176 smooth)  |  46.5 µs | 2.80 µs | 1.25 µs (37x) | 1.33 µs (35x) |
| lpb 31, 3 large primes, base bound 2^22, 32-92 bits (1500, 72 smooth)   | 4.02 ms  | 63.8 µs | 56.1 µs (72x) | 33.4 µs (120x) |
| lpb 36, 3 large primes, base bound 2^25, 38-108 bits (800, 27 smooth)   | 6.06 ms  |  157 µs |  133 µs (46x) | 85.3 µs (71x) |

No answer was wrong and no smooth cofactor was missed on these sets, in either
mode. On sets of smooth cofactors only, the mode giving up misses 0.25% to 0.30%
of them.

### What changed

Correctness:

- `tifa_factor` did no trial division above 26 bits, ran its p-1 and p+1
  pre-stages in the caller's mode, and left the multiplicities unset when an
  algorithm reported a complete factorization.
- SIQS seeded its polynomials from the clock, could use a prime dividing kN in
  'a', gave up on duplicate relations, kept or rejected the last positions of
  each negative-side block on the wrong sieve values, could run forever on
  unlucky parameters (it now detects a stalled collection and retries wider),
  read out of bounds with small factor bases, and wrapped around when
  interpolating decreasing table columns.
- The batch smoothness filter (CFRAC) dropped or retested candidates when its
  output filled up, and cleared integers still in use; ECM could step past the
  prime table; Fermat took its
  split modulo kN instead of n; p-1 failed when stage 1 caught every factor at
  once; `modinv_ui` shifted signed values.

Thread safety: no global state, a locked tuning cache, the caller's parameters
left untouched, per-call random curves, and a test running factorizations
concurrently (`test_threads`). Throughput scales to 14-15x on 18 threads.

New code:

- `factor_u64`: primality (deterministic, the Miller-Rabin bases run in
  lockstep) and factorization (Hart's method, then ECM three curves at a time)
  in machine words below 2^64.
- Montgomery arithmetic on 1 to 4 limbs for ECM, p-1 and p+1, up to 256 bits;
  ECM curves with torsion Z/2 x Z/8, a normalized stage 2, and bounds chosen by
  the size of the factor sought.
- `tifa_cofactor` (`cofactor.h`): the cofactorization of sieves, with a bounded
  ECM search and a strategy shared by the threads.
- Primality: `is_prime_u64` below 2^64, Baillie-PSW alone above.

SIQS: the standard Knuth-Schroeppel multiplier, polynomial switches without
divisions, sieve chunks sized to the L1 data cache, one sieve kernel for all
primes chosen by hits per chunk, a faster scan, native trial division of the
relations, a direct batch smoothness test for small factor bases, a cheaper
setup, and the defaults retuned (the large prime variation from 103 bits).

Parameter tuning (`tifa_tune`): searches both ways around the defaults with
pruning, tuning files interpolated between rows, `--max` always tuned, SIQS
tuned from 65 bits.

## 10 June 2011

TIFA has been registered in June 2011 at the French agency for software
protection ([APP](https://www.app.asso.fr/en)) with the Inter Deposit Digital
Number `IDDN.FR.001.220019.000.S.A.2011.000.31235`.
