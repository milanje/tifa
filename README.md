# TIFA &mdash; Tools for Integer FActorization

<img src="readme/tifacube.png" align="left" width="160" hspace="12" alt="TIFA logo">

`IDDN.FR.001.220019.000.S.A.2011.000.31235`

`TIFA` (_Tools for Integer FActorization_) is a small C library that implements
utilities and algorithms to perform integer factorizations, with an emphasis on
small to medium-size composites, say from 40 to 200 bits. It can be used on its
own, or inside a larger program such as a number field sieve, to factor its
cofactors one per thread.

This library was written around 2010 when I was working at
[LIX](https://www.lix.polytechnique.fr), the Computer Science Laboratory
of [Ecole Polytechnique](https://www.polytechnique.edu/en). The source code was
uploaded "as is" to [GitHub](https://github.com) in April 2020. It was
originally built with [SCons](https://scons.org); the build system has since
been migrated to [CMake](https://cmake.org).

<br clear="left">

## What's new

Version 0.2.0 (October 2026) reworks most of the library for speed,
correctness and thread safety. Against the code before, on one thread of an
Apple M5 Pro:

- `tifa_factor` is 6x to 100x faster from 32 to 110 bits, and 1.4x to 8x from
  128 to 200 bits. It now always returns complete factorizations, which it often
  failed to do above 64 bits.
- The new `tifa_cofactor` answers whether a sieve's cofactor splits into small
  enough primes 35x to 120x faster than with `tifa_factor` before.
- ECM is 14x to 20x faster, p-1 and p+1 4x to 6x, SIQS 1.4x to 1.9x, and the
  primality tests are thousands of times faster.
- The factorization functions can be called from several threads at once.

[readme/news.md](./readme/news.md) has all the measurements and the list of
changes.

## What's in there

- `tifa_factor()` (`tifa_factor.h`) factors an integer completely, choosing the
  methods by size: machine-word arithmetic below 2^64; above, trial division by
  the primes below 1000, a short ECM search, single p-1 and p+1 runs from 90
  bits, then SIQS.
- `tifa_cofactor()` (`cofactor.h`) decides whether an integer splits into primes
  below a bound, as the cofactorization step of a sieve asks (see
  [Cofactorization](#cofactorization) below).
- The algorithms themselves, each with its own header and parameters: SIQS
  (`siqs.h`), CFRAC (`cfrac.h`), ECM (`ecm.h`), Pollard's p-1 (`p_minus_1.h`),
  Williams' p+1 (`p_plus_1.h`), Pollard's rho (`rho.h`), SQUFOF (`squfof.h`),
  McKee's variant of Fermat's method (`fermat.h`) and trial division
  (`tdiv.h`).
- Below 2^64, deterministic primality testing and factorization in machine
  words (`factor_u64.h`); ECM, p-1 and p+1 in Montgomery arithmetic on one to
  four words, up to 256 bits (`mont_ecm.h`).

The [readme.txt](./readme/readme.txt) file in the `readme` directory has more
detail on the algorithms and the source layout; its build and installation
steps use CMake as well, matching the instructions below.

## Building

Requirements:

- a C99 compiler,
- [GMP](https://gmplib.org) (>= 6.0), the sole external dependency,
- POSIX threads,
- [CMake](https://cmake.org) >= 3.20.

Optional: [Doxygen](https://www.doxygen.nl) (+ LaTeX for the PDF) for the API
documentation.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build           # the test suite
cmake --install build            # honors -DCMAKE_INSTALL_PREFIX=...
```

The machine-word code (below 2^64) and the Montgomery ECM, p-1 and p+1 need a
compiler with 128-bit integers (GCC and Clang have them on 64-bit platforms);
without them, the library falls back on its GMP code.

If GMP lives in a non-standard location, pass `-DGMP_INCDIR=<dir>` and
`-DGMP_LIBDIR=<dir>`. Build options are set with `-D<OPTION>=<value>` (run
`cmake -LH build` to list them); for example `-DTIFA_BITSTRING_T=...`,
`-DTIFA_ALLOW_VERBOSE=ON`, `-DTIFA_ALLOW_TIMING=ON`.

The build produces both a static (`libtifa.a`) and a shared (`libtifa.so` /
`libtifa.dylib`) library, plus the command-line factorization programs and the
`tifa_tune` parameter-tuning tool (see [Parameter tuning](#parameter-tuning)
below). `tifa_factor` factors a number with the methods the `tifa_factor()`
function chooses:

```sh
./build/tools/tifa_factor 19409887069306060631113
```

The other programs run one algorithm each (`cfrac`, `ecm`, `fermat`, `siqs`,
`squfof`, `rho`, `p_minus_1`, `p_plus_1`, `tdiv`); given only the number, they
use the algorithm's default parameters. Installation also
provides a `tifa.pc` pkg-config file and a CMake package, so downstream projects
can consume TIFA with either:

```cmake
find_package(Tifa REQUIRED)
target_link_libraries(your_app PRIVATE Tifa::tifa)          # or Tifa::tifa_static
```

```sh
cc app.c $(pkg-config --cflags --libs tifa) -o app
```

Other targets: `cmake --build build --target doc` (Doxygen documentation, if
available) and `cpack --config build/CPackSourceConfig.cmake` (source tarball).

The [2011 user's guide](https://www.lix.polytechnique.fr/Labo/Jerome.Milan/tifa/doc/tifa_user.pdf)
(PDF, on the
[project's LIX webpage](https://www.lix.polytechnique.fr/Labo/Jerome.Milan/tifa/tifa.xhtml))
describes the library before the 2026 rework; the Doxygen documentation
(`cmake --build build --target doc`) describes the current one.

## Cofactorization

For the cofactors left by the sieve of the number field sieve,
`tifa_cofactor()` (`cofactor.h`) answers whether every prime factor is below
2^lpb, and stops as soon as the answer is known instead of factoring
completely. The bounds go into a strategy built once, which threads may share:

```c
cofactor_params_t params = {
    .lpb              = 31,         // accept prime factors below 2^31
    .fb_bound         = 1UL << 22,  // the cofactors have no smaller factor
    .max_large_primes = 3,
    .max_siqs_bits    = 62,         // above 2 lpb bits, ECM only: give up
                                    // rather than run SIQS (0: always exact)
};
cofactor_strategy_t* strategy = alloc_cofactor_strategy(&params);

if (tifa_cofactor(factors, multis, n, strategy) == COFACTOR_SMOOTH) {
    // factors and multis hold the prime factorization of n
}
free_cofactor_strategy(strategy);
```

Above `max_siqs_bits`, a composite that ECM does not split within the time SIQS
would take is given up (`COFACTOR_GAVE_UP`). With three large primes and
`max_siqs_bits = 2 lpb`, that takes 60% to 65% of the time of exact answers and
misses 0.25% to 0.30% of the smooth cofactors (see
[readme/news.md](./readme/news.md#cofactorization)).

## Parameter tuning

Each algorithm depends on numeric parameters (smoothness bounds `b1`/`b2`, curve
counts, sieve widths, factor-base sizes, ...) whose best values depend on the
bit-size of the number being factored. `set_<algo>_params_to_default(n, params)`
fills these; by default it uses built-in, size-indexed values. The `tifa_tune`
tool lets you measure better values on your own machine and feed them back to
the library.

There are three ways parameters reach the library, tried in order, with the
built-in defaults always as the final fallback:

1. **A runtime tuning file** — `<dir>/<algo>.tuning`, where `<dir>` is the one
   set by `tifa_tuning_set_dir()` (`tuning.h`, for embedders), else the
   `TIFA_TUNING_DIR` environment variable, else the directory compiled in with
   `-DTIFA_TUNING_DIR=<dir>` (default `<prefix>/share/tifa/tuning`). You can
   swap files without recompiling.
2. **Compiled-in ("baked") tables** — enabled with `-DTIFA_USE_BAKED_TUNING=ON`;
   consulted when no runtime file is found (see below).
3. **The built-in defaults** — used whenever neither of the above applies, so an
   untouched build behaves exactly as before this feature existed.

A tuning file lists rows of parameters for increasing bit-sizes. Between two
rows, most quantities (bounds, sizes, counts) are interpolated; the flags,
seeds, choices and the SIQS sieve width (whose cost jumps at each multiple of
the sieve chunk) are taken from the row below.

### Generating tuning files

```sh
cmake --build build --target tifa_tune
./build/tools/tune/tifa_tune --algo all --out ~/.tifa/tuning
```

`tifa_tune` builds test composites (products of similar-size primes), sweeps a
parameter grid around a size-dependent center for each bit-size, times each
run, and writes the combination minimizing the expected time to find a factor.
A combination is dropped as soon as its runs have taken longer in total than
the best one's. Useful options (`--help` lists
all): `--algo <names|all>` (a comma-separated list of algorithm names, or `all`;
`--help` lists the valid names), `--min`/`--max`/`--step` (bit-size range; `--max` is always tuned), `--trials`,
`--seed` (reproducible composites), `--timeout <sec>` (per-run wall-clock budget;
a run that exceeds it is killed in its child process and scored as a failure, so
a pathological parameter/input combination can never stall the tune — default
30 s, `0` disables), and `--quick` (a fast, coarse preset). More
test composites are used for small numbers (twice as many below 90 bits, four
times below 64), where runs are cheap and timing noise is otherwise larger. When
a given size takes more than a few seconds, a progress bar with an ETA is shown —
only on an interactive terminal, so redirected output stays plain. A
convenience target runs the quick preset into `<build>/tuning`:

```sh
cmake --build build --target tune-quick
```

Point the library at the generated files at runtime:

```sh
TIFA_TUNING_DIR=~/.tifa/tuning ./build/tools/siqs <number_to_factor>
```

### Baking tuned values into the library

To compile the tuned tables in (so no runtime file or environment variable is
needed), regenerate the committed baked tables and enable them:

```sh
# One-off: run a full tuning pass and regenerate lib/algo/{src,include}/tuning_baked.*
cmake -S . -B build -DTIFA_REGEN_TUNING=ON
cmake --build build --target regen-tuning

# Then build with the baked tables compiled in
cmake -S . -B build -DTIFA_USE_BAKED_TUNING=ON
cmake --build build
```

The committed `tuning_baked.c` ships empty (the accessor returns nothing, so the
built-in defaults are used); `regen-tuning` overwrites it with your machine's
measured values, which you may then commit. Tuning data is machine-specific, so
the shipped baseline is intentionally empty.

## Licence

The TIFA library is Copyright (C) 2011 CNRS, Ecole Polytechnique and INRIA.
The files added in 2026 are Copyright (C) 2026 Jerome Milan.

The TIFA library is free software; you can redistribute it and/or modify it
under the terms of the
[GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl.html)
as published by the Free Software Foundation; either
[version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
of the License, or (at your option) any later version.

The TIFA library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

The texts of the licenses are in the `readme` directory:
[LGPL 2.1](./readme/license_v2.1.txt), [LGPL 3](./readme/license_v3.txt) and
the [GPL 3](./readme/gpl_v3.txt) it refers to.

The TIFA "puzzle cube" logo was created by Jerome Milan and is
Copyright (C) 2011 [CNRS](https://www.cnrs.fr),
[INRIA](https://www.inria.fr/en) and
[Ecole Polytechnique](https://www.polytechnique.edu/en). All rights reserved.

Finally, note that TIFA has been registered in June 2011 at the French
agency for software protection ([APP](https://www.app.asso.fr/en)) with the
Inter Deposit Digital Number `IDDN.FR.001.220019.000.S.A.2011.000.31235`.
