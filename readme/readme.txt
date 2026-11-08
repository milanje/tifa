             -----------------------------------------------------
             TIFA LIBRARY - Tools for Integer FActorization
             -----------------------------------------------------

------------------------------------------------------------------------------
0. Contents
------------------------------------------------------------------------------

     1.  About TIFA
     2.  Legal
     3.  What is it?
     4.  Requirements
     5.  How to install
     6.  Where do things get installed?
     7.  What's in there?
     8.  Documentation
     9.  Quick start
    10.  Troubleshooting

------------------------------------------------------------------------------
1. About TIFA
------------------------------------------------------------------------------

  TIFA was written around 2010 at LIX, the Computer Science Laboratory of
  Ecole Polytechnique, uploaded "as is" to GitHub in April 2020, and
  completely revamped in the fall of 2026 for speed, correctness and thread
  safety: its SIQS is now 2.3 to 5 times as fast as the 2020 code from 100 to
  220 bits, and tifa_factor() 2.2 to 3.3 times from 100 to 180 bits. The
  file news.md in this directory lists the releases, and releases/ has the
  notes of each, with the details.

  TIFA stands for (T)ools for (I)nteger (FA)ctorization. Yes, this is a rather
  dull name but I figured it would still be better than MyCfracAndQSLibrary.
  As usual, your mileage may vary. Feel free to suggest alternatives.

  Constructive criticisms, suggestions and bug reports, are of course warmly
  welcomed. Send them to: milanje at gmail dot com.

------------------------------------------------------------------------------
2. Legal
------------------------------------------------------------------------------

  The TIFA library is Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
  The files added in 2026 are Copyright (C) 2026 Jerome Milan.

  The TIFA library is free software; you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by the
  Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The TIFA library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
  for more details.

  Copies of the GNU Lesser General Public License version 2.1 and version 3 are
  given in the files license_v2.1.txt and license_v3.txt in the same directory
  as this file.

  A copy of the GNU General Public License version 3 (referenced by the GNU
  Lesser General Public License version 3) is given in the file gpl_v3.txt in
  the same directory as this file.

  Finally, note that the TIFA library has been registered in June 2011 at the
  French agency for software protection (APP - https://www.app.asso.fr) with
  the Inter Deposit Digital Number:

      IDDN.FR.001.220019.000.S.A.2011.000.31235.

------------------------------------------------------------------------------
3. What is it?
------------------------------------------------------------------------------

  The TIFA library implements utilities and algorithms to perform integer
  factorizations, with an emphasis on small to medium-sized composites, say
  from 40 to 280 bits.

  The following factorization algorithms are implemented:

      - CFRAC   (Continued FRACtion factorization)
      - ECM     (Elliptic Curve Method)
      - Fermat  (McKee's "fast" variant of Fermat's algorithm)
      - SIQS    (Self-Initializing Quadratic Sieve)
      - SQUFOF  (SQUare FOrm Factorization)
      - Rho     (Pollard's rho method)
      - PM1     (Pollard's p-1 method)
      - PP1     (Williams's p+1 method)
      - TDIV    (trial division)

  On top of them, tifa_factor() factors an integer completely, choosing the
  method by size, and tifa_cofactor() decides whether an integer splits into
  primes below a bound, as the cofactorization step of a sieve asks. Below
  2^64, primality testing and factorization run in machine words.
  tifa_prove_prime() proves that an integer is prime, with a certificate
  that tifa_check_prime_cert() (or anyone) checks, that
  tifa_prime_cert_to_str() writes as text and tifa_prime_cert_from_str()
  reads back. tifa_batch_gcd() finds, among many integers (RSA moduli, say),
  those sharing a prime with another, from product and remainder trees.
  Their options (tifa_options.h) bound the time of a call, approximately,
  let a function interrupt it, run SIQS and the ECM search before it in
  several threads, for one large number, and set the size above which
  tifa_factor() searches composites with ECM alone, so that it always ends.

  Needless to say, this library is not intended to break any factorization
  records. :-)

------------------------------------------------------------------------------
4. Requirements
------------------------------------------------------------------------------

  First of all, this software relies on the GNU Multi Precision library (GMP)
  version 6.0 or higher, which is downloadable from:

      https://gmplib.org

  It will not build if GMP is not installed on your system.

  The library also needs POSIX threads (its tuning cache is shared between
  threads, and SIQS keeps its largest arrays per thread from one run to the
  next).

  Documentation is embedded in the code and extracted with the Doxygen tool,
  available from:

      https://www.doxygen.nl

  The TIFA code matches the C99 standard so you'll obviously need a
  C99-compliant compiler. TIFA was originally developed with GCC 3.3-4.2 and
  builds cleanly with modern GCC and Clang. The machine-word code (below 2^64),
  the Montgomery ECM, p-1 and p+1, and the values, smoothness test and square
  roots of SIQS's candidates need 128-bit integers, which GCC and Clang have
  on 64-bit platforms; without them, the library falls back on its slower GMP
  code.

  Lastly, you will need the CMake build tool version 3.20 or higher, available
  at:

      https://cmake.org

  to configure and build the TIFA library. (Earlier releases used SCons; the
  build system has since been migrated to CMake.)

  The Python package (see section 5) needs Python 3.10 or higher, and CMake
  3.26 or higher to be built from source.

------------------------------------------------------------------------------
5. How to install?
------------------------------------------------------------------------------

  TIFA uses the CMake build system. Building follows the standard
  out-of-source CMake workflow:

  1) Configure the build (into a "build" directory):

         prompt> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

     Set build parameters on this command line with -D<OPTION>=<value>. Run

         prompt> cmake -LH build

     to list them with their descriptions. Commonly used options include:

         -DCMAKE_INSTALL_PREFIX=<dir>  installation prefix (default /usr/local)
         -DGMP_INCDIR=<dir>            GMP headers, if in a non-standard place
         -DGMP_LIBDIR=<dir>            GMP library, if in a non-standard place
         -DTIFA_BITSTRING_T=<type>     native C type for a string of bits
         -DTIFA_ALLOW_VERBOSE=ON       turn on per-algorithm verbose output
         -DTIFA_ALLOW_TIMING=ON        turn on per-algorithm timing messages
         -DTIFA_TUNING_DIR=<dir>       default directory of the tuning files
                                       (<prefix>/share/tifa/tuning)
         -DTIFA_BUILD_BENCH=ON         build the benchmark programs of bench/
                                       (see README.md, not installed)

  2) Compile the library and its programs:

         prompt> cmake --build build

  3) Optionally, run the test suite:

         prompt> ctest --test-dir build

  4) Optionally, generate the documentation (requires Doxygen):

         prompt> cmake --build build --target doc

  5) Finally, install the library, its headers and the programs:

         prompt> cmake --install build

  A few other useful commands:

     cmake --build build --target clean            : remove build artifacts
     cpack --config build/CPackSourceConfig.cmake  : create a tifa-X.X.X.tar.gz
                                                     source archive

  Both a static (libtifa.a) and a shared (libtifa.so / libtifa.dylib) version
  of the library are built and installed. Installation also provides a
  pkg-config file (tifa.pc) and a CMake package config, so that downstream
  projects can use TIFA through either:

     find_package(Tifa REQUIRED)
     target_link_libraries(your_app PRIVATE Tifa::tifa)   # or Tifa::tifa_static

  or, with pkg-config:

     cc app.c $(pkg-config --cflags --libs tifa) -o app

  The Python package, tifa, is installed by pip rather than by CMake: its
  tifa.factor(n) is tifa_factor() (tifa.factor(n, proof=True) also proves
  the primes), tifa.batch_gcd(moduli) is tifa_batch_gcd(), and a
  tifa.Cofactorizer gives tifa_cofactor() for one integer or a whole list at
  a time. From version 0.4.6 on, with wheels for Linux and macOS (GMP in
  them):

         prompt> pip install tifa-factor

  or, from this tree, at its top (pip builds the library with CMake):

         prompt> pip install ".[test]"
         prompt> python -m pytest

  GMP_ROOT=<dir> in the environment takes GMP from <dir>/include and
  <dir>/lib. The file python/README.md has the details.

------------------------------------------------------------------------------
6. Where do things get installed?
------------------------------------------------------------------------------

  The installation layout follows the standard GNU directory conventions
  (CMake's GNUInstallDirs), rooted at CMAKE_INSTALL_PREFIX (default
  /usr/local):

    - The TIFA libraries (static libtifa.a and shared libtifa.so/.dylib) are
      installed in <prefix>/lib.

    - TIFA header files are installed in <prefix>/include/tifa. The sources
      include their headers by bare name, so consumers add that directory to
      their include path (handled automatically by the pkg-config file and the
      CMake package). Only the public API is installed: tifa.h and the
      headers it includes, factor_u64.h and tuning.h. The shared library
      exports nothing else (the declarations marked TIFA_API); the internal
      headers stay in the source tree.

    - Programs are installed in <prefix>/bin.

    - The pkg-config file is installed in <prefix>/lib/pkgconfig, and the
      CMake package in <prefix>/lib/cmake/Tifa.

    - The library looks for tuning files in <prefix>/share/tifa/tuning by
      default; installation does not create it (see "Parameter tuning" in
      section 7).

  Run "cmake -LH build" to see all installation-related variables.

------------------------------------------------------------------------------
7. What's in there?
------------------------------------------------------------------------------

  The source tree is organized as follow:

  tifa/

      lib/  ............... Files pertaining to the TIFA library in itself

          algo/ ........... Implementation of factorization algorithms
              include/
              src/
          data/ ........... Data files and global variables
              include/
              scripts/ .... Script (re-)generating the table of primes
              src/
          utils/ .......... Other functions of the library
              include/
              src/

      tools/ .............. Factorization programs using the TIFA library

          include/
          src/
          tune/ ........... tifa-tune: the parameter-tuning tool (see the
                            "Parameter tuning" section below)
          z2z8/ ........... Search for the ECM curves of z2z8_curves.c

      python/ ............. The Python package, tifa, and its tests
      tests/ .............. The test suite (run by ctest)
      bench/ .............. Benchmark programs and ab.sh, to compare builds
      cmake/ .............. CMake modules and package templates
      readme/ ............. This file, the licenses, news.md, todo.txt...
          releases/ ....... The notes of each release (unreleased.md: master)

  The table of primes is generated by a Python 3 script:

  tifa/lib/data/scripts/genprimes.py:
  ------------------------------------

    This script generates the files tifa/lib/data/include/first_primes.h and
    tifa/lib/data/src/first_primes.c by generating the smallest prime numbers.

    Type "python3 genprimes.py --help" in the tifa/lib/data/scripts/ directory
    for more information; its first lines give the command that writes the
    two files as they are.

  The ECM curves with torsion Z/2 x Z/8, tifa/lib/data/src/z2z8_curves.c, come
  from tifa/tools/z2z8/search.c and table.py: their first lines say how to run
  them.

  Parameter tuning (the tifa-tune tool):
  --------------------------------------

    Each algorithm depends on numeric parameters (smoothness bounds b1/b2,
    curve counts, sieve widths, factor-base sizes, ...) whose best values
    depend on the bit-size of the number to factor. The function
    set_<algo>_params_to_default(n, params) fills these in; by default it uses
    built-in, size-indexed values.

    The self-contained C tool tifa/tools/tune/ (target "tifa_tune", built with
    the rest of the programs) measures better values on the current machine and
    writes them to per-algorithm text files "<algo>.tuning". The library then
    consults these values, in order, always falling back to the built-in
    defaults:

      1. a runtime file <dir>/<algo>.tuning, where <dir> is the one set by
         tifa_tuning_set_dir() (tuning.h, for embedders), else the
         TIFA_TUNING_DIR environment variable (ignored by setuid and setgid
         programs), else the directory compiled in with -DTIFA_TUNING_DIR
         (default <prefix>/share/tifa/tuning; none in the Python package);
      2. compiled-in "baked" tables, if built with -DTIFA_USE_BAKED_TUNING=ON;
      3. the built-in defaults (an untouched build behaves as before), also
         for sizes more than 8 bits beyond the rows of a table;
         tifa_tuning_set_dir("") selects them alone, skipping 1 and 2.

    Between two rows of a tuning file, most quantities (bounds, sizes,
    counts) are interpolated; the flags, seeds, choices, the SIQS sieve
    width and its extra primes in 'a' are taken from the row below.

    Typical use:

        prompt> cmake --build build --target tifa_tune
        prompt> ./build/tools/tune/tifa-tune --algo all --out ~/.tifa/tuning
        prompt> TIFA_TUNING_DIR=~/.tifa/tuning ./build/tools/tifa-siqs <number>

    "tifa-tune --help" lists all options (algorithm, bit-size range or list,
    search, fit, trials, RNG seed, per-run time limit, and a fast "--quick"
    preset). The convenience target "cmake --build build --target tune-quick"
    runs the quick preset into <build>/tuning.

    From about 240 bits, a SIQS run takes from seconds to minutes, and the
    best parameters differ between processors (the built-in sieve half widths
    took up to 9% longer than twice as wide ones on an AMD EPYC 7763 from 264
    to 280 bits). To tune a few large sizes:

        prompt> ./build/tools/tune/tifa-tune --algo siqs \
                    --sizes 232,248,264,280 --search descent --trials 2 \
                    --fit 8 --out ~/.tifa/tuning

    "--search descent" starts from the defaults and tries each value of one
    parameter at a time, the others at their best yet, for at most two passes
    (14 combinations per pass for SIQS, which then also varies its sieve
    threshold, rather than the grid's 100). "--fit 8" writes a row every 8
    bits from the smallest tuned size to the largest, each a least-squares
    quadratic on the 5 nearest tuned sizes, kept within their values; the
    tuned rows are kept as comments at the top of the file.

    To bake tuned values into the library, regenerate the committed tables and
    enable them:

        prompt> cmake -S . -B build -DTIFA_REGEN_TUNING=ON
        prompt> cmake --build build --target regen-tuning
        prompt> cmake -S . -B build -DTIFA_USE_BAKED_TUNING=ON
        prompt> cmake --build build

    The committed lib/algo/src/tuning_baked.c ships empty (so the built-in
    defaults are used); regen-tuning overwrites it with the current machine's
    measured values. Tuning data is machine-specific, hence the empty baseline.

------------------------------------------------------------------------------
8. Documentation
------------------------------------------------------------------------------

  Library files are documented by embedded Doxygen comments which can be
  extracted and compiled as either HTML, PDF, RTF or XML files by running
  "cmake --build build --target doc" (requires Doxygen). The output formats are
  selected with the -DTIFA_GENERATE_DOXYGEN_* options; run "cmake -LH build"
  for the full list.

  The Doxygen documentation covers the library and the headers the programs
  share, not how to run the programs: section 9 shows that, and each program
  prints its usage when given a wrong number of arguments.

------------------------------------------------------------------------------
9. Quick start
------------------------------------------------------------------------------

  First, you'll have to compile and install the TIFA library and its
  associated programs as previously explained. In short, issue these commands
  in a shell:

      prompt> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
      prompt> cmake --build build
      prompt> cmake --install build

  To factor a number, use the tifa-factor program: it chooses the methods
  and their parameters, as the tifa_factor() function does.

      prompt> tifa-factor 19409887069306060631113

  With --prove, it also proves the prime factors and prints their
  certificates, in JSON, as tifa.check_certificate (Python) takes them:

      prompt> tifa-factor --prove 340282366920938463463374607431768211457

  Every program prints its version with --version. The other programs run
  one algorithm each: tifa-cfrac, tifa-ecm, tifa-fermat, tifa-siqs,
  tifa-squfof, tifa-rho, tifa-tdiv, tifa-p-minus-1 and tifa-p-plus-1 (the
  Continued FRACtion, Elliptic Curve Method, McKee-Fermat, Self-Initializing
  Quadratic Sieve, SQUare FOrm Factorization, Pollard's rho, trial division,
  Pollard's p-1 and Williams' p+1 algorithms). Given only the number, they
  use the algorithm's default parameters:

      prompt> tifa-siqs 23283795980376989165117

  They also take the algorithm's parameters (all but tifa-ecm's seed),
  before the number, in the order their usage message gives (run a program
  with -h or --help to see it). Run without argument, a program asks for
  the number to factor.

  All the factorization programs exit with status 0 for a complete
  factorization (or a prime, or 1), 1 for a partial one, 2 if no factor was
  found, and 3 for invalid input, so that "tifa-factor N && ..." reads as
  expected. With --prove, a factor not proven prime makes the factorization
  partial.

  To find, among many integers (RSA moduli, say), those that share a prime
  with another, use the tifa-batch-gcd program on files of them (or on its
  standard input): decimal integers separated by spaces or newlines, '#'
  starting a comment, or hexadecimal ones with --hex:

      prompt> tifa-batch-gcd moduli.txt

  For each integer sharing a prime with another, it prints its position
  (from 1) and its gcd with the product of the others: a factor of it, or
  the integer itself if all its primes are shared, as with a duplicate.
  With --all, it prints every gcd, in order. It exits with status 0 if some
  integers share primes, 2 if none do, and 3 for invalid input.

  To choose better parameters for your machine, see the tifa-tune tool
  ("Parameter tuning" in section 7).

------------------------------------------------------------------------------
10. Troubleshooting
------------------------------------------------------------------------------

  I can't even compile TIFA!
  --------------------------

    Check TIFA's requirements at the beginning of this file. In particular,
    you should have CMake v3.20 or newer and a C99 compiler installed on your
    system, together with (of course) the GMP library.

  Can TIFA run under OS X,Y,Z?
  ----------------------------

    It is built and tested on Linux (x86-64, with GCC and Clang) and macOS
    (arm64) on each push (.github/workflows/ci.yml). No particular effort
    was made to ensure compatibility with non-unix systems (such as
    Windows), but TIFA should be pretty easy to adapt to platforms that are
    POSIX compliant.
    That said, we do not plan to support non-unix operating systems, so in
    that case you're pretty much on your own.

  Your factoring programs yield no factor!
  ----------------------------------------

    For congruences of square methods, this can happen when the number of
    congruences of square obtained is too small. In this case, chances are that
    only trivial factors were found. CFRAC and SIQS then retry by themselves
    (with more relations and other multipliers for CFRAC; with more
    relations, then a larger factor base and sieve for SIQS). If a
    factorization still fails, run the program again with a larger number of
    relations (its nrelations parameter). In practice 32 should be enough
    most of the time.

  What? Do you mean that my program ran for hours for nothing and I have to
  -------------------------------------------------------------------------
  do it all over again?
  ---------------------

    Unfortunately yes: TIFA does not save its relations, to reload them
    before collecting more. That said, it focuses on small to medium-sized
    numbers, whose factorizations take seconds to minutes (on an Apple M5
    Pro, about 2 s at 220 bits, 3 to 6 s on x86-64, and 2 minutes at 280
    bits), and SIQS collects more relations by itself when the ones it has
    do not split the number.
