                -----------------------------------------------
                TIFA LIBRARY - development version - 2011-06-10
                -----------------------------------------------

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

  This is a development version of the TIFA library, dated June 10, 2011.

  TIFA stands for (T)ools for (I)nteger (FA)ctorization. Yes, this is a rather
  dull name but I figured it would still be better than MyCfracAndQSLibrary.
  As usual, your mileage may vary. Feel free to suggest alternatives.

  Constructive criticisms, suggestions and bug reports, are of course warmly
  welcomed. Send them to: milanje at gmail dot com.

------------------------------------------------------------------------------
2. Legal
------------------------------------------------------------------------------

  The TIFA library is Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.

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
  french agency for software protection (APP - http://app.legalis.net) with
  the Inter Deposit Digital Number:

      IDDN.FR.001.220019.000.S.A.2011.000.31235.

------------------------------------------------------------------------------
3. What is it?
------------------------------------------------------------------------------

  The TIFA library implements utilities and algorithms to perform integer
  factorizations, with an emphasis on small to medium-sized composites, say
  from 40 to 200 bits.

  For the time being, the following factorization algorithms have been
  implemented:

      - CFRAC   (Continued FRACtion factorization)
      - ECM     (Elliptic Curve Method)
      - Fermat  (McKee's "fast" variant of Fermat's algorithm)
      - SIQS    (Self-Initializing Quadratic Sieve)
      - SQUFOF  (SQUare FOrm Factorization)
      - Rho     (Pollard's rho method)
      - PM1     (Pollard's p-1 method)
      - PP1     (Williams's p+1 method)
      - TDIV    (trial division)

  Needless to say, this library is not intended to break any factorization
  records. :-)

------------------------------------------------------------------------------
4. Requirements
------------------------------------------------------------------------------

  First of all, this software relies on the GNU Multi Precision library (GMP)
  version 4.2 or (preferably) higher, which is downloadable from:

      http://gmplib.org

  It will not build if GMP is not installed on your system.

  To correctly use the scripts provided in the tifa/tools directory, you will
  also need the GMP Perl modules. They are available as part of the GMP
  distribution but not built by default, so if you did not compile the GMP
  library by yourself, chances are that these modules were not generated. In
  this case, you'll have to head over to http://gmplib.org to download
  a whole new distribution and generate the modules. Make sure that the
  version you download matches the one installed on your system, otherwise
  the modules could fail to compile.

  Needless to say, using the Perl scripts requires a distribution of Perl
  available from:

      http://www.perl.com/download.csp

  The code has been developed and tested using Perl 5.8.x.

  Additionally, the plotting Perl script in tifa/tools/scripts/ relies on
  Gnuplot to generate plots. Gnuplot is available from:

      http://www.gnuplot.info

  Documentation is embedded in the code and extracted with the Doxygen tool,
  available from:

      http://www.doxygen.org

  The TIFA code matches the C99 standard so you'll obviously need a
  C99-compliant compiler. TIFA was originally developed with GCC 3.3-4.2 and
  builds cleanly with modern GCC and Clang.

  Lastly, you will need the CMake build tool version 3.20 or higher, available
  at:

      https://cmake.org

  to configure and build the TIFA library. (Earlier releases used SCons; the
  build system has since been migrated to CMake.)

------------------------------------------------------------------------------
5. How to install?
------------------------------------------------------------------------------

  TIFA uses the CMake build system. Building follows the standard
  out-of-source CMake workflow:

  1) Configure the build (into a "build" directory):

         prompt> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

     Build parameters are set on this command line with -D<OPTION>=<value>;
     they are the successors of the old BuildOptions.py entries. Run

         prompt> cmake -LH build

     to list them with their descriptions. Commonly used options include:

         -DCMAKE_INSTALL_PREFIX=<dir>  installation prefix (default /usr/local)
         -DGMP_INCDIR=<dir>            GMP headers, if in a non-standard place
         -DGMP_LIBDIR=<dir>            GMP library, if in a non-standard place
         -DTIFA_BITSTRING_T=<type>     native C type for a string of bits
         -DTIFA_ALLOW_VERBOSE=ON       turn on per-algorithm verbose output
         -DTIFA_ALLOW_TIMING=ON        turn on per-algorithm timing messages

  2) Compile the library and its programs:

         prompt> cmake --build build

  3) Optionally, generate the documentation (requires Doxygen):

         prompt> cmake --build build --target doc

  4) Finally, install the library, programs, Perl modules/scripts and man
     pages:

         prompt> cmake --install build

  A few other useful commands:

     cmake --build build --target man              : (re)generate the man pages
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
      CMake package).

    - Programs and the Perl scripts are installed in <prefix>/bin.

    - Configuration files for the scripts go in <prefix>/share/tifa/conf.

    - Perl modules are installed in <prefix>/share/tifa/perl (override with
      -DTIFA_INSTALL_PERLMODDIR=<dir>).

    - Man pages are installed in <prefix>/share/man/man3.

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
              scripts/ .... Script to (re-)generates the data files
              src/
          utils/ .......... Other functions of the library
              include/
              src/

      tools/ .............. Factorization programs using the TIFA library

          include/
          scripts/ ........ Various perl scripts: wrappers for the C programs,
                            benchmarks scripts, plotting scripts (need gnuplot)

              Tifa/ ....... Perl modules used by the scripts and associated
                            man pages
          src/
          tune/ ........... tifa_tune: the parameter-tuning tool (see the
                            "Parameter tuning" section below)

  As you probably guessed, in addition to the library per se, this distribution
  includes a few Perl scripts (requiring Perl version 5.8 or higher) that I
  found helpful during the development stage.

  These scripts are:

  tifa/lib/data/scripts/genprimes.pl:
  ------------------------------------

    This script generates the files tifa/lib/data/include/first_primes.h and
    tifa/lib/data/src/first_primes.c by generating the smallest prime numbers.

    Type "./genprimes.pl --help" in the tifa/lib/data/scripts/ directory for
    more information.

  tifa/tools/scripts/factorize.pl:
  --------------------------------

    This is merely a Perl wrapper to make the use of TIFA's factoring programs
    more user friendly.

    Type "./factorize.pl --help" in the tifa/tools/scripts/ directory for more
    information.

  tifa/tools/scripts/benchmarker.pl:
  ----------------------------------

    This script creates another Perl script aimed at benchmarking TIFA's
    factoring programs. The generated benchmark script is not automatically
    executed and should therefore be launched manually.

    Type "./benchmarker.pl --help" in the tifa/tools/scripts/ directory
    for more information.

  tifa/tools/scripts/extractres.pl:
  ---------------------------------

    Extracts timing information from trace files to complete result files.

    Type "./extractres.pl --help" in the tifa/tools/scripts/ directory
    for more information.

  tifa/tools/scripts/plotmaker.pl:
  --------------------------------

    plotmaker.pl takes as input a result file from a script generated by
    tifa//tools/scripts/benchmarker.pl and creates several plots according to
    different parameters.

    For the time being, only gnuplot (http://www.gnuplot.info/) is supported.
    ROOT support (http://root.cern.ch/) is somewhat considered for future
    versions although you should definitely not hold your breath waiting
    for it.

    Type "./plotmaker.pl --help" in the tifa/tools/scripts/ directory for more
    information.

  These scripts rely on the following perl modules in the directory
  tifa/tools/scripts/Tifa:

  tifa/tools/scripts/Tifa/Bencher.pm:
  -----------------------------------

    Implements a general purpose benchmark framework and is used by the
    benchmarker.pl script to benchmark the various factorization programs for
    a wide selection of parameters values.

    Type "man Tifa::Bencher" for more information.

  tifa/tools/scripts/Tifa/DataDescriptor.pm:
  ------------------------------------------

    Writes and reads data description/entries in a file according to a
    TIFA specific text file format.

    Type "man Tifa::DataDescriptor" for more information.

  tifa/tools/scripts/Tifa/FormatConverter.pm:
  -------------------------------------------

    A mere wrapper for format conversion utilities.

    Type "man Tifa::FormatConverter" for more information.

  tifa/tools/scripts/Tifa/GnuPlotter.pm:
  --------------------------------------

    Generates 2D plots from 3D data with gnuplot.

    Type "man Tifa::GnuPlotter" for more information.

  tifa/tools/scripts/Tifa/NumberGenerator.pm:
  -------------------------------------------

    Generates prime or composite integers using GMP.

    Type "man Tifa::NumberGenerator" for more information.

  tifa/tools/scripts/Tifa/Program.pm:
  -----------------------------------

    Provides an abstraction of a 'launchable' and 'benchmarkable' command
    line program.

    Type "man Tifa::Program" for more information.

  tifa/tools/scripts/Tifa/ProgramRepository.pm:
  ---------------------------------------------

    A repository of all available (factoring) Tifa::Program's.

    Type "man Tifa::ProgramRepository" for more information.

  tifa/tools/scripts/Tifa/Rotor.pm:
  ---------------------------------

    Implements an odometer-like counter that can be used to generate all
    possible n-tuples from n sets of values.

    Type "man Tifa::Rotor" for more information.

  tifa/tools/scripts/Tifa/SimpleConfigReader.pm:
  ----------------------------------------------

    A truly minimalist configuration file reader used by TIFA's perl scripts.

    Type "man Tifa::SimpleConfigReader" for more information.

  Parameter tuning (the tifa_tune tool):
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

      1. a runtime file <dir>/<algo>.tuning, where <dir> comes from the
         TIFA_TUNING_DIR environment variable (or tifa_tuning_set_dir() for
         embedders);
      2. compiled-in "baked" tables, if built with -DTIFA_USE_BAKED_TUNING=ON;
      3. the built-in defaults (an untouched build behaves as before).

    Typical use:

        prompt> cmake --build build --target tifa_tune
        prompt> ./build/tools/tune/tifa_tune --algo all --out ~/.tifa/tuning
        prompt> TIFA_TUNING_DIR=~/.tifa/tuning ./build/tools/siqs <number>

    "tifa_tune --help" lists all options (algorithm, bit-size range, trials,
    RNG seed, and a fast "--quick" preset). The convenience target
    "cmake --build build --target tune-quick" runs the quick preset into
    <build>/tuning.

    To bake tuned values into the library, regenerate the committed tables and
    enable them:

        prompt> cmake -S . -B build -DTIFA_REGEN_TUNING=ON
        prompt> cmake --build build --target regen-tuning
        prompt> cmake -S . -B build -DTIFA_USE_BAKED_TUNING=ON
        prompt> cmake --build build

    The committed lib/algo/src/tuning_baked.c ships empty (so the built-in
    defaults are used); regen-tuning overwrites it with the current machine's
    measured values. Tuning data is machine-specific, hence the empty baseline.

    tifa_tune supersedes the older Perl benchmarking scripts for parameter
    selection specifically; those scripts remain useful for ad-hoc benchmarking
    and gnuplot visualization.

------------------------------------------------------------------------------
8. Documentation
------------------------------------------------------------------------------

  Library files are documented by embedded Doxygen comments which can be
  extracted and compiled as either HTML, PDF, RTF or XML files by running
  "cmake --build build --target doc" (requires Doxygen). The output formats are
  selected with the -DTIFA_GENERATE_DOXYGEN_* options; run "cmake -LH build"
  for the full list.

  The C programs in the directory tifa/tools don't currently have their
  documentation generated by Doxygen as it would not be particularly helpful.
  Documentation for these programs should probably be written independently and
  will hopefully come soon.

  Perl scripts included in this distribution are not documented per se,
  although they are quite extensively commented and have all a hopefully
  verbose enough --help option. They will also be the subject of doxygen
  independent documentation soon.

  Finally the Tifa::* perl modules are documented by man pages generated from
  embedded comments by pod2man.

------------------------------------------------------------------------------
9. Quick start
------------------------------------------------------------------------------

  First, you'll have to compile and install the TIFA library and its
  associated programs as previously explained. In short, issue these commands
  in a shell:

      prompt> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
      prompt> cmake --build build
      prompt> cmake --install build

  There are five programs that you will be likely to use: cfrac, ecm,
  fermat, siqs and squfof which as their names imply, factor integers using
  respectively the Continued FRACtion, Elliptic Curve Method, McKee-Fermat,
  Self-Initializing Quadratic Sieve and SQUare FOrm Factorization algorithms.
  However, it is recommended to use the provided wrapper script factorize.pl
  instead of invoking directly these programs. Indeed this perl wrapper
  provides a more friendly user interface and offers more options than the
  C programs.

  Using the factorize.pl script: factoring with CFRAC, SQUFOF, SIQS, etc.
  -----------------------------------------------------------------------

     This script wraps the cfrac, ecm, fermat, siqs and squfof programs.
     Typing "factorize.pl --help" on the command line gives the list of the
     available options.

     > Example: factor a number with CFRAC using the /script/ default parameter
                values

       factorize.pl --exe=cfrac 19409887069306060631113

     > Example: factor a number with SIQS using custom parameter values

       Alternatively, you can specify some custom parameter values. Parameters
       not defined will be set to the script's default values.

       factorize.pl --exe=siqs
                    --sieve_half_width=200000 \
                    --nprimes_in_factor_base=256 \
                    --nprimes_tdiv_smooth_nb=256 \
                    --nrelations=32 \
                    --linalg_method=0 \
                    --use_large_primes \
                    --nprimes_tdiv=256 \
                    23283795980376989165117

     > Example: factor a number using "optimal" parameters with CFRAC

       By default, the factoring programs determine by themselves the values of
       the parameters according to the size of the number to factor. These
       values are somewhat optimal provided that the size of the number is less
       than 200 bits.

       factorize.pl --exe=cfrac --use_defaults 31418716710282142372441

     > Example: factor with SIQS a list of numbers given in a file

       You can also input a whole list of numbers to factor by writing them
       on a file (comments beginning by # are allowed). In this case, you
       should also provide an output directory where the traces of the C
       program executions will be saved.

       factorize.pl --exe=siqs --use_defaults \
                    --number_file=numbers.txt --outdir=./traces

     > Example: use a configuration file

       Since it can quickly become annoying to type all the options directly on
       the command line, factorize.pl can optionally read them in a
       configuration file with a "<option_name> = <value>" kind of syntax.

       factorize.pl --conf=factorize.conf 40705262292383555756957

       Refer to the provided factorize.conf file for a detailed example of
       a configuration file.

  Using the benchmarker.pl script: Execute some timing benchmarks
  ---------------------------------------------------------------

     The benchmarker.pl script is used to generate timing benchmarks for TIFA's
     implementation of the factorization algorithms. More precisely it generates
     (but does not execute!) another Perl script which will (once manually
     launched) perform the benchmarks.

     The benchmarking script is generated in the following way:

         1. benchmarker.pl generates a whole range of numbers to factors
            according to some parameters described in a configuration file, or,
            optionally, simply reads these numbers from a file.

         2. benchmarker.pl then reads ranges for the factorization algorithm
            parameters and creates all possible combinations of these
            parameters' values.

         3. benchmarker.pl then creates a benchmarking Perl script which can
            then be used to perform the factorization on all the generated (or
            read) numbers, each number being factorized several times, once for
            each possible combination of the algorithm parameters' values.

     Type "benchmarker.pl --help" for more information.

     > Example:

       benchmarker.pl --conf=benchmarker.conf --outdir=results \
                      --macro=actual_bench_macro.pl

       Refer to the provided benchmarker.conf file for a detailed example of
       a configuration file.

   Using the plotmaker.pl script: Plot timing results
   --------------------------------------------------

      Once the benchmarks are completed, it is possible to plot timings as a
      function of time for all the combination of the algorithm's parameters
      thanks to the plotmaker.pl script.

      The plotmaker.pl script takes as input the result file produced by the
      benchmarking perl script generated by benchmarker.pl and produces plots
      in several formats (png, eps or pdf). For the time being, plotmaker.pl
      relies on gnuplot (http://www.gnuplot.info/) for the actual plotting
      engine. Due to the modular conception of this script, it can be extended
      to create graphics with other plotting programs.

      Type "plotmaker.pl --help" for more information.

      > Example:

        plotmaker.pl --program=gnuplot --outdir=plots --in=results.txt \
                     --format=pdf

------------------------------------------------------------------------------
10. Troubleshooting
------------------------------------------------------------------------------

  I can't even compile TIFA!
  --------------------------

    Check TIFA's requirements at the beginning of this file. In particular,
    you should have CMake v3.20 or newer and a C99 compiler installed on your
    system, together with (of course) the GMP library.

  The compilation fails: the linker complains about duplicate symbols!
  --------------------------------------------------------------------

    Let me guess: using Mac OS X? Apparently there is a "problem" with
    inlined functions with the linker... Ugly temporary workaround: comment
    the line

        #define __GMP_EXTERN_INLINE      extern __inline__

    in gmp.h. Hey, don't look at me this way! I told you this was ugly!

  The Perl scripts do not seem to work!
  -------------------------------------

    Please check:

      - that the directory containing these scripts is included in your
        PATH environment variable,

      - that TIFA's Perl modules can be found by your perl interpreter.
        Check and/or modify your PERL5LIB environment variable accordingly,

      - that the path to the C program is correct. Use the --exe option
        or edit the script's configuration files to set the full path to
        the factoring programs.

      - that you indeed have a not too old version of Perl installed on your
        system! :-)

  Can TIFA run under OS X,Y,Z?
  ----------------------------

    Maybe. But then again, maybe not. No particular effort was made to ensure
    compatibility with non-unix systems (such as Windows) but our bet is that
    TIFA should be pretty easy to adapt to non-unix platforms provided that
    they are POSIX compliant. That said, we do not plan to support non-unix
    operating systems anyway, so in that case you're pretty much on your own.

  Your factoring programs yield no factor!
  ----------------------------------------

    For congruences of square methods, this can happen when the number of
    congruences of square obtained is too small. In this case, chances are that
    only trivial factors were found. Try to run the program again with a larger
    number of relations using the --nrelations option. In practice 32 should be
    enough most of the time.

  What? Do you mean that my program ran for hours for nothing and I have to
  -------------------------------------------------------------------------
  do it all over again?
  ---------------------

    Unfortunately yes. The solution to this problem is obvious: save the data
    (essentially the matrix and relations) obtained during the course of the
    program and reload it before generating new relations. This will be solved
    in the next release... or the one after... That said, TIFA focuses on small
    to medium-sized numbers so factorizations should not take more than a few
    seconds...
