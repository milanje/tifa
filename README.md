# TIFA &mdash; Tools for Integer FActorization

`IDDN.FR.001.220019.000.S.A.2011.000.31235`

`TIFA` (_Tools for Integer FActorization_), is a small
C library that implements utilities and algorithms to perform integer
factorizations, with an emphasis on small to medium-size composites,
say from 40 to 200 bits.

This library was written around 2010 when I was working at
[LIX](https://www.lix.polytechnique.fr), the Computer Science Laboratory
of [Ecole Polytechnique](https://www.polytechnique.edu/en). It has been
mostly abandoned during the past few years and is probably not competitive
anymore.

The source code was uploaded "as is" to [GitHub](https://github.com) in
April 2020. It was originally built with [SCons](https://scons.org); the build
system has since been migrated to [CMake](https://cmake.org) (the old SCons
files have been removed).

The [readme.txt](./readme/readme.txt) file in the `readme` directory has more
detail on the algorithms, the source layout and the Perl tooling; its build and
installation steps use CMake as well, matching the instructions below.

## Building

Requirements:

- a C99 compiler,
- [GMP](https://gmplib.org) (>= 4.2), the sole external dependency,
- [CMake](https://cmake.org) >= 3.20.

Optional: [Doxygen](https://www.doxygen.nl) (+ LaTeX for the PDF) for the API
documentation, and Perl (with `pod2man`) plus [gnuplot](http://www.gnuplot.info)
for the benchmarking scripts.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build            # honors -DCMAKE_INSTALL_PREFIX=...
```

If GMP lives in a non-standard location, pass `-DGMP_INCDIR=<dir>` and
`-DGMP_LIBDIR=<dir>`. Build options mirror the old `BuildOptions.py` and are set
with `-D<OPTION>=<value>` (run `cmake -LH build` to list them); for example
`-DTIFA_BITSTRING_T=...`, `-DTIFA_ALLOW_VERBOSE=ON`, `-DTIFA_ALLOW_TIMING=ON`.

The build produces both a static (`libtifa.a`) and a shared (`libtifa.so` /
`libtifa.dylib`) library, plus the command-line factorization programs
(`cfrac`, `ecm`, `fermat`, `siqs`, `squfof`, `rho`, `tdiv`). Installation also
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

A [very crude documentation](https://www.lix.polytechnique.fr/Labo/Jerome.Milan/tifa/doc/tifa_user.pdf) is available as a PDF on the
[project's LIX webpage](https://www.lix.polytechnique.fr/Labo/Jerome.Milan/tifa/tifa.xhtml).

# Licence

The TIFA library is Copyright(C) 2011 CNRS, Ecole Polytechnique and INRIA.

The TIFA library is free software; you can redistribute it and/or modify it
under the terms of the
[GNU Lesser General Public License](http://www.gnu.org/licenses/lgpl.html)
as published by the Free Software Foundation; either
[version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
of the License, or (at your option) any later version.

The TIFA library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

The TIFA "puzzle cube" logo was created by Jerome Milan and is
Copyright (C) 2011 [CNRS](http://www.cnrs.fr/index.php),
[INRIA](http://en.inria.fr/) and
[Ecole Polytechnique](http://www.polytechnique.edu/en). All rights reserved.

Finally, note that TIFA has been registered in June 2011 at the european
agency for software protection ([APP](https://www.app.asso.fr/en)) with the
Inter Deposit Digital Number `IDDN.FR.001.220019.000.S.A.2011.000.31235`.
