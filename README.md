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
April 2020 to adapt the build process, which uses [SCons](https://scons.org),
to Python 3. To be honest, that crufty code shows its age and does need a good
clean up.

In the meantime, a proper [readme.txt](./readme/readme.txt) file can be found in the `readme` directory. This go over the requirements and the installation
steps.

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
