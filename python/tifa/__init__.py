# Copyright (C) 2026 Jerome Milan.
#
# This file is part of TIFA.
#
# TIFA is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""Fast factorization of integers of up to about 220 bits (66 digits).

>>> import tifa
>>> tifa.factor(360)
{2: 3, 3: 2, 5: 1}
>>> tifa.factor(2**64 + 1)
{274177: 1, 67280421310721: 1}
"""

from __future__ import annotations

import operator
from typing import Literal, Optional, SupportsIndex

from ._tifa import LIBRARY_VERSION, _factor

__all__ = ["FactorizationError", "factor"]
__version__ = LIBRARY_VERSION

Hint = Literal["balanced", "unbalanced"]

_HINTS = {None: 0, "balanced": 1, "unbalanced": 2}


class FactorizationError(RuntimeError):
    """The factorization of ``n`` could not be completed.

    ``factors`` holds what was found, as ``{factor: exponent}``: the factors
    multiply to ``n``, but some of them may not be prime.
    """

    def __init__(self, n: int, factors: dict[int, int]) -> None:
        super().__init__(f"could not factor {n} completely")
        self.n = n
        self.factors = factors


def factor(n: SupportsIndex, *, hint: Optional[Hint] = None) -> dict[int, int]:
    """Return the prime factorization of ``n`` as ``{prime: exponent}``.

    The primes come in increasing order, and ``factor(1)`` is ``{}``.

    >>> factor(360)
    {2: 3, 3: 2, 5: 1}
    >>> factor(97)
    {97: 1}

    ``n`` is a positive integer: an ``int``, or any integer type with
    ``__index__`` (NumPy's, gmpy2's ``mpz``). Up to about 220 bits, ``n`` is
    factored in seconds at most; beyond, the time keeps growing quickly.

    ``hint`` says what is known of ``n``, to save time: ``"balanced"`` if it
    is the product of two primes of about the same size (an RSA modulus),
    ``"unbalanced"`` if it has a prime factor of less than about a third of
    its size. A hint never changes the result, but a wrong one can cost much
    more time than none.

    The GIL is released while ``n`` is factored, so that threads factor in
    parallel::

        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor() as pool:
            results = list(pool.map(tifa.factor, numbers))

    Raises ``TypeError`` if ``n`` is not an integer, ``ValueError`` if it is
    not positive or ``hint`` is not one of the above, and
    ``FactorizationError`` in the rare event that ``n`` could not be factored
    completely.
    """
    n = operator.index(n)
    if n < 1:
        raise ValueError(f"n must be positive, not {n}")
    try:
        code = _HINTS[hint]
    except (KeyError, TypeError):
        raise ValueError(
            f"hint must be 'balanced', 'unbalanced' or None, not {hint!r}"
        ) from None
    if n == 1:
        return {}
    complete, pairs = _factor(n, code)
    factors: dict[int, int] = {}
    for p, e in sorted(pairs):
        factors[p] = factors.get(p, 0) + e
    if not complete:
        raise FactorizationError(n, factors)
    return factors
