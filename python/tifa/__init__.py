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
>>> tifa.factor(2**128 + 1, proof=True)
{59649589127497217: 1, 5704689200685129054721: 1}
>>> cert = tifa.prime_certificate(5704689200685129054721)
>>> tifa.check_certificate(5704689200685129054721, cert)
True
"""

from __future__ import annotations

import math
import operator
from typing import Iterable, Literal, Optional, Sequence, SupportsIndex

from ._tifa import LIBRARY_VERSION, _factor, _prime_certificate

__all__ = [
    "FactorizationError",
    "check_certificate",
    "factor",
    "prime_certificate",
]
__version__ = LIBRARY_VERSION

Hint = Literal["balanced", "unbalanced"]

Certificate = list[tuple[int, list[tuple[int, int]]]]

_HINTS = {None: 0, "balanced": 1, "unbalanced": 2}

_PROVEN, _NOT_PRIME = 0, 1

_SMALL = 1 << 64


class FactorizationError(RuntimeError):
    """The factorization of ``n`` could not be completed.

    ``factors`` holds what was found, as ``{factor: exponent}``: the factors
    multiply to ``n``, but some of them may not be prime.
    """

    def __init__(
        self, n: int, factors: dict[int, int], message: Optional[str] = None
    ) -> None:
        super().__init__(message or f"could not factor {n} completely")
        self.n = n
        self.factors = factors


def factor(
    n: SupportsIndex, *, hint: Optional[Hint] = None, proof: bool = False
) -> dict[int, int]:
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

    The primes are told from composites by a Baillie-PSW test, which no
    composite is known to pass. With ``proof=True``, those above 2^64 are
    also proven prime, as by ``prime_certificate``; below 2^64, the test is
    a proof already. This takes more time: up to a few times that of
    factoring ``n`` if it is a large prime, little otherwise.

    The GIL is released while ``n`` is factored, so that threads factor in
    parallel::

        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor() as pool:
            results = list(pool.map(tifa.factor, numbers))

    Raises ``TypeError`` if ``n`` is not an integer, ``ValueError`` if it is
    not positive or ``hint`` is not one of the above, and
    ``FactorizationError`` in the rare event that ``n`` could not be factored
    completely (or, with ``proof=True``, a factor not proven prime).
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
    if proof:
        for p in factors:
            if p >= _SMALL and _prime_certificate(p)[0] != _PROVEN:
                raise FactorizationError(
                    n, factors, f"could not prove that {p} is prime"
                )
    return factors


def prime_certificate(p: SupportsIndex) -> Certificate:
    """Prove that ``p`` is prime, and return the proof.

    The certificate is a list of steps ``(prime, [(q, a), ...])``, each
    proving a prime from primes ``q`` of ``prime - 1`` (below 2^64, or
    proven by an earlier step) and integers ``a``; the last step proves
    ``p``. A prime below 2^64 needs no step:

    >>> prime_certificate(97)
    []
    >>> prime_certificate(2**89 - 1)[-1][0] == 2**89 - 1
    True

    ``check_certificate`` checks it, in Python only (see its documentation
    for the theorems used). Up to about 220 bits, a proof takes a fraction
    of a second; beyond, it can take much longer, as ``p - 1`` has to be
    partly factored.

    Raises ``TypeError`` if ``p`` is not an integer, ``ValueError`` if it is
    not prime, and ``RuntimeError`` in the rare event that it could not be
    proven prime.
    """
    p = operator.index(p)
    status, steps = _prime_certificate(p)
    if status == _NOT_PRIME:
        raise ValueError(f"{p} is not prime")
    if status != _PROVEN:
        raise RuntimeError(f"could not prove that {p} is prime")
    return steps


def check_certificate(
    p: SupportsIndex, certificate: Iterable[tuple[int, Sequence[tuple[int, int]]]]
) -> bool:
    """Return whether ``certificate`` proves that ``p`` is prime.

    This is plain Python, which does not use TIFA's C code. Each step
    ``(prime, [(q_1, a_1), ...])`` holds if every ``q_i`` is below 2^64 and
    prime, or proven by an earlier step, if ``q_i`` divides ``prime - 1``,
    and if

        a_i^(prime - 1) = 1 (mod prime) and gcd(a_i^((prime - 1) / q_i) - 1,
        prime) = 1:

    then every prime factor of ``prime`` is 1 mod F, the part of
    ``prime - 1`` made of the ``q_i`` (Pocklington, 1914), and ``prime`` is
    prime if (F + 1)^2 > ``prime``, or if F^3 > ``prime`` and c1^2 - 4 c2
    is not a square, where ``prime`` = c2 F^2 + c1 F + 1 and 0 <= c1 < F
    (Brillhart, Lehmer and Selfridge, 1975). The last step must prove ``p``.
    The primes below 2^64 are those that are strong probable primes to the
    bases 2, 3, 5, ..., 37, the first 12 primes (Sorenson and Webster,
    2017).

    >>> check_certificate(2**89 - 1, prime_certificate(2**89 - 1))
    True
    >>> check_certificate(2**89 + 1, prime_certificate(2**89 - 1))
    False
    """
    try:
        p = operator.index(p)
        steps = [
            (
                operator.index(prime),
                [(operator.index(q), operator.index(a)) for q, a in pairs],
            )
            for prime, pairs in certificate
        ]
    except (TypeError, ValueError):
        return False
    if not steps:
        return p < _SMALL and _is_small_prime(p)
    if steps[-1][0] != p:
        return False
    proven: set[int] = set()
    for prime, pairs in steps:
        if not _check_step(prime, pairs, proven):
            return False
        proven.add(prime)
    return True


def _is_small_prime(n: int) -> bool:
    if n < 2:
        return False
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for b in bases:
        if n % b == 0:
            return n == b
    s = ((n - 1) & (1 - n)).bit_length() - 1
    d = (n - 1) >> s
    for b in bases:
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _check_step(p: int, pairs: list[tuple[int, int]], proven: set[int]) -> bool:
    if p < 3 or p % 2 == 0 or not pairs:
        return False
    rest = p - 1
    for q, a in pairs:
        if not (_is_small_prime(q) if q < _SMALL else q in proven):
            return False
        if rest % q != 0:
            return False
        while rest % q == 0:
            rest //= q
        x = pow(a, (p - 1) // q, p)
        if pow(x, q, p) != 1 or math.gcd(x - 1, p) != 1:
            return False
    f = (p - 1) // rest
    if (f + 1) ** 2 > p:
        return True
    if f**3 <= p:
        return False
    c2, c1 = divmod((p - 1) // f, f)
    d = c1 * c1 - 4 * c2
    return d < 0 or math.isqrt(d) ** 2 != d
