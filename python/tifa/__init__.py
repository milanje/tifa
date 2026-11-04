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
from typing import Any, Iterable, Literal, Optional, SupportsIndex, Union

from ._tifa import LIBRARY_VERSION, _factor, _prime_certificate

__all__ = [
    "FactorizationError",
    "check_certificate",
    "factor",
    "prime_certificate",
]
__version__ = LIBRARY_VERSION

Hint = Literal["balanced", "unbalanced"]

Pairs = list[tuple[int, int]]
Step = Union[tuple[int, Pairs], tuple[int, Pairs, int, Pairs]]
Certificate = list[Step]

_HINTS = {None: 0, "balanced": 1, "unbalanced": 2}

_PROVEN, _NOT_PRIME = 0, 1

_SMALL = 1 << 64
_WITNESSES = 1 << 32
_DISCRIMINANTS = 1 << 31


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
    a proof already. This takes 2% to 25% more time on random numbers of 100
    to 200 bits, but 40 ms on average for a prime of 224 bits, which
    ``factor`` alone answers in 0.1 ms.

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

    The certificate is a list of steps, the last one proving ``p``. A step
    ``(prime, [(q, a), ...])`` proves ``prime`` from primes ``q`` of
    ``prime - 1`` and integers ``a``; a step
    ``(prime, [(q, a), ...], D, [(r, P), ...])`` also uses primes ``r`` of
    ``prime + 1``, with Lucas sequences of discriminant ``D`` and parameters
    ``P``. Each ``q`` and ``r`` is below 2^64 or proven by an earlier step.
    A prime below 2^64 needs no step:

    >>> prime_certificate(97)
    []
    >>> prime_certificate(2**89 - 1)[-1][0] == 2**89 - 1
    True

    ``check_certificate`` checks it, in Python only (see its documentation
    for the theorems used). Up to 256 bits, a proof takes a fraction of a
    second on average, as parts of ``p - 1`` and ``p + 1`` have to be
    factored; beyond, it can take much longer.

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


def check_certificate(p: SupportsIndex, certificate: Iterable[Any]) -> bool:
    """Return whether ``certificate`` proves that ``p`` is prime.

    This is plain Python, which does not use TIFA's C code. In a step
    ``(prime, [(q_1, a_1), ...], D, [(r_1, P_1), ...])`` (or
    ``(prime, [(q_1, a_1), ...])``, with no ``r_j``), each ``q_i`` and
    ``r_j`` must be a prime below 2^64 or proven by an earlier step; the
    primes below 2^64 are those that are strong probable primes to the
    bases 2, 3, 5, ..., 37, the first 12 primes (Sorenson and Webster,
    2017). Let F1 be the part of ``prime - 1`` made of the ``q_i``, and F2
    that of ``prime + 1`` made of the ``r_j``. If

        a_i^(prime - 1) = 1 (mod prime) and
        gcd(a_i^((prime - 1) / q_i) - 1, prime) = 1

    for each ``i``, every prime factor of ``prime`` is 1 mod F1 (Pocklington,
    1914). If ``gcd(prime, D) = 1`` and, for each ``j``, Q_j =
    (P_j^2 - D) / 4 is an integer prime to ``prime``, and the Lucas sequence
    U of parameters P_j and Q_j has

        U_(prime + 1) = 0 (mod prime) and
        gcd(U_((prime + 1) / r_j), prime) = 1,

    every prime factor s of ``prime`` is (D / s) = +-1 mod F2 (Morrison,
    1975). Then ``prime`` is prime if (F1 + 1)^2 > ``prime``; or if F1^3 >
    ``prime`` and c1^2 - 4 c2 is not a square, where ``prime`` =
    c2 F1^2 + c1 F1 + 1 and 0 <= c1 < F1 (Brillhart, Lehmer and Selfridge,
    1975); or if G^2 > ``prime``, G = lcm(F1, F2), and c = ``prime`` mod G
    is 1, at least ``prime``, or no divisor of ``prime`` (a prime factor
    below the square root would be 1 or c mod G). The last step must prove
    ``p``, and the others primes below ``p``, one step each. As in C, the
    ``a_i`` and ``P_j`` are below 2^32, and D, not 0, below 2^31 in absolute
    value.

    >>> check_certificate(2**89 - 1, prime_certificate(2**89 - 1))
    True
    >>> check_certificate(2**89 + 1, prime_certificate(2**89 - 1))
    False
    """
    try:
        p = operator.index(p)
        steps = [_parse_step(step) for step in certificate]
    except (TypeError, ValueError):
        return False
    if not steps:
        return p < _SMALL and _is_small_prime(p)
    if steps[-1][0] != p or any(step[0] >= p for step in steps[:-1]):
        return False
    proven: set[int] = set()
    for step in steps:
        if step[0] in proven or not _check_step(*step, proven):
            return False
        proven.add(step[0])
    return True


def _parse_step(step: Any) -> tuple[int, Pairs, int, Pairs]:
    if len(step) == 2:
        prime, minus = step
        d, plus = 0, []
    elif len(step) == 4:
        prime, minus, d, plus = step
        d = operator.index(d)
        if d == 0 or abs(d) >= _DISCRIMINANTS:
            raise ValueError("D out of range")
    else:
        raise ValueError("not a step")

    def pairs(ps: Any) -> Pairs:
        out = [(operator.index(q), operator.index(a)) for q, a in ps]
        if not all(0 <= a < _WITNESSES for _, a in out):
            raise ValueError("witness out of range")
        return out

    return operator.index(prime), pairs(minus), d, pairs(plus)


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


def _part(x: int, primes: list[int]) -> Optional[int]:
    rest = x
    for q in primes:
        if q < 2 or rest % q != 0:
            return None
        while rest % q == 0:
            rest //= q
    return x // rest


def _half(x: int, n: int) -> int:
    x %= n
    return (x if x % 2 == 0 else x + n) // 2


def _lucas_u(k: int, p: int, q: int, d: int, n: int) -> int:
    u, v, qk = 1, p % n, q % n
    for bit in bin(k)[3:]:
        u, v, qk = u * v % n, (v * v - 2 * qk) % n, qk * qk % n
        if bit == "1":
            u, v = _half(p * u + v, n), _half(d * u + p * v, n)
            qk = qk * q % n
    return u


def _check_step(p: int, minus: Pairs, d: int, plus: Pairs, proven: set[int]) -> bool:
    if p < 3 or p % 2 == 0 or not (minus or plus):
        return False
    if plus and (d == 0 or math.gcd(p, d) != 1):
        return False
    f1 = _part(p - 1, [q for q, _ in minus])
    f2 = _part(p + 1, [r for r, _ in plus])
    if f1 is None or f2 is None:
        return False
    for q, _ in minus + plus:
        if not (_is_small_prime(q) if q < _SMALL else q in proven):
            return False
    for q, a in minus:
        x = pow(a, (p - 1) // q, p)
        if pow(x, q, p) != 1 or math.gcd(x - 1, p) != 1:
            return False
    for r, pp in plus:
        if (pp * pp - d) % 4 != 0:
            return False
        qq = (pp * pp - d) // 4
        if math.gcd(qq, p) != 1 or _lucas_u(p + 1, pp, qq, d, p) != 0:
            return False
        if math.gcd(_lucas_u((p + 1) // r, pp, qq, d, p), p) != 1:
            return False
    if (f1 + 1) ** 2 > p:
        return True
    if f1**3 > p:
        c2, c1 = divmod((p - 1) // f1, f1)
        disc = c1 * c1 - 4 * c2
        if disc < 0 or math.isqrt(disc) ** 2 != disc:
            return True
    g = math.lcm(f1, f2)
    if g * g > p:
        c = p % g
        return c == 1 or c >= p or p % c != 0
    return False
