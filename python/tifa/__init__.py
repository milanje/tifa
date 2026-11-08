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

"""Fast factorization of integers of up to about 280 bits (84 digits).

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
import numbers
import operator
import os
import threading
import time
from typing import Any, Iterable, Literal, Optional, SupportsIndex, Union

from ._tifa import (
    LIBRARY_VERSION,
    _batch_gcd,
    _cofactor,
    _cofactor_strategy,
    _factor,
    _prime_certificate,
    _set_tuning_dir,
)

__all__ = [
    "Cofactorizer",
    "FactorizationError",
    "FactorizationTimeout",
    "batch_gcd",
    "check_certificate",
    "factor",
    "prime_certificate",
    "set_tuning_dir",
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
        super().__init__(message or f"could not factor {_text(n)} completely")
        self.n = n
        self.factors = factors


class FactorizationTimeout(FactorizationError, TimeoutError):
    """The time limit of ``factor`` passed before the factorization was done.

    ``factors`` holds what was found, as for ``FactorizationError``.
    """


def _text(n: int) -> str:
    # str() refuses ints of more than sys.get_int_max_str_digits() digits.
    try:
        return str(n)
    except ValueError:
        return f"<{n.bit_length()}-bit integer>"


def _deadline(timeout: Optional[float]) -> Optional[float]:
    if timeout is None:
        return None
    if isinstance(timeout, bool) or not isinstance(timeout, numbers.Real):
        raise TypeError(f"timeout must be a number or None, not {timeout!r}")
    try:
        seconds = float(timeout)
    except OverflowError:
        seconds = math.inf if timeout > 0 else -math.inf
    if not seconds > 0:
        raise ValueError(f"timeout must be positive, not {timeout!r}")
    return time.monotonic() + seconds


def _threads(threads: int) -> int:
    if isinstance(threads, bool) or not isinstance(threads, numbers.Integral):
        raise TypeError(f"threads must be an integer, not {threads!r}")
    if threads < 1:
        raise ValueError(f"threads must be positive, not {threads!r}")
    return min(int(threads), 1024)


def _bits(name: str, bits: Optional[int]) -> int:
    # What the C functions take: 0 for the default.
    if bits is None:
        return 0
    if isinstance(bits, bool) or not isinstance(bits, numbers.Integral):
        raise TypeError(f"{name} must be an integer or None, not {bits!r}")
    if bits < 1:
        raise ValueError(f"{name} must be positive, not {bits!r}")
    return min(int(bits), 2**32 - 1)


def _interruptible() -> bool:
    # Python runs the signal handlers in the main thread only.
    return threading.current_thread() is threading.main_thread()


def _left(deadline: Optional[float]) -> float:
    # What the C functions take: 0 for no limit.
    if deadline is None:
        return 0.0
    return max(deadline - time.monotonic(), 1e-9)


def factor(
    n: SupportsIndex,
    *,
    hint: Optional[Hint] = None,
    proof: bool = False,
    timeout: Optional[float] = None,
    threads: int = 1,
    siqs_max_bits: Optional[int] = None,
    ecm_max_bits: Optional[int] = None,
) -> dict[int, int]:
    """Return the prime factorization of ``n`` as ``{prime: exponent}``.

    The primes come in increasing order, and ``factor(1)`` is ``{}``.

    >>> factor(360)
    {2: 3, 3: 2, 5: 1}
    >>> factor(97)
    {97: 1}

    ``n`` is a positive integer: an ``int``, or any integer type with
    ``__index__`` (NumPy's, gmpy2's ``mpz``). Up to about 220 bits, ``n`` is
    factored in seconds at most, and up to about 280 bits in minutes; beyond,
    the time keeps growing quickly. Composites of more than ``siqs_max_bits``
    bits (320 by default, about 40 minutes of SIQS) are searched by ECM
    alone, for factors of up to ``ecm_max_bits`` bits (85 by default, at
    most 100; each 5 bits more take about twice as long): ``factor`` always
    returns, but those it does not split raise ``FactorizationError``.

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

    With ``timeout``, a number of seconds, ``FactorizationTimeout`` is raised
    if ``n`` is not factored (and its factors proven, with ``proof=True``)
    within about that time: the algorithms check the time now and then, so
    that it can be exceeded, by a few tens of milliseconds at most (up to a
    second for numbers of 16000 bits, whose primality tests take longer):

    >>> factor(2**128 + 1, timeout=60)
    {59649589127497217: 1, 5704689200685129054721: 1}

    With ``threads`` above 1, SIQS collects its relations in up to that many
    threads, at most one per 100 primes of its factor base (none below about
    100 bits, 2 at 100 bits, 4 from about 120), and the ECM search before it
    shares its curves among as many. On an Apple M5 Pro, with 6 threads, a
    140-bit number takes 0.41 of the time of one thread, a 180-bit one 0.27,
    and a 240-bit one 0.22, for more CPU time in all (a quarter more at 140
    bits, 9% at 200). To factor many numbers, factoring several at once, one
    per thread, is better.

    In the main thread, Ctrl-C interrupts ``factor`` within about a
    twentieth of a second, raising ``KeyboardInterrupt`` (or whatever a
    signal handler raises). Python runs no signal handlers in other threads:
    there, ``timeout`` is the way to bound ``factor``.

    The GIL is released while ``n`` is factored, so that threads factor in
    parallel::

        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor() as pool:
            results = list(pool.map(tifa.factor, numbers))

    Raises ``TypeError`` if ``n``, ``threads``, ``siqs_max_bits`` or
    ``ecm_max_bits`` is not an integer or ``timeout`` not a number,
    ``ValueError`` if one of them is not positive or ``hint`` is not one of
    the above, ``FactorizationTimeout`` (a ``TimeoutError``) as said above,
    and ``FactorizationError`` if ``n`` could not be factored completely: a
    composite left above ``siqs_max_bits``, or, rarely, a failure (or, with
    ``proof=True``, a factor not proven prime).
    """
    n = operator.index(n)
    if n < 1:
        raise ValueError(f"n must be positive, not {_text(n)}")
    try:
        code = _HINTS[hint]
    except (KeyError, TypeError):
        raise ValueError(
            f"hint must be 'balanced', 'unbalanced' or None, not {hint!r}"
        ) from None
    deadline = _deadline(timeout)
    nthreads = _threads(threads)
    sizes = (_bits("siqs_max_bits", siqs_max_bits),
             _bits("ecm_max_bits", ecm_max_bits))
    if n == 1:
        return {}
    complete, pairs, late, too_large = _factor(
        n, code, _left(deadline), _interruptible(), nthreads, *sizes
    )
    factors: dict[int, int] = {}
    for p, e in sorted(pairs):
        factors[p] = factors.get(p, 0) + e
    if not complete:
        if late:
            raise FactorizationTimeout(
                n, factors, f"could not factor {_text(n)} in {timeout} s"
            )
        if too_large:
            raise FactorizationError(
                n,
                factors,
                f"could not factor {_text(n)} completely: composite factors"
                " too large for SIQS, without small factors for ECM",
            )
        raise FactorizationError(n, factors)
    if proof:
        for p in factors:
            if p < _SMALL:
                continue
            status, _, late = _prime_certificate(
                p, _left(deadline), _interruptible(), nthreads, *sizes
            )
            if late and status != _PROVEN:
                raise FactorizationTimeout(
                    n,
                    factors,
                    f"could not prove that {_text(p)} is prime in {timeout} s",
                )
            if status != _PROVEN:
                raise FactorizationError(
                    n, factors, f"could not prove that {_text(p)} is prime"
                )
    return factors


def set_tuning_dir(path: Union[str, bytes, os.PathLike, None]) -> None:
    """Take the parameters of the algorithms from the tuning files of ``path``.

    ``tifa-tune``, installed with the package, measures the parameters that
    suit the machine and writes them to files ``<algo>.tuning`` in a
    directory: ``tifa-tune --help`` tells how. Each file is read when its
    algorithm first needs it, and applies to the sizes it was tuned for, give
    or take 8 bits; the built-in parameters to the others.

    ``path`` holds for the whole process, all threads included, instead of the
    ``TIFA_TUNING_DIR`` environment variable. ``None`` goes back to that
    variable, and ``""`` to the built-in parameters alone.

    Raises ``TypeError`` if ``path`` is not a path or ``None``, and
    ``ValueError`` if it holds a null character.
    """
    _set_tuning_dir(None if path is None else os.fsencode(path))


def prime_certificate(
    p: SupportsIndex,
    *,
    timeout: Optional[float] = None,
    threads: int = 1,
    siqs_max_bits: Optional[int] = None,
    ecm_max_bits: Optional[int] = None,
) -> Certificate:
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
    factored; beyond, it can take much longer: with ``timeout``, a number of
    seconds, ``TimeoutError`` is raised if the proof is not done within that
    time. In the main thread, Ctrl-C interrupts it, as it does ``factor``.
    ``threads``, ``siqs_max_bits`` and ``ecm_max_bits`` are as for
    ``factor``, for the factorizations of the proof: composites above
    ``siqs_max_bits`` never go to SIQS, so that a proof always ends.

    Raises ``TypeError`` if ``p``, ``threads``, ``siqs_max_bits`` or
    ``ecm_max_bits`` is not an integer or ``timeout`` not a number,
    ``ValueError`` if ``p`` is not prime or one of the others not positive,
    ``TimeoutError`` as said above, and ``RuntimeError`` in the rare event
    that ``p`` could not be proven prime.
    """
    p = operator.index(p)
    deadline = _deadline(timeout)
    nthreads = _threads(threads)
    sizes = (_bits("siqs_max_bits", siqs_max_bits),
             _bits("ecm_max_bits", ecm_max_bits))
    status, steps, late = _prime_certificate(
        p, _left(deadline), _interruptible(), nthreads, *sizes
    )
    if late and status != _PROVEN:
        raise TimeoutError(
            f"could not prove that {_text(p)} is prime in {timeout} s"
        )
    if status == _NOT_PRIME:
        raise ValueError(f"{_text(p)} is not prime")
    if status != _PROVEN:
        raise RuntimeError(f"could not prove that {_text(p)} is prime")
    return steps


def batch_gcd(
    moduli: Iterable[SupportsIndex], *, timeout: Optional[float] = None
) -> list[int]:
    """Return the gcd of each of ``moduli`` with the product of the others.

    A gcd of 1 means that the integer shares no prime with the others; the
    integer itself, that all its primes are shared (as with a duplicate). Of
    RSA moduli, those whose gcd is neither 1 nor themselves give away a
    prime, so their factorization:

    >>> batch_gcd([15, 21, 22, 143])
    [3, 3, 11, 11]
    >>> batch_gcd([35, 35, 11])
    [35, 35, 1]

    The gcds come from Bernstein's product and remainder trees, as in
    Heninger et al.'s "Mining Your Ps and Qs" (2012): 30 s for 100,000 RSA
    moduli of 2048 bits on an Apple M5 Pro, the memory some
    log2(len(moduli)) times the size of their product (0.8 GB there). With
    ``timeout``, a number of seconds, ``TimeoutError`` is raised if the gcds
    are not done within that time; in the main thread, Ctrl-C interrupts it,
    as it does ``factor``.

    Raises ``TypeError`` if some modulus is not an integer or ``timeout`` not
    a number, ``ValueError`` if some modulus is not positive or ``timeout``
    not positive, and ``TimeoutError`` as said above.
    """
    ns = [operator.index(n) for n in moduli]
    for n in ns:
        if n < 1:
            raise ValueError(f"moduli must be positive, not {_text(n)}")
    deadline = _deadline(timeout)
    gcds, late = _batch_gcd(ns, _left(deadline), _interruptible())
    if late:
        raise TimeoutError(f"could not compute the gcds in {timeout} s")
    return gcds


class Cofactorizer:
    """Tells which integers split into primes below ``2**lpb``.

    This is the cofactorization step of a sieve (the number field sieve, or
    a quadratic sieve with large primes): of the integers left by the sieve,
    keep those whose prime factors are all below ``2**lpb``, and factor
    them. It stops as soon as the answer is known, much faster than
    ``factor``, and the work that depends only on the bounds is done once,
    here:

    >>> cofactorize = Cofactorizer(31, fb_bound=2**22, max_large_primes=3)
    >>> cofactorize(4294967311 * 4294967357)  # two primes above 2**32
    >>> cofactorize(2147483647 * 2147483629)
    {2147483629: 1, 2147483647: 1}
    >>> cofactorize.many([2147483647 * 2147483629, 4294967311 * 4294967357])
    [{2147483629: 1, 2147483647: 1}, None]

    A call returns ``{prime: exponent}``, the primes in increasing order, if
    every prime factor is below ``2**lpb`` and, with ``max_large_primes``,
    there are at most that many of them from ``fb_bound`` on (counted with
    their exponents); ``None`` otherwise. ``many`` answers for a whole list
    at once, without the per-call cost of Python (about a microsecond, when
    the C function takes 1 to 100 microseconds on the integers a sieve
    leaves).

    ``fb_bound`` says that the integers have no prime factor below it (the
    factor base bound of the sieve); below 1000 (0 if nothing is known), the
    primes below 1000 are divided out first. A wrong bound gives wrong
    answers.

    With ``max_siqs_bits``, a composite of more bits than that is never
    factored by SIQS: once ECM has searched it as long as SIQS would take,
    without finding a factor, it is given up, and the answer is ``None``
    even if it splits into small enough primes. With three large primes and
    ``max_siqs_bits = 2 * lpb``, this takes 60% to 65% of the time of exact
    answers and misses 0.25% to 0.30% of the smooth integers. 0, the
    default, gives exact answers.

    Calls release the GIL, and threads may share a ``Cofactorizer``. In the
    main thread, Ctrl-C interrupts them, as it does ``factor``.

    Raises ``TypeError`` if a bound is not an integer, and ``ValueError`` if
    ``lpb`` is not in [1, 63] or another bound is negative or too large.
    """

    __slots__ = ("_strategy", "_bounds")

    def __init__(
        self,
        lpb: SupportsIndex,
        *,
        fb_bound: SupportsIndex = 0,
        max_large_primes: SupportsIndex = 0,
        max_siqs_bits: SupportsIndex = 0,
    ) -> None:
        bounds = (
            operator.index(lpb),
            operator.index(fb_bound),
            operator.index(max_large_primes),
            operator.index(max_siqs_bits),
        )
        if not 1 <= bounds[0] <= 63:
            raise ValueError(f"lpb must be in [1, 63], not {bounds[0]}")
        if not 0 <= bounds[1] < 1 << 64:
            raise ValueError(f"fb_bound must be in [0, 2**64), not {bounds[1]}")
        for name, value in zip(("max_large_primes", "max_siqs_bits"), bounds[2:]):
            if not 0 <= value < 1 << 32:
                raise ValueError(f"{name} must be in [0, 2**32), not {value}")
        self._bounds = bounds
        self._strategy = _cofactor_strategy(*bounds)

    @property
    def lpb(self) -> int:
        return self._bounds[0]

    @property
    def fb_bound(self) -> int:
        return self._bounds[1]

    @property
    def max_large_primes(self) -> int:
        return self._bounds[2]

    @property
    def max_siqs_bits(self) -> int:
        return self._bounds[3]

    def __reduce__(self) -> tuple[Any, tuple[int, int, int, int]]:
        return (_cofactorizer, self._bounds)

    def __repr__(self) -> str:
        lpb, fb_bound, max_large_primes, max_siqs_bits = self._bounds
        return (
            f"Cofactorizer({lpb}, fb_bound={fb_bound}, "
            f"max_large_primes={max_large_primes}, max_siqs_bits={max_siqs_bits})"
        )

    def __call__(self, n: SupportsIndex) -> Optional[dict[int, int]]:
        """Return the factorization of ``n`` if it is smooth, else ``None``.

        Raises ``TypeError`` if ``n`` is not an integer, ``ValueError`` if it
        is not positive.
        """
        n = operator.index(n)
        if n < 1:
            raise ValueError(f"n must be positive, not {_text(n)}")
        return _cofactor(self._strategy, [n], _interruptible())[0]

    def many(self, ns: Iterable[SupportsIndex]) -> list[Optional[dict[int, int]]]:
        """Return what calls on each of ``ns`` would, in a list.

        Raises ``TypeError`` if some integer is not one, ``ValueError`` if
        some is not positive.
        """
        values = [operator.index(n) for n in ns]
        for n in values:
            if n < 1:
                raise ValueError(f"ns must be positive, not {_text(n)}")
        return _cofactor(self._strategy, values, _interruptible())


def _cofactorizer(
    lpb: int, fb_bound: int, max_large_primes: int, max_siqs_bits: int
) -> Cofactorizer:
    # For pickle (multiprocessing): the strategy is rebuilt from the bounds.
    return Cofactorizer(
        lpb,
        fb_bound=fb_bound,
        max_large_primes=max_large_primes,
        max_siqs_bits=max_siqs_bits,
    )


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
