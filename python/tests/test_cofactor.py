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

import math
import pickle
import random
from concurrent.futures import ThreadPoolExecutor

import pytest

import tifa


def is_prime(n):
    return tifa.factor(n) == {n: 1}


def random_prime(lo_bits, hi_bits, rng):
    while True:
        p = rng.randrange(1 << lo_bits, 1 << hi_bits) | 1
        if is_prime(p):
            return p


def expected(n, lpb, fb_bound, max_large_primes):
    factors = tifa.factor(n)
    large = sum(e for p, e in factors.items() if p >= fb_bound)
    if any(p >= 1 << lpb for p in factors):
        return None
    if max_large_primes and large > max_large_primes:
        return None
    return factors


def sieve_leftovers(count, lpb, fb_bits, max_bits, rng):
    # Products of primes from 2^fb_bits on, some of them below 2^lpb.
    numbers = []
    while len(numbers) < count:
        n = 1
        for _ in range(rng.randint(1, 4)):
            hi = lpb if rng.random() < 0.7 else lpb + 8
            n *= random_prime(fb_bits, hi, rng)
        if n.bit_length() <= max_bits:
            numbers.append(n)
    return numbers


@pytest.mark.parametrize(
    "lpb, fb_bits, max_large_primes, max_bits",
    [(20, 12, 0, 60), (28, 18, 2, 64), (31, 22, 3, 100), (36, 25, 3, 110)],
)
def test_exact_answers(lpb, fb_bits, max_large_primes, max_bits):
    rng = random.Random(lpb)
    cofactorize = tifa.Cofactorizer(
        lpb, fb_bound=1 << fb_bits, max_large_primes=max_large_primes
    )
    numbers = sieve_leftovers(60, lpb, fb_bits, max_bits, rng)
    answers = cofactorize.many(numbers)
    assert answers == [cofactorize(n) for n in numbers]
    smooth = 0
    for n, answer in zip(numbers, answers):
        assert answer == expected(n, lpb, 1 << fb_bits, max_large_primes)
        if answer is not None:
            assert list(answer) == sorted(answer)
            assert math.prod(p**e for p, e in answer.items()) == n
            smooth += 1
    assert 0 < smooth < len(numbers)


def test_small_primes_divided_out():
    cofactorize = tifa.Cofactorizer(20)
    assert cofactorize(2**10 * 3**5 * 997 * 1048573) == {
        2: 10,
        3: 5,
        997: 1,
        1048573: 1,
    }
    assert cofactorize(2 * 1048583) is None
    assert cofactorize(1) == {}


def test_giving_up_is_never_wrong():
    rng = random.Random(2)
    cofactorize = tifa.Cofactorizer(
        31, fb_bound=1 << 22, max_large_primes=3, max_siqs_bits=62
    )
    exact = tifa.Cofactorizer(31, fb_bound=1 << 22, max_large_primes=3)
    numbers = sieve_leftovers(60, 31, 22, 92, rng)
    for answer, right in zip(cofactorize.many(numbers), exact.many(numbers)):
        assert answer is None or answer == right


def test_many_empty_and_iterables():
    cofactorize = tifa.Cofactorizer(31)
    assert cofactorize.many([]) == []
    assert cofactorize.many(iter([2147483647, True])) == [
        {2147483647: 1},
        {},
    ]


def test_bounds():
    cofactorize = tifa.Cofactorizer(
        40, fb_bound=12345, max_large_primes=2, max_siqs_bits=80
    )
    assert (
        cofactorize.lpb,
        cofactorize.fb_bound,
        cofactorize.max_large_primes,
        cofactorize.max_siqs_bits,
    ) == (40, 12345, 2, 80)
    assert repr(cofactorize) == (
        "Cofactorizer(40, fb_bound=12345, max_large_primes=2, max_siqs_bits=80)"
    )
    with pytest.raises(AttributeError):
        cofactorize.lpb = 41


@pytest.mark.parametrize(
    "args, kwargs",
    [
        ((0,), {}),
        ((64,), {}),
        ((31,), {"fb_bound": -1}),
        ((31,), {"fb_bound": 1 << 64}),
        ((31,), {"max_large_primes": -1}),
        ((31,), {"max_siqs_bits": 1 << 32}),
    ],
)
def test_bad_bounds(args, kwargs):
    with pytest.raises(ValueError):
        tifa.Cofactorizer(*args, **kwargs)


def test_bad_integers():
    cofactorize = tifa.Cofactorizer(31)
    with pytest.raises(TypeError):
        tifa.Cofactorizer(31.0)
    with pytest.raises(TypeError):
        cofactorize(2.0)
    with pytest.raises(TypeError):
        cofactorize.many([3, "5"])
    for n in (0, -7):
        with pytest.raises(ValueError):
            cofactorize(n)
        with pytest.raises(ValueError):
            cofactorize.many([3, n])


def test_pickle():
    cofactorize = tifa.Cofactorizer(31, fb_bound=1 << 22, max_large_primes=3)
    copy = pickle.loads(pickle.dumps(cofactorize))
    assert repr(copy) == repr(cofactorize)
    assert copy(2147483647 * 2147483629) == {2147483629: 1, 2147483647: 1}


def test_threads_share_one():
    rng = random.Random(3)
    cofactorize = tifa.Cofactorizer(31, fb_bound=1 << 22, max_large_primes=3)
    numbers = sieve_leftovers(40, 31, 22, 92, rng)
    alone = cofactorize.many(numbers)
    with ThreadPoolExecutor(4) as pool:
        assert list(pool.map(cofactorize, numbers)) == alone
        chunks = [numbers[i : i + 10] for i in range(0, len(numbers), 10)]
        assert sum(pool.map(cofactorize.many, chunks), []) == alone
