import doctest
import math
import random
from concurrent.futures import ThreadPoolExecutor

import pytest

import tifa


def product(factors):
    return math.prod(p**e for p, e in factors.items())


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x in (1, n - 1) or a % n == 0:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def random_prime(bits, rng):
    while True:
        p = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(p):
            return p


def check(n, factors):
    assert product(factors) == n
    assert list(factors) == sorted(factors)
    assert all(is_prime(p) and e >= 1 for p, e in factors.items())


@pytest.mark.parametrize(
    "n, expected",
    [
        (1, {}),
        (2, {2: 1}),
        (4, {2: 2}),
        (360, {2: 3, 3: 2, 5: 1}),
        (97, {97: 1}),
        (2**61 - 1, {2**61 - 1: 1}),
        (2**64, {2: 64}),
        (2**64 + 1, {274177: 1, 67280421310721: 1}),
        (3**100, {3: 100}),
        (10**30, {2: 30, 5: 30}),
    ],
)
def test_known(n, expected):
    assert tifa.factor(n) == expected


def test_small_numbers():
    for n in range(1, 3000):
        check(n, tifa.factor(n))


@pytest.mark.parametrize("bits", [40, 64, 65, 100, 128, 160])
def test_semiprimes(bits):
    rng = random.Random(bits)
    for _ in range(4):
        p, q = random_prime(bits // 2, rng), random_prime(bits - bits // 2, rng)
        assert tifa.factor(p * q) == dict(sorted({p: 1, q: 1}.items()))


def test_mixed_composites():
    rng = random.Random(7)
    for _ in range(30):
        n = rng.getrandbits(rng.randrange(2, 140)) | 1
        check(n, tifa.factor(n))


def test_large_prime_powers():
    p = random_prime(70, random.Random(3))
    assert tifa.factor(p**3 * 12) == {2: 2, 3: 1, p: 3}


@pytest.mark.parametrize("hint", ["balanced", "unbalanced", None])
def test_hints(hint):
    rng = random.Random(11)
    p, q = random_prime(50, rng), random_prime(70, rng)
    assert tifa.factor(p * q, hint=hint) == {p: 1, q: 1}


def test_index_types():
    class Index:
        def __index__(self):
            return 360

    assert tifa.factor(Index()) == {2: 3, 3: 2, 5: 1}
    assert tifa.factor(True) == {}


@pytest.mark.parametrize("n", [0, -1, -12])
def test_not_positive(n):
    with pytest.raises(ValueError):
        tifa.factor(n)


@pytest.mark.parametrize("n", [12.0, "12", None])
def test_not_an_integer(n):
    with pytest.raises(TypeError):
        tifa.factor(n)


@pytest.mark.parametrize("hint", ["Balanced", "rsa", 1, ""])
def test_bad_hint(hint):
    with pytest.raises(ValueError):
        tifa.factor(15, hint=hint)


def test_hint_is_keyword_only():
    with pytest.raises(TypeError):
        tifa.factor(15, "balanced")


def test_threads():
    rng = random.Random(5)
    numbers = [random_prime(45, rng) * random_prime(45, rng) for _ in range(32)]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(tifa.factor, numbers))
    for n, factors in zip(numbers, results):
        check(n, factors)


def test_version():
    assert tifa.__version__.count(".") == 2


def test_docstrings():
    assert doctest.testmod(tifa).failed == 0
