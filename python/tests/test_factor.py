import doctest
import math
import os
import random
import signal
import subprocess
import sys
import time
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


def test_giving_up_raises_factorization_error(tmp_path):
    # A tuning file with which SIQS gives up (its factor base is too small).
    n = (2**61 - 1) * (2**89 - 1)
    (tmp_path / "siqs.tuning").write_text(f"{n.bit_length()} 100 1 1 1 0 0 0\n")
    code = (
        "import tifa\n"
        "try:\n"
        f"    tifa.factor({n})\n"
        "except tifa.FactorizationError as e:\n"
        f"    print(e.factors == {{{n}: 1}})\n"
    )
    env = {**os.environ, "TIFA_TUNING_DIR": str(tmp_path)}
    run = subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True
    )
    assert run.stdout.strip() == "True", run.stderr


def test_timeout():
    rng = random.Random(220)
    n = random_prime(110, rng) * random_prime(110, rng)
    start = time.monotonic()
    with pytest.raises(tifa.FactorizationTimeout) as raised:
        tifa.factor(n, timeout=0.05)
    assert time.monotonic() - start < 2
    assert isinstance(raised.value, TimeoutError)
    assert isinstance(raised.value, tifa.FactorizationError)
    assert product(raised.value.factors) == n
    assert tifa.factor(2**64 + 1, timeout=10) == {274177: 1, 67280421310721: 1}
    assert tifa.factor(2**64 + 1, timeout=math.inf) == {274177: 1, 67280421310721: 1}


def test_small_factors_of_large_numbers():
    # ECM finds them up to 512 bits; SIQS would take minutes, or forever.
    rng = random.Random(512)
    for sizes in ([45, 235], [40, 50, 300], [30, 40, 45, 385]):
        primes = [random_prime(bits, rng) for bits in sizes]
        expected = {p: 1 for p in primes}
        assert tifa.factor(math.prod(primes), timeout=60) == expected


def test_timeout_proof():
    p = random_prime(400, random.Random(400))
    with pytest.raises(TimeoutError):
        tifa.prime_certificate(p, timeout=0.001)
    with pytest.raises(tifa.FactorizationTimeout):
        tifa.factor(p * 3, proof=True, timeout=0.001)


@pytest.mark.parametrize("timeout", [0, -1, math.nan])
def test_timeout_must_be_positive(timeout):
    with pytest.raises(ValueError):
        tifa.factor(15, timeout=timeout)
    with pytest.raises(ValueError):
        tifa.prime_certificate(97, timeout=timeout)


@pytest.mark.parametrize("timeout", ["1", True, 1j])
def test_timeout_must_be_a_number(timeout):
    with pytest.raises(TypeError):
        tifa.factor(15, timeout=timeout)


class Alarm(Exception):
    pass


def test_signal_handler_interrupts():
    def ring(signum, frame):
        raise Alarm

    rng = random.Random(241)
    n = random_prime(120, rng) * random_prime(120, rng)
    old = signal.signal(signal.SIGALRM, ring)
    try:
        signal.setitimer(signal.ITIMER_REAL, 0.1)
        start = time.monotonic()
        with pytest.raises(Alarm):
            tifa.factor(n)
        assert time.monotonic() - start < 2
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)
    assert tifa.factor(2**64 + 1) == {274177: 1, 67280421310721: 1}


def test_ctrl_c_interrupts():
    rng = random.Random(240)
    n = random_prime(120, rng) * random_prime(120, rng)
    code = (
        "import tifa\n"
        "print('ready', flush=True)\n"
        "try:\n"
        f"    tifa.factor({n})\n"
        "except KeyboardInterrupt:\n"
        "    print('interrupted', flush=True)\n"
    )
    child = subprocess.Popen(
        [sys.executable, "-c", code], stdout=subprocess.PIPE, text=True
    )
    try:
        assert child.stdout.readline().strip() == "ready"
        time.sleep(0.3)
        start = time.monotonic()
        child.send_signal(signal.SIGINT)
        out, _ = child.communicate(timeout=30)
        assert time.monotonic() - start < 5
        assert out.strip() == "interrupted"
    finally:
        child.kill()
        child.wait()
