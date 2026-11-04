import json
import random

import pytest

import tifa

from test_factor import is_prime, random_prime

M89 = 2**89 - 1
M127 = 2**127 - 1


@pytest.mark.parametrize("p", [2, 3, 37, 41, 2**61 - 1, 2**64 - 59])
def test_small_primes_need_no_step(p):
    assert tifa.prime_certificate(p) == []
    assert tifa.check_certificate(p, [])


# Strong pseudoprimes to the bases up to 7, and up to 31.
@pytest.mark.parametrize("n", [0, 1, 4, 2**64 - 1, 3215031751, 3825123056546413051])
def test_small_composites(n):
    with pytest.raises(ValueError):
        tifa.prime_certificate(n)
    assert not tifa.check_certificate(n, [])


@pytest.mark.parametrize("n", [-7, 2**64 + 1, 2**128 + 1, M89 * M127])
def test_composites(n):
    with pytest.raises(ValueError):
        tifa.prime_certificate(n)


@pytest.mark.parametrize("n", [12.0, "12", None])
def test_not_an_integer(n):
    with pytest.raises(TypeError):
        tifa.prime_certificate(n)


@pytest.mark.parametrize("p", [2**64 + 13, M89, M127])
def test_known_primes(p):
    cert = tifa.prime_certificate(p)
    assert cert[-1][0] == p
    assert tifa.check_certificate(p, cert)


@pytest.mark.parametrize("bits", [65, 80, 100, 128, 160, 190, 220])
def test_random_primes(bits):
    rng = random.Random(bits)
    for _ in range(4):
        p = random_prime(bits, rng)
        cert = tifa.prime_certificate(p)
        assert tifa.check_certificate(p, cert)
        for step in cert:
            prime = step[0]
            assert is_prime(prime)
            assert all((prime - 1) % q == 0 for q, _ in step[1])
            if len(step) == 4:
                assert all((prime + 1) % r == 0 for r, _ in step[3])


def certificate_where(bits, condition):
    rng = random.Random(bits)
    while True:
        p = random_prime(bits, rng)
        cert = tifa.prime_certificate(p)
        if condition(cert):
            return p, cert


def altered(cert, step, part, pair, q=None, a=None):
    steps = [[x if i % 2 == 0 else list(x) for i, x in enumerate(s)] for s in cert]
    q0, a0 = steps[step][part][pair]
    steps[step][part][pair] = (q0 if q is None else q, a0 if a is None else a)
    return steps


def test_altered_certificates():
    p, cert = certificate_where(160, lambda c: len(c) >= 2)
    assert tifa.check_certificate(p, cert)
    q = cert[-1][1][0][0]
    assert not tifa.check_certificate(p, altered(cert, -1, 1, 0, a=1))
    assert not tifa.check_certificate(p, altered(cert, -1, 1, 0, a=0))
    assert not tifa.check_certificate(p, altered(cert, -1, 1, 0, q=q + 2))
    assert not tifa.check_certificate(p, altered(cert, 0, 1, 0, a=1))
    assert not tifa.check_certificate(p, cert[1:])
    assert not tifa.check_certificate(p, cert[:-1])
    assert not tifa.check_certificate(p, cert[::-1])
    assert not tifa.check_certificate(p + 2, cert)


def jacobi_is_minus_one(d, p):
    return pow(d % p, (p - 1) // 2, p) == p - 1


def test_steps_with_p_plus_1():
    p, cert = certificate_where(200, lambda c: len(c[-1]) == 4)
    assert tifa.check_certificate(p, cert)
    _, minus, d, plus = cert[-1]
    assert jacobi_is_minus_one(d, p)
    r, pp = plus[0]
    assert not tifa.check_certificate(p, altered(cert, -1, 3, 0, a=pp + 1))
    assert not tifa.check_certificate(p, altered(cert, -1, 3, 0, q=r + 2))
    assert not tifa.check_certificate(p, cert[:-1] + [(p, minus, 0, plus)])
    assert not tifa.check_certificate(p, cert[:-1] + [(p, minus, d + 1, plus)])
    assert not tifa.check_certificate(p, cert[:-1] + [(p, minus, d, [])])


def test_p_plus_1_alone():
    # 2^89 - 1, from p + 1 = 2^89 only: the D accepted are non-residues.
    accepted = [
        (d, pp)
        for a in range(5, 60, 2)
        for d in [a if (a // 2) % 2 == 0 else -a]
        for pp in range(1, 40, 2)
        if tifa.check_certificate(M89, [(M89, [], d, [(2, pp)])])
    ]
    assert accepted
    assert all(jacobi_is_minus_one(d, M89) for d, _ in accepted)
    d, pp = accepted[0]
    assert not tifa.check_certificate(M89 + 2, [(M89 + 2, [], d, [(2, pp)])])
    assert not tifa.check_certificate(M89, [(M89, [], d, [(2, pp + 1)])])


def test_certificate_of_a_composite():
    n = 2**64 + 1
    assert not tifa.check_certificate(n, [(n, [(2, 3)])])
    assert not tifa.check_certificate(n, [(n, [(2, 3), (3, 2)])])


@pytest.mark.parametrize(
    "cert",
    [
        None,
        12,
        [12],
        [("x", [])],
        [(M89, None)],
        [(M89, [(2,)])],
        [(M89, [])],
        [(M89, [], 5)],
        [(M89, [], 0, [(2, 1)])],
    ],
)
def test_malformed_certificates(cert):
    assert not tifa.check_certificate(M89, cert)


def test_json_round_trip():
    p, cert = certificate_where(200, lambda c: len(c[-1]) == 4)
    assert tifa.check_certificate(p, json.loads(json.dumps(cert)))


def test_factor_with_proof():
    rng = random.Random(9)
    for bits in (65, 80):
        p, q = random_prime(bits, rng), random_prime(bits + 10, rng)
        assert tifa.factor(p * q * 12, proof=True) == {2: 2, 3: 1, p: 1, q: 1}
    p = random_prime(200, rng)
    assert tifa.factor(p, proof=True) == {p: 1}
    assert tifa.factor(1, proof=True) == {}
    assert tifa.factor(2**128 + 1, proof=True) == tifa.factor(2**128 + 1)


def test_proof_is_keyword_only():
    with pytest.raises(TypeError):
        tifa.factor(15, None, True)
