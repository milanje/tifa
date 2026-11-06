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
import random

import pytest

import tifa


def direct(moduli):
    return [
        math.gcd(n, math.prod(moduli[:i] + moduli[i + 1 :]))
        for i, n in enumerate(moduli)
    ]


def random_prime(bits, rng):
    while True:
        p = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if pow(2, p - 1, p) == 1 and pow(3, p - 1, p) == 1:
            return p


def test_shared_primes():
    rng = random.Random(2012)
    primes = [random_prime(256, rng) for _ in range(60)]
    moduli = [primes[i] * primes[i + 1] for i in range(0, 60, 2)]
    moduli += [primes[0] * random_prime(256, rng), moduli[5]]
    gcds = tifa.batch_gcd(moduli)
    assert gcds == direct(moduli)
    assert gcds[0] == primes[0] and gcds[5] == moduli[5]
    assert gcds.count(1) == len(moduli) - 4


@pytest.mark.parametrize("moduli", [[], [7], [1, 1], [6, 10, 15], range(1, 40)])
def test_small(moduli):
    assert tifa.batch_gcd(moduli) == direct(list(moduli))


def test_not_positive():
    with pytest.raises(ValueError):
        tifa.batch_gcd([15, 0])
    with pytest.raises(ValueError):
        tifa.batch_gcd([15, -21])


def test_not_integers():
    with pytest.raises(TypeError):
        tifa.batch_gcd([15, 2.0])
    with pytest.raises(TypeError):
        tifa.batch_gcd(15)


def test_timeout():
    rng = random.Random(1)
    moduli = [rng.getrandbits(2048) | 1 for _ in range(2000)]
    with pytest.raises(TimeoutError):
        tifa.batch_gcd(moduli, timeout=0.001)
    with pytest.raises(ValueError):
        tifa.batch_gcd(moduli, timeout=0)
