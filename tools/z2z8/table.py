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
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.

"""Builds z2z8_curves (lib/data/src/z2z8_curves.c) from search.c's output.

Keeps, for each curve, the first point of infinite order with the fewest
bits in its coordinates; checks that 16 divides the group order modulo 40
primes; prints the 27 curves with the fewest bits (ties in search order) as
C rows: d, y0 (numerators and denominators).

Usage: python3 table.py < candidates.txt
"""
import sys
from fractions import Fraction as F
from math import isqrt


def add(P, Q, d):
    (x1, y1), (x2, y2) = P, Q
    t = d * x1 * x2 * y1 * y2
    if t in (1, -1):
        return None
    return ((x1 * y2 + y1 * x2) / (1 + t), (y1 * y2 - x1 * x2) / (1 - t))


def is_torsion(P, d):
    # Torsion points of Z/2 x Z/8 have order at most 8.
    Q = P
    for _ in range(8):
        if Q is None or Q == (0, 1):
            return True
        Q = add(Q, P, d)
    return Q is None or Q == (0, 1)


def order(d, p):
    # Montgomery form B.v^2 = u^3 + A.u^2 + u, A = 2(1 + d)/(1 - d),
    # B = 4/(1 - d): #E = p + 1 + sum over u of (B.f(u) / p).
    dm = d.numerator * pow(d.denominator, -1, p) % p
    inv = pow((1 - dm) % p, -1, p)
    A, B = 2 * (1 + dm) * inv % p, 4 * inv % p
    s = 0
    for u in range(p):
        r = B * (u * u * u + A * u * u + u) % p
        s += 0 if r == 0 else (1 if pow(r, (p - 1) // 2, p) == 1 else -1)
    return p + 1 + s


best = {}
for line in sys.stdin:
    _, d, x = line.split()
    d, x = F(d), F(x)
    y2 = (1 - x * x) / (1 - d * x * x)
    yn, yd = isqrt(y2.numerator), isqrt(y2.denominator)
    if yn * yn != y2.numerator or yd * yd != y2.denominator:
        continue
    y = F(yn, yd)
    if y == 0 or is_torsion((x, y), d):
        continue
    height = max(max(abs(v.numerator), v.denominator) for v in (d, x, y)).bit_length()
    if d not in best or height < best[d][0]:
        best[d] = (height, y)

primes = [p for p in range(1009, 1700) if all(p % q for q in range(2, 42))][:40]
curves = sorted(best.items(), key=lambda kv: kv[1][0])[:27]
for d, (_, y) in curves:
    for p in primes:
        if d.numerator % p and d.denominator % p and (d.numerator - d.denominator) % p:
            assert order(d, p) % 16 == 0, (d, p)
    print("    { %dULL, %dULL, %dULL, %dULL }," % (d.numerator, d.denominator, y.numerator, y.denominator))
