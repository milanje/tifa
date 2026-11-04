# tifa

Fast factorization of integers of up to about 220 bits (66 digits), with
the [TIFA](https://github.com/milanje/tifa) C library: its quadratic sieve,
elliptic curves and smaller methods, chosen by the size of the number.

```python
>>> import tifa
>>> tifa.factor(360)
{2: 3, 3: 2, 5: 1}
>>> tifa.factor(2**64 + 1)
{274177: 1, 67280421310721: 1}
>>> tifa.factor(1)
{}
```

`factor(n)` returns the prime factorization of a positive integer as
`{prime: exponent}`, primes in increasing order, like SymPy's `factorint`.
`n` can be an `int` or any integer type with `__index__` (NumPy's, gmpy2's).

## Speed

One thread, balanced semiprimes (the hardest numbers of a size), on an AMD
EPYC 9V45:

| bits | 100 | 140 | 180 | 220 |
|---|---:|---:|---:|---:|
| `tifa.factor` | 0.7 ms | 9.6 ms | 0.18 s | 3.0 s |
| PARI/GP 2.15 `factorint` | 1.9 ms | 39 ms | 0.55 s | 9.0 s |

Numbers with small factors take much less. From 180 to 220 bits the time
doubles about every 10 bits, and keeps growing beyond.

## Hints

When something is known of `n`, a hint saves time:

```python
tifa.factor(n, hint="balanced")    # two primes of about the same size
tifa.factor(n, hint="unbalanced")  # a prime factor under a third of its size
```

A hint never changes the result, but a wrong one can cost much more time
than none.

## Threads

`factor` releases the GIL, so that threads factor in parallel:

```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor() as pool:
    results = list(pool.map(tifa.factor, numbers))
```

## Errors

`TypeError` if `n` is not an integer, `ValueError` if it is not positive
or the hint is unknown. In the rare event that `n` can not be factored
completely (or, with `proof=True`, a factor not proven prime),
`tifa.FactorizationError` is raised, its `factors` attribute holding what
was found (factors that multiply to `n`, some of them not prime).

## Proofs

The primes are told from composites by a Baillie-PSW test, which no
composite is known to pass, but which proves nothing above 2^64.
`proof=True` also proves them prime:

```python
tifa.factor(n, proof=True)
```

`prime_certificate(p)` returns the proof, which anyone can check:
`check_certificate` does, in plain Python.

```python
>>> cert = tifa.prime_certificate(2**89 - 1)
>>> cert
[(618970019642690137449562111, [(2, 3), (3, 3), (5, 3), (17, 3), (23, 3), (89, 2), (353, 3), (397, 3)])]
>>> tifa.check_certificate(2**89 - 1, cert)
True
```

Each step `(p, [(q, a), ...])` proves `p` prime from primes `q` of `p - 1`
and integers `a` (Pocklington's theorem, and Brillhart, Lehmer and
Selfridge's); a step `(p, [(q, a), ...], D, [(r, P), ...])` also uses
primes `r` of `p + 1`, with Lucas sequences of discriminant `D` and
parameters `P` (Morrison's theorem). Each `q` and `r` is below 2^64 or
proven by an earlier step; a prime below 2^64 needs no step. The
certificate is plain data: `json.dumps` saves it.

On an Apple M5 Pro, a proof takes 0.5 ms on average at 128 bits, 5.5 ms at
192 bits, 53 ms at 224 bits and 66 ms at 256 bits (up to about a second
when both `p - 1` and `p + 1` are hard to factor). On random numbers of 100
to 200 bits, `factor(n, proof=True)` takes 2% to 25% more time than
`factor(n)`. Without `proof=True`, nothing changes.

## Installing

From PyPI:

```sh
pip install tifa-factor
```

Wheels exist for Linux (x86-64 and arm64, glibc and musl) and macOS (arm64
and x86-64, 11.0 and later), for Python 3.10 and later; GMP comes in them.
Elsewhere, pip builds the package from source.

From a checkout of [TIFA](https://github.com/milanje/tifa), at the top of
the tree:

```sh
pip install ".[test]"    # or "pip install ." without pytest
python -m pytest         # the tests
```

A build from source needs a C compiler, CMake 3.26 or later (pip fetches it
if the system's is older) and GMP with its headers (`libgmp-dev` on
Debian and Ubuntu, `gmp-devel` on Fedora, `brew install gmp` on macOS). A
GMP installed elsewhere than in the usual places is found through
`GMP_ROOT`, the directory holding its `include` and `lib`:

```sh
GMP_ROOT=/opt/gmp pip install .
```

## Licence

LGPL 2.1 or later, as TIFA. The wheels contain GMP 6.3.0
([gmplib.org](https://gmplib.org)), under the LGPL 3 (readme/license_v3.txt),
linked statically: `pip install --no-binary tifa-factor tifa-factor` builds
the package against the GMP of the system instead.
