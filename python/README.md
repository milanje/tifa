# tifa

Fast factorization of integers of up to about 280 bits (84 digits), with
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

Numbers with small factors take much less. From 180 bits the time doubles
about every 10 to 12 bits: on an Apple M5 Pro, about 9 s at 240 bits, 24 s
at 256 and 2 minutes at 280.

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

## Time limits

```python
try:
    factors = tifa.factor(n, timeout=60)
except tifa.FactorizationTimeout as e:
    print("not done in a minute:", e.factors)
```

`FactorizationTimeout` is a `TimeoutError`, and a `FactorizationError`:
its `factors` hold what was found, which multiply to `n`.
`prime_certificate(p, timeout=60)` raises `TimeoutError`. The limit is
approximate: the algorithms check the time now and then (SIQS every 16
chunks of its sieve, ECM after each curve), so that a call can end a little
after it, by a millisecond from 240 to 280 bits on an Apple M5 Pro, a few
tens at most.

In the main thread, Ctrl-C interrupts `factor` and `prime_certificate`
within about a twentieth of a second, raising `KeyboardInterrupt` (or what
another signal handler raises). Python runs no signal handlers in other
threads: there, a time limit is the way to bound them.

## Tuning

The parameters of the algorithms were tuned on an Apple M5 Pro and recent
x86-64 processors; others may want their own, above 240 bits most. The
package installs `tifa-tune`, which measures them on the machine and writes
them to files that `set_tuning_dir` points the package at:

```sh
tifa-tune --algo siqs --sizes 232,248,264,280 --search descent --trials 2 \
    --fit 8 --out /path/to/tuning
```

```python
tifa.set_tuning_dir("/path/to/tuning")
```

Each file applies to the sizes it was tuned for, give or take 8 bits.
`set_tuning_dir` holds for the whole process; until it is called, the
directory is that of the `TIFA_TUNING_DIR` environment variable, if set.
`set_tuning_dir("")` keeps the built-in parameters alone. `tifa-tune --help`
lists its options, which
[TIFA's README](https://github.com/milanje/tifa#parameter-tuning) explains;
the example above takes a few hours.

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

## Shared primes

`batch_gcd(moduli)` returns the gcd of each integer with the product of the
others: 1 when it shares no prime with them. Among RSA moduli, a gcd other
than 1 and the modulus itself is a prime factor of it, found without
factoring anything:

```python
>>> tifa.batch_gcd([15, 21, 22, 143])
[3, 3, 11, 11]
```

Product and remainder trees (Heninger et al., "Mining Your Ps and Qs",
2012) take 30 s for 100,000 moduli of 2048 bits on an Apple M5 Pro, and
0.8 GB; `timeout` and Ctrl-C work as for `factor`.

## Cofactorization

The cofactorization step of a sieve (the number field sieve, or a
quadratic sieve with large primes) asks which of the integers it leaves
split into primes below a bound, 2^lpb. A `Cofactorizer` answers, much
faster than `factor`, as it stops as soon as the answer is known:

```python
>>> cofactorize = tifa.Cofactorizer(31, fb_bound=2**22, max_large_primes=3)
>>> cofactorize(2147483647 * 2147483629)
{2147483629: 1, 2147483647: 1}
>>> print(cofactorize(4294967311 * 4294967357))  # primes above 2^31
None
```

`cofactorize.many(numbers)` returns the answers for a whole list, without
the cost of a Python call for each. `fb_bound` says that the integers have
no prime factor below it (the factor base bound); `max_large_primes` limits
the number of prime factors from `fb_bound` on; with `max_siqs_bits`, a
composite of more bits than that is given up (`None`) once ECM has searched
it as long as SIQS would take, which saves 35% to 40% of the time and
misses 0.25% to 0.30% of the smooth integers with three large primes and
`max_siqs_bits = 2 * lpb`. Calls release the GIL, and threads may share a
`Cofactorizer`; it pickles, for processes. On an Apple M5 Pro, on integers
of 33 to 64 bits without prime factors below 2^23 (lpb 32, two large
primes), `many` takes 1.53 us per integer, 3% more than the C function, and
calls one by one 1.92 us.

## Command line

`python -m tifa` prints the prime factors of the integers of its command
line, or else of its standard input, as GNU `factor` does:

```console
$ python -m tifa 360 0x10001
360: 2 2 2 3 3 5
65537: 65537
$ python -m tifa --exponents 360
360: 2^3 3^2 5
```

`--proof`, `--hint` and `--timeout` are those of `factor`; an integer that
is invalid or not factored is reported on the standard error, and makes the
exit status 1. The package also installs the programs of the C library
`tifa-factor`, `tifa-batch-gcd` (see [Shared primes](#shared-primes)) and
`tifa-tune` (see [Tuning](#tuning)).

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
