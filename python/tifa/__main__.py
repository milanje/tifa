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

"""python -m tifa: factor the integers of the command line or of the input."""

from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable, Iterator, Optional

import tifa


def _words(lines: Iterable[str]) -> Iterator[str]:
    for line in lines:
        yield from line.split()


def _line(n: int, factors: dict[int, int], exponents: bool) -> str:
    if exponents:
        parts = [f"{p}^{e}" if e > 1 else str(p) for p, e in factors.items()]
    else:
        parts = [str(p) for p, e in factors.items() for _ in range(e)]
    return f"{n}: {' '.join(parts)}".rstrip()


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tifa",
        description=(
            "Print the prime factors of each integer, as GNU factor does: "
            "those of the command line, or else those of the standard input, "
            "separated by spaces or newlines."
        ),
        epilog=(
            "Exit status: 0 if every integer was factored, 1 if some could "
            "not be (invalid, not factored in time, or not proven)."
        ),
    )
    parser.add_argument("numbers", nargs="*", metavar="n", help="a positive integer")
    parser.add_argument(
        "-e", "--exponents", action="store_true",
        help="print p^e rather than p repeated e times",
    )
    parser.add_argument(
        "--proof", action="store_true",
        help="also prove the primes (see tifa.prime_certificate)",
    )
    parser.add_argument(
        "--hint", choices=["balanced", "unbalanced"],
        help="what is known of the integers (see tifa.factor)",
    )
    parser.add_argument(
        "--timeout", type=float, metavar="SECONDS",
        help="give up on an integer after about this long",
    )
    parser.add_argument(
        "--version", action="version", version=f"tifa {tifa.__version__}"
    )
    args = parser.parse_args(argv)
    if args.timeout is not None and not args.timeout > 0:
        parser.error("--timeout must be positive")

    words = args.numbers if args.numbers else _words(sys.stdin)
    status = 0
    for word in words:
        try:
            try:
                n = int(word)
            except ValueError:
                n = int(word, 0)  # 0x..., 0o..., 0b...
            if n < 1:
                raise ValueError
        except ValueError:
            print(f"python -m tifa: '{word}' is not a positive integer", file=sys.stderr)
            status = 1
            continue
        try:
            factors = tifa.factor(
                n, hint=args.hint, proof=args.proof, timeout=args.timeout
            )
        except tifa.FactorizationError as error:
            print(f"python -m tifa: {error}", file=sys.stderr)
            status = 1
            continue
        print(_line(n, factors, args.exponents), flush=True)
    return status


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(130)
    except BrokenPipeError:
        # A reader that went away (python -m tifa ... | head): no traceback,
        # and none either when Python flushes stdout at exit.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(1)
