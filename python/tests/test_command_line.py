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

import shutil
import subprocess
import sys
import sysconfig

import pytest

import tifa


def run(*args, stdin=None):
    return subprocess.run(
        [sys.executable, "-m", "tifa", *args],
        input=stdin, capture_output=True, text=True, timeout=120,
    )


def test_numbers_of_the_command_line():
    out = run("360", "97", "1", "0x10001", "007")
    assert out.returncode == 0, out.stderr
    assert out.stdout.splitlines() == [
        "360: 2 2 2 3 3 5",
        "97: 97",
        "1:",
        "65537: 65537",
        "7: 7",
    ]


def test_numbers_of_the_input_and_exponents():
    out = run("-e", stdin="360\n18446744073709551617  12\n")
    assert out.returncode == 0, out.stderr
    assert out.stdout.splitlines() == [
        "360: 2^3 3^2 5",
        "18446744073709551617: 274177 67280421310721",
        "12: 2^2 3",
    ]


def test_proof_and_hint():
    n = (2**89 - 1) * (2**61 - 1)
    out = run("--proof", "--hint", "balanced", str(n))
    assert out.returncode == 0, out.stderr
    assert out.stdout == f"{n}: {2**61 - 1} {2**89 - 1}\n"


def test_invalid_numbers_go_on():
    out = run("12", "abc", "-5", "15")
    assert out.returncode == 1
    assert out.stdout.splitlines() == ["12: 2 2 3", "15: 3 5"]
    assert "'abc' is not a positive integer" in out.stderr
    assert "'-5' is not a positive integer" in out.stderr


def test_timeout():
    n = (2**107 - 1) * (2**127 - 1)
    out = run("--timeout", "0.01", "--hint", "balanced", str(n), "6")
    assert out.returncode == 1
    assert out.stdout == "6: 2 3\n"
    assert "could not factor" in out.stderr
    assert run("--timeout", "0", "6").returncode == 2


def test_version():
    out = run("--version")
    assert out.returncode == 0
    assert out.stdout.strip() == f"tifa {tifa.__version__}"


@pytest.mark.parametrize("program", ["tifa-factor", "tifa-batch-gcd", "tifa-tune"])
def test_programs_are_installed(program):
    path = shutil.which(program, path=sysconfig.get_path("scripts"))
    assert path is not None
    out = subprocess.run([path, "--version"], capture_output=True, text=True)
    assert out.returncode == 0
    assert out.stdout.strip() == f"{program} (TIFA) {tifa.__version__}"


def test_installed_programs_work():
    scripts = sysconfig.get_path("scripts")
    factor = shutil.which("tifa-factor", path=scripts)
    out = subprocess.run([factor, "1000000016000000063"], capture_output=True, text=True)
    assert out.returncode == 0
    assert "1000000007" in out.stdout and "1000000009" in out.stdout
    gcd = shutil.which("tifa-batch-gcd", path=scripts)
    out = subprocess.run([gcd], input="15 21 22 143\n", capture_output=True, text=True)
    assert out.returncode == 0
    assert out.stdout.splitlines() == ["1 3", "2 3", "3 11", "4 11"]
