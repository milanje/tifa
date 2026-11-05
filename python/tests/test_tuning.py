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
import sysconfig

import pytest

import tifa

P, Q = 2**61 - 1, 2**89 - 1


@pytest.fixture(autouse=True)
def restore_tuning_dir():
    yield
    tifa.set_tuning_dir(None)


def test_set_tuning_dir(tmp_path):
    # A table with which SIQS gives up (its factor base is too small).
    n = P * Q
    (tmp_path / "siqs.tuning").write_text(f"{n.bit_length()} 100 1 1 1 0 0 0\n")
    tifa.set_tuning_dir(tmp_path)
    with pytest.raises(tifa.FactorizationError):
        tifa.factor(n)
    tifa.set_tuning_dir("")
    assert tifa.factor(n) == {P: 1, Q: 1}


def test_set_tuning_dir_rejects_non_paths():
    with pytest.raises(TypeError):
        tifa.set_tuning_dir(3)
    with pytest.raises(ValueError):
        tifa.set_tuning_dir("a\0b")


def test_tifa_tune_is_installed(tmp_path):
    tune = shutil.which("tifa-tune", path=sysconfig.get_path("scripts"))
    assert tune is not None
    run = subprocess.run(
        [tune, "--algo", "siqs", "--sizes", "71", "--trials", "1",
         "--out", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert run.returncode == 0, run.stderr
    assert (tmp_path / "siqs.tuning").is_file()
    tifa.set_tuning_dir(tmp_path)
    p, q = 2**35 - 31, 2**36 - 5
    assert tifa.factor(p * q) == {p: 1, q: 1}
