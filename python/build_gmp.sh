#!/usr/bin/env bash
#
# Builds GMP into $GMP_ROOT for the wheels (cibuildwheel's before-all, see
# pyproject.toml): a static, position-independent library linked into the
# extension module, with the code of every x86-64 processor (--enable-fat)
# chosen at run time.
#
set -euo pipefail

version=6.3.0
sha256=a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
prefix=${GMP_ROOT:?GMP_ROOT must name the installation directory}

if ! command -v m4 > /dev/null; then
    dnf -y install m4 || yum -y install m4 || apk add --no-cache m4
fi

work=$(mktemp -d)
cd "$work"
# gmplib.org blocks GitHub's runners.
curl -fsSL --retry 3 -o gmp.tar.xz "https://ftpmirror.gnu.org/gmp/gmp-$version.tar.xz"
if command -v sha256sum > /dev/null; then
    echo "$sha256  gmp.tar.xz" | sha256sum -c -
else
    echo "$sha256  gmp.tar.xz" | shasum -a 256 -c -
fi
tar xf gmp.tar.xz
cd "gmp-$version"

fat=
case "$(uname -m)" in
    x86_64 | amd64) fat=--enable-fat ;;
esac
# GMP 6.3.0's configure tests are not C23, GCC 15's default.
./configure CC="${CC:-cc} -std=gnu17" --prefix="$prefix" \
    --enable-static --disable-shared --with-pic $fat
jobs=$(getconf _NPROCESSORS_ONLN 2> /dev/null || sysctl -n hw.ncpu)
make -j"$jobs"
make -j"$jobs" check
make install
