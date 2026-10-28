#!/bin/sh
#
# ab.sh <base> <new> [<args>...]
#
# Runs the commands <base> and <new> (each a program and its own arguments,
# split on spaces) with the arguments <args>, alternately over four rounds,
# each order twice, and prints the mean time per run of each and the ratio
# new/base. Both must print "...: <t> us per run", as the benchmark programs
# do. For example:
#
#   bench/ab.sh old/bench/siqs-bench new/bench/siqs-bench 70 100 3
#   bench/ab.sh build/bench/siqs-bench \
#       "build/bench/siqs-bench --threshold-delta -1" 120 30 3
#
if [ $# -lt 2 ]; then
    echo "usage: ab.sh <base> <new> [<args>...]" >&2
    exit 2
fi
base=$1
new=$2
shift 2

time_of() {
    t=$("$@" | sed -n 's/.*: \([0-9.]*\) us per run.*/\1/p')
    if [ -z "$t" ]; then
        echo "ab.sh: no time per run from: $*" >&2
        exit 1
    fi
    echo "$t"
}

tb=0
tn=0
for r in 1 2 3 4; do
    if [ $((r % 2)) -eq 1 ]; then
        b=$(time_of $base "$@") || exit 1
        n=$(time_of $new "$@") || exit 1
    else
        n=$(time_of $new "$@") || exit 1
        b=$(time_of $base "$@") || exit 1
    fi
    tb=$(awk -v s="$tb" -v t="$b" 'BEGIN { print s + t }')
    tn=$(awk -v s="$tn" -v t="$n" 'BEGIN { print s + t }')
done
awk -v b="$tb" -v n="$tn" -v a="$*" 'BEGIN {
    printf "%-24s base %10.1f us  new %10.1f us  ratio %.3f\n", a, b / 4, n / 4, n / b
}'
