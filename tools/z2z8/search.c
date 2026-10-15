// Candidate curves and points for z2z8_curves in lib/utils/src/factor_u64.c.
//
// Edwards curves x^2 + y^2 = 1 + d.x^2.y^2 with torsion Z/2 x Z/8: a point
// (x8, x8) of order 8 needs d = (2.x8^2 - 1)/x8^4, and d is a square (full
// 2-torsion) for x8 = (k^2 - 2k + 2)/(k^2 - 2), k = p/q. For each such d,
// prints up to 3 points with x = a/b and y^2 = (1 - x^2)/(1 - d.x^2) a
// square; table.py drops the torsion points.
//
// Usage: cc -O2 search.c -lm -o search && ./search 30 3000 > candidates.txt

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef __int128 i128;

static uint64_t gcd(uint64_t a, uint64_t b) {
    while (b != 0) {
        const uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int is_square(const u128 x) {
    u128 y = (u128) sqrtl((long double) x);
    while (y * y > x) {
        y--;
    }
    while ((y + 1) * (y + 1) <= x) {
        y++;
    }
    return y * y == x;
}

int main(int argc, char** argv) {
    const int P = atoi(argv[1]);
    const int A = atoi(argv[2]);

    for (int p = -P; p <= P; p++) {
        for (int q = 1; q <= P; q++) {
            if (gcd((uint64_t) abs(p), (uint64_t) q) != 1) {
                continue;
            }
            const i128 xn = (i128) p * p - 2 * (i128) p * q + 2 * (i128) q * q;
            const i128 xd = (i128) p * p - 2 * (i128) q * q;
            if ((xn == 0) || (xd == 0) || (xn == xd) || (xn == -xd)) {
                continue;
            }
            i128 dn = (2 * xn * xn - xd * xd) * xd * xd;
            i128 dd = xn * xn * xn * xn;
            if ((dn <= 0) || (dn > (i128) UINT64_MAX) || (dd > (i128) UINT64_MAX)) {
                continue;
            }
            const uint64_t g = gcd((uint64_t) dn, (uint64_t) dd);
            dn /= g;
            dd /= g;
            if (dn == dd) {
                continue;
            }
            int npoints = 0;
            for (int b = 1; (b <= A) && (npoints < 3); b++) {
                for (int a = 1; (a <= A) && (npoints < 3); a++) {
                    if ((a == b) || (gcd((uint64_t) a, (uint64_t) b) != 1)) {
                        continue;
                    }
                    //
                    // y^2 = N / M: a square iff N.M is (both of one sign).
                    //
                    const i128 N = ((i128) b * b - (i128) a * a) * dd;
                    const i128 M = dd * b * b - dn * (i128) a * a;
                    if ((N == 0) || (M == 0) || ((N < 0) != (M < 0))) {
                        continue;
                    }
                    const u128 an = (u128) (N < 0 ? -N : N);
                    const u128 am = (u128) (M < 0 ? -M : M);
                    if (am > (u128) -1 / an) {
                        continue;
                    }
                    if (is_square(an * am)) {
                        printf("%d/%d %llu/%llu %d/%d\n", p, q,
                               (unsigned long long) dn,
                               (unsigned long long) dd, a, b);
                        npoints++;
                    }
                }
            }
        }
    }
    return 0;
}
