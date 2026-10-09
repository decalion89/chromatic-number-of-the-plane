"""Referee E: torsion characters with kappa > 2/7 for the one-prime group <i, rho> and the two-prime group
<i, rho, sigma>, exact integer arithmetic.

A character of order M (prime to the primes inverted) is c = w/M, w in Z[i]/M; its values on the group H are
Re(conj(w) gamma)/M mod 1, gamma in H (reduced mod M).  kappa(w) = min over the H-orbit of conj(w) of
min(||Re z/M||, ||Im z/M||) (i is in H).  We compute the largest H-invariant subset of
{z : ||Re z/M||, ||Im z/M|| > 2/7} by pruning, then the exact kappa of each surviving orbit of exact order M.

Usage: python3 ref_torsion.py one|two MMAX
"""
import sys
from math import gcd
from fractions import Fraction as Fr
import numpy as np


def inv(a, M):
    return pow(a, -1, M)


def run(which, MMAX):
    res = []
    for M in range(2, MMAX + 1):
        if M % 5 == 0:
            continue
        if which == "two" and M % 13 == 0:
            continue
        gens = []
        i5 = inv(5, M)
        gens.append(((3 * i5) % M, (4 * i5) % M))   # rho = (3+4i)/5
        if which == "two":
            i13 = inv(13, M)
            gens.append(((5 * i13) % M, (12 * i13) % M))  # sigma = (5+12i)/13
        gens.append((0, 1))  # i
        u = np.arange(M, dtype=np.int64)[:, None] * np.ones((1, M), dtype=np.int64)
        v = np.ones((M, 1), dtype=np.int64) * np.arange(M, dtype=np.int64)[None, :]
        du = np.minimum(u, M - u)
        dv = np.minimum(v, M - v)
        alive = (7 * du > 2 * M) & (7 * dv > 2 * M)
        # images under multiplication by each generator (a+bi)(u+vi) = (au - bv) + (av + bu) i
        imgs = []
        for (a, b) in gens:
            imgs.append((np.mod(a * u - b * v, M), np.mod(a * v + b * u, M)))
        while True:
            new = alive.copy()
            for (iu, iv) in imgs:
                new &= alive[iu, iv]
            if (new == alive).all():
                break
            alive = new
        if not alive.any():
            continue
        # orbits of exact order M
        seen = set()
        pts = list(zip(*np.nonzero(alive)))
        for (a0, b0) in pts:
            a0, b0 = int(a0), int(b0)
            if (a0, b0) in seen:
                continue
            orb = {(a0, b0)}
            stack = [(a0, b0)]
            while stack:
                (x, y) = stack.pop()
                for (a, b) in gens:
                    z = ((a * x - b * y) % M, (a * y + b * x) % M)
                    if z not in orb:
                        orb.add(z)
                        stack.append(z)
            seen |= orb
            if gcd(gcd(a0, b0), M) != 1:
                continue  # lower order: counted at its own order
            kap = min(Fr(min(x, M - x, y, M - y), M) for (x, y) in orb)
            res.append((M, kap, len(orb), (a0, b0)))
    return res


def main():
    which, MMAX = sys.argv[1], int(sys.argv[2])
    res = run(which, MMAX)
    orders = sorted(set(r[0] for r in res))
    print("group %s, orders up to %d: exact orders with an orbit of kappa > 2/7: %s" % (which, MMAX, orders))
    best = {}
    for M, kap, size, rep in res:
        if M not in best or kap > best[M][0]:
            best[M] = (kap, size, rep)
    for M in orders:
        kap, size, rep = best[M]
        print("  M=%d: max kappa %s (%.6f), orbit size %d, representative w=%d+%di" % (M, kap, float(kap), size, rep[0], rep[1]))
    if res:
        mk = max(r[1] for r in res)
        print("largest kappa over all these orbits: %s" % mk)


if __name__ == "__main__":
    main()
