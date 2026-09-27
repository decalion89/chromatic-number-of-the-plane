"""The exact condition under which the two blocking tests agree.

M is the group the edge vectors generate, M_sat its saturation in Z^d.  From
0 -> M -> M_sat -> T -> 0 with T = M_sat/M finite, restriction
Hom(Z^d, Z/n) -> Hom(M, Z/n) is onto exactly when Ext^1(T, Z/n) = T/nT
vanishes -- that is, when no invariant factor of M shares a prime with n.
And M has an invariant factor divisible by p exactly when the rank of the
direction matrix drops mod p.  So:

    the Z^d test and the module test agree at n  <=>  rank_p = rank_Q
    for every prime p dividing n.

Two Gaussian eliminations, no SAT.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

t0 = time.time()


def rank_mod(rows, p):
    A = [[x % p for x in r] for r in rows]
    d, rk, row = len(A[0]), 0, 0
    for c in range(d):
        piv = next((i for i in range(row, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        inv = pow(A[row][c], -1, p)
        A[row] = [(x * inv) % p for x in A[row]]
        for i in range(len(A)):
            if i != row and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[row])]
        row += 1
        rk += 1
    return rk


def rank_q(rows):
    A = [[Fr(x) for x in r] for r in rows]
    d, rk, row = len(A[0]), 0, 0
    for c in range(d):
        piv = next((i for i in range(row, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        f = A[row][c]
        A[row] = [x / f for x in A[row]]
        for i in range(len(A)):
            if i != row and A[i][c]:
                g = A[i][c]
                A[i] = [x - g * y for x, y in zip(A[i], A[row])]
        row += 1
        rk += 1
    return rk


def check(rows, name):
    rq = rank_q(rows)
    drops = []
    for p in (2, 3, 5):
        if rank_mod(rows, p) < rq:
            drops.append(p)
    verdict = ("every verdict stands" if not drops
               else f"AT RISK at n divisible by {drops}")
    print(f"{name}: {len(rows)} vectors, dim {len(rows[0])}, rank_Q {rq}; "
          f"rank drops mod {drops or 'nothing'} -- {verdict}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return drops


from hn.homcol import denominator_29_directions, edge_vectors
check(denominator_29_directions(3), "denominator-29 set")
from hn.degrey import build_G
from hn.graph import build_graph
check(edge_vectors(build_graph(build_G(as_graph=False))), "de Grey's G")
