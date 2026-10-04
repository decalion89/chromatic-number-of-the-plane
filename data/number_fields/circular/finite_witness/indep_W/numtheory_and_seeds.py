#!/usr/bin/env python3
"""Arithmetic facts used in Section 10 / README, and the seed graphs (referee code).
Usage: python3 numtheory_and_seeds.py REPO_ROOT"""
import json, gzip, os, sys, itertools
from pysat.solvers import Solver

REPO = sys.argv[1]
FW = os.path.join(REPO, 'data/number_fields/circular/finite_witness')


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_in_Q7(d):
    # d squarefree integer: sqrt d in Q_7 iff 7 does not divide d and d is a nonzero square mod 7 (Hensel, p odd)
    return d % 7 != 0 and legendre(d, 7) == 1


for d in (11, 191, 455):
    a_holds = d % 4 != 3          # (a) for real quadratic fields, as stated in the paper's introduction
    b_holds = d % 3 != 2          # (b)
    l7 = legendre(d, 7)
    beh = 'ramifies' if d % 7 == 0 else ('splits' if l7 == 1 else 'inert')
    # residue degree 1 above 7  <=>  split or ramified (disc of Q(sqrt d) is 4d for d = 3 mod 4)
    print(f"d={d}: d mod 12 = {d % 12}; (a) holds: {a_holds}; (b) holds: {b_holds}; 7 {beh} in Q(sqrt{d}) "
          f"(Legendre ({d}/7) = {l7}); condition (7): {beh != 'inert'}; sqrt{d} in Q_7: {sqrt_in_Q7(d)}; "
          f"Prop. C1 (residue field F_7) applies: {beh != 'inert'}")
# Hensel lift check for sqrt 11 and sqrt 191 modulo 7^10 (explicit)
for d in (11, 191):
    r = next(t for t in range(7) if (t * t - d) % 7 == 0)
    mod = 7
    for _ in range(9):
        mod *= 7
        # Newton step
        r = (r - (r * r - d) * pow(2 * r, -1, mod)) % mod
    print(f"  sqrt{d} mod 7^10 = {r}; check (r^2 - d) mod 7^10 = {(r * r - d) % mod}")
print("  455 = 5*7*13: v_7(455) = 1 is odd, so 455 is not a square in Q_7; Q_7(sqrt455) is ramified of degree 2 "
      "(residue degree 1)")

# seed graphs: sizes, denominators, 3-colourability, how many seed points survive in the witness
for d, W in ((11, 'witness_q11'), (191, 'witness_q191'), (455, 'witness_q455')):
    q = json.load(open(os.path.join(REPO, 'data/quadratic_planes', f'q{d}.json')))
    P = [tuple(p) for p in q['points']]
    D = q['D']
    n = len(P)
    E = [(i, j) for i in range(n) for j in range(i + 1, n)
         if (lambda A, B, C, Ee: A * A + d * B * B + C * C + d * Ee * Ee == D * D and A * B + C * Ee == 0)
         (*(P[j][k] - P[i][k] for k in range(4)))]
    s = Solver(name='glucose4')
    x = lambda v, c: 3 * v + c + 1
    for v in range(n):
        s.add_clause([x(v, c) for c in range(3)])
    for (i, j) in E:
        for c in range(3):
            s.add_clause([-x(i, c), -x(j, c)])
    three = s.solve()
    w = json.load(gzip.open(os.path.join(FW, W + '.json.gz'), 'rt'))
    inter = len(set(P) & set(map(tuple, w['points'])))
    print(f"seed q{d}.json: {n} points, D = {D}, {len(E)} unit pairs (file lists {len(q['edges'])}), "
          f"3-colourable: {three}; seed points in {W}: {inter} of {n}; witness D = {w['denominator']}")
