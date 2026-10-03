# Lemma W on Z: G(Z, D) 3-colourable  <=>  exists alpha with ||alpha d|| >= 1/3 for all d in D  (kappa(D) >= 1/3).
# kappa test exact (endpoints of the strips); 3-colourability: SAT on the segment [0, N] (UNSAT proves chi >= 4).
import itertools
from fractions import Fraction as Fr
from pysat.solvers import Solver
def kappa_ok(D):
    cands = set()
    for d in D:
        for k in range(d):
            cands.add(Fr(3 * k + 1, 3 * d)); cands.add(Fr(3 * k + 2, 3 * d))
    for a in cands:
        if all(Fr(1, 3) <= (a * d) % 1 <= Fr(2, 3) for d in D):
            return True
    return False
def seg3(D, N):
    s = Solver(name="cadical195"); v = lambda x, c: 3 * x + c + 1
    for x in range(N):
        s.add_clause([v(x, c) for c in range(3)])
        for d in D:
            if x + d < N:
                for c in range(3): s.add_clause([-v(x, c), -v(x + d, c)])
    return s.solve()
bad = 0; tot = 0; four = []
for D in itertools.combinations(range(1, 25), 3):
    from math import gcd
    if gcd(gcd(D[0], D[1]), D[2]) != 1: continue
    tot += 1
    k = kappa_ok(D); s = seg3(D, 12 * max(D) + 60)
    if not k: four.append(D)
    if k != s:
        bad += 1; print("MISMATCH", D, "kappa>=1/3:", k, "segment 3-colourable:", s)
print("triples", tot, "mismatches", bad, "4-chromatic:", len(four))
# Zhu's classification: chi = 4 iff D = {1, 2, 3n} or D = {x, y, x + y} with x != y mod 3 (gcd 1)
def zhu4(D):
    a, b, c = D
    if a == 1 and b == 2 and c % 3 == 0: return True
    if a + b == c and (a - b) % 3 != 0: return True
    return False
print("disagreements with Zhu's list:", [D for D in four if not zhu4(D)][:10], [D for D in itertools.combinations(range(1,25),3) if zhu4(D) and D not in four and gcd(gcd(*D[:2]),D[2])==1][:10])
