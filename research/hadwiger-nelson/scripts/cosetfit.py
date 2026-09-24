"""How close is a grown graph's colouring to a coset colouring, and where are the defects?
For every admissible psi, the best relabelling pi of the 5 colours maximises #{x : pi(c(x)) = psi(x)};
report the best fit, then where the mismatched vertices sit (distance from the pair)."""
import sys, json, itertools, math
from fractions import Fraction as Fr
from math import gcd
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
d = json.load(open(sys.argv[1])); col = d.get("colouring")
assert col, "checkpoint has no colouring"
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]; A, Bi = d["A"], d["B"]
n = len(V)
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
raw = [tuple((p - V[A]).x.c) + tuple((p - V[A]).y.c) for p in V]
den = 1
for v in Ev + raw:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
B = echelon(E); r = len(B)
Cu = [coords(B, e) for e in E]
co = np.array([[int(t) % 5 for t in coords(B, tuple(int(Fr(x) * den) for x in v))] for v in raw], dtype=np.int64)
adm = [psi for psi in itertools.product(range(5), repeat=r) if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in Cu)]
colv = np.array(col, dtype=np.int64)
best = (-1, None, None)
perms = list(itertools.permutations(range(5)))
for psi in adm:
    pv = (co @ np.array(psi, dtype=np.int64)) % 5
    T = np.zeros((5, 5), dtype=np.int64)
    np.add.at(T, (colv, pv), 1)
    sc, pi = max((sum(T[c][p[c]] for c in range(5)), p) for p in perms)
    if sc > best[0]: best = (sc, psi, pi)
sc, psi, pi = best
pv = (co @ np.array(psi, dtype=np.int64)) % 5
match = np.array([pi[colv[i]] == pv[i] for i in range(n)])
print(f"{sys.argv[1]}: n={n}, rank {r}, {len(adm)} admissible psi; best coset fit {sc}/{n} = {sc/n:.3f} (psi {psi})")
fx = np.array([p.fx for p in V]); fy = np.array([p.fy for p in V])
dist = np.minimum(np.hypot(fx - fx[A], fy - fy[A]), np.hypot(fx - fx[Bi], fy - fy[Bi]))
for lo, hi in ((0, 0.5), (0.5, 1.01), (1.01, 1.5), (1.5, 2.01), (2.01, 3), (3, 4), (4, 6), (6, 99)):
    sel = (dist >= lo) & (dist < hi)
    if sel.sum(): print(f"  distance {lo:>4}-{hi:<4}: {sel.sum():5d} vertices, coset-matching {match[sel].mean():.3f}")
print(f"  c(a) = c(b): {col[A] == col[Bi]}")
