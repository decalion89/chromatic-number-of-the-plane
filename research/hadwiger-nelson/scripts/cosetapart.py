"""Which pairs can EVER be forced APART: those no coset colouring keeps alike.

(u, v) survives iff psi(u - v) != 0 for every admissible psi -- the difference
behaves like an edge vector for every coset colouring.  Tabulate, by distance
class, how many pairs of the graph survive.
"""
import sys, json, time, itertools
from fractions import Fraction as Fr
from math import gcd
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
NAME = sys.argv[1]; t0 = time.time()
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n
raw = [tuple(p.x.c) + tuple(p.y.c) for p in P]
den = 1
for v in raw:
    for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
Z = [tuple(int(Fr(q) * den) for q in v) for v in raw]
Ev = sorted({tuple(a - b for a, b in zip(Z[j], Z[i])) for i, j in G.edges()})
half = list({max(v, tuple(-x for x in v)) for v in Ev})
B = echelon(half); r = len(B)
C = np.array([coords(B, v) for v in half]) % 5
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
ok = np.ones(len(PS), bool)
for c in C: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
hx = np.array([p.fx for p in P]); hy = np.array([p.fy for p in P])
I, J = np.triu_indices(n, 1) if n <= 1000 else (None, None)
D2 = (hx[I] - hx[J]) ** 2 + (hy[I] - hy[J]) ** 2
keep = (D2 < 9.0) & (np.abs(D2 - 1) > 1e-9)
tot = Counter(); surv = Counter()
for i, j in zip(I[keep], J[keep]):
    w = tuple(a - b for a, b in zip(Z[j], Z[i]))
    try: wc = np.array(coords(B, w)) % 5
    except AssertionError: continue
    cl = round(float((hx[i]-hx[j])**2 + (hy[i]-hy[j])**2), 6)
    tot[cl] += 1
    if ((ADM @ wc) % 5 != 0).all(): surv[cl] += 1
print(f"{NAME}: {len(ADM)} admissible psi; distance classes with pairs no coset colouring keeps alike:")
for cl, m in sorted(surv.items(), key=lambda kv: -kv[1])[:25]:
    print(f"   d^2={cl:.6f}  {m} of {tot[cl]} pairs")
print(f"  total: {sum(surv.values())} of {sum(tot.values())} pairs   [{time.time()-t0:.0f}s]")
