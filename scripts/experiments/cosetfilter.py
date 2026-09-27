"""Which pairs can EVER be forced together: the coset colourings decide first.

An integral multiquadratic graph always has coset 5-colourings c = psi(p) with
psi : M -> Z/5 nonzero on every edge vector.  Such a colouring gives u and v the
same colour iff psi(u - v) = 0.  So (u, v) can be monochromatic in EVERY
5-colouring only if psi(u - v) = 0 for EVERY admissible psi -- and growth along
the existing directions keeps every admissible psi admissible.  Enumerate all
psi over the module's rank and filter the five-closable pairs.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, itertools
from fractions import Fraction as Fr
from math import gcd
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
ROOT = HN_DIR
NAME = sys.argv[1]
t0 = time.time()
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n
# ambient integer coordinates with one common denominator for points AND edges
raw = [tuple(p.x.c) + tuple(p.y.c) for p in P]
den = 1
for v in raw:
    for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
Z = [tuple(int(Fr(q) * den) for q in v) for v in raw]
Evec = sorted({tuple(a - b for a, b in zip(Z[j], Z[i])) for i, j in G.edges()})
half = list({max(v, tuple(-x for x in v)) for v in Evec})
B = echelon(half)
C = np.array([coords(B, v) for v in half]) % 5
r = len(B)
print(f"{NAME}: n={n}, {len(half)} directions, module rank {r}; enumerating 5^{r} functionals   [{time.time()-t0:.0f}s]", flush=True)
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
ok = np.ones(len(PS), bool)
for c in C:
    ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
print(f"  admissible coset colourings: {len(ADM)} of {len(PS)}   [{time.time()-t0:.0f}s]", flush=True)
def rank5(M):
    M = [[int(x) for x in r] for r in M % 5]; rk = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % 5), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], 3, 5)
        M[rk] = [(x * inv) % 5 for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % 5 for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
print(f"  admissible functionals span a space of dimension {rank5(ADM[:400])} of {r}: common kernel is "
      f"{'exactly 5M' if rank5(ADM[:400]) == r else 'larger than 5M'}", flush=True)
from collections import defaultdict
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 5.5), int(hy[i] // 5.5))].append(i) if False else None
hx = np.array([p.fx for p in P]); hy = np.array([p.fy for p in P])
targets = {Fr(5, 9): "cos 1/10", Fr(25): "cos 49/50", Fr(4, 3): "(unclosable)", Fr(3): "Moser 5/6 (integral)"}
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 5.5), int(hy[i] // 5.5))].append(i)
near = []
for i in range(n):
    cx, cy = int(hx[i] // 5.5), int(hy[i] // 5.5)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            near.extend((i, j) for j in cells.get((cx + dx, cy + dy), ()) if j > i)
NI = np.array([a for a, b in near]); NJ = np.array([b for a, b in near])
ND = (hx[NI] - hx[NJ])**2 + (hy[NI] - hy[NJ])**2
keep5 = []
for T, why in targets.items():
    msk = np.abs(ND - float(T)) < 1e-9
    sel = list(zip(NI[msk].tolist(), NJ[msk].tolist()))
    surv = 0; tot = 0
    for i, j in sel:
        w = tuple(a - b for a, b in zip(Z[j], Z[i]))
        try: wc = np.array(coords(B, w)) % 5
        except AssertionError: continue       # different components
        tot += 1
        if not ((ADM @ wc) % 5).any():
            surv += 1; keep5.append((int(i), int(j), str(T)))
    print(f"  d^2={T} [{why}]: {tot} pairs; {surv} with psi(u-v)=0 for EVERY admissible psi   [{time.time()-t0:.0f}s]", flush=True)

json.dump({"graph": NAME, "unsplittable_by_cosets": keep5}, open(f"{ROOT}/data/coset_unsplit_{NAME}", "w"))
print(f"  written data/coset_unsplit_{NAME}: {len(keep5)} pairs", flush=True)
