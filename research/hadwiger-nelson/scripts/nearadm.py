"""Blocked module: coset colourings fail only on the few unit vectors lying in 5M.  Look for
'nearly admissible' psi (zero exactly on those) and repair them with a small second coordinate
chi : M -> Z/m that is nonzero there: any proper 5-colouring of Cay(Z/5 x Z/m, (psi, chi)(U))
would be a periodic colouring of the whole blocked module -- and would end the plain search."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, time
from fractions import Fraction as Fr
import numpy as np
from pysat.solvers import Solver
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
ROOT = HN_DIR
name = sys.argv[1]; MS = [int(x) for x in sys.argv[2].split(",")]
d = json.load(open(f"{ROOT}/data/{name}")); F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P)); B = echelon(E); r = len(B)
C = [coords(B, v) for v in E]
from math import gcd
from functools import reduce
content = [abs(reduce(gcd, c)) for c in C]
special = [k for k, cc in enumerate(content) if cc % 5 == 0]
print(f"{name}: {len(C)} directions, rank {r}; directions in 5M: {len(special)}", flush=True)
C5 = np.array([[x % 5 for x in c] for c in C], dtype=np.int64)
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
V = (PS @ C5.T) % 5
zero = (V == 0)
ns = np.array([k for k in range(len(C)) if k not in special])
near = np.nonzero(~np.any(zero[:, ns], axis=1))[0]
print(f"  psi nonzero on every direction outside 5M: {len(near)}", flush=True)
t0 = time.time()
for m in MS:
    tried = colourable = 0
    Cm = np.array([[x % m for x in c] for c in C], dtype=np.int64)
    CH = np.array(list(itertools.product(range(m), repeat=r)), dtype=np.int64) if m ** r <= 2_000_000 else None
    W = (CH @ Cm.T) % m
    goodchi = np.nonzero(np.all(W[:, special] != 0, axis=1))[0]
    print(f"  m={m}: chi nonzero on the {len(special)} directions in 5M: {len(goodchi)} of {len(CH)}", flush=True)
    if len(near) == 0 or len(goodchi) == 0: continue
    rng = np.random.default_rng(3)
    for trial in range(3000):
        i = int(rng.choice(near)); j = int(rng.choice(goodchi))
        S = set()
        for k in range(len(C)):
            a, b = int(V[i][k]), int(W[j][k])
            if (a, b) == (0, 0): break
            S.add((a, b)); S.add(((-a) % 5, (-b) % m))
        else:
            tried += 1
            elts = [(a, b) for a in range(5) for b in range(m)]; idx = {e: t for t, e in enumerate(elts)}
            X = lambda v, c: 1 + v * 5 + c
            cl = [[X(v, c) for c in range(5)] for v in range(len(elts))]
            for v, (a, b) in enumerate(elts):
                for (s1, s2) in S:
                    w = idx[((a + s1) % 5, (b + s2) % m)]
                    if v < w:
                        for c in range(5): cl.append([-X(v, c), -X(w, c)])
            cl.append([X(0, 0)])
            s = Solver(name="cd19", bootstrap_with=cl)
            if s.solve():
                colourable += 1
                if colourable <= 3: print(f"    COLOURABLE: psi #{i}, chi #{j}, |S| = {len(S)}", flush=True)
            s.delete()
    print(f"  m={m}: {tried} loopless (psi, chi) tried, {colourable} 5-colourable   [{time.time()-t0:.0f}s]", flush=True)
