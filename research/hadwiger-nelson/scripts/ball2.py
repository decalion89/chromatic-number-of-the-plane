"""The complete combinatorial 2-ball of five_rho7: O + u + v for all unit vectors u, v.
Colour it by tabu search and measure the best coset fit R.  Does rigidity show up in a
graph that contains EVERY point within two unit steps of a vertex?"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, subprocess, time
from fractions import Fraction as Fr
from math import gcd
import numpy as np
sys.path.insert(0, HN_DIR)
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
t0 = time.time()
d = json.load(open(HN_DIR + "/data/five_rho7.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
U = list({key(P[b] - P[a]): P[b] - P[a] for a, b in g.edges()}.values()) + list({key(P[a] - P[b]): P[a] - P[b] for a, b in g.edges()}.values())
U = list({key(u): u for u in U}.values())
O = P[0] - P[0]
pts = {key(O): O}
for u in U: pts.setdefault(key(u), u)
for u in U:
    for v in U:
        w = u + v; pts.setdefault(key(w), w)
V = list(pts.values()); n = len(V)
vk = {key(p): i for i, p in enumerate(V)}
Uf = [(u.fx, u.fy) for u in U]
EA, EB = [], []
for i, p in enumerate(V):
    for ux, uy in Uf:
        j = vk.get((round(p.fx + ux, 9), round(p.fy + uy, 9)))
        if j is not None and j > i: EA.append(i); EB.append(j)
print(f"2-ball: {n} points, {len(EA)} edges, mean degree {2*len(EA)/n:.1f}   [{time.time()-t0:.0f}s]", flush=True)
Ev = sorted({tuple(u.x.c) + tuple(u.y.c) for u in U})
raw = [tuple(p.x.c) + tuple(p.y.c) for p in V]
den = 1
for v in Ev + raw:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
B = echelon(E); r = len(B); Cu = [coords(B, e) for e in E]
co = np.array([[int(t) % 5 for t in coords(B, tuple(int(Fr(x) * den) for x in v))] for v in raw], dtype=np.int64)
adm = [psi for psi in itertools.product(range(5), repeat=r) if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in Cu)]
PV = [(co @ np.array(psi, dtype=np.int64)) % 5 for psi in adm]
perms = list(itertools.permutations(range(5)))
TABU = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "tabucol")
for sd in (1, 2, 3):
    inp = f"{n} {len(EA)} 5 200000000 {sd}\n" + "\n".join(f"{a} {b}" for a, b in zip(EA, EB)) + "\n" + "\n".join(["-1"] * n) + "\n"
    out = subprocess.run([TABU], input=inp, capture_output=True, text=True).stdout.split("\n")
    if not out[0].startswith("OK"): print(f"  seed {sd}: tabu {out[0]}   [{time.time()-t0:.0f}s]", flush=True); continue
    colv = np.array([int(x) for x in out[1:1 + n]], dtype=np.int64)
    best = 0
    for pv in PV:
        T = np.zeros((5, 5), dtype=np.int64); np.add.at(T, (colv, pv), 1)
        best = max(best, max(sum(T[c][p[c]] for c in range(5)) for p in perms))
    print(f"  seed {sd}: tabu {out[0]}; best coset fit R = {best}/{n} = {best/n:.3f}   [{time.time()-t0:.0f}s]", flush=True)
