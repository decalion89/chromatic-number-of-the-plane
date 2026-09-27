"""Rigidity at finite scale: colour a graph by tabu search (several seeds) and measure the best
coset fit R = max_{psi, relabelling} #{x : pi(c(x)) = psi(x)} / n.  R = 1 for a coset colouring,
about 0.2 for an unstructured one.  Units are taken from the module file given (e.g. five_rho7.json)."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, subprocess
from fractions import Fraction as Fr
from math import gcd
import numpy as np
sys.path.insert(0, HN_DIR)
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
d = json.load(open(sys.argv[1])); seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]
g = build_graph(V); EA, EB = zip(*g.edges()); n = len(V)
Ev = sorted({tuple((V[b] - V[a]).x.c) + tuple((V[b] - V[a]).y.c) for a, b in g.edges()})
raw = [tuple((p - V[0]).x.c) + tuple((p - V[0]).y.c) for p in V]
den = 1
for v in Ev + raw:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
B = echelon(E); r = len(B); Cu = [coords(B, e) for e in E]
co = np.array([[int(t) % 5 for t in coords(B, tuple(int(Fr(x) * den) for x in v))] for v in raw], dtype=np.int64)
adm = [psi for psi in itertools.product(range(5), repeat=r) if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in Cu)]
PV = [(co @ np.array(psi, dtype=np.int64)) % 5 for psi in adm]
perms = list(itertools.permutations(range(5)))
print(f"{sys.argv[1]}: n={n} m={len(EA)} rank {r}, {len(adm)} admissible psi", flush=True)
for sd in range(1, seeds + 1):
    inp = f"{n} {len(EA)} 5 50000000 {sd}\n" + "\n".join(f"{a} {b}" for a, b in zip(EA, EB)) + "\n" + "\n".join(["-1"] * n) + "\n"
    out = subprocess.run([__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "tabucol")], input=inp, capture_output=True, text=True).stdout.split("\n")
    if not out[0].startswith("OK"): print(f"  seed {sd}: tabu {out[0]}"); continue
    colv = np.array([int(x) for x in out[1:1 + n]], dtype=np.int64)
    best = 0
    for pv in PV:
        T = np.zeros((5, 5), dtype=np.int64); np.add.at(T, (colv, pv), 1)
        best = max(best, max(sum(T[c][p[c]] for c in range(5)) for p in perms))
    print(f"  seed {sd}: tabu {out[0]}; best coset fit R = {best}/{n} = {best/n:.3f}", flush=True)
