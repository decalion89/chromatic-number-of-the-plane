"""Independent check of a claimed obstruction, before a word is said about it.

Three kinds of claim, all rebuilt from the saved exact coordinates:
  udg    the unit-distance graph on the points is not 5-colourable        (chi(R^2) >= 6)
  apart  the unit-distance graph plus c(A) = c(B) is not 5-colourable     (A, B forced APART)
  same   the unit-distance graph plus c(A) != c(B) is not 5-colourable    (A, B forced SAME: a spindle
         rotation about A that moves B by 1 then gives chi(R^2) >= 6 directly, when |AB| >= 1/2)
  two    the {1, d}-graph (edges at distance 1 and d) is not 5-colourable (chi(R^2, {1, d}) >= 6)
  two-apart, two-same   as apart / same, but on the {1, d}-graph (a gadget or forced pair for the
         cascade: with a unit-distance gadget for d they become unit-distance statements)
The kinds name the conclusion, as grow_lean.py's MODE does: MODE=apart growth is checked with kind apart,
MODE=same growth with kind same.

Every edge is re-derived from the exact field arithmetic, never from floats: a float KD-tree only
proposes pairs, and each is kept only if |p - q|^2 equals 1 (or d^2) exactly. No colour is pinned, so a
wrong pin cannot manufacture UNSAT. Three unrelated pysat solvers must agree, and if kissat and
drat-trim are given, kissat's DRAT proof is checked too.

usage: verify_pair.py <graph.json> udg|apart|same|two [--d2 a/b] [--kissat PATH --drat-trim PATH]
"""
import sys, json, time, argparse, subprocess, os
from fractions import Fraction as Fr
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from hn.field import Field
from hn.geometry import Point
from pysat.solvers import Solver

ap = argparse.ArgumentParser()
ap.add_argument("graph"); ap.add_argument("kind", choices=["udg", "apart", "same", "two", "two-apart", "two-same"])
ap.add_argument("--d2", default=None); ap.add_argument("--kissat", default=None); ap.add_argument("--drat-trim", default=None)
a = ap.parse_args()
t0 = time.time(); K = 5
d = json.load(open(a.graph)); F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(p, q) for p, q in x]), F.element([Fr(p, q) for p, q in y])) for x, y in d["points"]]
n = len(V)
keys = set((tuple(p.x.c), tuple(p.y.c)) for p in V)
print(f"{a.graph}: {n} points, {len(keys)} distinct   [{time.time()-t0:.0f}s]", flush=True)
assert len(keys) == n, "repeated points"
pts = np.array([[p.fx, p.fy] for p in V]); tree = cKDTree(pts)


def exact_pairs(dist2):
    r = float(dist2) ** 0.5
    pr = tree.query_pairs(r + 1e-6, output_type="ndarray")
    dd = ((pts[pr[:, 0]] - pts[pr[:, 1]]) ** 2).sum(1)
    pr = pr[np.abs(dd - float(dist2)) < 1e-6]
    target = F.element([Fr(dist2)] + [Fr(0)] * (F.dim - 1))
    return [(int(i), int(j)) for i, j in pr if V[i].dist2(V[j]) == target]


E = exact_pairs(Fr(1))
print(f"  exact unit edges: {len(E)}", flush=True)
if a.kind.startswith("two"):
    D2 = Fr(a.d2 or d.get("dist2"))
    E2 = exact_pairs(D2)
    print(f"  exact edges at d^2 = {D2}: {len(E2)}", flush=True)
    E = E + E2
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for i, j in E:
    for c in range(K): cnf.append([-X(i, c), -X(j, c)])
if a.kind in ("apart", "two-apart"):
    A, B = d["A"], d["B"]
    print(f"  pair A={A}, B={B}: d^2 = {V[A].dist2(V[B])}; imposing c(A) = c(B) (UNSAT = forced apart)", flush=True)
    for c in range(K): cnf += [[-X(A, c), X(B, c)], [X(A, c), -X(B, c)]]
if a.kind in ("same", "two-same"):
    A, B = d["A"], d["B"]
    print(f"  pair A={A}, B={B}: d^2 = {V[A].dist2(V[B])}; imposing c(A) != c(B) (UNSAT = forced same)", flush=True)
    for c in range(K): cnf.append([-X(A, c), -X(B, c)])
verdicts = {}
for name in ("cd19", "g4", "m22"):
    s = Solver(name=name, bootstrap_with=cnf); r = s.solve(); s.delete()
    verdicts[name] = r
    print(f"  {name}: {'SAT (colourable)' if r else 'UNSAT'}   [{time.time()-t0:.0f}s]", flush=True)
if a.kissat:
    path = a.graph.replace(".json", f".verify_{a.kind}.cnf"); proof = path + ".drat"
    nv = n * K
    open(path, "w").write(f"p cnf {nv} {len(cnf)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cnf))
    out = subprocess.run([a.kissat, "-n", "--no-binary", path, proof], capture_output=True, text=True).stdout
    st = next((l for l in out.split("\n") if l.startswith("s ")), "s UNKNOWN")
    print(f"  kissat: {st}   [{time.time()-t0:.0f}s]", flush=True)
    if "UNSAT" in st and a.drat_trim:
        dt = subprocess.run([a.drat_trim, path, proof], capture_output=True, text=True).stdout
        print("  drat-trim:", next((l for l in dt.split("\n") if "VERIFIED" in l or "NOT" in l), dt[-300:]), flush=True)
if len(set(verdicts.values())) > 1:
    print("DISAGREEMENT between solvers -- not a result")
elif not any(verdicts.values()):
    print(f"ALL SOLVERS: UNSAT ({a.kind}) -- still to be checked by hand before any claim")
else:
    print(f"colourable ({a.kind}): no obstruction")
