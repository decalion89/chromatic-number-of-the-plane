"""Build, verify (3 pysat solvers + kissat/DRAT) and shrink a lattice witness: points of Z[w]/sqrt(-3)
(written as Eisenstein integers a + b w, real distance^2 = N/3), edges at the forbidden norms."""
import sys, json, subprocess
from pysat.solvers import Solver
K = 5
norms = [int(x) for x in sys.argv[1].split(",")]; R = int(sys.argv[2]); out = sys.argv[3]
def norm(a, b): return a * a - a * b + b * b
pts = [(a, b) for a in range(-2 * R, 2 * R + 1) for b in range(-2 * R, 2 * R + 1) if norm(a, b) <= R * R]
def edges_of(P):
    E = []
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            if norm(P[i][0] - P[j][0], P[i][1] - P[j][1]) in norms: E.append((i, j))
    return E
def sat(P, E, name="cd19"):
    X = lambda v, c: 1 + v * K + c
    s = Solver(name=name)
    for v in range(len(P)): s.add_clause([X(v, c) for c in range(K)])
    for i, j in E:
        for c in range(K): s.add_clause([-X(i, c), -X(j, c)])
    r = s.solve(); s.delete(); return r
E = edges_of(pts)
print(f"norms {norms} (d^2 = {[f'{n}/3' for n in norms]}), R={R}: {len(pts)} points, {len(E)} edges")
print("  solvers:", {nm: sat(pts, E, nm) for nm in ("cd19", "g4", "m22")})
# shrink: drop vertices while still not 5-colourable
P = list(pts)
changed = True
while changed:
    changed = False
    for v in list(P):
        Q = [p for p in P if p != v]
        if not sat(Q, edges_of(Q)):
            P = Q; changed = True
E = edges_of(P)
print(f"  vertex-critical core: {len(P)} points, {len(E)} edges; solvers:", {nm: sat(P, E, nm) for nm in ("cd19", "g4", "m22")})
X = lambda v, c: 1 + v * K + c
cls = [[X(v, c) for c in range(K)] for v in range(len(P))] + [[-X(i, c), -X(j, c)] for i, j in E for c in range(K)]
cnf = out.replace(".json", ".cnf")
open(cnf, "w").write(f"p cnf {len(P) * K} {len(cls)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cls))
import os
KS = os.environ.get("KISSAT", "kissat")
DT = os.environ.get("DRAT_TRIM", "drat-trim")
ko = subprocess.run([KS, "-n", "--no-binary", cnf, cnf + ".drat"], capture_output=True, text=True).stdout
print("  kissat:", next(l for l in ko.split("\n") if l.startswith("s ")))
dt = subprocess.run([DT, cnf, cnf + ".drat"], capture_output=True, text=True).stdout
print("  drat-trim:", next((l for l in dt.split("\n") if "VERIFIED" in l), "NOT VERIFIED"))
json.dump({"lattice": "Z[w]/sqrt(-3): point (a, b) is (a + b w)/sqrt(-3), squared distance N(a + b w)/3",
           "forbidden_norms": norms, "forbidden_d2": [f"{n}/3" for n in norms], "points": P, "edges": E}, open(out, "w"))
