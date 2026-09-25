"""Is a graph, with its unit edges plus the edges at one or more Galois orbits of distances, 5-colourable?

A Galois automorphism of the CM field maps a unit-distance gadget at d^2 to one at sigma(d^2)
(notes/worker_jobs.md), so a witness that uses a whole orbit needs a single gadget. Orbits (d^2):
  14_2: (14 -+ 2 sqrt33)/3   d = 0.9149, 2.9149   (P(same) ~ 0.08 in the L16 seed)
  14_5: 14/3 -+ 5 sqrt33/9   d = 1.2146, 2.8032
  9_1:  (9 -+ sqrt33)/6      d = 0.7366, 1.5676   (two-step, P(same) ~ 0.21)
  43:   4/3                  d = 2/sqrt3          (P(same) ~ 0.08)
  none: unit edges only
Join several with '+', e.g. 14_2+43. Writes the CNF next to the scratch directory and runs tabu, then kissat.
usage: orbit_witness_test.py graph.json orbits ktime [kissat] [workdir]"""
import sys, json, time, subprocess
from fractions import Fraction as Fr
import numpy as np
from scipy.spatial import cKDTree
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from hn.field import Field
from hn.geometry import Point
path, oname, KT = sys.argv[1], sys.argv[2], int(sys.argv[3])
KISSAT = sys.argv[4] if len(sys.argv) > 4 else "kissat"
WORK = sys.argv[5] if len(sys.argv) > 5 else "."
d = json.load(open(path)); F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
n = len(P); xy = np.array([[p.fx, p.fy] for p in P]); tree = cKDTree(xy)
s33 = F.sqrt(33); q = lambda a, b=1: F.rational(Fr(a, b))
ORB = {"14_2": [q(14, 3) - s33 * q(2, 3), q(14, 3) + s33 * q(2, 3)],
       "14_5": [q(14, 3) - s33 * q(5, 9), q(14, 3) + s33 * q(5, 9)],
       "9_1": [q(3, 2) - s33 * q(1, 6), q(3, 2) + s33 * q(1, 6)],
       "43": [q(4, 3)], "none": []}
def pairs_at(t):
    r = float(t) ** .5
    pr = tree.query_pairs(r + 1e-6, output_type="ndarray")
    dd = ((xy[pr[:, 0]] - xy[pr[:, 1]]) ** 2).sum(1); pr = pr[np.abs(dd - r * r) < 1e-6]
    return [(int(i), int(j)) for i, j in pr if P[i].dist2(P[j]) == t]
E = pairs_at(q(1)); print(f"{path}: {n} points, {len(E)} unit edges", flush=True)
for part in oname.split("+"):
    for t in ORB[part]:
        ex = pairs_at(t); print(f"  + {len(ex)} edges at d^2 = {t} (d = {float(t)**.5:.4f})", flush=True); E += ex
t0 = time.time()
inp = f"{n} {len(E)} 5 30000000 7\n" + "".join(f"{i} {j}\n" for i, j in E) + " ".join(["-1"] * n) + "\n"
out = subprocess.run([os.path.join(os.path.dirname(os.path.abspath(__file__)), "tabucol")], input=inp, capture_output=True, text=True).stdout.split()
print(f"  tabu: {out[0]} {out[1]}  [{time.time()-t0:.0f}s]", flush=True)
if out[0] != "OK":
    cnf = os.path.join(WORK, os.path.basename(path).replace(".json", "") + f"_{oname}.cnf"); X = lambda v, c: 1 + v * 5 + c
    with open(cnf, "w") as fh:
        cl = [[X(v, c) for c in range(5)] for v in range(n)] + [[-X(i, c), -X(j, c)] for i, j in E for c in range(5)]
        fh.write(f"p cnf {n*5} {len(cl)}\n" + "".join(" ".join(map(str, c_)) + " 0\n" for c_ in cl))
    ko = subprocess.run([KISSAT, f"--time={KT}", "-n", cnf], capture_output=True, text=True).stdout
    print("  kissat:", next((l for l in ko.split("\n") if l.startswith("s ")), "s UNKNOWN"), f"[{time.time()-t0:.0f}s]", flush=True)
