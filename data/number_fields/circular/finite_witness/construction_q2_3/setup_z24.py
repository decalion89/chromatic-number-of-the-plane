"""setup_z24.py SPEC OUT.json: the unit vectors mu_24 * {g^l} of Q(sqrt2, sqrt3)^2 (as elements of Q(zeta_24), integer
coordinates over the power basis 1, zeta, ..., zeta^7 times a common denominator D), one per +- pair, and an
LLL-reduced Z-basis of their relation lattice (PARI matkerint + qflll).  SPEC as in kappa_z24.py (e.g. s2:2)."""
import sys, json, subprocess, os
sys.argv = [sys.argv[0]] + sys.argv[1].split() + [sys.argv[2]]
OUT = sys.argv.pop()
import importlib.util
from fractions import Fraction as Fr
from math import lcm
# reuse the field code of kappa_z24.py without running its MILP
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kappa_z24.py")).read()
src = src.split("den = lcm(")[0]
g = {}
exec(compile(src, "kappa_z24_head", "exec"), g)
U = g["U"]
D = lcm(*[x.denominator for u in U for x in u])
Ui = [[int(x * D) for x in u] for u in U]
n = len(Ui)
gp = "M = [" + ";".join(",".join(str(Ui[j][k]) for j in range(n)) for k in range(8)) + "];\n"
gp += "K = matkerint(M); K = K * qflll(K); f = \"kz.txt\"; for (j = 1, #K, write1(f, Vec(K[,j]))); quit;\n"
open("kz.gp", "w").write(gp)
if os.path.exists("kz.txt"):
    os.remove(os.path.abspath("kz.txt"))
subprocess.run(["gp", "-q", "kz.gp"], check=True)
rels = [json.loads(l) for l in open("kz.txt").read().replace("][", "]\n[").split("\n") if l.strip()]
for r in rels:
    assert all(sum(r[j] * Ui[j][k] for j in range(n)) == 0 for k in range(8))
json.dump({"field": "Q(zeta24) power basis", "D": D, "units": Ui, "relations": rels}, open(OUT, "w"))
print(f"{n} vectors, D = {D}, {len(rels)} basis relations, l1 norms {sorted(sum(abs(x) for x in r) for r in rels)}")
