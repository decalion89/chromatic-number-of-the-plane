"""Build the twisted coset colouring EXPLICITLY on a saved graph and check it.

Finds psi (admissible mod 5), t = -2 psi(e), and an exact rational functional
w on the ambient coordinates with <w,u> >= 0 on D_t and <w,e> > 0; then colours
every vertex by c(p) = psi(p) + t * floor((phi(p) - phi0) / L) with the level
boundary placed between a and b, and checks every edge exactly and c(a) = c(b).
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, time
from fractions import Fraction as Fr
from math import gcd, floor
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, HN_DIR)
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
path = sys.argv[1]
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
A, Bv = d["A"], d["B"]
g = build_graph(P); n = g.n; Eg = list(g.edges())
raw = [tuple(p.x.c) + tuple(p.y.c) for p in P]
den = 1
for v in raw:
    for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
Z = [tuple(int(Fr(q) * den) for q in v) for v in raw]
diff = lambda i, j: tuple(x - y for x, y in zip(Z[j], Z[i]))
dirs = sorted({max(v, tuple(-x for x in v)) for v in (diff(i, j) for i, j in Eg)})
Bm = echelon(dirs); r = len(Bm)
Cm = [coords(Bm, v) for v in dirs]
Ci = np.array([[x % 5 for x in c] for c in Cm], dtype=np.int64)
ev = diff(A, Bv); e2 = ev                       # b - a = 2e
e = tuple(x // 2 for x in e2); assert all(x % 2 == 0 for x in e2)
ke = dirs.index(max(e, tuple(-x for x in e))); sgn = 1 if dirs[ke] == e else -1
print(f"{path}: n={n}, {len(dirs)} directions, rank {r}; pair direction index {ke}", flush=True)
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64) if 5 ** r <= 400000 else np.random.default_rng(0).integers(0, 5, size=(400000, r))
ok = np.ones(len(PS), bool)
for c in Ci: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
Uv = dirs + [tuple(-x for x in v) for v in dirs]
Ui = np.vstack([Ci, (-Ci) % 5]); VAL = (ADM @ Ui.T) % 5
amb = np.array([[float(x) for x in v] for v in Uv]); scale = np.abs(amb).max(axis=1, keepdims=True); ambn = amb / scale
ef = np.array([float(x) for x in e]); efn = ef / np.abs(ef).max()
found = None
for i in range(len(ADM)):
    s = int(VAL[i][ke]) if sgn == 1 else (-int(VAL[i][ke])) % 5      # psi(e)
    t = (-2 * s) % 5
    D = ambn[VAL[i] == t]
    res = linprog(np.zeros(len(ef)), A_ub=np.vstack([-D, -efn[None, :]]), b_ub=np.concatenate([np.zeros(len(D)), [-1.0]]),
                  bounds=[(-1000, 1000)] * len(ef), method="highs")
    if res.status != 0: continue
    w = [Fr(x).limit_denominator(1000) for x in res.x]
    # exact re-check on the integer vectors
    Dint = [Uv[j] for j in range(len(Uv)) if VAL[i][j] == t]
    if all(sum(a * b for a, b in zip(w, u)) >= 0 for u in Dint) and sum(a * b for a, b in zip(w, e)) > 0:
        found = (ADM[i], t, w); break
if not found:
    print("  no exact twisted colouring found"); sys.exit()
psi, t, w = found
phi = lambda z: sum(a * b for a, b in zip(w, z))
L = max(max(abs(phi(u)) for u in Uv), 2 * phi(e))
ps = [sum(int(x) * int(y) for x, y in zip(psi, coords(Bm, diff(A, j)))) % 5 for j in range(n)]
ph = [phi(diff(A, j)) for j in range(n)]
phi0 = (ph[A] + ph[Bv]) / 2
col = [(ps[j] + t * floor((ph[j] - phi0) / L)) % 5 for j in range(n)]
bad = sum(1 for i, j in Eg if col[i] == col[j])
print(f"  psi={tuple(int(x) for x in psi)}, t={t}, L={float(L):.3f}: improper edges {bad} of {len(Eg)}; "
      f"c(a)={col[A]} c(b)={col[Bv]} -> {'KEPT ALIKE: no gadget at this pair on these directions' if col[A] == col[Bv] and bad == 0 else 'check failed'}", flush=True)
