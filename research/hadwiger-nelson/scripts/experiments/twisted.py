"""Twisted coset colourings: a coset colouring plus a colour shift at every
level of an integer linear functional on the edge module.

    c(p) = psi(p) + t * floor(phi(p) / L)      (mod 5)

psi : M -> Z/5 admissible (nonzero on every unit vector), phi : M -> Z linear,
L >= max |phi(u)|.  A unit step moves the level by 0 or +-1, so c is proper on
EVERY graph with these edge directions exactly when the class
D_t = {u : psi(u) = t} lies in the half-space phi >= 0.  It is not periodic and
not a coset colouring.  At the pair (a, a + 2e), with the level boundary placed
between them, the difference is 2 psi(e) + t sign(phi(e)); so the pair is kept
ALIKE when t = -2 psi(e) and phi(e) > 0 -- a refutation of any distance-2 gadget
along e -- and a 5e pair is SPLIT by any t with phi(e) != 0.  By Farkas the
half-space exists iff -e is not in the cone of D_t: one small LP each.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, itertools, time
from fractions import Fraction as Fr
from math import gcd
import numpy as np
from scipy.optimize import linprog
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(os.path.join(HN_DIR, "scripts", "gate.py")).read().split("def gate(g, label):")[0])
ROOT = HN_DIR
src = sys.argv[1]
if src.startswith("units:"):
    E = [tuple(v) for v in json.load(open(src[6:]))["directions"]]
elif src == "ei":
    sols = [(a,b,c,d) for a in range(-7,8) for b in range(-4,5) for c in range(-12,13) for d in range(-3,4)
            if 3*a*a + 11*b*b + c*c + 33*d*d == 144 and a*b == -c*d]
    E = list({max(v, tuple(-x for x in v)) for v in sols})
else:
    d = json.load(open(src if src.startswith("/") else f"{ROOT}/data/{src}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    E = edge_vectors(build_graph(P))
B = echelon(E); r = len(B)
Cm = [coords(B, v) for v in E]
Ci = np.array([[x % 5 for x in c] for c in Cm], dtype=np.int64)          # module coordinates mod 5
# the half-space test lives in M (x) R = the span of the ambient vectors: use those, normalised
amb = np.array([[float(x) for x in v] for v in E])
amb = amb / np.abs(amb).max(axis=1, keepdims=True)
C = amb; dim = C.shape[1]
U = np.vstack([C, -C]); Ui = np.vstack([Ci, (-Ci) % 5])
t0 = time.time()
if 5 ** r <= 400000:
    PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
else:
    rng = np.random.default_rng(0); PS = rng.integers(0, 5, size=(400000, r))
ok = np.ones(len(PS), bool)
for c in Ci: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
VAL = (ADM @ Ui.T) % 5
print(f"{src}: {len(E)} directions, rank {r}, {len(ADM)} admissible psi   [{time.time()-t0:.0f}s]", flush=True)
def halfspace(D, e):
    """w with <w,u> >= 0 on D and <w,e> >= 1 ?"""
    A_ub = np.vstack([-D, -e[None, :]]) if len(D) else -e[None, :]
    b_ub = np.concatenate([np.zeros(len(D)), [-1.0]])
    res = linprog(np.zeros(dim), A_ub=A_ub, b_ub=b_ub, bounds=[(-1000, 1000)] * dim, method="highs")
    return res.status == 0, (res.x if res.status == 0 else None)
MODE = sys.argv[2] if len(sys.argv) > 2 else "apart"
dirs = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else list(range(len(E)))
nE = len(E)
for k in dirs:
    hit = None; tried = 0
    for kk in (k, k + nE):                    # both orientations of the direction
        e = U[kk]
        for i in range(len(ADM)):
            s = int(VAL[i][kk])
            ts = [(-2 * s) % 5] if MODE == "apart" else [1, 2, 3, 4]
            for t in ts:
                D = U[VAL[i] == t]
                tried += 1
                feas, w = halfspace(D, e)
                if feas: hit = (tuple(int(x) for x in ADM[i]), t, "+" if kk == k else "-"); break
            if hit or tried > 8000: break
        if hit: break
    print(f"  direction {k}: {'REFUTED by a twisted colouring ' + str(hit) if hit else 'no twisted colouring found (' + str(tried) + ' LPs)'}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
