"""A module rich in unit vectors: does the twisted obstruction die?

Units of K = Q(sqrt-3, sqrt-11) integral at 5, as planar vectors with
coordinates in Q(sqrt3, sqrt11): products of w6 = (1+sqrt-3)/2,
mu = (5+sqrt-11)/6, k27 = (5+8 sqrt-11)/27, k46 = (35+9 sqrt-11)/46,
r7 = (1+4 sqrt-3)/7, r13 = (-11+4 sqrt-3)/13.  Their Z-span has rank <= 8 but
many unit directions; with enough of them every class psi = t should
positively span, and no twisted coset colouring survives.
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
t0 = time.time()
F = Field((3, 11)); r3, r11 = F.sqrt(3), F.sqrt(11)
q = lambda a, b=1: F.rational(Fr(a, b))
GEN = {"w6": (q(1, 2), r3 * q(1, 2)), "mu": (q(5, 6), r11 * q(1, 6)), "k27": (q(5, 27), r11 * q(8, 27)),
       "k46": (q(35, 46), r11 * q(9, 46)), "r7": (q(1, 7), r3 * q(4, 7)), "r13": (q(-11, 13), r3 * q(4, 13))}
def mul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def inv(a): return (a[0], -a[1])
def pw(a, k):
    r = (q(1), q(0)); b = a if k >= 0 else inv(a)
    for _ in range(abs(k)): r = mul(r, b)
    return r
names = sys.argv[1].split(",") if len(sys.argv) > 1 else ["mu", "k27", "k46", "r7", "r13"]
EXP = int(sys.argv[2]) if len(sys.argv) > 2 else 1
units = {}
for a in range(6):
    for ks in itertools.product(range(-EXP, EXP + 1), repeat=len(names)):
        z = pw(GEN["w6"], a)
        for nm, k in zip(names, ks): z = mul(z, pw(GEN[nm], k))
        v = tuple(z[0].c) + tuple(z[1].c)
        units[v] = z
den = 1
for v in units:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in units)})
B = echelon(E); r = len(B)
Cm = [coords(B, v) for v in E]
Ci = np.array([[x % 5 for x in c] for c in Cm], dtype=np.int64)
print(f"generators {names} exponents <= {EXP}: {len(E)} directions, module rank {r}, den {den}   [{time.time()-t0:.0f}s]", flush=True)
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
ok = np.ones(len(PS), bool)
for c in Ci: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
print(f"  admissible coset colourings: {len(ADM)}   [{time.time()-t0:.0f}s]", flush=True)
if not len(ADM): sys.exit()
amb = np.array([[float(x) for x in v] for v in E]); amb /= np.abs(amb).max(axis=1, keepdims=True)
U = np.vstack([amb, -amb]); Ui = np.vstack([Ci, (-Ci) % 5]); VAL = (ADM @ Ui.T) % 5
dim = amb.shape[1]
def halfspace(D, e):
    res = linprog(np.zeros(dim), A_ub=np.vstack([-D, -e[None, :]]), b_ub=np.concatenate([np.zeros(len(D)), [-1.0]]),
                  bounds=[(-1000, 1000)] * dim, method="highs")
    return res.status == 0
# fast: for each class (psi, t) decide first whether its dual cone is {0}
def dual_nonzero(D):
    for i in range(dim):
        for sgn in (1, -1):
            c = np.zeros(dim); c[i] = -sgn                      # maximise sgn * w_i
            res = linprog(c, A_ub=-D, b_ub=np.zeros(len(D)), bounds=[(-1, 1)] * dim, method="highs")
            if res.status == 0 and -res.fun > 1e-9: return True
    return False
live = {}
for i in range(len(ADM)):
    for t in (1, 2, 3, 4):
        D = U[VAL[i] == t]
        if dual_nonzero(D): live[(i, t)] = D
print(f"  classes (psi, t) whose dual cone is not zero: {len(live)} of {4 * len(ADM)}   [{time.time()-t0:.0f}s]", flush=True)
surv_apart = []; surv_pair = []
nE = len(E)
for k in range(nE):
    # both orientations: the pair (a, a+2e) read from either end is the same pair
    ap = False; pr = False
    for kk in (k, k + nE):
        e = U[kk]
        ap = ap or any(halfspace(live[(i, (-2 * VAL[i][kk]) % 5)], e) for i in range(len(ADM)) if (i, (-2 * VAL[i][kk]) % 5) in live)
        pr = pr or any(halfspace(D, e) for D in live.values())
    if not ap: surv_apart.append(k)
    if not pr: surv_pair.append(k)
print(f"  directions NOT refuted by any twisted colouring: apart {len(surv_apart)} of {len(E)}, pair {len(surv_pair)} of {len(E)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"names": names, "EXP": EXP, "den": den, "directions": [list(v) for v in E], "surv_apart": surv_apart, "surv_pair": surv_pair},
          open(f"/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/richmod_{'_'.join(names)}_{EXP}.json", "w"))
