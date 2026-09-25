"""The triple relaxation of the circular gate on the lambda-closure M' = M1 + lambda M1.

A circular colouring phi' of M' (with the units U1, lambda U1 and w_e = 5(1-lambda)e) restricts to
phi = phi'|M1 and phi~ = phi' o lambda, and needs
    phi in F(M1),  phi~ in F(M1),  chi = 5(phi - phi~) in F(M1).
This MILP drops the compatibility of phi and phi~ on M1 cap lambda M1, so it is a relaxation:
INFEASIBLE proves that the lambda-closure has no circular 5-colouring.
Variables: phi, phi~ in [0,1]^r; one integer per direction for each of the three families; slack s.

usage: python3 circtriple.py <M1.json> <time limit s> [k]"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, os, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
src = open(os.path.join(HN_DIR, "scripts", "circgate.py")).read()
pre = src[:src.index("# optional: also require the multiples")]
M1 = sys.argv[1]; TL = float(sys.argv[2]); KC = int(sys.argv[3]) if len(sys.argv) > 3 else 5
sys.argv = ["circgate.py", M1, str(TL)]
exec(pre)
reps, seenr = [], set()
for i in range(len(CU)):
    k_ = tuple(CU[i]); nk_ = tuple(-CU[i])
    if nk_ in seenr or k_ in seenr: continue
    seenr.add(k_); reps.append(i)
C = CU[reps].astype(float); m = len(reps)
lo = 1.0 / KC; hi = 1.0 - 1.0 / KC
# variable layout: phi (r) | phit (r) | n1 (m) | n2 (m) | n3 (m) | s
nv = 2 * r + 3 * m + 1
def block(coef_phi, coef_phit, which):
    A = np.zeros((m, nv))
    A[:, :r] = coef_phi; A[:, r:2 * r] = coef_phit
    A[np.arange(m), 2 * r + which * m + np.arange(m)] = -1.0
    return A
rows_lo, rows_hi = [], []
for which, (cp, ct) in enumerate([(C, 0 * C), (0 * C, C), (5 * C, -5 * C)]):
    A = block(cp, ct, which)
    A1 = A.copy(); A1[:, -1] = -1.0; rows_lo.append(A1)      # val - n - s >= lo
    A2 = A.copy(); A2[:, -1] = 1.0; rows_hi.append(A2)       # val - n + s <= hi
cons = [LinearConstraint(np.vstack(rows_lo), lb=lo, ub=np.inf), LinearConstraint(np.vstack(rows_hi), lb=-np.inf, ub=hi)]
bl = lambda M_: (np.floor(np.minimum(M_, 0).sum(axis=1)) - 1, np.ceil(np.maximum(M_, 0).sum(axis=1)) + 1)
l1, u1 = bl(C); l3, u3 = bl(np.hstack([5 * C, -5 * C]))
lb = np.concatenate([np.zeros(2 * r), l1, l1, l3, [0.0]]); ub = np.concatenate([np.ones(2 * r), u1, u1, u3, [0.5]])
integ = np.concatenate([np.zeros(2 * r), np.ones(3 * m), [0]])
cobj = np.zeros(nv); cobj[-1] = -1.0
print(f"{os.path.basename(M1)}: {m} directions, rank {r}; triple MILP with {nv} variables, window [1/{KC}, {KC-1}/{KC}]   [{time.time()-t0:.0f}s]", flush=True)
res = milp(cobj, constraints=cons, integrality=integ, bounds=Bounds(lb, ub), options={"time_limit": TL})
print(f"  status {res.status}: {res.message}", flush=True)
if res.x is not None:
    print(f"  slack {res.x[-1]:.6f}; phi = {np.round(res.x[:r], 6)}; phi~ = {np.round(res.x[r:2*r], 6)}")
print(f"  [{time.time()-t0:.0f}s]")
