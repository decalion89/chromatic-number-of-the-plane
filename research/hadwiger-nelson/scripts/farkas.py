"""Exact certificates that NO twisted colouring keeps a given 2e pair alike.

For each admissible psi and each orientation of e, with t = -2 psi(e), a
twisted colouring needs a functional >= 0 on the class D_t and > 0 on e.  By
Farkas it does not exist iff -e is a nonnegative combination of D_t.  Find the
combination by LP, rationalise it, and check the identity EXACTLY on the
integer vectors -- so the negative claim rests on arithmetic, not on the solver.
"""
import sys, json, itertools, time
from fractions import Fraction as Fr
from math import gcd
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/{sys.argv[1]}")); dirs = [int(x) for x in sys.argv[2].split(",")]
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P)); B = echelon(E); r = len(B)
Cm = [coords(B, v) for v in E]
Ci = np.array([[x % 5 for x in c] for c in Cm], dtype=np.int64)
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
ok = np.ones(len(PS), bool)
for c in Ci: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
Uint = [list(v) for v in E] + [[-x for x in v] for v in E]
Ui = np.vstack([Ci, (-Ci) % 5]); VAL = (ADM @ Ui.T) % 5
Uf = np.array([[float(x) for x in v] for v in Uint]); sc = np.abs(Uf).max(axis=1)
nE = len(E); t0 = time.time()
for k in dirs:
    certified = 0; failed = 0
    for kk in (k, k + nE):
        e = np.array(Uint[kk], dtype=float)
        for i in range(len(ADM)):
            t = (-2 * int(VAL[i][kk])) % 5
            idx = [j for j in range(2 * nE) if VAL[i][j] == t]
            # find lambda >= 0 with sum lambda_j u_j = -e  (scaled columns for conditioning)
            A = (Uf[idx] / sc[idx, None]).T
            res = linprog(np.ones(len(idx)), A_eq=A, b_eq=-e, bounds=[(0, None)] * len(idx), method="highs")
            if res.status != 0: failed += 1; continue
            lam = [Fr(x / sc[j]).limit_denominator(10 ** 9) for x, j in zip(res.x, idx)]
            # exact check: small positive support re-solved exactly would be ideal; here test the identity
            tot = [sum(l * Uint[j][c] for l, j in zip(lam, idx)) + Uint[kk][c] for c in range(len(Uint[0]))]
            if all(x == 0 for x in tot) and all(l >= 0 for l in lam): certified += 1
            else:
                # re-solve exactly on the LP's support by rational Gaussian elimination
                supp = [j for x, j in zip(res.x, idx) if x > 1e-12]
                import sympy
                Mx = sympy.Matrix([[Uint[j][c] for j in supp] for c in range(len(Uint[0]))])
                bv = sympy.Matrix([-Uint[kk][c] for c in range(len(Uint[0]))])
                try:
                    sol, params = Mx.gauss_jordan_solve(bv)
                    sol = sol.subs({p: 0 for p in params})
                    if all(s >= 0 for s in sol): certified += 1
                    else: failed += 1
                except ValueError:
                    failed += 1
    print(f"  direction {k}: exact Farkas certificates {certified} of {2 * len(ADM)} (psi, orientation) cases; uncertified {failed}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
