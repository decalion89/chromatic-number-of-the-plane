"""Exact certificates that NO twisted colouring exists at all on a module.

A twisted colouring psi + t floor(phi/L) needs a real functional phi != 0 with
phi >= 0 on the class D_t = {units u : psi(u) = t}.  By Stiemke's lemma there is
none iff D_t spans and has a strictly positive linear dependency
sum_u lambda_u u = 0, lambda_u > 0.  One LP per (psi, t) finds lambda >= 1;
all but a basis of the lambdas are rounded to rationals and the basis
coefficients are then solved EXACTLY -- the certificate is the exact identity
with every coefficient positive, checked in integer arithmetic.
Then: every coarse colouring psi + f (f constant on cells much larger than a
unit step) is a coset colouring, since a facet with jump delta would need
D_{-delta} in a half-space.
"""
import sys, json, itertools, time
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/{sys.argv[1]}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P)); B = echelon(E); r = len(B)
C = [coords(B, v) for v in E]                    # integer coordinates in a basis of M (huge: used mod 5 only)
Uint = [list(v) for v in E] + [[-x for x in v] for v in E]   # small integer FIELD coordinates for the geometry
DIM = len(Uint[0])
Ci = np.array([[x % 5 for x in c] for c in C], dtype=np.int64)
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
ok = np.ones(len(PS), bool)
for c in Ci: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
Ui = np.vstack([Ci, (-Ci) % 5]); VAL = (ADM @ Ui.T) % 5
Uf = np.array([[float(x) for x in v] for v in Uint]); sc = np.abs(Uf).max(axis=1)
def exact_rank(vecs):
    M = [[Fr(x) for x in v] for v in vecs]; rk = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
def solve_exact(Acols, b):
    """solve sum_k x_k Acols[k] = b exactly (Acols: r linearly independent integer vectors)."""
    n = len(Acols); M = [[Fr(Acols[k][i]) for k in range(n)] + [Fr(b[i])] for i in range(len(b))]
    row = 0; piv_cols = []
    for c in range(n):
        p = next((i for i in range(row, len(M)) if M[i][c] != 0), None)
        if p is None: return None
        M[row], M[p] = M[p], M[row]
        for i in range(len(M)):
            if i != row and M[i][c] != 0:
                f = M[i][c] / M[row][c]; M[i] = [a - f * bb for a, bb in zip(M[i], M[row])]
        piv_cols.append(c); row += 1
    if any(M[i][n] != 0 for i in range(row, len(M))): return None
    return [M[i][n] / M[i][i] for i in range(n)]
t0 = time.time(); cert = 0; fail = []
for i in range(len(ADM)):
    for t in range(1, 5):
        idx = [j for j in range(len(Uint)) if VAL[i][j] == t]
        if exact_rank([Uint[j] for j in idx]) < r: fail.append((i, t, "rank")); continue
        A = (Uf[idx] / sc[idx, None]).T
        res = linprog(np.zeros(len(idx)), A_eq=A, b_eq=np.zeros(DIM), bounds=[(1.0, None)] * len(idx), method="highs")
        if res.status != 0: fail.append((i, t, "lp")); continue
        lam = res.x / sc[idx]                             # coefficients on the integer vectors
        # choose a basis among idx (greedy, exact), round the rest, solve the basis exactly
        basis = []
        for j in sorted(range(len(idx)), key=lambda q: -lam[q]):
            if exact_rank([Uint[idx[q]] for q in basis + [j]]) > len(basis): basis.append(j)
            if len(basis) == r: break
        rest = [q for q in range(len(idx)) if q not in basis]
        lr = {q: Fr(float(lam[q])).limit_denominator(10 ** 6) for q in rest}
        rhs = [-sum(lr[q] * Uint[idx[q]][c] for q in rest) for c in range(DIM)]
        lb = solve_exact([Uint[idx[q]] for q in basis], rhs)
        if lb is not None and all(x > 0 for x in lb) and all(x > 0 for x in lr.values()):
            full = {**lr, **{q: x for q, x in zip(basis, lb)}}
            assert all(sum(full[q] * Uint[idx[q]][c] for q in range(len(idx))) == 0 for c in range(DIM))
            cert += 1
        else: fail.append((i, t, "rounding"))
    if (i + 1) % 120 == 0: print(f"  {i+1}/{len(ADM)} psi: {cert} exact certificates, {len(fail)} failures   [{time.time()-t0:.0f}s]", flush=True)
print(f"{sys.argv[1]}: rank {r}, {len(E)} directions, {len(ADM)} admissible psi -> exact Stiemke certificates {cert} of {4*len(ADM)} (psi, t); failures {fail[:5]}   [{time.time()-t0:.0f}s]")
