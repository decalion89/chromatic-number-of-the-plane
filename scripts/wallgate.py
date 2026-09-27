"""The wall gate: 1-dimensional colourings c(x) = G(phi(x), psi(x)).

phi an integer functional on M, psi : M -> Z/5.  Any proper 5-colouring G of
Cay(Z x Z/5, {(phi(u), psi(u))}) pulls back to M.  Twisted colourings are the
staircases G(n, a) = a + t*floor(n/L); here G is arbitrary -- a wall of any
thickness, or any 1-D pattern.  phi runs over the functionals with small values
on every unit (|phi(u)| <= B), found by Fincke-Pohst on sum_u phi(u)^2; psi
over all of Z/5-hom (admissible or not: a unit with psi(u) = 0 only needs
phi(u) != 0).  For each (phi, psi) a window n in [-L, L] is solved with the
witness at the centre:
  apart : G(0,0) = G(2 phi(e), 2 psi(e)) for some direction e
  pair  : G(0,0) != G(5 phi(e), 0)       for some direction e
SAT on a long window is evidence, not a colouring of M; a hit is re-searched
as a periodic colouring (Z/N x Z/5) before it counts.

usage: wallgate.py <module.json> <apart|pair> <B> [L]
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, time
from fractions import Fraction as Fr
import numpy as np
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
ROOT = HN_DIR
name, MODE, BND = sys.argv[1], sys.argv[2], int(sys.argv[3]); L = int(sys.argv[4]) if len(sys.argv) > 4 else 40
d = json.load(open(f"{ROOT}/data/{name}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P))
B = echelon(E); C = np.array([coords(B, v) for v in E], dtype=object); r = len(B); nd = len(C)
t0 = time.time()
# ---- functionals with |phi(u)| <= BND on every unit.  The echelon coordinates are huge, so
# LLL-reduce the value lattice {(phi(u))_u : phi in Z^r} first (exact, sympy), then Fincke-Pohst.
src = open(HN_DIR + "/scripts/allunits.py").read()
exec(src[src.index("def lll_exact"):src.index("Gr, Ur = lll_exact(G)")])
G0 = [[Fr(sum(int(C[k][i]) * int(C[k][j]) for k in range(nd))) for j in range(r)] for i in range(r)]
_, Tm = lll_exact(G0)                                                  # rows: new functionals in old coordinates
Bv = [[sum(Tm[i][j] * int(C[u][j]) for j in range(r)) for u in range(nd)] for i in range(r)]
G = np.array([[sum(a * b for a, b in zip(Bv[i], Bv[j])) for j in range(r)] for i in range(r)], dtype=float)
R = np.linalg.cholesky(G).T                                            # upper: k^T G k = |R k|^2
rad = nd * BND * BND + 1e-6
found = []
def fp(k, kv, partial):
    if k < 0:
        if any(kv):
            vals = [sum(kv[i] * Bv[i][u] for i in range(r)) for u in range(nd)]
            if max(abs(v) for v in vals) <= BND:
                phi = [sum(kv[i] * Tm[i][j] for i in range(r)) for j in range(r)]
                found.append(tuple(phi))
        return
    s = sum(R[k][j] * kv[j] for j in range(k + 1, r))
    rem = rad - partial
    if rem < 0: return
    w = (rem ** 0.5) / abs(R[k][k]); c0 = -s / R[k][k]
    for v in range(int(np.ceil(c0 - w - 1e-9)), int(np.floor(c0 + w + 1e-9)) + 1):
        kv[k] = v
        fp(k - 1, kv, partial + (R[k][k] * v + s) ** 2)
    kv[k] = 0
fp(r - 1, [0] * r, 0.0)
phis = []
for p_ in found:
    if tuple(-x for x in p_) not in phis: phis.append(p_)
print(f"{name}: {nd} directions, rank {r}; reduced value-vector norms {[sum(x*x for x in b) for b in Bv]}, max |value| {[max(abs(x) for x in b) for b in Bv]}; "
      f"{len(phis)} functionals (up to sign) with |phi(u)| <= {BND} on every unit   [{time.time()-t0:.0f}s]", flush=True)
allpsi = list(itertools.product(range(5), repeat=r))
ADMONLY = __import__('os').environ.get('ADMONLY', '1') == '1'
Cn = np.array([[int(x) % 5 for x in row] for row in C], dtype=np.int64)          # psi needs C mod 5 only
X = lambda v, c: 1 + v * 5 + c
alive = set(range(nd)); hits = []
for pi, phi in enumerate(phis):
    fv = np.array([sum(int(c) * p for c, p in zip(row, phi)) for row in C], dtype=np.int64)   # phi(u) per direction
    PS = (np.array(allpsi, dtype=np.int64) @ Cn.T) % 5   # psi(u) per psi, per direction
    ok = ~np.any((PS == 0) & (fv[None, :] == 0), axis=1)  # no unit on (0, 0)
    if ADMONLY: ok &= ~np.any(PS == 0, axis=1)
    cand = np.nonzero(ok)[0]
    width = 2 * L + 1; N = 5 * width
    vid = lambda n, a: (n + L) * 5 + a
    for ci in cand:
        psv = PS[ci]
        S = set()
        for k in range(nd):
            S.add((int(fv[k]), int(psv[k]))); S.add((-int(fv[k]), int(-psv[k]) % 5))
        cl = [[X(v, c) for c in range(5)] for v in range(N)]
        for n in range(-L, L + 1):
            for a in range(5):
                for (s, t) in S:
                    m = n + s
                    if -L <= m <= L:
                        i, j = vid(n, a), vid(m, (a + t) % 5)
                        if i < j:
                            for c in range(5): cl.append([-X(i, c), -X(j, c)])
        cl.append([X(vid(0, 0), 0)])
        sol = Solver(name="cd19", bootstrap_with=cl)
        top = N * 5 + 1; sel = []
        for k in alive:
            if MODE == "apart":
                n2, a2 = 2 * int(fv[k]), (2 * int(psv[k])) % 5
                if abs(n2) > L: continue
                z = top; top += 1; sol.add_clause([-z, X(vid(n2, a2), 0)]); sel.append((z, k))
            else:
                n5 = 5 * int(fv[k])
                if abs(n5) > L: continue
                z = top; top += 1; sol.add_clause([-z, -X(vid(n5, 0), 0)]); sel.append((z, k))
        if not sel: sol.delete(); continue
        sol.add_clause([z for z, _ in sel])
        if sol.solve():
            m = sol.get_model()
            col = lambda v: next(c for c in range(5) if m[X(v, c) - 1] > 0)
            hit = [k for z, k in sel if m[z - 1] > 0]
            hits.append((phi, tuple(int(x) for x in allpsi[ci]), hit))
            print(f"  phi #{pi} {phi}, psi {allpsi[ci]}: window colouring refutes {len(hit)} directions (e.g. {hit[:5]})", flush=True)
            alive -= set(hit)
        sol.delete()
    print(f"  phi #{pi} {phi}: max|phi(u)| = {int(max(abs(fv)))}, {len(cand)} psi loopless; directions left {len(alive)}   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"hits": hits}, open(__import__("os").environ.get("WALLOUT", f"wall_{MODE}_{BND}.json"), "w"))
print(f"  done: {len(hits)} window hits; directions unrefuted {len(alive)} of {nd}   [{time.time()-t0:.0f}s]", flush=True)
