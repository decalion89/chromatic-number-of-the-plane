"""The lattice N spanned by the translations that a grown colouring nearly respects, and the
quotient M / N it would factor through: Smith form, the images of the unit vectors, and whether
the quotient Cayley graph (made finite by a torus in the free part) is 5-colourable.

usage: TH=0.8 python3 quot2d.py <checkpoint.json> [torus sizes ...]"""
import sys, os, json, itertools, time
import numpy as np
import sympy
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "periodfind.py")).read()
src = src[:src.index("# ---- the lattice spanned by the best translations")]
sys.argv = [sys.argv[0], sys.argv[1], "300"] + sys.argv[2:]
TORI = [int(x) for x in sys.argv[3:]] or [10, 15, 20, 30]
exec(src)
TH = float(os.environ.get("TH", "0.8"))
name2t = {nm: t for t, nm in cands.items()}
good = [name2t[nm] for rate, same, cnt, nm in res if same >= TH]
Mx = sympy.Matrix([list(t) for t in good])
# Smith form with transforms: rows of Mx generate N inside Z^r (row convention)
def snf(A):
    A = sympy.Matrix(A); m, n_ = A.shape; Lm = sympy.eye(m); Rm = sympy.eye(n_)
    k = 0
    while k < min(m, n_):
        best = None
        for i in range(k, m):
            for j in range(k, n_):
                if A[i, j] != 0 and (best is None or abs(A[i, j]) < abs(A[best[0], best[1]])): best = (i, j)
        if best is None: break
        i, j = best
        A.row_swap(k, i); Lm.row_swap(k, i); A.col_swap(k, j); Rm.col_swap(k, j)
        done = True
        for i in range(k + 1, m):
            q = A[i, k] // A[k, k]
            if q: A[i, :] = A[i, :] - q * A[k, :]; Lm[i, :] = Lm[i, :] - q * Lm[k, :]
            if A[i, k] != 0: done = False
        for j in range(k + 1, n_):
            q = A[k, j] // A[k, k]
            if q: A[:, j] = A[:, j] - q * A[:, k]; Rm[:, j] = Rm[:, j] - q * Rm[:, k]
            if A[k, j] != 0: done = False
        if done: k += 1
    return A, Lm, Rm
D, Lm, Rm = snf(Mx)
dg = [abs(int(D[i, i])) if i < min(D.shape) else 0 for i in range(r)]
print(f"  N: {len(good)} translations (agreement >= {TH}); Smith invariants of Z^{r}/N: {dg}")
# class of x: (x Rm)_i mod d_i, with d_i = 0 meaning a free Z coordinate
Rn = [[int(Rm[i, j]) for j in range(r)] for i in range(r)]
def cls(x):
    y = [sum(int(x[i]) * Rn[i][j] for i in range(r)) for j in range(r)]
    return tuple(y[j] % dg[j] if dg[j] else y[j] for j in range(r))
free = [j for j in range(r) if dg[j] == 0]; tors = [j for j in range(r) if dg[j] > 1]
print(f"  M/N = Z^{len(free)} x " + " x ".join(f"Z/{dg[j]}" for j in tors))
ucls = [cls(c) for c in cu]
img = {}
for i, c in enumerate(ucls):
    key = tuple(c[j] for j in free + tors); img.setdefault(key, []).append(i)
zero = tuple([0] * (len(free) + len(tors)))
print(f"  the {nu} unit vectors map to {len(img)} classes; into zero: {len(img.get(zero, []))}")
for k_, v in sorted(img.items(), key=lambda kv: -len(kv[1]))[:12]:
    print(f"    {k_}: {len(v)} units")
# how well does the colouring factor through M/N?  majority colour per class
pc = [cls(c) for c in cp]
from collections import Counter, defaultdict
byc = defaultdict(Counter)
for x, k_ in enumerate(pc): byc[tuple(k_[j] for j in free + tors)][int(col[x])] += 1
maj = sum(c.most_common(1)[0][1] for c in byc.values())
print(f"  {len(byc)} classes occupied; majority colour covers {maj}/{n} = {maj / n:.3f} of the points")
if img.get(zero): print("  a unit vector lies in N: no colouring factors through M/N"); sys.exit()
# tori: free part mod T, torsion as is
from pysat.solvers import Solver
dims_t = [dg[j] for j in tors]
for T in TORI:
    dims = [T] * len(free) + dims_t
    Ng = int(np.prod(dims))
    if Ng > 400000: print(f"  torus {T}: {Ng} vertices, skipped"); continue
    st = [1] * len(dims)
    for i in range(len(dims) - 2, -1, -1): st[i] = st[i + 1] * dims[i + 1]
    el = np.array(list(itertools.product(*[range(x) for x in dims])), dtype=np.int64)
    ia = el @ np.array(st)
    dirs = set()
    for k_ in img:
        v = tuple(int(a) % d for a, d in zip(k_, dims))
        if not any(v): dirs = None; break
        if tuple((-a) % d for a, d in zip(v, dims)) not in dirs: dirs.add(v)
    if dirs is None: print(f"  torus {T}: a unit wraps to zero"); continue
    EP = []
    for v in dirs:
        ib = ((el + np.array(v)) % np.array(dims)) @ np.array(st)
        lo, hi = np.minimum(ia, ib), np.maximum(ia, ib); EP.append(lo * Ng + hi)
    EP = np.unique(np.concatenate(EP)); EA, EB = EP // Ng, EP % Ng
    with Solver(name="cadical195") as s:
        for v in range(Ng): s.add_clause([v * 5 + k + 1 for k in range(5)])
        for a, b in zip(EA.tolist(), EB.tolist()):
            for k in range(5): s.add_clause([-(a * 5 + k + 1), -(b * 5 + k + 1)])
        s.add_clause([1])
        ok = s.solve()
    print(f"  torus {T}: {Ng} vertices, {len(dirs)} directions, {len(EA)} edges: 5-colourable {ok}   [{time.time()-t0:.0f}s]", flush=True)
    if ok: break
