"""Does a grown colouring hide a period?  For many small lattice vectors t (differences and sums
of unit vectors, small multiples of units) measure how often c(x + t) = pi(c(x)) for the best
relabelling pi, over all x with both ends in the graph.  A colouring that is the restriction of a
periodic colouring of the whole module shows rate 1 on its period lattice.

usage: python3 periodfind.py <checkpoint.json> [min pairs]"""
import sys, json, os, itertools, time
from fractions import Fraction as Fr
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
t0 = time.time()
d = json.load(open(sys.argv[1])); col = np.array(d["colouring"], dtype=np.int64)
MINP = int(sys.argv[2]) if len(sys.argv) > 2 else 300
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]; A = d["A"]
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
raw = [tuple((p - V[A]).x.c) + tuple((p - V[A]).y.c) for p in V]
den = 1
for v in Ev + raw:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
B = echelon(E); r = len(B)
iv = lambda v: tuple(int(t) for t in coords(B, tuple(int(Fr(x) * den) for x in v)))
cu = np.array([iv(v) for v in Ev], dtype=object); cp = [iv(v) for v in raw]
idx = {c: i for i, c in enumerate(cp)}
n = len(cp); nu = len(cu)
print(f"{os.path.basename(sys.argv[1])}: {n} points, {nu} units, rank {r}   [{time.time()-t0:.0f}s]", flush=True)
cands = {}
for i in range(nu):
    for k in (2, 3, 4, 5, 6):
        cands.setdefault(tuple(k * cu[i]), f"{k}u{i}")
    for j in range(i + 1, nu):
        cands.setdefault(tuple(cu[i] - cu[j]), f"u{i}-u{j}"); cands.setdefault(tuple(cu[i] + cu[j]), f"u{i}+u{j}")
uset = {tuple(c) for c in cu}
cands = {t: nm for t, nm in cands.items() if any(t) and t not in uset}
print(f"  {len(cands)} candidate translations   [{time.time()-t0:.0f}s]", flush=True)
perms = list(itertools.permutations(range(5)))
# linear hash mod a 61-bit prime: h(x + t) = h(x) + h(t), so every lookup is one vectorised search
PR = (1 << 61) - 1
rr = [int(x) for x in np.random.default_rng(99).integers(1, 1 << 60, size=r)]
hh = lambda v: sum(int(a) * b for a, b in zip(v, rr)) % PR
H = np.array([hh(c) for c in cp], dtype=np.int64); order = np.argsort(H); Hs = H[order]
res = []
for t, nm in cands.items():
    q = (H + hh(t)) % PR
    pos = np.searchsorted(Hs, q); pos[pos >= n] = 0
    hit = Hs[pos] == q
    cnt = int(hit.sum())
    if cnt < MINP: continue
    T = np.zeros((5, 5), dtype=np.int64)
    np.add.at(T, (col[np.nonzero(hit)[0]], col[order[pos[hit]]]), 1)
    best = max(sum(T[a][p[a]] for a in range(5)) for p in perms)
    same = np.trace(T)
    res.append((best / cnt, same / cnt, cnt, nm))
res.sort(reverse=True)
print(f"  {len(res)} translations with >= {MINP} pairs   [{time.time()-t0:.0f}s]")
for rate, same, cnt, nm in res[:25]:
    print(f"    {nm:>12}: best-relabelled agreement {rate:.3f}, identical colour {same:.3f}  ({cnt} pairs)")
rates = np.array([x[0] for x in res])
print(f"  agreement quantiles: median {np.median(rates):.3f}, 90% {np.quantile(rates, 0.9):.3f}, max {rates.max():.3f}")
# ---- the lattice spanned by the best translations: its index, and whether it swallows a unit
TH = float(os.environ.get("TH", "0.8"))
import sympy
good = [t for (rate, same, cnt, nm), t in zip(res, [None] * len(res)) if False]
name2t = {nm: t for t, nm in cands.items()}
good = [name2t[nm] for rate, same, cnt, nm in res if same >= TH]
print(f"  {len(good)} translations with identical-colour agreement >= {TH}")
if good:
    Mx = sympy.Matrix([list(t) for t in good])
    rk_ = Mx.rank()
    print(f"  they span a lattice of rank {rk_}")
    if rk_ == r:
        from sympy.matrices.normalforms import hermite_normal_form
        Hm = hermite_normal_form(Mx.T)                    # columns generate the lattice
        idx_ = abs(Hm[:, :r].det()) if Hm.shape[1] >= r else None
        print(f"  index of that lattice in the edge module: {idx_}  (= {sympy.factorint(idx_) if idx_ else None})")
        Hr = Hm[:, :r]
        inL = 0
        for c in cu:
            sol = Hr.LUsolve(sympy.Matrix(list(c)))
            inL += all(x.is_integer for x in sol)
        print(f"  unit vectors lying in it: {inL} of {nu}")
