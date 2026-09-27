"""Combinatorial r-balls of five_rho7 in exact lattice coordinates, and what their colourings
do on the unit circles at the centre.

B_r = { u_1 + ... + u_s : s <= r } for the module's unit vectors, built in the integer
coordinates of the edge module (rank 8) so that r = 3 (about a million points) is cheap.
The ball is peeled to its 5-core with the centre's closed neighbourhood protected (a vertex
of degree <= 4 can always be coloured last, so peeling never changes a colouring question
about the protected vertices), coloured by tabu search, and then for every centre x whose
whole neighbourhood lies in the ball the pairs x + u, x + v with u - v in 5M are compared
with the other non-adjacent pairs: in a coset colouring the first are always alike.

usage: python3 ballr.py <r> [seeds] [out.json]"""
import sys, json, os, time, subprocess, itertools
from fractions import Fraction as Fr
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
t0 = time.time()
R = int(sys.argv[1]); SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
OUT = sys.argv[3] if len(sys.argv) > 3 else None
d = json.load(open(os.path.join(os.path.dirname(HERE), "data", "five_rho7.json")))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
U = {}
for a, b in g.edges():
    for w in (P[b] - P[a], P[a] - P[b]): U.setdefault(key(w), w)
U = list(U.values()); nu = len(U)
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
den = 1
for v in Ev:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
B = echelon(E); rk = len(B)
C0 = [[int(t) for t in coords(B, tuple(int(Fr(x) * den) for x in v))] for v in Ev]      # huge
# LLL-reduce the lattice under the trace form (every unit has the same length there),
# so that the units get small integer coordinates in the reduced basis
src = open(os.path.join(HERE, "allunits.py")).read()
exec(src[src.index("def lll_exact"):src.index("Gr, Ur = lll_exact(G)")])
Bv = [list(b) for b in B]
wt = list(F._prod) * 2                        # the trace form: coefficient of sqrt(q) weighs q
G0 = [[Fr(sum(w * a * b for w, a, b in zip(wt, Bv[i], Bv[j]))) for j in range(rk)] for i in range(rk)]
_, T = lll_exact(G0)                          # new basis rows = T * old basis rows
# coordinates in the new basis: c' = c T^{-1}
import sympy
Ti = sympy.Matrix(T).inv()
CU = np.array([[int(x) for x in (sympy.Matrix([c]) * Ti)] for c in C0], dtype=np.int64)
print(f"  unit coordinates after LLL: max |c| = {np.abs(CU).max()}", flush=True)
assert CU.shape == (nu, rk)
# hash an 8-vector of integers to one int64 key (random odd multipliers, wrapping arithmetic);
# every lookup below is re-checked coordinate by coordinate, so a collision cannot go unnoticed
W = np.random.default_rng(12345).integers(1, 1 << 62, size=rk, dtype=np.int64) * 2 + 1
def pack(A):
    with np.errstate(over="ignore"): return (A * W).sum(axis=1)
pts = np.zeros((1, rk), dtype=np.int64)
layer = pts
for s in range(R):
    nxt = (layer[:, None, :] + CU[None, :, :]).reshape(-1, rk)
    allp = np.concatenate([pts, nxt]); k = pack(allp)
    _, first = np.unique(k, return_index=True)
    new = allp[np.sort(first)]
    layer = new[len(pts):] if len(new) > len(pts) else new[:0]
    pts = new
    print(f"  B_{s + 1}: {len(pts)} points   [{time.time() - t0:.0f}s]", flush=True)
n = len(pts); K = pack(pts); order = np.argsort(K); Ks = K[order]
# edges: for each unit direction u (one of each +-pair) look up p + u
EA, EB = [], []
done = set()
for i in range(nu):
    neg = tuple(-CU[i])
    if neg in done: continue
    done.add(tuple(CU[i]))
    Q = pts + CU[i]; q = pack(Q); pos = np.searchsorted(Ks, q); pos[pos >= n] = 0
    hit = Ks[pos] == q
    assert (pts[order[pos[hit]]] == Q[hit]).all(), "hash collision"

    EA.append(np.nonzero(hit)[0]); EB.append(order[pos[hit]])
EA = np.concatenate(EA); EB = np.concatenate(EB)
deg = np.bincount(EA, minlength=n) + np.bincount(EB, minlength=n)
print(f"B_{R}: {n} points, {len(EA)} edges, mean degree {2 * len(EA) / n:.1f}, full-degree vertices {(deg == nu).sum()}   [{time.time() - t0:.0f}s]", flush=True)
# 5-core, protecting the centre and its neighbours
prot = np.zeros(n, dtype=bool); prot[np.nonzero((np.abs(pts).sum(axis=1) == 0))[0]] = True
c0 = np.nonzero(prot)[0][0]
prot[EB[EA == c0]] = True; prot[EA[EB == c0]] = True
alive = np.ones(n, dtype=bool)
while True:
    dg = np.bincount(EA[alive[EA] & alive[EB]], minlength=n) + np.bincount(EB[alive[EA] & alive[EB]], minlength=n)
    kill = alive & ~prot & (dg <= 4)
    if not kill.any(): break
    alive &= ~kill
keep = np.nonzero(alive)[0]; remap = -np.ones(n, dtype=np.int64); remap[keep] = np.arange(len(keep))
m = alive[EA] & alive[EB]; ea, eb = remap[EA[m]], remap[EB[m]]
print(f"5-core: {len(keep)} points, {len(ea)} edges, mean degree {2 * len(ea) / len(keep):.1f}   [{time.time() - t0:.0f}s]", flush=True)
cp = pts[keep]; nk = len(keep)
kk = pack(cp); ordk = np.argsort(kk); kks = kk[ordk]
def look(Q):
    q = pack(Q); pos = np.searchsorted(kks, q); pos[pos >= nk] = 0
    r = np.where(kks[pos] == q, ordk[pos], -1)
    assert (cp[r[r >= 0]] == Q[r >= 0]).all(), "hash collision"
    return r
# unit-pair classes: adjacent / 5M / other
cuset = {tuple(r) for r in CU}
pair5 = np.zeros((nu, nu), dtype=np.int8)      # 0 skip, 1 5M, 2 other
for i in range(nu):
    for j in range(i + 1, nu):
        dv = CU[i] - CU[j]
        if tuple(dv) in cuset or not dv.any(): continue
        pair5[i, j] = 1 if (dv % 5 == 0).all() else 2
# centres: core vertices that keep at least THR of their unit neighbours inside the core
THR = int(os.environ.get("THR", "120"))
dgk = np.bincount(ea, minlength=nk) + np.bincount(eb, minlength=nk)
full = np.nonzero(dgk >= THR)[0]
nbr = np.stack([look(cp[full] + CU[i]) for i in range(nu)], axis=1)   # len(full) x nu
print(f"  {len(full)} centres keep >= {THR} of their {nu} neighbours in the core", flush=True)
TAB = os.path.join(HERE, "tabucol")
if not os.path.exists(TAB): subprocess.run(["gcc", "-O2", "-o", TAB, os.path.join(HERE, "tabucol.c")], check=True)
I5 = np.argwhere(pair5 == 1); IO = np.argwhere(pair5 == 2)
res = []
for sd in range(1, SEEDS + 1):
    inp = f"{nk} {len(ea)} 5 2000000000 {sd}\n" + "\n".join(f"{a} {b}" for a, b in zip(ea.tolist(), eb.tolist())) + "\n" + "\n".join(["-1"] * nk) + "\n"
    out = subprocess.run([TAB], input=inp, capture_output=True, text=True).stdout.split("\n")
    if not out[0].startswith("OK"):
        print(f"  seed {sd}: tabu {out[0]}   [{time.time() - t0:.0f}s]", flush=True); continue
    colv = np.array([int(x) for x in out[1:1 + nk]], dtype=np.int64)
    a5 = n5 = ao = no = 0
    for row in nbr:
        cn = np.where(row >= 0, colv[row], -1)
        ok5 = (cn[I5[:, 0]] >= 0) & (cn[I5[:, 1]] >= 0); oko = (cn[IO[:, 0]] >= 0) & (cn[IO[:, 1]] >= 0)
        e5 = cn[I5[ok5, 0]] == cn[I5[ok5, 1]]; eo = cn[IO[oko, 0]] == cn[IO[oko, 1]]
        a5 += e5.sum(); n5 += len(e5); ao += eo.sum(); no += len(eo)
    print(f"  seed {sd}: tabu {out[0]};  on the {len(full)} full centres: 5M pairs alike {a5}/{n5} = {a5 / max(n5, 1):.3f},"
          f" other pairs {ao}/{no} = {ao / max(no, 1):.3f}   [{time.time() - t0:.0f}s]", flush=True)
    res.append(colv.tolist())
if OUT and res:
    json.dump({"r": R, "core_points": cp.tolist(), "edges": [ea.tolist(), eb.tolist()], "colourings": res}, open(OUT, "w"))
