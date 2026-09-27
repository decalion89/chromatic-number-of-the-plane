"""Unwrap the UNSAT core of an ideal torus.

(module set-up as in idealquot.py)

A module whose unit vectors come in hexagons is a Z[omega]-module, so for alpha = a + b*omega the
sublattice alpha*M is a period lattice of index N(alpha)^(rank/2).  A 5-colouring of the finite
Cayley graph Cay(M / alpha M, U) pulls back to an alpha*M-periodic colouring of the whole
unit-distance graph on M.  When N(alpha) is prime to 5 such a colouring is never a coset
colouring (alpha acts invertibly on M/5M, so a coset colouring is never alpha*M-periodic):
finding one would show that the module has colourings the rigidity theorems do not describe.

usage: python3 idealquot.py <module.json> <a> <b> [time limit, s]"""
import sys, json, os, time, itertools
from fractions import Fraction as Fr
from math import gcd
import numpy as np
import sympy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
t0 = time.time()
path, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
TL = float(sys.argv[4]) if len(sys.argv) > 4 else 3600
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(p, q) for p, q in xy[0]]), F.element([Fr(p, q) for p, q in xy[1]]))
if "units" in d:
    U = [mk(xy) for xy in d["units"]]
else:
    P = [mk(xy) for xy in d["points"]]; g = build_graph(P); Ud = {}
    for i, j in g.edges():
        for w in (P[j] - P[i], P[i] - P[j]): Ud.setdefault((round(w.fx, 9), round(w.fy, 9)), w)
    U = list(Ud.values())
# close the unit set under the 60-degree turn (exactly): the edge directions of a finite graph
# need not be all the unit vectors of the module they span
cc = [Fr(0)] * F.dim; cc[1] = Fr(1); S3 = F.element(cc)          # sqrt(3): bitmask 1 -> gens[0] = 3
assert F.gens[0] == 3
turn = lambda p: Point(p.x * Fr(1, 2) - p.y * S3 * Fr(1, 2), p.x * S3 * Fr(1, 2) + p.y * Fr(1, 2))
fkey = lambda x, y: (round(x, 7), round(y, 7))
n0 = len(U); seen = {fkey(u.fx, u.fy) for u in U}
for u in list(U):
    w = u
    for _ in range(5):
        w = turn(w)
        if fkey(w.fx, w.fy) not in seen: seen.add(fkey(w.fx, w.fy)); U.append(w)
nu = len(U)
uidx = {fkey(u.fx, u.fy): i for i, u in enumerate(U)}
rot = [uidx.get(fkey(turn(u).fx, turn(u).fy)) for u in U]
assert all(r is not None for r in rot)
print(f"  {n0} edge directions, {nu} after closing under the 60-degree turn", flush=True)
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
den = 1
for v in Ev:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev[:n0])})
B = echelon(E); rk = len(B)
C0 = [coords(B, tuple(int(Fr(x) * den) for x in v)) for v in Ev]
assert all(Fr(t).denominator == 1 for c in C0 for t in c), "a turned unit is not in the module"
C0 = [[int(t) for t in c] for c in C0]
src = open(os.path.join(HERE, "allunits.py")).read()
exec(src[src.index("def lll_exact"):src.index("Gr, Ur = lll_exact(G)")])
wt = list(F._prod) * 2
G0 = [[Fr(sum(w * p * q for w, p, q in zip(wt, B[i], B[j]))) for j in range(rk)] for i in range(rk)]
_, T = lll_exact(G0)
Ti = sympy.Matrix(T).inv()
CU = sympy.Matrix([[int(x) for x in (sympy.Matrix([c]) * Ti)] for c in C0])       # nu x rk
# omega as an integer matrix: CU[rot[i]] = CU[i] * Om
basis = []
for i in range(nu):
    if sympy.Matrix([CU.row(j) for j in basis + [i]]).rank() > len(basis): basis.append(i)
    if len(basis) == rk: break
Bm = sympy.Matrix([CU.row(i) for i in basis]); Rm = sympy.Matrix([CU.row(rot[i]) for i in basis])
Om = Bm.inv() * Rm
assert all(x == int(x) for x in Om), "omega does not act integrally"
assert CU * Om == sympy.Matrix([CU.row(rot[i]) for i in range(nu)]), "omega matrix fails on some unit"
Al = a * sympy.eye(rk) + b * Om                     # rows of Al generate alpha*M (row convention)
N = abs(int(Al.det()))
print(f"{os.path.basename(path)}: {nu} units, rank {rk}; alpha = {a} + {b} omega, N(alpha) = {a*a + a*b + b*b}, |M / alpha M| = {N}   [{time.time()-t0:.0f}s]", flush=True)
# reduce with the Smith form instead: M / alpha M = sum Z/d_i, coordinates via a unimodular change
from sympy.matrices.normalforms import smith_normal_form
def snf_with_transform(A):
    # returns (D, Lm, Rm) with Lm * A * Rm = D diagonal, Lm, Rm unimodular (simple elimination)
    A = sympy.Matrix(A); m, n = A.shape
    Lm = sympy.eye(m); Rm = sympy.eye(n)
    for k in range(min(m, n)):
        while True:
            # move the smallest nonzero entry of the remaining block to (k, k)
            best = None
            for i in range(k, m):
                for j in range(k, n):
                    if A[i, j] != 0 and (best is None or abs(A[i, j]) < abs(A[best[0], best[1]])): best = (i, j)
            if best is None: return A, Lm, Rm
            i, j = best
            A.row_swap(k, i); Lm.row_swap(k, i); A.col_swap(k, j); Rm.col_swap(k, j)
            done = True
            for i in range(k + 1, m):
                q = A[i, k] // A[k, k]
                if q: A[i, :] = A[i, :] - q * A[k, :]; Lm[i, :] = Lm[i, :] - q * Lm[k, :]
                if A[i, k] != 0: done = False
            for j in range(k + 1, n):
                q = A[k, j] // A[k, k]
                if q: A[:, j] = A[:, j] - q * A[:, k]; Rm[:, j] = Rm[:, j] - q * Rm[:, k]
                if A[k, j] != 0: done = False
            if done:
                # divisibility fix so that d_k | d_{k+1} is not needed for reduction
                break
    return A, Lm, Rm
D, Lm, Rmx = snf_with_transform(Al)                   # Lm * Al * Rmx = D
dg = [abs(int(D[i, i])) for i in range(rk)]
assert int(np.prod(dg)) == N
# x in Z^rk (row vector) is in alpha*M iff x = y Al for integer y iff x Rmx = y Lm^{-1} D ... use:
# alpha*M = Z^rk Al; x Rmx = (y Lm^{-1}) D, so class of x = (x Rmx)_i mod d_i
Rn = np.array([[int(Rmx[i, j]) for j in range(rk)] for i in range(rk)], dtype=object)
CUn = np.array([[int(CU[i, j]) for j in range(rk)] for i in range(nu)], dtype=object)
cls = lambda x: tuple(int(v) % dd for v, dd in zip(np.dot(x, Rn), dg))
ucl = [cls(CUn[i]) for i in range(nu)]
zero = tuple([0] * rk)
loops = sum(1 for c in ucl if c == zero)
print(f"  invariants {sorted(dg)}; unit vectors in alpha*M: {loops}   [{time.time()-t0:.0f}s]", flush=True)
if loops:
    print("  a unit vector lies in alpha*M: no periodic colouring through this ideal"); sys.exit()
# enumerate the group Z/d_1 x ... and build the Cayley graph
dims = dg
strides = [1] * rk
for i in range(rk - 2, -1, -1): strides[i] = strides[i + 1] * dims[i + 1]
idx = lambda c: sum(ci * s for ci, s in zip(c, strides))
dirs = {}
for c in ucl:
    neg = tuple((-x) % dd for x, dd in zip(c, dims))
    if neg not in dirs: dirs[c] = 1
dirs = list(dirs)
el = np.array(list(itertools.product(*[range(x) for x in dims])), dtype=np.int64)
st_ = np.array(strides, dtype=np.int64); dm_ = np.array(dims, dtype=np.int64)
ia_ = el @ st_
EP = []
for c in dirs:
    ib_ = ((el + np.array(c, dtype=np.int64)) % dm_) @ st_
    lo, hi = np.minimum(ia_, ib_), np.maximum(ia_, ib_)
    keep_ = lo != hi
    EP.append(lo[keep_] * N + hi[keep_])
EP = np.unique(np.concatenate(EP)); EA_, EB_ = EP // N, EP % N; del EP
edges = None
print(f"  Cayley graph: {N} vertices, {len(dirs)} directions, {len(EA_)} edges   [{time.time()-t0:.0f}s]", flush=True)

# ---------------------------------------------------------------------------------------------
# The torus T = Cay(M / alpha M, U) has no 5-colouring.  Extract an UNSAT core (vertex selectors,
# failed assumptions, then deletion passes), and ask whether the core can be LIFTED to M: lifts
# lam(v) in v + alpha Z^r with lam(b) - lam(a) a genuine unit vector on every core edge.  A lift
# is an honest finite unit-distance graph isomorphic to the core, hence not 5-colourable.
# A core that only exists because it wraps around the torus cannot be lifted.
from pysat.solvers import Solver
K = 5
var = lambda v, k: v * K + k + 1
Nt = N
sel = lambda v: Nt * K + v + 1
adj = [[] for _ in range(Nt)]
for a_, b_ in zip(EA_.tolist(), EB_.tolist()): adj[a_].append(b_); adj[b_].append(a_)
def core_of(verts, rounds=1):
    verts = sorted(verts)
    with Solver(name="cadical195") as s:
        vs = set(verts)
        for v in verts:
            s.add_clause([-sel(v)] + [var(v, k) for k in range(K)])
        for a_, b_ in zip(EA_.tolist(), EB_.tolist()):
            if a_ in vs and b_ in vs:
                for k in range(K): s.add_clause([-var(a_, k), -var(b_, k)])
        ok = s.solve(assumptions=[sel(v) for v in verts])
        if ok: return None
        core = sorted(abs(l) - Nt * K - 1 for l in s.get_core())
    return core
t1 = time.time()
core = core_of(range(Nt))
print(f"  failed-assumption core: {len(core)} of {Nt} vertices   [{time.time()-t0:.0f}s]", flush=True)
for it in range(int(os.environ.get("SHRINK", "3"))):
    c2 = core_of(core)
    if c2 is None or len(c2) >= len(core): break
    core = c2
    print(f"    re-solve on the core: {len(core)} vertices   [{time.time()-t0:.0f}s]", flush=True)
# lift: BFS over the core's induced subgraph; each torus edge class may come from several units
cs = set(core)
units_of_class = {}
CUi = np.array([[int(CU[i, j]) for j in range(rk)] for i in range(nu)], dtype=object)
for i, c in enumerate(ucl): units_of_class.setdefault(tuple(c), []).append(i)
elv = el                                     # element index -> class vector
def cls_of_idx(p): return tuple(int(x) for x in elv[p])
def diff_class(a_, b_):
    return tuple(int((x - y) % dd) for x, y, dd in zip(elv[b_], elv[a_], dims))
lift = {}; tree_bad = 0; edges_core = 0; bad = 0; multi = 0
from collections import deque
for root in core:
    if root in lift: continue
    lift[root] = np.zeros(rk, dtype=object)       # arbitrary base lift per component (translation)
    dq = deque([root])
    while dq:
        a_ = dq.popleft()
        for b_ in adj[a_]:
            if b_ not in cs or b_ in lift: continue
            reps = units_of_class.get(diff_class(a_, b_), [])
            if not reps: tree_bad += 1; continue
            if len(reps) > 1: multi += 1
            lift[b_] = lift[a_] + CUi[reps[0]]
            dq.append(b_)
uset = {tuple(int(x) for x in r) for r in CUi}
for a_, b_ in zip(EA_.tolist(), EB_.tolist()):
    if a_ in cs and b_ in cs:
        edges_core += 1
        dv = tuple(int(x) for x in (lift[b_] - lift[a_]))
        if dv not in uset: bad += 1
print(f"  core: {len(core)} vertices, {edges_core} edges; tree steps with several unit representatives: {multi}")
print(f"  greedy lift: {edges_core - bad} of {edges_core} core edges become genuine unit steps, {bad} do not")
json.dump({"module": os.path.basename(path), "alpha": [a, b], "core": core}, open(os.environ.get("OUT", "/tmp/torecore.json"), "w"))
