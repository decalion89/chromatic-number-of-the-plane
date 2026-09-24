"""Colourings periodic under ker(psi) cap alpha*M, and whether they can split a 5e pair.

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
sys.path = [p for p in sys.path if "scratchpad" not in p]
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
# Colourings periodic under L = ker(psi) cap alpha*M.  Since 5 and N(alpha) are coprime,
# M / L = Z/5 x M/alpha M, and a colouring is ANY function f(psi(x), x mod alpha M) that is proper:
# f(s + psi(u), y + u) != f(s, y).  The coset colourings f(s, y) = s are among them; the twisted
# family f = s + k(y) is a small part.  Query: can such a colouring split a 5e pair,
# f(0, 0) != f(0, 5e)?  (5e has psi-value 0.)  SAT would be a non-coset colouring of the whole
# module that splits a 5e pair; UNSAT rules out this whole periodic family for that pair.
Cm5 = [[int(t) % 5 for t in c] for c in C0]
adm = [psi for psi in itertools.product(range(5), repeat=rk) if all(sum(p * c for p, c in zip(psi, cu)) % 5 for cu in Cm5)]
PSIS = [int(x) for x in os.environ.get("PSIS", os.environ.get("PSI", "0")).split(",")]
EUS = [int(x) for x in os.environ.get("EUS", "0").split(",")]
K = 5
var = lambda v, k: v * K + k + 1
KISSAT = os.environ["KISSAT"]; CNF = os.environ.get("CNF", "/tmp/iq2.cnf")
uclA = np.array(ucl, dtype=np.int64)
import subprocess
summary = []
for PSI in PSIS:
    psi = adm[PSI % len(adm)]
    pv = np.array([sum(p * c for p, c in zip(psi, cu)) % 5 for cu in Cm5], dtype=np.int64)
    NV = 5 * N
    seen = set(); EP = []
    yy = ia_
    for i in range(nu):
        key = (int(pv[i]), tuple(uclA[i].tolist()))
        neg = (int((-pv[i]) % 5), tuple(((-uclA[i]) % dm_).tolist()))
        if key in seen or neg in seen: continue
        seen.add(key)
        yb = ((el + uclA[i]) % dm_) @ st_
        for s_ in range(5):
            a_ = s_ * N + yy; b_ = ((s_ + int(pv[i])) % 5) * N + yb
            lo, hi = np.minimum(a_, b_), np.maximum(a_, b_)
            EP.append(lo * NV + hi)
    EP = np.unique(np.concatenate(EP)); EA2, EB2 = EP // NV, EP % NV; del EP
    nb0 = int(EB2[EA2 == 0][0])
    print(f"  psi #{PSI % len(adm)} = {psi}: {NV} vertices, {len(seen)} directions, {len(EA2)} edges   [{time.time()-t0:.0f}s]", flush=True)
    with open(CNF, "w") as fh:
        fh.write(f"p cnf {NV * K} {NV + len(EA2) * K + 2}\n")
        vv = np.arange(NV, dtype=np.int64) * K
        np.savetxt(fh, np.stack([vv + k + 1 for k in range(K)] + [np.zeros(NV, dtype=np.int64)], axis=1), fmt="%d")
        for k in range(K):
            np.savetxt(fh, np.stack([-(EA2 * K + k + 1), -(EB2 * K + k + 1), np.zeros(len(EA2), dtype=np.int64)], axis=1), fmt="%d")
        fh.write(f"{var(0, 0)} 0\n{var(nb0, 1)} 0\n")
    if os.environ.get("LOCAL"):
        # is the pair query already UNSAT on the torus ball of radius LOCAL around 0 and 5e?
        from pysat.solvers import Solver
        R_ = int(os.environ["LOCAL"])
        nbrs = [[] for _ in range(NV)]
        for a_, b_ in zip(EA2.tolist(), EB2.tolist()): nbrs[a_].append(b_); nbrs[b_].append(a_)
        for i in EUS:
            y5 = int(((5 * uclA[i]) % dm_) @ st_)
            ball = {0, y5}; front = {0, y5}
            for _ in range(R_):
                front = {w for v_ in front for w in nbrs[v_]} - ball; ball |= front
            bs = sorted(ball); bset = set(bs)
            with Solver(name="cadical195") as sv:
                for v_ in bs: sv.add_clause([var(v_, k) for k in range(K)])
                for a_, b_ in zip(EA2.tolist(), EB2.tolist()):
                    if a_ in bset and b_ in bset:
                        for k in range(K): sv.add_clause([-var(a_, k), -var(b_, k)])
                sv.add_clause([var(0, 0)]); sv.add_clause([-var(y5, 0)])
                ok = sv.solve()
            print(f"    unit {i}: radius-{R_} torus ball around 0 and 5e ({len(bs)} vertices): split possible = {ok}", flush=True)
        continue
    if os.environ.get("WRITEONLY"):
        # one CNF per unit with the query clause and an exact header, for proof/core extraction
        for i in EUS:
            y5 = int(((5 * uclA[i]) % dm_) @ st_)
            qp = f"{CNF}.u{i}.cnf"
            with open(qp, "w") as fo, open(CNF) as fi:
                hdr = fi.readline().split(); fo.write(f"p cnf {hdr[2]} {int(hdr[3]) + 1}\n")
                for line in fi: fo.write(line)
                fo.write(f"-{var(y5, 0)} 0\n")
            print(f"    wrote {qp}; vertex of 5e = {y5}; NV = {NV}", flush=True)
        continue
    for i in EUS:
        y5 = int(((5 * uclA[i]) % dm_) @ st_)        # 5e: psi-coordinate 0, alpha-part 5e
        q = subprocess.run(f"cat {CNF} - | {KISSAT} --relaxed -n --time={int(TL)}", shell=True, input=f"-{var(y5, 0)} 0\n", capture_output=True, text=True).stdout
        v = [l for l in q.splitlines() if l.startswith("s ")]
        verdict = v[0] if v else f"no verdict within {int(TL)} s"
        summary.append((PSI % len(adm), i, verdict))
        print(f"    unit {i}: split (0, 5e): {verdict}   [{time.time()-t0:.0f}s]", flush=True)
print("  summary: " + ", ".join(f"{a}/{b}:{'U' if 'UNSAT' in c else ('S' if 'SATIS' in c else '?')}" for a, b, c in summary))
