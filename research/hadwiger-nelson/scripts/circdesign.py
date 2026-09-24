"""Design a module with no circular 5-colouring, by cutting planes.

Loop:  (1) MILP for a character phi of M with every unit in [1/5, 4/5] (scripts/circgate.py's model);
       (2) if one exists, look for unit vectors w = rho * u (u a unit of M, rho a rotation of the field)
           that NO extension of phi to M + Zw can colour: w has order N modulo M, phi(Nw) is fixed,
           and all N candidate values (phi(Nw) + j)/N fall outside [1/5, 4/5];
       (3) add the cheapest such w (smallest N), close under the 60-degree turn, repeat.
Stops when the MILP is infeasible (F(M) empty: no coset and no circular colouring) or when no
candidate cuts phi.  The module stays of rank 8 (every rotation lies in the field); it only gets finer.

Note: phi is expressed in this script's own basis (echelon of ALL units after the 60-degree
closure), which is not circgate.py's basis (echelon of the listed units only).  A single new unit of
order N >= 2 modulo M can never cut phi (two candidate values 1/N apart cannot both miss a window of
width 3/5); cuts come from hidden units (order 1) or from several new units added jointly.

usage: python3 circdesign.py <start.json> <rounds> <milp seconds> [out.json]"""
import sys, json, os, time, itertools
from fractions import Fraction as Fr
from math import gcd, floor
import numpy as np
import sympy
from scipy.optimize import milp, LinearConstraint, Bounds
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
t0 = time.time()
path = sys.argv[1]; ROUNDS = int(sys.argv[2]); TL = float(sys.argv[3]); OUT = sys.argv[4] if len(sys.argv) > 4 else None
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(p, q) for p, q in xy[0]]), F.element([Fr(p, q) for p, q in xy[1]]))
if "units" in d: U = [mk(xy) for xy in d["units"]]
else:
    P = [mk(xy) for xy in d["points"]]; g = build_graph(P); Ud = {}
    for i, j in g.edges():
        for w in (P[j] - P[i], P[i] - P[j]): Ud.setdefault((round(w.fx, 9), round(w.fy, 9)), w)
    U = list(Ud.values())
def fe(**kw):
    c = [Fr(0)] * F.dim
    for m in range(F.dim):
        if F._prod[m] in kw.get("rad", {}): c[m] = Fr(kw["rad"][F._prod[m]])
    return F.element(c)
def rot(cos_c, sin_rad):                      # rotation with cos = rational, sin = a * sqrt(q)
    cth = F.rational(cos_c); sth = fe(rad={q: a for q, a in sin_rad.items()})
    return lambda p: Point(p.x * cth - p.y * sth, p.x * sth + p.y * cth)
ROT = {
    "omega": rot(Fr(1, 2), {3: Fr(1, 2)}),
    "rho7": rot(Fr(1, 7), {3: Fr(4, 7)}), "rho7^-1": rot(Fr(1, 7), {3: Fr(-4, 7)}),
    "kappa": rot(Fr(-11, 14), {3: Fr(5, 14)}), "kappa^-1": rot(Fr(-11, 14), {3: Fr(-5, 14)}),
    "sigma": rot(Fr(5, 6), {11: Fr(1, 6)}), "sigma^-1": rot(Fr(5, 6), {11: Fr(-1, 6)}),
    "lambda": rot(Fr(49, 50), {11: Fr(3, 50)}), "lambda^-1": rot(Fr(49, 50), {11: Fr(-3, 50)}),
}
if 247 in F._prod:
    ROT["tau247"] = rot(Fr(-3, 16), {247: Fr(1, 16)}); ROT["tau247^-1"] = rot(Fr(-3, 16), {247: Fr(-1, 16)})
key = lambda p: (round(p.fx, 7), round(p.fy, 7))
def close60(U):
    seen = {key(u) for u in U}; out = list(U)
    for u in list(U):
        w = u
        for _ in range(5):
            w = ROT["omega"](w)
            if key(w) not in seen: seen.add(key(w)); out.append(w)
    return out
src = open(os.path.join(HERE, "allunits.py")).read()
exec(src[src.index("def lll_exact"):src.index("Gr, Ur = lll_exact(G)")])
wt = list(F._prod) * 2
def setup(U):
    Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
    den = 1
    for v in Ev:
        for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
    B = echelon(E); r = len(B)
    G0 = [[Fr(sum(w_ * p * q for w_, p, q in zip(wt, B[i], B[j]))) for j in range(r)] for i in range(r)]
    _, T = lll_exact(G0); Ti = sympy.Matrix(T).inv()
    Bm = sympy.Matrix([[Fr(x) for x in b] for b in B])            # r x 2D, rows = echelon basis
    piv = list(Bm.T.rref()[1])[:r] if False else None
    # choose r independent columns once
    cols = []
    for j in range(Bm.shape[1]):
        if sympy.Matrix.hstack(*([Bm[:, k] for k in cols] + [Bm[:, j]])).rank() > len(cols): cols.append(j)
        if len(cols) == r: break
    Bsq = Bm[:, cols]; Bsq_inv = Bsq.inv()
    def lc(p):                                  # rational LLL coordinates of any vector of M (x) Q
        v = [Fr(x) * den for x in (tuple(p.x.c) + tuple(p.y.c))]
        c = sympy.Matrix([[v[j] for j in cols]]) * Bsq_inv          # c * Bm = v on the pivot columns
        if list(c * Bm) != [sympy.Rational(x.numerator, x.denominator) for x in v]: return None
        return [Fr(int(x.p), int(x.q)) for x in (c * Ti)]
    CU = [lc(u) for u in U]
    assert all(x.denominator == 1 for c in CU for x in c)
    return lc, np.array([[int(x) for x in c] for c in CU], dtype=np.int64), r
def solve(CU, r):
    reps, seenr = [], set()
    for i in range(len(CU)):
        k_ = tuple(CU[i]); nk_ = tuple(-CU[i])
        if nk_ in seenr or k_ in seenr: continue
        seenr.add(k_); reps.append(i)
    C = CU[reps].astype(float); m = len(reps)
    lo_n = np.floor(np.minimum(C, 0).sum(axis=1)) - 1; hi_n = np.ceil(np.maximum(C, 0).sum(axis=1)) + 1
    nv = r + m + 1; cobj = np.zeros(nv); cobj[-1] = -1.0
    A1 = np.zeros((m, nv)); A1[:, :r] = C; A1[np.arange(m), r + np.arange(m)] = -1.0; A1[:, -1] = -1.0
    A2 = np.zeros((m, nv)); A2[:, :r] = C; A2[np.arange(m), r + np.arange(m)] = -1.0; A2[:, -1] = 1.0
    cons = [LinearConstraint(A1, lb=0.2, ub=np.inf), LinearConstraint(A2, lb=-np.inf, ub=0.8)]
    lb = np.concatenate([np.zeros(r), lo_n, [0.0]]); ub = np.concatenate([np.ones(r), hi_n, [0.5]])
    integ = np.concatenate([np.zeros(r), np.ones(m), [0]])
    res = milp(cobj, constraints=cons, integrality=integ, bounds=Bounds(lb, ub), options={"time_limit": TL})
    return res, m
fr = lambda q: q - floor(q)
def slack_of(v): return min(v - Fr(1, 5), Fr(4, 5) - v)
U = close60(U)
history = []
for rd in range(ROUNDS):
    lc, CU, r = setup(U)
    res, m = solve(CU, r)
    if res.status == 2 or res.x is None:
        print(f"round {rd}: {len(U)} units ({m} directions): MILP {res.message}   [{time.time()-t0:.0f}s]", flush=True)
        if res.status == 2: print("  ==> NO circular 5-colouring (F(M) empty)", flush=True)
        break
    phi = [Fr(x).limit_denominator(10 ** 6) for x in res.x[:r]]
    sl = min(slack_of(fr(sum(Fr(int(c)) * p for c, p in zip(cu, phi)))) for cu in CU)
    print(f"round {rd}: {len(U)} units ({m} directions), rank {r}: phi found, slack {float(res.x[-1]):.5f} (exact after rounding {float(sl):.5f})   [{time.time()-t0:.0f}s]", flush=True)
    if sl < 0:
        print("  rounded phi lost properness; stopping"); break
    # candidate cuts
    seen = {key(u) for u in U}; cands = {}
    for name, R_ in ROT.items():
        if name == "omega": continue
        for u in U:
            w = R_(u)
            if key(w) in seen or key(w) in cands: continue
            c = lc(w)
            if c is None: continue
            N = 1
            for x in c: N = N * x.denominator // gcd(N, x.denominator)
            base = sum(N * x * p for x, p in zip(c, phi))          # phi(N w), a fixed value mod 1
            best = max(slack_of(fr((base + j) / N)) for j in range(N))
            if best < 0: cands[key(w)] = (N, name, w, float(best))
    if not cands:
        print("  no rotation of a unit cuts this phi; stopping"); break
    Nmin = min(v[0] for v in cands.values())
    cut = [v for v in cands.values() if v[0] == Nmin]
    print(f"  {len(cands)} cutting units; adding the {len(cut)} of least order N = {Nmin} (e.g. {cut[0][1]})", flush=True)
    history.append({"round": rd, "units": len(U), "phi": [str(x) for x in phi], "cut_order": Nmin, "added": len(cut)})
    U = close60(U + [v[2] for v in cut])
if OUT:
    json.dump({"field_generators": list(F.gens), "units": [[[[x.numerator, x.denominator] for x in u.x.c], [[x.numerator, x.denominator] for x in u.y.c]] for u in U],
               "points": [[[[0, 1]] * F.dim, [[0, 1]] * F.dim]], "A": 0, "B": 0, "history": history}, open(OUT, "w"))
    print(f"  wrote {OUT}")
