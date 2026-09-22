"""Three copies about one centre, which needs only a palette of two.

Everyone glues two copies, because the spindle is a two-copy device.  Three is
available and cheaper in what it demands.

If |s - q|^2 = 1/3 then q, rot120_s(q) and rot240_s(q) are an equilateral
triangle of side exactly 1:  2 * (1/sqrt3) * sin(60 deg) = 1.  So in

    U  =  H  u  rot120_s(H)  u  rot240_s(H)

the three images of every such q are mutually adjacent, and U is k-colourable
only if three k-colourings of H, agreeing wherever the copies overlap, give q
three DIFFERENT colours.  Refusing that needs no forced-equal pair -- it is
enough that q has at most two colours available once the pivot's is fixed.
A cap of two where the two-copy route needs a cap of one.

Four copies cannot be used: four points pairwise one apart do not exist in the
plane.  So three is the whole of the extra room, and rot120 costs sqrt3, which
every carrier here already has.
"""
import sys, time, random, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

WHICH, K = sys.argv[1], int(sys.argv[2])
LO = int(sys.argv[3]); HI = int(sys.argv[4])

if WHICH in ("Z", "H"):
    F = Field((3, 11, 247))
    base = build_Sa(F)
    r1 = rotation_joining(Fr(1), F).about(base[25])
    seen, H = set(), []
    for p in base:
        for q in (p, r1(p)):
            if q not in seen: seen.add(q); H.append(q)
    if WHICH == "H":
        pts = H
    else:
        r2 = rotation_joining(Fr(64, 9), F).about(H[157])
        seen, pts = set(), []
        for p in H:
            for q in (p, r2(p)):
                if q not in seen: seen.add(q); pts.append(q)
else:
    F = Field((3, 5, 7, 11))
    pts = {"Sa": build_Sa, "Y": build_Y,
           "G": lambda f: build_G(f, as_graph=False)}[WHICH](F)
n0 = len(pts)
S = set(pts)
half = F.rational(Fr(-1, 2))
r120 = rotation_joining(Fr(1, 3), F)          # cos -1/2, sin sqrt3/2
print(f"carrier {WHICH}: {n0} points, k={K}, centres [{LO},{HI})", flush=True)

t0 = time.time()
best = (0, None)
for si in range(LO, min(HI, n0)):
    s = pts[si]
    # the triangle only exists if some point sits at squared distance 1/3
    ring = sum(1 for p in pts if (p - s).norm2() == Fr(1, 3))
    if ring == 0:
        continue
    a = r120.about(s)
    b = a.about if False else None
    seen, U = set(), []
    for p in pts:
        q1 = a(p); q2 = a(q1)
        for z in (p, q1, q2):
            if z not in seen: seen.add(z); U.append(z)
    g = build_graph(U)
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    sv = Solver(name="m22", bootstrap_with=cnf)
    ok = sv.solve()
    sv.delete()
    ov = 3 * n0 - n
    if not ok:
        print(f"  !!! v{si}: ring={ring} n={n} overlap={ov} "
              f"NOT {K}-COLOURABLE  <<<<<<<<<<<<<<<<<<", flush=True)
        pickle.dump((WHICH, K, si), open(
            f"/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/TRIP_{WHICH}_{K}_{si}.pkl", "wb"))
        continue
    if ov > best[0]:
        best = (ov, si)
    if si % 20 == 0:
        print(f"  v{si}: ring={ring} n={n} overlap={ov} {K}-colourable "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"done; best overlap {best[0]} at v{best[1]}  [{time.time()-t0:.0f}s]", flush=True)
