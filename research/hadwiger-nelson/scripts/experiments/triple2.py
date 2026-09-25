"""Three copies about one centre, ranked by overlap, on a 5-chromatic carrier.

If |s-q|^2 = 1/3 then q and its two images under the 120- and 240-degree
rotations about s form an equilateral triangle of side exactly 1.  So

    U  =  H  u  rot120_s(H)  u  rot240_s(H)

is k-colourable only if three k-colourings of H that agree everywhere the
copies overlap give q three different colours.  Overlap is the whole point:
with the copies nearly disjoint the three colourings are nearly independent
and any three available colours will do, while at high overlap they are almost
the same colouring and cannot spread.

Four copies are not available -- four points pairwise one apart do not exist
in the plane -- so this is the whole of the extra room over the usual pair.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

WHICH, K, TOP = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
if WHICH in ("Z", "H"):
    F = Field((3, 11, 247))
    base = build_Sa(F)
    r1 = rotation_joining(Fr(1), F).about(base[25])
    seen, H = set(), []
    for p in base:
        for q in (p, r1(p)):
            if q not in seen: seen.add(q); H.append(q)
    if WHICH == "H": pts = H
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
n0 = len(pts); S = set(pts)
r120 = rotation_joining(Fr(1, 3), F)
print(f"carrier {WHICH}: {n0} points, k={K}", flush=True)

t0 = time.time()
ranked = []
for si, s in enumerate(pts):
    if not any((p - s).norm2() == Fr(1, 3) for p in pts): continue
    a = r120.about(s)
    seen = set()
    for p in pts:
        q1 = a(p); seen.add(p); seen.add(q1); seen.add(a(q1))
    if len(seen) == n0: continue
    ranked.append((3 * n0 - len(seen), si, len(seen)))
ranked.sort(key=lambda x: -x[0])
print(f"  {len(ranked)} usable centres; best overlaps "
      f"{[r[0] for r in ranked[:8]]}   [{time.time()-t0:.0f}s]", flush=True)

for ov, si, nn in ranked[:TOP]:
    a = r120.about(pts[si])
    seen, U = set(), []
    for p in pts:
        q1 = a(p); q2 = a(q1)
        for z in (p, q1, q2):
            if z not in seen: seen.add(z); U.append(z)
    g = build_graph(U)
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(g.n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf); ok = s.solve(); s.delete()
    m = sum(len(x) for x in g.adj) // 2
    print(f"  v{si}: overlap={ov:<6} n={g.n:<6} m={m:<7} deg={2.0*m/g.n:.2f} "
          f"{K}-colourable={ok}" + ("" if ok else "   <<<<<<<<<<<<<<<<<<<<"),
          flush=True)
    if not ok:
        pickle.dump((WHICH, K, si), open(
            f"/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/TRIPHIT_{WHICH}_{K}_{si}.pkl", "wb"))
        break
