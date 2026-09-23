"""Which angular separations actually occur inside a candidate neighbourhood?

The rotation stack about p is the one construction whose hypothesis survives
the rotation: if a 5-colouring of G' = union of rho^i G has |c(N_{G'}(p))| <= 2
then each copy's restriction is an ESCAPE colouring of that copy, because
N_{rho^i G}(p) = rho^i(N_G(p)) is a subset of N_{G'}(p).  So every copy
contributes its own escape conclusion, and the conclusions chain.

A conclusion "u, v monochromatic, both in N(p)" is an angular relation: u and v
lie on the unit circle about p at separation theta, and the stack links a point
to the point theta further round.  Starting from relations at separations S the
colour classes are cosets of <S>, and the contradiction -- two monochromatic
points a unit apart -- appears exactly when 60 degrees lies in <S>.

Which separations are available?  A separation inside <rho> is a rational
multiple of pi, and a monochromatic pair needs a rational squared distance in
the relation, so by Niven's theorem theta is 60, 90, 120 or 180 degrees:

    60 deg   d = 1       ADJACENT -- never monochromatic, so unusable
    90 deg   d = sqrt2   <90> = {0,90,180,270}, no 60          on its own: no
   120 deg   d = sqrt3   <120> = {0,120,240}, no 60            on its own: no
   180 deg   d = 2       <180> = {0,180}, no 60                on its own: no

    <120, 180> = <60>    CLOSES
    <90, 120>  = <30>    CLOSES
    <90, 180>  = <90>    no

The escape analysis already returns a 120-degree relation, uniformly, at 25 of
the 30 richest candidate points of the 803-graph.  So the entire remaining
requirement is one more relation, at 90 or at 180 degrees:

    every escape colouring must also have a monochromatic pair at
    distance sqrt2, or a monochromatic ANTIPODAL pair through p.

That is a statement about escape colourings, not an unconditional forced pair,
so it does not meet the wall the certificate pricing just found.  The first
question is whether the separations are even present: a relation at 180 degrees
needs N(p) to contain an antipodal pair in the first place.
"""
import sys, time, json, math
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
for name in ("five_247_c.json", "five_247.json", "five_tuned_1_1.json",
             "five_dense_2.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    cells = defaultdict(list)
    for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
    # candidate points: intersections of two unit circles, collected by float key
    cand = defaultdict(int)
    for i in range(n):
        cx, cy = int(hx[i] // 2), int(hy[i] // 2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cells.get((cx + dx, cy + dy), ()):
                    if j <= i: continue
                    ex, ey = hx[j] - hx[i], hy[j] - hy[i]
                    dd = ex * ex + ey * ey
                    if dd <= 1e-12 or dd >= 4.0: continue
                    h = math.sqrt(max(0.0, 1.0 - dd / 4.0)); r = math.sqrt(dd)
                    mx, my = (hx[i] + hx[j]) / 2, (hy[i] + hy[j]) / 2
                    ux, uy = -ey / r, ex / r
                    for s in (+1, -1):
                        cand[(round(mx + s * h * ux, 9), round(my + s * h * uy, 9))] += 1
    rich = sorted(cand.items(), key=lambda kv: -kv[1])[:4000]
    tally = Counter(); withsep = Counter(); best = None
    for (kx, ky), _ in rich:
        nb = [i for i in range(n) if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
        if len(nb) < 8: continue
        seps = Counter()
        for u, v in combinations(nb, 2):
            dd = (hx[u]-hx[v])**2 + (hy[u]-hy[v])**2
            for lbl, val in (("60", 1.0), ("90", 2.0), ("120", 3.0), ("180", 4.0)):
                if abs(dd - val) < 1e-9: seps[lbl] += 1; break
            else: seps["other"] += 1
        tally[len(nb)] += 1
        for lbl in ("90", "120", "180"):
            if seps[lbl]: withsep[lbl] += 1
        if seps["180"] or seps["90"]:
            if best is None or len(nb) > best[0]:
                best = (len(nb), kx, ky, dict(seps))
    tot = sum(tally.values())
    print(f"\n{name}  n={n}: {tot} candidate points with |N| >= 8   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"   with a  90 deg pair (sqrt2)    : {withsep['90']}", flush=True)
    print(f"   with a 120 deg pair (sqrt3)    : {withsep['120']}", flush=True)
    print(f"   with a 180 deg pair (antipodal): {withsep['180']}", flush=True)
    if best: print(f"   richest such point: |N|={best[0]} seps={best[3]}", flush=True)
