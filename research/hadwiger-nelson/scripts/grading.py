"""The 60-degree grading: why a rotation stack about the target can never close.

Inside N(p) every edge is a 60-degree step, so along any path the angle moves
by +-60 per step and the parity of the path length equals the parity of the
angle index (angle difference divided by 60).  Two consequences, and they are
not about any particular graph:

  * within a connected component of N(p), a 2-colouring IS the angle parity.
    Separation 120 deg (index 2, even) is therefore monochromatic automatically;
    separation 60 (index 1) and 180 (index 3) are forced APART automatically.
  * across components nothing is related at all, and the 2-colourings of N(p)
    are exactly the 10 * 2^c patterns: a pair of colours, and a parity offset
    per component.

So the only relation the escape hypothesis can ever hand a rotation stack is a
120-degree one -- which is exactly what the escape analysis returns, uniformly,
at 25 of the 30 richest candidate points.  And <120> = {0, 120, 240} does not
contain 60, so the stack's colour classes are the two inscribed triangles and
no two of them are ever a unit apart.  The stack reproduces the hypothesis
instead of contradicting it.

The one opening the argument leaves is an antipodal pair whose two points lie
in DIFFERENT components: there the parities are independent, so the pair can be
monochromatic, and <120, 180> = <60> would close.  That is not a forcing, only
a freedom -- but it is the thing to count, so count it.

This script checks the grading claim rather than assuming it, and counts the
opening.
"""
import sys, time, json, math
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict, Counter, deque
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
for name in ("five_247_c.json", "five_tuned_1_1.json", "five_dense_2.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    cells = defaultdict(list)
    for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
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
    rich = sorted(cand.items(), key=lambda kv: -kv[1])[:3000]
    viol = 0; pts = 0; comps = Counter()
    anti_same = 0; anti_split = 0; split_points = []
    for (kx, ky), _ in rich:
        nb = [i for i in range(n) if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
        if len(nb) < 8: continue
        pts += 1
        S = set(nb)
        nadj = {u: [v for v in nb if v != u and
                    abs((hx[u]-hx[v])**2 + (hy[u]-hy[v])**2 - 1.0) < 1e-9] for u in nb}
        # bipartition by BFS, and angle index relative to the component root
        part = {}; root = {}; ang = {q: math.atan2(hy[q]-ky, hx[q]-kx) for q in nb}
        cid = 0
        for u in nb:
            if u in part: continue
            part[u] = 0; root[u] = cid; q = deque([u])
            while q:
                w = q.popleft()
                for z in nadj[w]:
                    if z not in part:
                        part[z] = 1 - part[w]; root[z] = cid; q.append(z)
            cid += 1
        comps[cid] += 1
        # the grading claim: same component => parity of colour == parity of index
        for u, v in combinations(nb, 2):
            if root[u] != root[v]: continue
            da = (ang[v] - ang[u]) % (2 * math.pi)
            idx = da / (math.pi / 3)
            if abs(idx - round(idx)) > 1e-6:
                viol += 1; continue                    # not a multiple of 60 deg
            if (round(idx) % 2) != ((part[u] + part[v]) % 2):
                viol += 1
        for u, v in combinations(nb, 2):
            dd = (hx[u]-hx[v])**2 + (hy[u]-hy[v])**2
            if abs(dd - 4.0) < 1e-9:
                if root[u] == root[v]: anti_same += 1
                else:
                    anti_split += 1
                    if len(split_points) < 3: split_points.append((kx, ky, len(nb)))
    print(f"\n{name}  n={n}: {pts} candidate points with |N| >= 8   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"   grading violations (index parity != bipartition): {viol}", flush=True)
    print(f"   components per neighbourhood: {dict(sorted(comps.items()))}", flush=True)
    print(f"   antipodal pairs inside one component (forced APART): {anti_same}", flush=True)
    print(f"   antipodal pairs across components (FREE, the opening): {anti_split}",
          flush=True)
    for s in split_points: print(f"      e.g. |N|={s[2]} at ({s[0]:.6f},{s[1]:.6f})", flush=True)
