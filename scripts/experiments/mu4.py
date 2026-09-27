"""mu at FOUR colours, where the structure exists -- calibration, and a long shot.

Every mu measured at five colours in this project is 2, its floor.  The
instrument reaches its ceiling on the Moser spindle, so the floor is real and
not an artefact.  What has never been measured is mu at FOUR, on the carriers
that are rigid at four -- and that is where the interesting distribution
should live, because free@4 = 0.00 % on Sa means its vertices really are
pinned, which is the local shadow of what a blocked point needs.

Two reasons to run it.

CALIBRATION.  If mu_4 reaches 3 often on Sa and never 4, that is what "close"
looks like one level down, and it prices the gap at five honestly.  If it never
even reaches 3, the five-colour floor of 2 is not a five-colour problem at all.

A LONG SHOT.  If mu_4 = 4 anywhere on Sa, then Sa + p has no 4-colouring: a
5-chromatic unit-distance graph on 398 vertices, against a record of 509.  That
is almost certainly not there -- de Grey would have found it, and his
construction needed 1581 points precisely because no single point does it --
but the check costs seconds and the payoff is not small.

Candidates are generated in-field, as images of the vertex set under rotations
the field already holds, and scored with a float filter in front of the exact
comparison.  mu itself comes from hn.blocked.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa, build_Y
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from hn.blocked import MuSolver

t0 = time.time()
F = Field((3, 5, 7, 11))
rot60 = _rot60(F)
Sa = build_Sa(F)

def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for r in [Point(z.x, -z.y) for z in list(out)]:
        if r not in out: out.append(r)
    return out

def grow(b):
    seen, U = set(Sa), list(Sa)
    for w in orb(Sa[b]):
        rot = rotation_joining(Fr(1), F).about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    return U

CASES = [("Sa", list(Sa)), ("Y", list(build_Y(F))), ("Sa[25]", grow(25))]
half = F.rational(Fr(1, 2))
r30 = Rotation(F.sqrt(3) * half, half)
one = F.rational(Fr(1))

for name, U in CASES:
    g = build_graph(U); n = g.n
    m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {name}  n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    S = set(g.vertices)
    deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
    cand = set()
    for c in deg_order[:min(60, n)]:
        for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
            for p in g.vertices:
                z = rot(p)
                if z not in S:
                    cand.add(z)
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
    gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
    scored, seen = [], set()
    for z in cand:
        zx, zy = z.fx, z.fy
        cx, cy = int(zx // 1), int(zy // 1)
        nb = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in cells.get((cx + dx, cy + dy), ()):
                    ex, ey = gx[i] - zx, gy[i] - zy
                    if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                       (g.vertices[i] - z).norm2() == one:
                        nb.append(i)
        if len(nb) >= 4:
            t = tuple(sorted(nb))
            if t not in seen:
                seen.add(t); scored.append((len(nb), t, z))
    scored.sort(key=lambda u: -u[0])
    print(f"    {len(cand)} candidates, {len(scored)} distinct neighbourhoods "
          f"of size >= 4, largest {scored[0][0] if scored else 0}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not scored:
        continue
    hist = Counter()
    with MuSolver(g, k=4, budget=4_000_000) as ms:
        if not ms.colourable:
            print("    (not 4-colourable -- skipped)", flush=True); continue
        for k, nb, z in scored[:600]:
            v = ms.mu(nb)
            hist[v] += 1
            if v is not None and v >= 4:
                print(f"    *** BLOCKED AT FOUR, |N|={k}: {name} + p is "
                      f"5-chromatic on {n+1} vertices ***", flush=True)
                print(f"        p = ({z.x}, {z.y})", flush=True)
            elif v == 3:
                print(f"    |N|={k}: mu_4 = 3  (one short)", flush=True)
    print(f"    mu_4 over the richest {min(600, len(scored))}: {dict(hist)}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
