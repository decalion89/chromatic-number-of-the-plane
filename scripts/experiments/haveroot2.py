"""Does sigma manufacture sqrt2 -- and with it, square centres, natively?

The square closure wants a vertex h that is the midpoint of two graph points
sqrt2 apart.  Triangular-lattice carriers have no such pair at all, since 2 is
not Loeschian, and translating a copy by a vector of that length couples badly:
three translates of the 803-graph share only 169 edges between 2409 vertices,
whichever in-field offset of length 1/sqrt2 is used, because a translate by a
non-lattice vector lands nowhere near the lattice.

But sigma is not a translation.  It multiplies by (1 - sqrt11 i)/6, whose norm
is 1/3, so the points it makes carry sqrt11 in a way the carrier never did --
and the squared distances it creates need not be Loeschian at all.  So ask
directly, of the graph it builds: are there pairs at d^2 = 2, and at d^2 = 1/2,
and does any pair at 2 have its own midpoint present?
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = HN_DIR
t0 = time.time()
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 1
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
r11 = F.sqrt(11); sixth = F.rational(Fr(1, 6))
def sig(p, c, conj):
    s = -r11 if conj else r11
    dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + (dx - s * dy) * sixth, c.y + (s * dx + dy) * sixth)
pts = {}
for x, y in d["points"]:
    q = Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
    pts[(round(q.fx, 9), round(q.fy, 9))] = q
for _ in range(ROUNDS):
    G0 = build_graph(list(pts.values()))
    for c, u in G0.edges():
        for (a, b) in ((c, u), (u, c)):
            for conj in (False, True):
                p = sig(G0.vertices[b], G0.vertices[a], conj)
                k = (round(p.fx, 9), round(p.fy, 9))
                if k not in pts: pts[k] = p
G = build_graph(list(pts.values())); n = G.n
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
key = {(round(hx[i], 9), round(hy[i], 9)): i for i in range(n)}
print(f"sigma^{ROUNDS} of {NAME}: n={n}   [{time.time()-t0:.0f}s]", flush=True)
two = F.rational(2); half = F.rational(Fr(1, 2))
at2 = []; at_half = 0
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if abs(dd - 2.0) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == two:
            at2.append((i, j))
        elif abs(dd - 0.5) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == half:
            at_half += 1
print(f"  pairs at d^2 = 2   : {len(at2)}", flush=True)
print(f"  pairs at d^2 = 1/2 : {at_half}   [{time.time()-t0:.0f}s]", flush=True)
squares = []
for i, j in at2:
    mx = (hx[i] + hx[j]) / 2; my = (hy[i] + hy[j]) / 2
    h = key.get((round(mx, 9), round(my, 9)))
    if h is not None: squares.append((h, i, j))
print(f"  square centres present (midpoint in the graph): {len(squares)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
if squares:
    adj = defaultdict(set)
    for x, y in G.edges():
        adj[x].add(y); adj[y].add(x)
    squares.sort(key=lambda t: -(len(adj[t[0]]) + len(adj[t[1]]) + len(adj[t[2]])))
    print(f"    richest: hub degree {len(adj[squares[0][0]])}, "
          f"diagonal degrees {len(adj[squares[0][1]])},{len(adj[squares[0][2]])}",
          flush=True)
    json.dump({"source": NAME, "rounds": ROUNDS,
               "squares": [[a, b, c] for a, b, c in squares[:500]]},
              open(f"{ROOT}/data/squares_{NAME}", "w"))
