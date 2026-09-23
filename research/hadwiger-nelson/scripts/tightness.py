"""A robust tightness measure, and a test of the one lever found so far.

Solver time is a noisy shadow of how few 5-colourings a graph has: the same
4159-point hexagon graph took 452 s on one run and 53 s on another.  So measure
CONFLICTS instead of seconds, over several solves each perturbed by pinning a
random vertex to a random colour -- different starting points for the search,
same question -- and report the median.  Colour symmetry makes every such pin
harmless: some 5-colouring puts any vertex in any colour.

The hypothesis under test is the one lever this session found.  Six copies of
the vertex-critical 803 on a hexagon, translated, were about as tight as the 803;
the same six copies TURNED by multiples of 60 degrees were an order of magnitude
tighter.  If orientation diversity is what matters, then turning by an angle with
no alignment to the lattice at all should be tighter still -- and the Moser angle,
cos = 5/6 and sin = sqrt11/6, lies in the field already.

Built the same way each time: the densest vertex of each copy placed on one
hexagon vertex, the copy turned about that vertex by i times the chosen angle,
and the unit circle about the centre pruned to the hexagon.
"""
import sys, time, json, random, statistics
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(P)
adj0 = defaultdict(set)
for x, y in g0.edges(): adj0[x].add(y); adj0[y].add(x)
v0 = max(range(g0.n), key=lambda v: len(adj0[v]))
r3, r11 = F.sqrt(3), F.sqrt(11)
half, sixth = F.rational(Fr(1, 2)), F.rational(Fr(1, 6))
ANGLES = {"translate": (F.rational(1), F.rational(0)),
          "60deg": (half, r3 * half),
          "moser": (F.rational(Fr(5, 6)), r11 * sixth)}

def hexgraph(mode):
    c, s = ANGLES[mode]
    O = Point(F.rational(0), F.rational(0))
    hexa = []; x, y = F.rational(1), F.rational(0)
    for _ in range(6):
        hexa.append(Point(x, y)); x, y = x * half - y * r3 * half, x * r3 * half + y * half
    pts = {(0.0, 0.0): O}
    for hv in hexa: pts[(round(hv.fx, 9), round(hv.fy, 9))] = hv
    bx, by = g0.vertices[v0].x, g0.vertices[v0].y
    for i, hv in enumerate(hexa):
        cc, ss = F.rational(1), F.rational(0)
        for _ in range(i): cc, ss = cc * c - ss * s, cc * s + ss * c
        for p in P:
            dx, dy = p.x - bx, p.y - by
            r = Point(hv.x + dx * cc - dy * ss, hv.y + dx * ss + dy * cc)
            pts[(round(r.fx, 9), round(r.fy, 9))] = r
    hk = {(round(hv.fx, 9), round(hv.fy, 9)) for hv in hexa}
    keep = [p for k, p in pts.items()
            if not (abs(p.fx*p.fx + p.fy*p.fy - 1.0) < 1e-9 and k not in hk)]
    return build_graph(keep)

def tightness(G, runs=3, seed=0):
    n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K): cnf.append([-X(a, c), -X(b, c)])
    rng = random.Random(seed); out = []
    for _ in range(runs):
        s = Solver(name="cd19", bootstrap_with=cnf)
        v = rng.randrange(n); c = rng.randrange(K)
        t1 = time.time(); ok = s.solve(assumptions=[X(v, c)])
        st = s.accum_stats(); s.delete()
        out.append((st.get("conflicts", 0), time.time() - t1, ok))
    return out

which = sys.argv[1:] or ["base", "translate", "60deg", "moser"]
for mode in which:
    G = g0 if mode == "base" else hexgraph(mode)
    res = tightness(G)
    conf = [r[0] for r in res]; secs = [r[1] for r in res]
    print(f"  {mode:<10s} n={G.n:<5d} edges={sum(1 for _ in G.edges()):<6d} "
          f"conflicts {conf}  median {statistics.median(conf):>9.0f}  "
          f"seconds median {statistics.median(secs):6.1f}  "
          f"colourable {[r[2] for r in res]}   [{time.time()-t0:.0f}s]", flush=True)
