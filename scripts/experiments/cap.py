"""de Grey's actual mechanism, tested on the tightest graph: the palette CAP.

This session measured mu -- the MINIMUM number of colours a set is forced to
show -- which is the blocking direction.  What made four go to five was the
opposite quantity.  The project identified it earlier, exactly: de Grey's Y
forces its ring's antipodal pair because, after the bite, the six points of Sa's
radius-2 ring can show AT MOST TWO colours in any 4-colouring.  A cap, not a
floor:

    cap(S) = max over proper colourings of |c(S)|

and a ring whose cap is small is nearly monochromatic in every colouring, which
is where forced pairs come from.  On de Grey's own G at five colours the census
found no capped ring at any distance other than 1 -- and a unit ring's cap of 4
is only the statement that its centre gets a colour.

The tight hexagon graph is new, and a cap is exactly what tightness should
produce: in a uniquely colourable graph every set has one colour set, so every
cap falls to the number of colours it actually shows.

Screen cheaply first: Kempe swaps from one colouring generate thousands of
colourings for free, and any ring seen showing all five in one of them is not
capped.  Survivors get the real question, one SAT call each on a warm solver:
five clauses, "colour c appears somewhere on the ring", guarded by a selector.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, math, random
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "tight_hexagon_4159.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n; E = list(G.edges())
adj = [[] for _ in range(n)]
for a, b in E: adj[a].append(b); adj[b].append(a)
adjs = [set(a) for a in adj]
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K): cnf.append([-X(v, a), -X(v, b)])
for a, b in E:
    for c in range(K): cnf.append([-X(a, c), -X(b, c)])
s = Solver(name="cd19", bootstrap_with=cnf)
assert s.solve()
m = s.get_model()
col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
print(f"{NAME}: n={n}; first colouring   [{time.time()-t0:.0f}s]", flush=True)

# rings: for every vertex as centre, group the other vertices by exact squared
# distance; keep classes of six or more that are not the unit ring
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 3), int(hy[i] // 3))].append(i)
rings = []
for c0 in range(n):
    cx, cy = int(hx[c0] // 3), int(hy[c0] // 3)
    byd = defaultdict(list)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j == c0: continue
                dd = (hx[j]-hx[c0])**2 + (hy[j]-hy[c0])**2
                if dd < 9.0 and abs(dd - 1.0) > 1e-9:
                    byd[round(dd, 7)].append(j)
    for dd, pts in byd.items():
        if len(pts) >= 6: rings.append((c0, dd, pts))
print(f"  {len(rings)} non-unit rings of six or more points, radius < 3   "
      f"[{time.time()-t0:.0f}s]", flush=True)

random.seed(3)
alive = [r for r in rings if len({col[v] for v in r[2]}) < K]
def kempe(col):
    v = random.randrange(n); a = col[v]
    b = random.choice([c for c in range(K) if c != a])
    comp = {v}; q = deque([v])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in comp and col[w] in (a, b): comp.add(w); q.append(w)
    for u in comp: col[u] = b if col[u] == a else a
for t in range(1, 30001):
    kempe(col)
    if t % 500 == 0:
        alive = [r for r in alive if len({col[v] for v in r[2]}) < K]
    if t % 10000 == 0:
        print(f"    {t} Kempe swaps: {len(alive)} rings never seen with all five"
              f"   [{time.time()-t0:.0f}s]", flush=True)
alive = [r for r in alive if len({col[v] for v in r[2]}) < K]
print(f"  Kempe-screened: {len(alive)} of {len(rings)} rings still candidates   "
      f"[{time.time()-t0:.0f}s]", flush=True)

# Largest rings first: a six-point ring can fail to show five colours by
# accident in every colouring a walk happens to visit, but a ring of twelve or
# more that never shows five is genuinely squeezed -- and every new colouring
# the solver returns is used at once to discard every ring it shows uncapped.
alive.sort(key=lambda r: -len(r[2]))
MINPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
alive = [r for r in alive if len(r[2]) >= MINPTS]
print(f"  of those, {len(alive)} have {MINPTS} or more points", flush=True)
nsel = n * K + 1; capped = []; calls = 0
while alive:
    c0, dd, pts = alive[0]
    sel = nsel; nsel += 1
    for c in range(K):
        s.add_clause([-sel] + [X(v, c) for v in pts])
    r = s.solve(assumptions=[sel])
    c2 = None
    if r:
        mm = s.get_model()
        c2 = [next(c for c in range(K) if mm[X(v, c) - 1] > 0) for v in range(n)]
    s.add_clause([-sel])
    if not r:
        capped.append((c0, dd, len(pts)))
        print(f"  *** CAPPED: centre {c0}, d^2 = {dd:.6f}, {len(pts)} points "
              f"never show all five ***   [{time.time()-t0:.0f}s]", flush=True)
        json.dump({"graph": NAME, "capped": capped},
                  open(f"{ROOT}/data/tight_caps.json", "w"))
        alive = alive[1:]
    else:
        for _ in range(200): kempe(c2)
        alive = [q for q in alive[1:] if len({c2[v] for v in q[2]}) < K]
    calls += 1
    if calls % 10 == 0:
        print(f"    {calls} solver calls: {len(alive)} rings left   "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\n  capped non-unit rings: {len(capped)}   [{time.time()-t0:.0f}s]",
      flush=True)
