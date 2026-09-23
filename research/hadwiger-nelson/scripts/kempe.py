"""Kempe chains: thousands of 5-colourings for the price of one solve.

The forced-pair scan on the tight graph stalled for the reason it was built to
exploit: every new colouring of a graph that close to the edge costs the solver
minutes, and eliminating 782 776 candidates needs dozens of them.

But one colouring already contains the others.  Pick two colours a and b, take a
connected component of the subgraph those two colours induce, and swap a with b
inside it: every edge inside the component still joins a to b, every edge
leaving it goes to a third colour, so the result is again a proper 5-colouring.
That is a Kempe swap, it costs one breadth-first search, and it splits exactly
the monochromatic pairs that straddle the component's boundary.

So walk the Kempe graph: thousands of random swaps from the first colouring,
discarding every candidate that any visited colouring splits.  A pair that
survives is monochromatic across the whole Kempe class of that colouring --
which is not yet a proof, since other classes may exist, so every survivor then
gets the real test: one assumption call asking for a colouring that splits it.
Correctness never rests on the walk; the walk only decides what is worth asking.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "tight_hexagon_4159.json"
SWAPS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n
E = list(G.edges())
adj = [[] for _ in range(n)]
for a, b in E: adj[a].append(b); adj[b].append(a)
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
print(f"{NAME}: n={n} edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            cnf.append([-X(v, a), -X(v, b)])
for a, b in E:
    for c in range(K):
        cnf.append([-X(a, c), -X(b, c)])
s = Solver(name="cd19", bootstrap_with=cnf)
assert s.solve(), "not 5-colourable"
def colouring():
    m = s.get_model()
    return [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
col = colouring()
print(f"  first colouring   [{time.time()-t0:.0f}s]", flush=True)
assert all(col[a] != col[b] for a, b in E)
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
adjset = [set(a) for a in adj]
cand = []
for i in range(n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i or j in adjset[i]: continue
                dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
                if 1e-12 < dd < 4.0 - 1e-9 and col[i] == col[j]:
                    cand.append((i, j))
print(f"  {len(cand)} spindle-able candidates   [{time.time()-t0:.0f}s]", flush=True)

random.seed(1)
def kempe(col):
    v = random.randrange(n)
    a = col[v]; b = random.choice([c for c in range(K) if c != a])
    comp = {v}; q = deque([v])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in comp and col[w] in (a, b):
                comp.add(w); q.append(w)
    for u in comp:
        col[u] = b if col[u] == a else a
    return len(comp)
sizes = []
for t in range(1, SWAPS + 1):
    sizes.append(kempe(col))
    if t % 200 == 0:
        cand = [(i, j) for i, j in cand if col[i] == col[j]]
    if t % 4000 == 0:
        assert all(col[a] != col[b] for a, b in E), "a Kempe swap broke properness"
        print(f"    {t} swaps (mean component {sum(sizes[-4000:])/4000:.1f}): "
              f"{len(cand)} candidates left   [{time.time()-t0:.0f}s]", flush=True)
cand = [(i, j) for i, j in cand if col[i] == col[j]]
print(f"  after {SWAPS} Kempe swaps: {len(cand)} candidates survive   "
      f"[{time.time()-t0:.0f}s]", flush=True)

nsel = n * K + 1; forced = []; calls = 0
while cand:
    i, j = cand[0]
    sel = nsel; nsel += 1
    for c in range(K):
        s.add_clause([-sel, -X(i, c), -X(j, c)])
    r = s.solve(assumptions=[sel]); calls += 1
    # READ THE MODEL BEFORE TOUCHING THE SOLVER.  add_clause moves CaDiCaL out
    # of its satisfied state, and a first version retired the selector first and
    # then asked for the colouring -- a fatal API error, raised on the very first
    # call that came back satisfiable.
    c2 = colouring() if r else None
    s.add_clause([-sel])
    if not r:
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        forced.append((i, j, dd))
        print(f"  *** FORCED TOGETHER AT FIVE: ({i},{j}) d^2 = {dd:.6f} ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        json.dump({"graph": NAME, "forced": forced},
                  open(f"{ROOT}/data/kempe_forced.json", "w"))
        cand = cand[1:]; continue
    for _ in range(300): kempe(c2)
    cand = [(a, b) for a, b in cand if c2[a] == c2[b]]
    print(f"    solver call {calls}: {len(cand)} left   [{time.time()-t0:.0f}s]",
          flush=True)
print(f"\n  forced pairs at five: {len(forced)}   [{time.time()-t0:.0f}s]", flush=True)
