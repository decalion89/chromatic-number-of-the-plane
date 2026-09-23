"""Split every candidate at once: one solve, and a core if it fails.

Eliminating candidate forced pairs one at a time slows to a crawl as the
colourings the warm solver returns start to resemble one another -- 54 965,
11 999, 8 444, 7 794, 7 514.  Ask the opposite question instead, for all of them
together: give each candidate pair (u,v) a selector meaning "u and v differ",
assume every selector, and solve.

  * SATISFIABLE: one colouring splits every candidate, so none is forced, and
    the whole class is settled by a single call.
  * UNSATISFIABLE: the core is a set K of pairs of which at least one is
    monochromatic in EVERY 5-colouring -- a disjunctive forcing statement over
    specific spindle-able pairs, derived rather than hunted.  |K| = 1 is a
    forced pair outright; a small K is the narrow disjunction the multispindle
    needs.

The pigeonhole trap from the deficiency episode is waiting here: thousands of
split constraints can complete a K6, which needs six colours whatever the graph
does.  So candidates are taken ONE DISTANCE CLASS AT A TIME.  A two-distance set
in the plane has at most five points, so the {1, d} graph can never contain K6,
and neither can anything selected from it; any core found is about the graph.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "tight_hexagon_4159.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n; E = list(G.edges())
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(b); adj[b].add(a)
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for a, b in E:
    for c in range(K): base.append([-X(a, c), -X(b, c)])
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
classes = defaultdict(list)
for i in range(n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i or j in adj[i]: continue
                dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
                if 1e-12 < dd < 4.0 - 1e-9:
                    classes[round(dd, 7)].append((i, j))
ranked = sorted(classes.items(), key=lambda kv: -len(kv[1]))
print(f"{NAME}: n={n}; {len(classes)} spindle-able distance classes, richest "
      f"{[len(v) for _, v in ranked[:6]]}   [{time.time()-t0:.0f}s]", flush=True)
TOP = int(sys.argv[2]) if len(sys.argv) > 2 else 12
for d2, pairs in ranked[:TOP]:
    S = lambda t: n * K + 1 + t
    cnf = list(base)
    for t, (i, j) in enumerate(pairs):
        for c in range(K): cnf.append([-S(t), -X(i, c), -X(j, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    ass = [S(t) for t in range(len(pairs))]
    s.conf_budget(20_000_000)
    r = s.solve_limited(assumptions=ass)
    if r is None:
        print(f"  d^2={d2:<11.6f} {len(pairs):>6d} pairs: budget out   "
              f"[{time.time()-t0:.0f}s]", flush=True); s.delete(); continue
    if r:
        print(f"  d^2={d2:<11.6f} {len(pairs):>6d} pairs: one colouring splits "
              f"them all -- none forced   [{time.time()-t0:.0f}s]", flush=True)
        s.delete(); continue
    core = sorted(s.get_core() or ass)
    i = 0
    while i < len(core):
        trial = core[:i] + core[i+1:]
        s.conf_budget(20_000_000)
        if trial and s.solve_limited(assumptions=trial) is False:
            c2 = set(s.get_core() or trial)
            core = [a for a in trial if a in c2] or trial; i = 0
        else:
            i += 1
    K_pairs = [pairs[a - (n * K + 1)] for a in core]
    tag = "*** A FORCED PAIR ***" if len(K_pairs) == 1 else "a disjunction"
    print(f"  d^2={d2:<11.6f} {len(pairs):>6d} pairs: UNSAT, minimal core "
          f"{len(K_pairs)} -- {tag}   [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"graph": NAME, "d2": d2, "core": K_pairs},
              open(f"{ROOT}/data/splitall_{round(d2,6)}.json", "w"))
    s.delete()
