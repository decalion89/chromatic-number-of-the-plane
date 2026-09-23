"""The glue, one rung higher.

At four colours the ladder was:  Sa is 4-chromatic and has NO pair forced to
share a colour in every 4-colouring; glue Sa to a rotated copy of itself about
a vertex and the union HAS forced pairs.  The glue manufactured the forcing out
of a carrier that had none.  Spindle one of them and the result refuses four.

The exact analogue one rung up is what is needed for six, and it has never been
run: take a 5-chromatic graph, glue it to a rotated copy of itself about a
vertex, and ask for a pair forced equal in every FIVE-colouring.  If the glue
manufactures forcing at five the way it did at four, one rotation closes the
problem.

The carrier is the smallest 5-chromatic graph in the project, 803 points in
Q(sqrt3, sqrt11, sqrt247).  Size is the whole point: G has 1582 vertices and a
single 5-colouring of the 7141-point symmetric graph costs cadical 322 s, so a
filter needing hundreds of them is only affordable on a small object.  At 1600
vertices and five colours the instance is tight (it sits at its chromatic
number), which is the regime where minisat with randomised phases and blocking
beats cadical.

Two things this gets right that earlier passes got wrong:

  * at-most-one clauses are included, so a colour read off a model is the
    vertex's colour and not the smallest of three true variables;
  * the glue is chosen once and held fixed, so the filter measures the union
    and not a re-choice of the best overlap.
"""
import sys, time, json, random
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 247))
d = json.load(open(ROOT + "/data/five_247_c.json"))
H = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
gH = build_graph(H)
print(f"carrier n={gH.n} m={sum(len(a) for a in gH.adj)//2}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

r60 = _rot60(F)
ROTS = [("rot60", r60), ("rot120", Rotation(r60.cos * r60.cos - r60.sin * r60.sin,
                                            r60.cos * r60.sin + r60.sin * r60.cos)),
        ("rot180", Rotation(F.rational(Fr(-1)), F.rational(Fr(0))))]

# Rank the glues by overlap.  A glue that shares many points ties the two
# copies together at many places at once, which is what forces agreement; a
# glue sharing one point is a hinge and provably cannot raise anything.
S = set(H)
deg = sorted(range(gH.n), key=lambda v: -len(gH.adj[v]))
best = []
for tag, r in ROTS:
    for w in deg[:240]:
        rot = r.about(H[w])
        ov = sum(1 for p in H if rot(p) in S)
        best.append((ov, tag, w))
best.sort(reverse=True)
print(f"  best glues: {[(o, t, w) for o, t, w in best[:8]]}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

ov, tag, w = best[0]
rot = dict(ROTS)[tag].about(H[w])
seen, V = set(H), list(H)
for p in H:
    q = rot(p)
    if q not in seen:
        seen.add(q); V.append(q)
g = build_graph(V)
n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"  glue {tag} about vertex {w}: overlap {ov}, union n={n} m={m} "
      f"deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]", flush=True)

EDGES = list(g.edges())
def base_clauses(K):
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):                      # at most one -- so a model reads
        for a in range(K):                  # back as an honest colouring
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in EDGES:
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    return cl, X

for K in (4, 5):
    cl, X = base_clauses(K)
    t1 = time.time()
    s = Solver(name="cd19", bootstrap_with=cl); ok = s.solve(); s.delete()
    print(f"  {K}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if K == 4 and ok:
        print("  the union lost the carrier -- abort", flush=True); sys.exit()

# ---- the forced-pair filter, at FIVE ----------------------------------------
K = 5
cl, X = base_clauses(K)
tri = g.find_clique(K - 1) or []
for i, x in enumerate(tri):                 # kill the colour permutations
    cl.append([X(x, i)])
print(f"  pinned a clique of {len(tri)}", flush=True)

rnd = random.Random(20260923)
parts = [list(range(n))]
def npairs(ps): return sum(len(p) * (len(p) - 1) // 2 for p in ps)
extra = []
BLOCK = 40
for it in range(400):
    s = Solver(name="m22", bootstrap_with=cl)
    for e in extra:
        s.add_clause(e)
    s.set_phases([rnd.choice([1, -1]) * (1 + v) for v in range(n * K)])
    if not s.solve():
        print(f"  round {it}: exhausted (blocking made it UNSAT)", flush=True)
        s.delete(); break
    mod = s.get_model(); s.delete()
    pos = set(l for l in mod if l > 0)
    col = [0] * n
    for v in range(n):
        for c in range(K):
            if X(v, c) in pos:
                col[v] = c; break
    out = []
    for p in parts:
        b = {}
        for v in p:
            b.setdefault(col[v], []).append(v)
        out.extend(x for x in b.values() if len(x) > 1)
    parts = out
    sample = rnd.sample(range(n), BLOCK)
    extra.append([-X(v, col[v]) for v in sample])
    if it % 10 == 0 or npairs(parts) < 200:
        print(f"  round {it}: {npairs(parts)} candidate pairs in "
              f"{len(parts)} classes   [{time.time()-t0:.0f}s]", flush=True)
    if not parts:
        print("  *** every pair separated: NO forced pair at five ***",
              flush=True); break

cand = [(p[i], p[j]) for p in parts for i in range(len(p))
        for j in range(i + 1, len(p))]
print(f"  {len(cand)} candidates survive   [{time.time()-t0:.0f}s]", flush=True)

# A candidate is only a candidate.  Ask the solver directly: forbid the two
# from sharing a colour and see whether five colours are still possible.
conf = []
for (u, v) in cand[:400]:
    s = Solver(name="cd19", bootstrap_with=cl)
    for c in range(K):
        s.add_clause([-X(u, c), -X(v, c)])
    if not s.solve():
        conf.append((u, v))
        dd = (g.vertices[u] - g.vertices[v]).norm2()
        print(f"  *** FORCED AT FIVE: {u},{v}  d^2 = {dd} ***", flush=True)
    s.delete()
print(f"  confirmed {len(conf)} forced pairs at five   [{time.time()-t0:.0f}s]",
      flush=True)
if conf:
    json.dump({"field_generators": list(F.gens), "n": n, "m": m,
               "glue": [tag, w], "forced_at_five": conf,
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g.vertices]},
              open(ROOT + "/data/forced_at_five.json", "w"))
    print("  *** written data/forced_at_five.json ***", flush=True)
