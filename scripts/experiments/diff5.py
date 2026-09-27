"""Forced DIFFERENT, the dual target, never searched.

Every filter here has hunted pairs forced to AGREE, because that is what the
spindle consumes.  The bipartite-neighbourhood theorem points at the other
relation.

  If free@5 = 0 then every degree-4 vertex has its four neighbours taking four
  distinct colours in every 5-colouring.  Inside a neighbourhood the induced
  graph has maximum degree 2 and no 4-cycle -- four steps of +-60 degrees
  cannot sum to 360 -- so at most 3 of those 6 pairs are edges.  At least 3
  must be forced APART without being adjacent.

And the relation closes the problem on its own.  Five points pairwise forced
apart use all five colours in every 5-colouring; a sixth forced apart from all
five has nothing left.  So a "rainbow 5-set plus one" is chi(R^2) >= 6, with no
rotation and no spindle -- a completely different shape of certificate from
de Grey's.

Note what cannot be found inside a 5-COLOURABLE graph: six pairwise forced
apart, because that is exactly a graph with no 5-colouring.  Five is the most
that can exist, and finding five is the construction step; the sixth point is
what a rotation would then have to supply.

The filter is the mirror of the forced-equal one and just as cheap.  Sample
proper colourings; a pair that AGREES in any one of them is proved not forced
apart, so the survivors are certificates the other way round.  Each survivor is
then confirmed directly: add the clauses saying the two share a colour and ask
whether five colours remain possible.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
K = 5
# The two small ones only: at 10 % attrition per round the filter needs
# hundreds of rounds to empty, which is seconds at 800 points and an hour
# at 2000.
for name in ("five_247_c.json", "five_247.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {name}  n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if not s.solve():
        print("    refuses five!", flush=True); s.delete(); continue
    rng = random.Random(97)
    # A pair survives only if it differs in EVERY sampled colouring.  Start
    # from the non-edges that differ in the first, then intersect.
    keep = None
    rounds = 0
    while rounds < 400 and (keep is None or keep):
        rounds += 1
        ass = []
        if keep is not None:
            picked = []
            for v in rng.sample(range(n), 60):
                if all(u not in g.adj[v] for u in picked):
                    picked.append(v)
                if len(picked) == 6:
                    break
            ass = [X(v, rng.randrange(K)) for v in picked]
        s.conf_budget(600_000)
        if s.solve_limited(assumptions=ass) is not True:
            continue
        pos = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
        by = defaultdict(list)
        for v in range(n):
            by[col[v]].append(v)
        same = set()
        for vs in by.values():
            for i, x in enumerate(vs):
                for y in vs[i+1:]:
                    same.add((x, y))
        if keep is None:
            keep = set()
            for x in range(n):
                for y in range(x + 1, n):
                    if col[x] != col[y] and y not in g.adj[x]:
                        keep.add((x, y))
        else:
            keep -= same
        if rounds % 25 == 0 or not keep:
            print(f"    round {rounds}: {len(keep)} non-adjacent pairs still "
                  f"always apart   [{time.time()-t0:.0f}s]", flush=True)
        if not keep:
            break
    if not keep:
        print("    none -- every non-adjacent pair agrees somewhere",
              flush=True); s.delete(); continue
    conf = []
    # a RANDOM sample: sorting puts every low index first, which is not a
    # sample of the survivors but a sample of one corner of the graph
    probe = rng.sample(sorted(keep), min(600, len(keep)))
    for (x, y) in probe:
        s.conf_budget(2_000_000)
        r = s.solve_limited(assumptions=[X(x, 0), X(y, 0)])
        if r is False:
            conf.append((x, y, str((g.vertices[x] - g.vertices[y]).norm2())))
            print(f"    *** FORCED APART: {x},{y} "
                  f"d^2={(g.vertices[x]-g.vertices[y]).norm2()} ***", flush=True)
    s.delete()
    print(f"    {len(conf)} confirmed forced-apart non-adjacent pairs"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if conf:
        json.dump({"graph": name, "forced_apart": conf},
                  open(f"{ROOT}/data/forced_apart_{name}", "w"))
