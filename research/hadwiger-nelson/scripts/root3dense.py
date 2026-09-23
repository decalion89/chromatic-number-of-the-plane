"""Is any pair at distance sqrt3 forced APART at five colours?  The named target.

The escape analysis returned a uniform answer: of the thirty richest candidate
points in the 803-graph, twenty-five need exactly ONE pair forced apart, and
every one of those pairs sits at squared distance 3.  Any of them will do.  So
the whole first step -- mu_5 from 2 to 3 -- reduces to a single question:

        can two points at distance sqrt3 be forced to differ
        in every proper 5-colouring of some unit-distance graph?

It is cheap to ask.  A pair is forced apart exactly when the graph plus the
clause "these two share a colour" has no 5-colouring, so each candidate is one
assumption call, and only the pairs at d^2 = 3 need asking -- a few thousand
instead of a few million.

Asked of the bigger graphs, not just the one the requirement came from: a pair
that is free in the 803-point graph may be pinned inside a larger one, and the
requirement is about the pair, not about which graph it lives in.

The rhombus is what makes sqrt3 interesting.  Two points that far apart have
exactly two common unit-distance neighbours, and those two are adjacent to each
other -- so at THREE colours the configuration forces the pair EQUAL.  At five
it forces nothing, and what is wanted is the opposite.
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
t0 = time.time()
K = 5
# The densest 5-chromatic objects, where N(u) u N(v) is largest.
#
# Forcing u and v apart is exactly the contracted graph refusing five, and
# the contracted vertex is adjacent to N(u) u N(v).  At mean degree 16.8 that
# union runs to about thirty points, against the twelve to seventeen a single
# unit circle offers -- so a sqrt3 pair is a far larger target than a blocked
# point, and the place to aim it is where the degree is highest.
for name in ("five_dense_2.json", "five_twotune_small.json"):
    try:
        d = json.load(open(f"{ROOT}/data/{name}"))
    except FileNotFoundError:
        continue
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    print(f"\n  {name}  n={n}   [{time.time()-t0:.0f}s]", flush=True)
    three = F.rational(Fr(3))
    # pairs at squared distance 3, found with a float filter and confirmed
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 2), int(p.fy // 2))].append(i)
    gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
    pairs = []
    for i in range(n):
        cx, cy = int(gx[i] // 2), int(gy[i] // 2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cells.get((cx + dx, cy + dy), ()):
                    if j <= i:
                        continue
                    ex, ey = gx[j] - gx[i], gy[j] - gy[i]
                    if abs(ex * ex + ey * ey - 3.0) < 1e-9 and \
                       (g.vertices[i] - g.vertices[j]).norm2() == three:
                        pairs.append((i, j))
    print(f"    {len(pairs)} pairs at d^2 = 3   [{time.time()-t0:.0f}s]",
          flush=True)
    if not pairs:
        continue
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
        print("    refuses five", flush=True); s.delete(); continue
    # One solve per pair is the wrong shape.  On the 1139-point graph each
    # call is 17 ms and 3138 of them cost a minute; on the 2041-point graph
    # the base colouring alone is 50 s, every assumption call restarts cadical,
    # and 8256 of them is two days.
    #
    # Filter first, and the filter CERTIFIES: a pair that shares a colour in
    # one exhibited proper colouring is PROVED not forced apart, so a handful
    # of colourings eliminate almost everything and only the survivors are
    # worth a solve.  Diversity comes from pinning a few independent vertices
    # to random colours, since cadical ignores set_phases.
    import random
    rng = random.Random(3)
    alive = set(pairs)
    for rnd in range(25):
        ass = []
        if rnd:
            picked = []
            for v in rng.sample(range(n), 80):
                if all(u not in g.adj[v] for u in picked):
                    picked.append(v)
                if len(picked) == 8:
                    break
            ass = [X(v, rng.randrange(K)) for v in picked]
        s.conf_budget(2_000_000)
        if s.solve_limited(assumptions=ass) is not True:
            continue
        pos = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
        alive = {(a, b) for (a, b) in alive if col[a] != col[b]}
        print(f"    round {rnd}: {len(alive)} pairs still never sharing"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        if not alive:
            break
    forced, tested, budget_out = [], 0, 0
    for (a, b) in sorted(alive):
        tested += 1
        s.conf_budget(4_000_000)
        r = s.solve_limited(assumptions=[X(a, 0), X(b, 0)])
        if r is False:
            forced.append((a, b))
            print(f"    *** FORCED APART AT FIVE: {a},{b} at d^2 = 3 ***",
                  flush=True)
        elif r is None:
            budget_out += 1
    s.delete()
    print(f"    {len(pairs)} pairs, {tested} survived the filter, "
          f"{len(forced)} forced apart, "
          f"{budget_out} out of budget   [{time.time()-t0:.0f}s]", flush=True)
    if forced:
        json.dump({"graph": name, "forced_apart_at_root3": forced},
                  open(f"{ROOT}/data/root3_forced_{name}", "w"))
