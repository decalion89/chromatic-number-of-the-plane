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
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
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
for name in ("five_247_c.json", "five_247.json", "five_tuned_1_1.json",
             "five_twotune_small.json"):
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
    # forced apart <=> "they share colour 0" is unsatisfiable, by symmetry
    forced, tested, budget_out = [], 0, 0
    for (a, b) in pairs:
        tested += 1
        s.conf_budget(2_000_000)
        r = s.solve_limited(assumptions=[X(a, 0), X(b, 0)])
        if r is False:
            forced.append((a, b))
            print(f"    *** FORCED APART AT FIVE: {a},{b} at d^2 = 3 ***",
                  flush=True)
        elif r is None:
            budget_out += 1
    s.delete()
    print(f"    {tested} tested, {len(forced)} forced apart, "
          f"{budget_out} out of budget   [{time.time()-t0:.0f}s]", flush=True)
    if forced:
        json.dump({"graph": name, "forced_apart_at_root3": forced},
                  open(f"{ROOT}/data/root3_forced_{name}", "w"))
