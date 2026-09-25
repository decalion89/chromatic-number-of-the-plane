"""free@5 against mean degree, averaged instead of sampled once.

free@k counts the vertices with a spare colour in ONE proper colouring, and
different colourings differ -- the same 803-point graph came back 11.46 % in
one pass and 14.82 % in another.  Comparing a single measurement to

    free@k ~ (k-1) (1 - 1/(k-1))^deg

therefore carries a few points of noise.  Averaging over a dozen colourings,
made diverse by blocking a random sample of the previous one and randomising
the solver's phases, removes most of it and lets the ratio to the law be read
as a real number rather than an impression.

The point of the curve is the exponent.  If the 5-chromatic graphs track the
law, local rigidity at five is unreachable, because the law needs mean degree
about 35 for free@5 to fall below 1/n and a planar unit-distance graph has
O(n^{1/3}) mean degree.  If they beat it by the factor Sa beats it at four --
twenty -- then degree 19 suffices and the tuned chain already reaches 16.8.

At-most-one clauses are included here, and they matter: a colour is READ off
each model, and without them a vertex with three true colour variables reads
as the smallest, which quietly biases every vertex towards colour 0 and
inflates free@k.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, time, random, statistics
from fractions import Fraction as Fr
sys.path.insert(0, HN_DIR)
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
# The tuned graphs first: they are the ones at the new degrees, and the
# question they answer is the one that matters.
FILES = ["five_tuned_1_1.json", "five_tuned_4_1.json", "five_tuned_1_3.json",
         "five_247_c.json", "five_247_b.json", "five_247.json"]
K = 5
def law(deg):
    return (K - 1) * ((K - 2) / (K - 1)) ** deg

rows = []
for name in FILES:
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
    deg = 2.0 * m / n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    # Diversity has to come from assumptions, not phases.  cadical ignores
    # set_phases, and minisat -- which honours them -- stalls on these loose
    # instances, so neither of the usual knobs works.  Pinning twenty random
    # vertices to random colours for each round does work: it is a different
    # sub-problem every time, cadical solves each one fast, and the resulting
    # colourings are genuinely unrelated rather than one colouring nudged.
    # A pinning that happens to be unsatisfiable is simply redrawn.
    # The first call is UNBUDGETED and its answer is read for what it is.
    # A budgeted call returns None when it runs out and False when the
    # instance is unsatisfiable, and treating both as "try again" would throw
    # away the only result that matters: a 5-chromatic graph that turns out to
    # refuse five is chi(R^2) >= 6.
    s = Solver(name="cd19", bootstrap_with=cnf)
    t1 = time.time()
    if not s.solve():
        print(f"  *** {name}: REFUSES FIVE  n={n} deg={deg:.2f} "
              f"[{time.time()-t1:.0f}s] ***", flush=True)
        json.dump({"field_generators": list(F.gens), "n": n, "m": m,
                   "source": name,
                   "points": [[[[c.numerator, c.denominator] for c in q.x.c],
                               [[c.numerator, c.denominator] for c in q.y.c]]
                              for q in g.vertices]},
                  open(f"{ROOT}/data/six_candidate.json", "w"))
        print("  *** written data/six_candidate.json ***", flush=True)
        s.delete(); continue
    print(f"    {name}: base 5-colouring in {time.time()-t1:.0f}s",
          flush=True)
    rng = random.Random(41)
    # Twenty random vertices pinned to random colours is too much: many such
    # pinnings are unsatisfiable or merely hard, and the retries dominate.
    # Six independent vertices and a conflict budget keep every try cheap --
    # a try that runs over budget is abandoned, not ground out.
    vals, tries = [], 0
    while len(vals) < 12 and tries < 30:
        tries += 1
        ass = []
        if vals:
            picked = []
            for v in rng.sample(range(n), 60):
                if all(u not in g.adj[v] for u in picked):
                    picked.append(v)
                if len(picked) == 6:
                    break
            ass = [X(v, rng.randrange(K)) for v in picked]
        s.conf_budget(400_000)
        if s.solve_limited(assumptions=ass) is not True:
            continue
        pos = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
        vals.append(sum(1 for v in range(n)
                        if len(set(col[u] for u in g.adj[v])) < K - 1) / n)
    s.delete()
    mu = statistics.mean(vals)
    sd = statistics.pstdev(vals) if len(vals) > 1 else 0.0
    rows.append((deg, name, n, mu, sd, len(vals)))
    print(f"  {name:22s} n={n:5d} deg={deg:5.2f}  free@5 = "
          f"{100*mu:6.3f} % +- {100*sd:.3f}  over {len(vals)} colourings "
          f"(law {100*law(deg):6.3f} %, {law(deg)/max(mu,1e-9):5.2f}x better)"
          f"   [{time.time()-t0:.0f}s]", flush=True)

print("\n  by degree:", flush=True)
for deg, name, n, mu, sd, k in sorted(rows):
    print(f"    deg {deg:5.2f}  free@5 {100*mu:6.3f} %  law {100*law(deg):6.3f} %"
          f"  ratio {law(deg)/max(mu,1e-9):5.2f}   {name}", flush=True)
