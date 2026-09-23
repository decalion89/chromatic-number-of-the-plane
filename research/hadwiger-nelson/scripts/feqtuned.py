"""The forced-equal filter at five on the tuned graphs, run to exhaustion.

The tuned chain produced 5-chromatic graphs at degrees the project had never
reached -- 13.35 and 13.82 against a previous ceiling of 13.78 only at 4273
points -- and none of them has been asked the question that matters: is any
pair forced to share a colour in every 5-colouring?

Attrition is fast in this direction.  Two random vertices agree with
probability about 1/5, so the survivor count falls by a factor of five per
colouring and a couple of million pairs are gone in ten rounds, against the
hundred and forty the forced-APART filter needed.  So the filter is run until
it empties rather than for a fixed budget, and every pair it drops is dropped
because an exhibited proper colouring separates it -- a certificate, not a
sample.

Anything that survives is confirmed directly: forbid the pair from sharing any
colour and ask whether five colours remain possible.  A no is chi(R^2) >= 6 by
one rotation, since the spindle angle of any distance above 1/2 exists.
"""
import sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, required_radical
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
K = 5
for name in ("five_tuned_1_1.json", "five_tuned_4_1.json", "five_tuned_1_3.json"):
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
    t1 = time.time()
    if not s.solve():
        print(f"    *** REFUSES FIVE ***   [{time.time()-t1:.0f}s]", flush=True)
        json.dump({"field_generators": list(F.gens), "n": n, "m": m,
                   "source": name,
                   "points": [[[[c.numerator, c.denominator] for c in q.x.c],
                               [[c.numerator, c.denominator] for c in q.y.c]]
                              for q in g.vertices]},
                  open(f"{ROOT}/data/six_candidate.json", "w"))
        s.delete(); continue
    print(f"    base colouring in {time.time()-t1:.0f}s", flush=True)
    rng = random.Random(1729)
    parts = [list(range(n))]
    rounds = 0
    def npairs(ps):
        return sum(len(p) * (len(p) - 1) // 2 for p in ps)
    while parts and rounds < 60:
        rounds += 1
        ass = []
        if rounds > 1:
            picked = []
            for v in rng.sample(range(n), 80):
                if all(u not in g.adj[v] for u in picked):
                    picked.append(v)
                if len(picked) == 8:
                    break
            ass = [X(v, rng.randrange(K)) for v in picked]
        s.conf_budget(3_000_000)
        if s.solve_limited(assumptions=ass) is not True:
            continue
        pos = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
        out = []
        for p in parts:
            b = defaultdict(list)
            for v in p:
                b[col[v]].append(v)
            out.extend(x for x in b.values() if len(x) > 1)
        parts = out
        print(f"    round {rounds}: {npairs(parts)} pairs still always equal"
              f"   [{time.time()-t0:.0f}s]", flush=True)
    cand = [(p[i], p[j]) for p in parts for i in range(len(p))
            for j in range(i + 1, len(p))]
    conf = []
    for (x, y) in cand[:800]:
        s.conf_budget(4_000_000)
        if s.solve_limited(assumptions=[X(x, 0), X(y, 1)]) is False:
            dd = (g.vertices[x] - g.vertices[y]).norm2()
            conf.append((x, y, str(dd)))
            print(f"    *** FORCED AT FIVE: {x},{y} d^2={dd} "
                  f"radical {required_radical(dd)} ***", flush=True)
    s.delete()
    print(f"    {len(cand)} candidates, {len(conf)} confirmed forced"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if conf:
        json.dump({"graph": name, "forced_at_five": conf},
                  open(f"{ROOT}/data/forced_five_{name}", "w"))
