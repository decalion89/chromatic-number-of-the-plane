"""Six copies, and only ONE endpoint has to sit at the special radius.

The three-copy gadget demanded |v - q1| = |v - q2| = 1/sqrt3.  That is two
constraints, and the exhaustive searches returned zero.  Six copies demand
only one.

Take rotations about v by the six angles

        {0, 2pi/3, 4pi/3}  u  {beta, beta + 2pi/3, beta + 4pi/3}

for any beta, and let U be the union of the six rotated copies of H.  Two
copies clash on option i exactly when their angle difference is +- alpha_i,
where alpha_i = 2 arcsin(1/(2 r_i)) is the spindle angle of that option's
radius.  So

  option 1 at r1 = 1/sqrt3   has alpha_1 = 2pi/3, and its clash graph is the
                             TWO TRIANGLES (each coset is a triangle):
                             independence 2
  option 2 at any r2 > 1/2   has alpha_2 = beta, and its clash graph is the
                             perfect MATCHING pairing each coset member with
                             its beta-translate: independence 3

Two options can therefore cover at most 2 + 3 = 5 of the six copies, and every
copy must be covered.  Contradiction: U has no k-colouring.

Nothing has to be avoided.  Extra coincidences among the angles only add edges
to the clash graphs, which only shrinks the independent sets.  If beta = 2pi/3
the two cosets merge and it degenerates to the three-copy version, which also
works.  The only requirement left is r2 > 1/2, so that the spindle angle beta
exists at all.

Why one endpoint can never be freed.  Each clash graph has maximum degree 2, so
its independence number is at least the number of vertices over 3, with
equality only for disjoint triangles; a union of paths has independence at
least half.  Two path-unions give at least s, and the pigeonhole never bites.
So at least one clash graph must contain an odd cycle, i.e. its alpha must be a
rational multiple of 2pi.  The triangle needs alpha = 2pi/3 and radius
1/sqrt3, which Q(sqrt3) holds; the next odd cycle is the pentagon, which needs
sin 72 = sqrt(10 + 2 sqrt5)/4, a nested radical no multiquadratic field
contains.  1/sqrt3 is the only special radius available here.

The search is the conditional forced-pair filter: fix v and q1, assume
c(v) = 0 and c(q1) != 0, sample proper colourings under those assumptions, and
keep the vertices that take colour 0 in every one of them.  Those are the only
possible q2, and each survivor is confirmed by one more solve.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()

def load(name, gens):
    F = Field(gens)
    d = json.load(open(f"{ROOT}/data/{name}"))
    return F, [Point(F.element([Fr(a, b) for a, b in x]),
                     F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]

CASES = [("Sa", Field((3, 5, 7, 11)), list(build_Sa(Field((3, 5, 7, 11)))), 4)]
F1, Z = load("five_247_c.json", (3, 11, 247))
CASES.append(("five_247_c", F1, Z, 5))
F2, Zb = load("five_247.json", (3, 11, 247))
CASES.append(("five_247", F2, Zb, 5))

MODELS, BLOCK = 12, 25
for name, F, U, K in CASES:
    THIRD = F.rational(Fr(1, 3))   # r2 needs only float(d^2) > 1/4,
    # the condition for the spindle angle beta to exist at all
    g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
    d2 = [[None] * n for _ in range(n)]
    at13 = defaultdict(list)
    for i in range(n):
        pi = g.vertices[i]
        for j in range(i + 1, n):
            e = (pi - g.vertices[j]).norm2()
            d2[i][j] = d2[j][i] = e
            if e == THIRD:
                at13[i].append(j); at13[j].append(i)
    pairs = [(v, q) for v, qs in at13.items() for q in qs]
    print(f"  {name:14s} n={n} m={m} K={K}: {len(pairs)} ordered (v,q1) at "
          f"d^2=1/3   [{time.time()-t0:.0f}s]", flush=True)
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for c in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, c)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    rnd = random.Random(7)
    hits, forced = [], []
    for idx, (v, q1) in enumerate(pairs):
        s = Solver(name="cd19", bootstrap_with=cnf)
        ass = [X(v, 0), -X(q1, 0)]
        if not s.solve(assumptions=ass):
            forced.append((v, q1))
            print(f"    *** {v},{q1} is a FORCED PAIR at {K} ***", flush=True)
            s.delete(); continue
        cand = None
        for _ in range(MODELS):
            if not s.solve(assumptions=ass):
                break
            pos = set(l for l in s.get_model() if l > 0)
            col = [next(c for c in range(K) if X(u, c) in pos) for u in range(n)]
            here = set(u for u in range(n) if col[u] == 0)
            cand = here if cand is None else (cand & here)
            if len(cand) <= 1:
                break
            s.add_clause([-X(u, col[u]) for u in rnd.sample(range(n), BLOCK)])
        cand = {u for u in (cand or set())
                if u != v and u != q1 and d2[v][u] is not None and float(d2[v][u]) > 0.25}
        for u in sorted(cand):
            if not s.solve(assumptions=[X(v, 0), -X(q1, 0), -X(u, 0)]):
                hits.append((v, q1, u, str(d2[v][u])))
                print(f"    *** DISJUNCTION at {K}: v={v} q1={q1} (1/3) "
                      f"q2={u} (d^2={d2[v][u]}) ***", flush=True)
        s.delete()
        if idx % 200 == 0:
            print(f"      ..{idx}/{len(pairs)}  {len(hits)} so far"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
    print(f"    {name}: {len(forced)} forced pairs, {len(hits)} disjunctions "
          f"over {len(pairs)} pivots   [{time.time()-t0:.0f}s]", flush=True)
    if hits:
        json.dump({"graph": name, "K": K, "hits": hits[:200]},
                  open(f"{ROOT}/data/disjunctions_{name}.json", "w"))
        print(f"    written data/disjunctions_{name}.json", flush=True)
