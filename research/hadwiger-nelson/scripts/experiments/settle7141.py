"""Settle the 2199 candidates on the symmetric graph, one at a time and soundly.

The group-expanded filter cut 25 493 370 pairs of the 7141-point C6-invariant
graph down to 2199 candidates at five colours.  They were never settled: the
one-call certificate -- is G plus all 2199 candidate edges still 5-colourable?
-- is a very hard SAT instance and never returned.  Each candidate on its own
is easy, and any single one that survives is the whole problem.

Soundness first.  The base CNF pins a triangle, which is free for DECIDING
colourability but destroys colour symmetry, so "x = 0 and y = 1 is impossible"
no longer means x and y must agree.  The sound test is the direct one: forbid
them from sharing ANY colour and ask whether five colours remain possible.
That is five clauses, not an assumption -- so each pair gets a selector s_k,
the five clauses are gated behind it once, and the test is solve(assumptions =
[s_k]).  Unasserted selectors let the solver switch the other pairs off, so one
warm solver settles them all, keeping every learned clause.

The group halves the work twice over.  Forcing is equivariant, so a pair is
forced exactly when its whole orbit is, and only orbit representatives need a
call.

An UNSAT is chi(R^2) >= 6: the pair is forced equal in every 5-colouring, and
the rotation by 2 arcsin(1/(2d)) about one endpoint carries the other to
distance exactly 1 from itself.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, _rot60, required_radical
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
SCR = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
t0 = time.time()
d = json.load(open(f"{ROOT}/data/five_symmetric.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n; K = 5
cand = [tuple(x) for x in json.load(open(f"{SCR}/cand7141.json"))]
print(f"n={n} m={sum(len(a) for a in g.adj)//2}, {len(cand)} candidates",
      flush=True)

pos = {p: i for i, p in enumerate(g.vertices)}
rot = _rot60(F)
perm = None
if all(rot(p) in pos for p in g.vertices):
    perm = [pos[rot(p)] for p in g.vertices]
    print("  C6 acts; deduplicating by orbit", flush=True)
reps, seen = [], set()
for a, b in cand:
    key = (min(a, b), max(a, b))
    if key in seen:
        continue
    reps.append(key)
    if perm:
        x, y = a, b
        for _ in range(6):
            x, y = perm[x], perm[y]
            seen.add((min(x, y), max(x, y)))
    else:
        seen.add(key)
print(f"  {len(reps)} orbit representatives   [{time.time()-t0:.0f}s]",
      flush=True)

X = lambda v, c: 1 + v * K + c
SEL = lambda k: 1 + n * K + k
# No at-most-one clauses here.  They are needed when a COLOUR is read off a
# model, and nothing is read off these: the question is decision only.  The
# edge clauses already make adjacent colour SETS disjoint, so a satisfying
# assignment yields a proper colouring by picking any true colour per vertex,
# and the pair test -- forbid a and b from sharing any colour -- is sound the
# same way.  Dropping them removes 71 410 clauses and the search they induce.
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
for i, v in enumerate(g.find_clique(3) or []):
    cnf.append([X(v, i)])
    for c in range(K):
        if c != i:
            cnf.append([-X(v, c)])
for k, (a, b) in enumerate(reps):
    for c in range(K):
        cnf.append([-SEL(k), -X(a, c), -X(b, c)])
print(f"  cnf built: {len(cnf)} clauses   [{time.time()-t0:.0f}s]", flush=True)

# Every selector must be assumed OFF except the one under test.  Left free
# they are unconstrained variables the solver may set TRUE, which switches on
# all 376 "must differ" constraints at once and turns an easy instance into a
# very hard one -- and cadical ignores set_phases, so the default cannot be
# steered.  Assumptions are the only reliable switch.
OFF = [-SEL(k) for k in range(len(reps))]
s = Solver(name="cd19", bootstrap_with=cnf)
if not s.solve(assumptions=OFF):
    print("  *** the graph itself refuses five ***", flush=True)
    sys.exit()
print(f"  base 5-colouring found   [{time.time()-t0:.0f}s]", flush=True)

forced = []
for k, (a, b) in enumerate(reps):
    t1 = time.time()
    r = s.solve(assumptions=[SEL(k)] + OFF[:k] + OFF[k+1:])
    if not r:
        dd = (g.vertices[a] - g.vertices[b]).norm2()
        forced.append((a, b, str(dd)))
        print(f"  *** FORCED AT FIVE: {a},{b}  d^2 = {dd}  "
              f"radical {required_radical(dd)} ***", flush=True)
        json.dump({"n": n, "forced_at_five": forced},
                  open(f"{ROOT}/data/forced_at_five_7141.json", "w"))
    if k % 20 == 0:
        print(f"    ..{k}/{len(reps)}  {len(forced)} forced  "
              f"(last call {time.time()-t1:.1f}s)   [{time.time()-t0:.0f}s]",
              flush=True)
s.delete()
print(f"  DONE: {len(forced)} forced pairs at five of {len(reps)} "
      f"representatives   [{time.time()-t0:.0f}s]", flush=True)
if not forced:
    print("  every candidate separated -- no forcing at five here", flush=True)
