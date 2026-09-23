"""The smallest RAINBOW SET: a set of vertices that must show all five colours.

mu has only ever been measured on unit circles here, where the induced graph is
bipartite and the answer is 2.  But mu is defined for any set, and for a general
set S of a 5-chromatic graph C,

    S is a rainbow set  <=>  no 5-colouring of C leaves a colour off S
                        <=>  mu_5(C, S) = 5

and S = everything always works, because C is 5-chromatic.  So the question is
the MINIMUM size, and it is one UNSAT core: assume colour 0 absent from every
vertex, read the core, shrink.

Why the minimum matters.  Adjoin to C a new point h placed with no unit-distance
neighbour at all.  Then in every 5-colouring c(h) is some colour, that colour
appears on S, and so

    every 5-colouring of C + h gives  c(h) = c(v)  for some v in S

which is the consumable shape -- a statement about ONE point, chainable by
rotations about it, with all the partners it can name staying on their own
circles about h where the clash geometry is computable.  And here h is a FREE
parameter, and so are the rotation angles; the only thing not free is |S|, which
is the width of the disjunction and the number of copies it will cost.

A rainbow set lying ON a unit circle would be a blocked point outright, and the
grading forbids nothing of the kind -- chi(N(p)) = 2 bounds the LOCAL structure,
not mu.  The hub route needs no such thing: S may sit anywhere in C.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAMES = sys.argv[1:] or ["five_247_c.json"]
for NAME in NAMES:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    E = list(g.edges())
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in E:
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    ass = [-X(v, 0) for v in range(n)]
    ok = s.solve(assumptions=ass)
    print(f"\n{NAME} n={n}: colour 0 absent everywhere -> {ok} "
          f"(False is the graph being 5-chromatic)   [{time.time()-t0:.0f}s]",
          flush=True)
    if ok:
        s.delete(); continue
    core = set(s.get_core() or [])
    print(f"  first core: {len(core)}   [{time.time()-t0:.0f}s]", flush=True)
    for _ in range(60):
        a2 = sorted(core)
        if s.solve(assumptions=a2): break
        c2 = set(s.get_core() or [])
        if len(c2) >= len(core): break
        core = c2
        print(f"    core -> {len(core)}   [{time.time()-t0:.0f}s]", flush=True)
    cur = sorted(core); i = 0
    while i < len(cur):
        trial = cur[:i] + cur[i+1:]
        if trial and not s.solve(assumptions=trial):
            c2 = set(s.get_core() or trial)
            cur = [a for a in trial if a in c2] or trial
            i = 0
        else:
            i += 1
        if i and i % 40 == 0:
            print(f"      {len(cur)} left, probe {i}   [{time.time()-t0:.0f}s]",
                  flush=True)
    S = sorted(-a // K if False else (-a - 1) // K for a in cur)
    print(f"  MINIMAL RAINBOW SET: {len(S)} vertices   [{time.time()-t0:.0f}s]",
          flush=True)
    print(f"    {S[:40]}{' ...' if len(S) > 40 else ''}", flush=True)
    sub = [(u, v) for u, v in E if u in set(S) and v in set(S)]
    print(f"    unit edges inside it: {len(sub)}", flush=True)
    json.dump({"graph": NAME, "rainbow_set": S, "size": len(S)},
              open(f"{ROOT}/data/rainbow_{NAME}", "w"))
    s.delete()
