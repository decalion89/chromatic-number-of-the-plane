"""How much graph does ambient forcing need?  Measure the support.

Sa has rho = 5 at four colours, with a witness of five points that induce a
handful of edges and are only 3-chromatic.  The forcing comes from the other
392 vertices -- but how many of them?

Give every vertex OUTSIDE the witness a selector, assume them all, and ask for
a proper 4-colouring leaving colour 0 off the witness.  It is UNSAT, and the
solver's core names a subset of the ambient vertices that already suffices.
Iterating shrinks it.  What comes back is the AMBIENT SUPPORT: the smallest
piece of Sa that makes those five points force all four colours.
\nIf the support is small the mechanism is local and transplantable; if it is
most of the graph, ambient forcing is as expensive as the graph itself, and
the five-colour version would need something the size of de Grey's G to make
63 points force.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from hn import degrey
from hn.graph import build_graph
from rhotight import run
from pysat.solvers import Solver

g = build_graph(degrey.build_Sa())
k = 4
print(f"Sa: {g.n} vertices, {g.m} edges", flush=True)

# the witness, by decision at the known value
S = None
for budget in (5,):
    v = run(g, k, budget, rounds=100000, report=10 ** 9)
print("  (rho = 5 established earlier; recovering a witness)", flush=True)

# recover one explicitly: grow a family, take the greedy cover that forces
def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])


def greedy_cover(fam, n):
    if not fam:
        return []
    sets = [set(cl) for cl in fam]
    live = set(range(len(fam)))
    cnt = [0] * n
    for i in live:
        for u in sets[i]:
            cnt[u] += 1
    out = []
    while live:
        u = max(range(n), key=lambda z: cnt[z])
        out.append(u)
        for i in [i for i in live if u in sets[i]]:
            live.discard(i)
            for z in sets[i]:
                cnt[z] -= 1
    return out


colour = Solver(name="cd19", bootstrap_with=cls)
fam = []
for rnd in range(100000):
    C = greedy_cover(fam, g.n)
    if not colour.solve(assumptions=[-x(v, 0) for v in C]):
        S = sorted(C)
        print(f"  witness of {len(S)}: {S}  (round {rnd})", flush=True)
        break
    m = set(colour.get_model())
    fam.append([v for v in range(g.n) if x(v, 0) in m])
colour.delete()

# now the support: selectors on everything outside S
NV = g.n * k
sel = lambda v: NV + 1 + v
base = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        base.append([-x(a, c), -x(b, c)])
for v in S:
    base.append([-x(v, 0)])          # the witness must avoid colour 0
outside = [v for v in range(g.n) if v not in set(S)]
# a selector switches a vertex's "needs a colour" clause on
hard = [cl for cl in base if len(cl) != k or cl[0] % k != 1]
enc = []
for v in range(g.n):
    if v in set(S):
        enc.append([x(v, c) for c in range(k)])
    else:
        enc.append([-sel(v)] + [x(v, c) for c in range(k)])
form = enc + [cl for cl in base if len(cl) == 2 or len(cl) == 1]

cur, t0 = outside, time.time()
for rnd in range(10):
    with Solver(name="cd19", bootstrap_with=form) as s:
        if s.solve(assumptions=[sel(v) for v in cur]):
            print(f"  round {rnd}: SAT -- {len(cur)} ambient vertices are not "
                  f"enough", flush=True)
            break
        core = sorted({abs(l) - NV - 1 for l in s.get_core()})
    print(f"  round {rnd}: ambient support {len(cur)} -> {len(core)}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if len(core) >= len(cur):
        break
    cur = core
print(f"\nambient support: {len(cur)} of {g.n - len(S)} outside vertices "
      f"({100*len(cur)/g.n:.0f} per cent of Sa) make {len(S)} points force "
      f"all four colours")
