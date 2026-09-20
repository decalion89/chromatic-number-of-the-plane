"""The ambient support, by deletion rather than by core.

A single UNSAT core came back the full size of its input, which says the
solver did not shrink it -- not that nothing can be dropped.  Deletion settles
it: take the witness S, remove one ambient vertex at a time, and keep the
removal whenever S still forces all four colours in what is left.

What survives is a genuinely minimal ambient support, and its size is the
price of ambient forcing: five points forcing four colours, propped up by how
much graph?
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.graph import build_graph
from pysat.solvers import Solver

g = build_graph(degrey.build_Sa())
k = 4
S = [49, 91, 210, 215, 330]
print(f"Sa: {g.n} vertices; witness {S}", flush=True)


def forces(keep):
    """Does S force all k colours in the subgraph induced on `keep`?"""
    idx = {v: i for i, v in enumerate(keep)}
    n = len(keep)

    def x(i, c):
        return 1 + i * k + c

    cls = [[x(i, c) for c in range(k)] for i in range(n)]
    for a, b in g.edges():
        if a in idx and b in idx:
            for c in range(k):
                cls.append([-x(idx[a], c), -x(idx[b], c)])
    for c in range(k):
        with Solver(name="cd19",
                    bootstrap_with=cls + [[-x(idx[v], c)] for v in S]) as s:
            if s.solve():
                return False
    return True


keep = list(range(g.n))
assert forces(keep), "the witness must force in Sa itself"
print("  confirmed forcing in the whole graph", flush=True)

rng = random.Random(3)
order = [v for v in range(g.n) if v not in set(S)]
rng.shuffle(order)
t0, dropped = time.time(), 0
kset = set(keep)
for v in order:
    trial = sorted(kset - {v})
    if forces(trial):
        kset = set(trial)
        dropped += 1
        if dropped % 25 == 0:
            print(f"  dropped {dropped}, {len(kset)} left  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\nminimal ambient support: {len(kset)} vertices of {g.n} "
      f"({100*len(kset)/g.n:.0f} per cent), of which {len(S)} are the witness")
print(f"  so {len(kset) - len(S)} ambient vertices make {len(S)} points force "
      f"all four colours")
