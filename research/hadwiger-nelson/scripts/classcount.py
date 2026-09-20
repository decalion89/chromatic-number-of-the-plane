"""How many distinct colour classes does a graph actually admit?

rho is the minimum set meeting every realisable colour class, so it is small
when the classes are few or when they overlap heavily.  Counting them is
therefore the structural quantity behind ambient forcing, and it is directly
measurable: ask for a colour-0 class, forbid exactly that class, ask again.

Sa reaches rho = 5 at four colours.  If its classes number in the hundreds
that explains it outright, and if de Grey's G has vastly more at five colours
then the difficulty there is counted rather than guessed.

Forbidding a class means forbidding the assignment that produced it, not the
colouring -- so each round removes one class from the family rather than one
colouring, and the count is of CLASSES.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver
from hn import degrey
from hn.graph import build_graph


def count_classes(g, k, cap=200000, report=2000):
    n = g.n

    def x(v, c):
        return 1 + v * k + c

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    s = Solver(name="cd19", bootstrap_with=cls)
    seen, t0 = set(), time.time()
    try:
        for rnd in range(cap):
            if not s.solve():
                return len(seen), True, time.time() - t0
            m = set(s.get_model())
            C = frozenset(v for v in range(n) if x(v, 0) in m)
            seen.add(C)
            # forbid exactly this class: not every member takes colour 0 while
            # every non-member avoids it
            s.add_clause([-x(v, 0) for v in C] + [x(v, 0) for v in range(n)
                                                  if v not in C])
            if rnd and rnd % report == 0:
                print(f"    {len(seen)} classes  [{time.time()-t0:.0f}s]",
                      flush=True)
    finally:
        s.delete()
    return len(seen), False, time.time() - t0


for name, pts, k in (("Sa", degrey.build_Sa(), 4),
                     ("de Grey S", degrey.build_S(), 3)):
    g = build_graph(pts)
    print(f"{name} at k = {k}: {g.n} vertices", flush=True)
    n, done, el = count_classes(g, k, cap=40000)
    print(f"  {'exactly' if done else 'at least'} {n} distinct colour classes"
          f"  [{el:.0f}s]\n", flush=True)
