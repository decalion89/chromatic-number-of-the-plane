"""Shrink the hitting set GREEDILY, which is what actually forces diversity.

The previous shrink only dropped vertices that no class needed uniquely, so it
returned sets of 40 to 56 where a single vertex sometimes hit the whole
family -- the minimum read 1 at fifteen hundred and three thousand classes.
A set that loose barely constrains the escaping colouring, which is the same
disease as before in a milder form.

Greedy set cover is far closer to the minimum: take the vertex hitting the
most classes, remove those classes, repeat.  The result is a handful of
vertices rather than dozens, and a colouring that must escape a handful is
pushed much further from the ones already seen.

Same decision, same cheap side: if no 63-set hits the family, rho > 63.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
from pysat.solvers import Solver
from hn import degrey


def greedy_cover(fam, n):
    """A near-minimum hitting set: most-covering vertex first."""
    if not fam:
        return []
    remaining = list(range(len(fam)))
    sets = [set(cl) for cl in fam]
    chosen, live = [], set(remaining)
    count = [0] * n
    for i in live:
        for v in sets[i]:
            count[v] += 1
    while live:
        v = max(range(n), key=lambda u: count[u])
        chosen.append(v)
        gone = [i for i in live if v in sets[i]]
        for i in gone:
            live.discard(i)
            for u in sets[i]:
                count[u] -= 1
    return chosen


g = degrey.build_G()
n, k, BUD = g.n, 5, 63


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
colour = Solver(name="cd19", bootstrap_with=cls)

print(f"G: {n} vertices; greedy cover drives the loop, budget {BUD}",
      flush=True)
fam, t0, best = [], time.time(), 0
try:
    for rnd in range(500000):
        S = greedy_cover(fam, n)
        if len(S) > BUD:
            print(f"  round {rnd}: greedy needs {len(S)} > {BUD} -- checking "
                  f"the exact minimum", flush=True)
            w = WCNF()
            for cl in fam:
                w.append([v + 1 for v in cl])
            for v in range(n):
                w.append([-(v + 1)], weight=1)
            with RC2(w) as rc2:
                mod = rc2.compute()
            lo = sum(1 for v in range(n) if mod[v] > 0)
            print(f"    exact minimum {lo}  =>  rho >= {lo}"
                  + ("   *** OVER 63 ***" if lo > BUD else "  (greedy "
                     "overshot; continuing)"), flush=True)
            if lo > BUD:
                break
            S = [v for v in range(n) if mod[v] > 0]
        if not colour.solve(assumptions=[-x(v, 0) for v in S]):
            print(f"  round {rnd}: S of {len(S)} is FORCING -- "
                  f"rho <= {len(S)}  [{time.time()-t0:.0f}s]", flush=True)
            break
        mm = set(colour.get_model())
        fam.append([v for v in range(n) if x(v, 0) in mm])
        if len(S) > best:
            best = len(S)
            print(f"  round {rnd}: |F| = {len(fam)}, greedy cover {len(S)}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
finally:
    colour.delete()
print(f"largest greedy cover seen: {best}")
