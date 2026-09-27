"""Track the LOWER bound on rho as the class family grows.

The decision loop at budget 63 keeps finding 63-sets, which says only that the
family collected so far can still be hit by 63 vertices.  The informative
number is the minimum that hits it: that minimum is a valid lower bound on
rho, since the family is a subfamily of all realisable classes, and it climbs
as the family grows.  When it passes 63, de Grey's G is ruled out for a core
of three; where it stalls is how close the graph really is.

So the loop is the same, but each report solves the hitting problem to
optimality with MaxSAT instead of merely satisfying a budget.  Escaping
classes are still generated against a cheap 63-bounded set, because that is
what makes the next class informative.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.card import CardEnc, EncType
from pysat.examples.rc2 import RC2
from pysat.formula import IDPool, WCNF
from pysat.solvers import Solver
from hn import degrey

g = degrey.build_G()
n, k, BUD = g.n, 5, 63


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
colour = Solver(name="cd19", bootstrap_with=cls)

pool = IDPool(start_from=n + 1)
card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=BUD, vpool=pool,
                      encoding=EncType.seqcounter)
hit = Solver(name="cd19", bootstrap_with=card.clauses)

print(f"G: {n} vertices; tracking the lower bound on rho(G,5)", flush=True)
fam, t0, best = [], time.time(), 0
try:
    for rnd in range(500000):
        if not hit.solve():
            print(f"  round {rnd}: no {BUD}-set hits {len(fam)} classes -- "
                  f"rho > {BUD}  [{time.time()-t0:.0f}s]", flush=True)
            break
        m = hit.get_model()
        S = [v for v in range(n) if m[v] > 0]
        if not colour.solve(assumptions=[-x(v, 0) for v in S]):
            print(f"  round {rnd}: S of {len(S)} is FORCING -- rho <= {BUD}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            break
        mm = set(colour.get_model())
        c0 = [v for v in range(n) if x(v, 0) in mm]
        fam.append(c0)
        hit.add_clause([v + 1 for v in c0])
        if rnd and rnd % 2000 == 0:
            w = WCNF()
            for cl in fam:
                w.append([v + 1 for v in cl])
            for v in range(n):
                w.append([-(v + 1)], weight=1)
            with RC2(w) as rc2:
                mod = rc2.compute()
            lo = sum(1 for v in range(n) if mod[v] > 0)
            best = max(best, lo)
            print(f"  round {rnd}: |F| = {len(fam)}, minimum hitting set "
                  f"{lo}  =>  rho >= {lo}"
                  + ("   *** OVER 63 ***" if lo > BUD else "")
                  + f"  [{time.time()-t0:.0f}s]", flush=True)
            if lo > BUD:
                break
finally:
    colour.delete()
    hit.delete()
print(f"best lower bound reached: rho(G,5) >= {best}")
