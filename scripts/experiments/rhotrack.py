"""The fixed loop, reporting the TRUE lower bound rather than a minimal set.

The shrink-to-minimal step makes the family diverse, which is what the loop
needs, but the size it reports is a MINIMAL hitting set and fluctuates: 46,
56, 49, 47 on de Grey's G without trending.  The number that means something
is the MINIMUM hitting set of the family collected, because the family is a
subfamily of all realisable colour classes and so

    minimum hitting set of F  <=  rho.

It is a valid lower bound at every moment, it only rises, and when it passes
63 the graph is ruled out for a core of three.  Computed by MaxSAT every few
thousand rounds -- too slow to do every round, which is why the loop itself
still runs on the cheap minimal shrink.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from pysat.card import CardEnc, EncType
from pysat.examples.rc2 import RC2
from pysat.formula import IDPool, WCNF
from pysat.solvers import Solver
from hn import degrey
from rhotight import shrink_hitting

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

print(f"G: {n} vertices; lower bound on rho(G,5), budget {BUD}", flush=True)
fam, t0, best = [], time.time(), 0
try:
    for rnd in range(500000):
        if not hit.solve():
            print(f"  round {rnd}: no {BUD}-set hits {len(fam)} classes -- "
                  f"rho > {BUD}  [{time.time()-t0:.0f}s]", flush=True)
            break
        m = hit.get_model()
        S = shrink_hitting([v for v in range(n) if m[v] > 0], fam)
        if not colour.solve(assumptions=[-x(v, 0) for v in S]):
            print(f"  round {rnd}: S of {len(S)} is FORCING -- "
                  f"rho <= {len(S)}  [{time.time()-t0:.0f}s]", flush=True)
            break
        mm = set(colour.get_model())
        c0 = [v for v in range(n) if x(v, 0) in mm]
        fam.append(c0)
        hit.add_clause([v + 1 for v in c0])
        if rnd and rnd % 1500 == 0:
            w = WCNF()
            for cl in fam:
                w.append([v + 1 for v in cl])
            for v in range(n):
                w.append([-(v + 1)], weight=1)
            with RC2(w) as rc2:
                mod = rc2.compute()
            lo = sum(1 for v in range(n) if mod[v] > 0)
            best = max(best, lo)
            print(f"  round {rnd}: |F| = {len(fam)}, MINIMUM hitting set "
                  f"{lo}  =>  rho >= {lo}"
                  + ("   *** OVER 63 ***" if lo > BUD else "")
                  + f"  [{time.time()-t0:.0f}s]", flush=True)
            if lo > BUD:
                break
finally:
    colour.delete()
    hit.delete()
print(f"best lower bound: rho(G,5) >= {best}")
