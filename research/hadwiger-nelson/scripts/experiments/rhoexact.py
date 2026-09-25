"""rho exactly, by implicit hitting set -- and a growing lower bound meanwhile.

S is rainbow-forcing exactly when it meets every colour class of every proper
k-colouring, so rho is a minimum hitting set over that family.  The family is
exponential, but it never has to be written down:

  1. keep a finite family F of classes seen so far;
  2. compute its minimum hitting set S -- since F is a subfamily, |S| <= rho,
     so every round prints a valid LOWER BOUND on rho;
  3. ask the solver for a proper colouring in which some class misses S;
  4. if it finds one, that class is a new member of F -- go to 1;
  5. if it cannot, S hits every class there is, so S is forcing and rho = |S|.

Random sampling does not work here: blocking one colouring at a time returns
near-identical ones and the bound saturates at 5 within a second.  Asking for
a colouring that DEFEATS the current S is what makes each new class
informative.

The expensive UNSAT comes once, at the end, on a set of size rho rather than
on all 1581 vertices.  And if it never arrives, the bound still climbs -- past
63 would rule the graph out for a core of three without any hard proof.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
from pysat.solvers import Solver


def rho_exact(graph, k, target=None, rounds=100000, report=25):
    n = graph.n

    def x(v, c):
        return 1 + v * k + c

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in graph.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    solver = Solver(name="cd19", bootstrap_with=cls)

    fam, t0, best = [], time.time(), 0
    try:
        for rnd in range(rounds):
            if fam:
                w = WCNF()
                for cl in fam:
                    w.append([v + 1 for v in cl])
                for v in range(n):
                    w.append([-(v + 1)], weight=1)
                with RC2(w) as rc2:
                    model = rc2.compute()
                S = [v for v in range(n) if model[v] > 0]
            else:
                S = []
            if len(S) > best:
                best = len(S)
                if target and best > target:
                    print(f"  round {rnd}: rho >= {best} > {target} -- "
                          f"RULED OUT  [{time.time()-t0:.0f}s]", flush=True)
                    return best, False
            if rnd % report == 0:
                print(f"  round {rnd}: |F| = {len(fam)}, rho >= {len(S)}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
            # a colouring whose colour-0 class misses S
            if not solver.solve(assumptions=[-x(v, 0) for v in S]):
                print(f"  round {rnd}: no such colouring -- S is forcing, "
                      f"rho = {len(S)}  [{time.time()-t0:.0f}s]", flush=True)
                return len(S), True
            m = set(solver.get_model())
            fam.append(frozenset(v for v in range(n) if x(v, 0) in m))
    finally:
        solver.delete()
    return best, False


if __name__ == "__main__":
    from hn import degrey
    from hn.graph import build_graph

    print("Sa at four colours, where greedy deletion found rho = 7:",
          flush=True)
    r, exact = rho_exact(build_graph(degrey.build_Sa()), 4, report=10)
    print(f"  rho = {r} (exact={exact})\n", flush=True)

    print("de Grey's G at five colours, target 63:", flush=True)
    r, exact = rho_exact(degrey.build_G(), 5, target=63, report=25)
    print(f"  rho {'=' if exact else '>='} {r}")
