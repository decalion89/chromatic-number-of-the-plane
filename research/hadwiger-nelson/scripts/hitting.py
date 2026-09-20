"""rho is a hitting-set number, and that gives a cheap necessary test.

S is rainbow-forcing exactly when it meets every colour class of every proper
k-colouring: if some class misses S, that colour is absent from S.  So

    rho = the minimum set meeting every REALISABLE colour class,

a hitting set over a family that is exponential but samplable -- each proper
colouring contributes k disjoint classes, and each is one fast SAT call.

The minimum hitting set of a SUBfamily is at most rho, so a sample gives a
LOWER bound, which is the useful direction here: if the sample already needs
more than 63 vertices to hit, the graph cannot host a core of three and is
ruled out without any hard UNSAT at all.

Computed exactly by MaxSAT over the sampled classes.  Sampling is randomised
by permuting the variable order between calls, and blocking clauses stop the
solver returning the same colouring twice.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
from pysat.solvers import Solver


def sample_classes(graph, k, rounds=60, seed=0):
    """Colour classes from many different proper k-colourings."""
    n = graph.n

    def x(v, c):
        return 1 + v * k + c

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in graph.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])

    rng = random.Random(seed)
    fam, s = [], Solver(name="cd19", bootstrap_with=cls)
    try:
        for _ in range(rounds):
            lits = [x(v, rng.randrange(k)) for v in rng.sample(range(n), n)]
            if not s.solve(assumptions=[]):
                break
            m = set(s.get_model())
            colour = [next(c for c in range(k) if x(v, c) in m)
                      for v in range(n)]
            for c in range(k):
                fam.append(frozenset(v for v in range(n) if colour[v] == c))
            # forbid this exact colouring so the next call differs
            s.add_clause([-x(v, colour[v]) for v in range(n)])
    finally:
        s.delete()
    return fam


def min_hitting_set(fam, n):
    w = WCNF()
    for cl in fam:
        if not cl:
            return None            # an empty class cannot be hit
        w.append([v + 1 for v in cl])
    for v in range(n):
        w.append([-(v + 1)], weight=1)
    with RC2(w) as rc2:
        model = rc2.compute()
    return sum(1 for v in range(n) if model[v] > 0)


if __name__ == "__main__":
    from hn import degrey
    from hn.graph import build_graph

    print("Sa at four colours first, where rho is 7:", flush=True)
    sa = build_graph(degrey.build_Sa())
    t = time.time()
    fam = sample_classes(sa, 4, rounds=40)
    lo = min_hitting_set(fam, sa.n)
    print(f"  {len(fam)} sampled classes -> hitting set {lo}  "
          f"(a lower bound on rho = 7)  [{time.time()-t:.0f}s]", flush=True)

    print("\nde Grey's G at five colours, where the target is 63:", flush=True)
    g = degrey.build_G()
    for rounds in (20, 60, 150, 400):
        t = time.time()
        fam = sample_classes(g, 5, rounds=rounds)
        lo = min_hitting_set(fam, g.n)
        print(f"  {len(fam)} classes from {rounds} colourings -> "
              f"rho(G,5) >= {lo}"
              + ("   *** ALREADY OVER 63 ***" if lo > 63 else "")
              + f"  [{time.time()-t:.0f}s]", flush=True)
        if lo > 63:
            break
