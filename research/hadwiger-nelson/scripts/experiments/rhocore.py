"""Minimise a rainbow-forcing set with UNSAT cores instead of by deletion.

rho is the least size of a set using all k colours in every proper
k-colouring.  Testing one set is cheap in one direction and expensive in the
other -- "not forcing" is SAT, "forcing" is UNSAT -- so greedy deletion pays
the expensive price once per vertex.  Cores pay it once per ROUND.

Give each vertex a selector s_v with the clause (not s_v or not x_{v,0}), and
assume s_v for every v in S.  The formula then says: a proper k-colouring
exists in which NO vertex of S takes colour 0.  If that is UNSAT then every
colouring puts colour 0 somewhere on S, and since nothing in the encoding
distinguishes the colours, permuting c to 0 shows the same for every colour.
So S is rainbow-forcing -- and the solver's UNSAT CORE is a subset of S that
is already forcing, handed over for free.

Iterating to a fixpoint gives a small forcing set in a handful of hard solves
rather than thousands.  No symmetry breaking anywhere: the colour-permutation
argument above is what makes one colour enough, and pinning would destroy it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.solvers import Solver


def shrink(graph, k, rounds=12, verbose=True):
    n = graph.n
    NV = n * k

    def x(v, c):
        return 1 + v * k + c

    def sel(v):
        return NV + 1 + v

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in graph.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    for v in range(n):
        cls.append([-sel(v), -x(v, 0)])

    S = list(range(n))
    t0 = time.time()
    for rnd in range(rounds):
        with Solver(name="g4", bootstrap_with=cls) as s:
            sat = s.solve(assumptions=[sel(v) for v in S])
            if sat:
                if verbose:
                    print(f"    round {rnd+1}: SAT -- {len(S)} vertices are "
                          f"NOT forcing, stop", flush=True)
                return S, False
            core = sorted({abs(l) - NV - 1 for l in s.get_core()})
        if verbose:
            print(f"    round {rnd+1}: {len(S)} -> {len(core)}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
        if len(core) >= len(S):
            return core, True
        S = core
    return S, True


if __name__ == "__main__":
    from hn import degrey
    from hn.graph import build_graph

    print("validation -- Sa at four colours, where rho is known to be 7:",
          flush=True)
    sa = build_graph(degrey.build_Sa())
    S, ok = shrink(sa, 4)
    print(f"  forcing set of {len(S)} vertices (forcing={ok})\n", flush=True)

    print("de Grey's G at five colours, the target being 63:", flush=True)
    g = degrey.build_G()
    S, ok = shrink(g, 5)
    print(f"  forcing set of {len(S)} of {g.n} vertices (forcing={ok})")
    print(f"  rho(G,5) <= {len(S)}" if ok else "  no forcing set reached")
