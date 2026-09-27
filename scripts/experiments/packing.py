"""A lower bound on rho from DISJOINT colour classes, using only SAT calls.

A forcing set must meet every realisable colour class.  Classes that are
pairwise disjoint therefore each demand their own vertex, so

    rho  >=  the largest number of pairwise disjoint realisable classes.

That is a packing bound, and every step of computing it is a cheap SAT call:
keep the union U of the classes chosen so far, and ask for a proper colouring
whose colour-0 class avoids U entirely.  Satisfiable means a new disjoint
class; unsatisfiable means the packing is maximal for this greedy path.

Two things fall out.  The size of the smallest class found bounds how many can
possibly be packed -- classes of 300 vertices cannot give more than five on
1581 -- so the same run measures whether the bound can even reach 63.  And
greedily preferring SMALL classes makes the packing longer, so the search asks
for a colouring whose colour-0 class avoids U and is as small as the solver
will make it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.solvers import Solver


def pack(graph, k, target=None, cap=400):
    n = graph.n

    def x(v, c):
        return 1 + v * k + c

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in graph.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    s = Solver(name="cd19", bootstrap_with=cls)
    used, sizes, t0 = set(), [], time.time()
    try:
        for i in range(cap):
            if not s.solve(assumptions=[-x(v, 0) for v in sorted(used)]):
                print(f"  no colour class avoids the {len(used)} used "
                      f"vertices: packing is maximal at {len(sizes)}",
                      flush=True)
                break
            m = set(s.get_model())
            C = [v for v in range(n) if x(v, 0) in m]
            # shrink it: any subset of an independent set stays independent,
            # but it must stay a whole COLOUR CLASS, so keep it as found
            used |= set(C)
            sizes.append(len(C))
            if len(sizes) % 10 == 0 or len(sizes) < 5:
                print(f"  {len(sizes)} disjoint classes, sizes "
                      f"{min(sizes)}..{max(sizes)}, {len(used)} of {n} "
                      f"vertices used  [{time.time()-t0:.0f}s]", flush=True)
            if target and len(sizes) > target:
                print(f"  *** rho >= {len(sizes)} > {target} ***", flush=True)
                break
    finally:
        s.delete()
    return len(sizes), sizes, len(used)


if __name__ == "__main__":
    from hn import degrey
    from hn.graph import build_graph

    print("Sa at four colours (rho = 7), as a check on the bound's strength:",
          flush=True)
    m, sizes, used = pack(build_graph(degrey.build_Sa()), 4)
    print(f"  packing {m}, class sizes {sorted(sizes)[:8]}...\n", flush=True)

    print("de Grey's G at five colours, target 63:", flush=True)
    g = degrey.build_G()
    m, sizes, used = pack(g, 5, target=63)
    print(f"  packing {m} => rho(G,5) >= {m}")
    print(f"  smallest class found: {min(sizes) if sizes else '-'}; "
          f"{used} of {g.n} vertices consumed")
    if sizes:
        print(f"  ceiling on any packing: {g.n // min(sizes)} classes, so this "
              f"bound {'can' if g.n // min(sizes) > 63 else 'CANNOT'} reach 63")
