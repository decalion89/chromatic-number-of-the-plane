"""rho with symmetry breaking that is actually valid, which is the delicate part.

The unbroken instance is slow because colours 1..4 are interchangeable. Pinning
a triangle to 1, 2, 3 removes that -- but the colour-permutation argument the
whole method rests on then needs care, because it is only allowed to permute
colours OTHER than 0.

    Given a proper colouring c avoiding colour 0 on S, the triangle carries
    three distinct colours. If none of them is 0, permuting {1,2,3,4} sends
    them to 1, 2, 3 while preserving "avoids 0 on S", and the pinned formula
    is contradicted. If one of them IS 0, no permutation fixes it and the
    argument fails.

So the pinning is sound exactly when the triangle is inside the set being
tested. The fix is one line: always report the triangle's three vertices as
part of the core. Then a colouring avoiding 0 on the reported set avoids 0 on
the triangle too, the permutation exists, and the conclusion holds. It costs
at most three vertices.

Run alongside the unbroken version rather than replacing it, so the two can be
compared -- a symmetry break that is subtly wrong would show up as a smaller
answer than the honest one.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver


def shrink(graph, k, rounds=20):
    n, NV = graph.n, graph.n * k

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

    tri = graph.find_clique(3)
    if tri:
        for i, u in enumerate(tri):
            cls.append([x(u, i + 1)])          # colours 1, 2, 3
        print(f"  triangle {tri} pinned to colours 1,2,3 "
              f"(and always kept in the core)", flush=True)
    tri = tri or []

    S, t0 = list(range(n)), time.time()
    for rnd in range(rounds):
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve(assumptions=[sel(v) for v in S]):
                print(f"    round {rnd+1}: SAT -- {len(S)} not forcing, stop",
                      flush=True)
                return S, False
            core = sorted(set(abs(l) - NV - 1 for l in s.get_core()) | set(tri))
        print(f"    round {rnd+1}: {len(S)} -> {len(core)}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if len(core) >= len(S):
            return core, True
        S = core
    return S, True


if __name__ == "__main__":
    from hn import degrey
    from hn.graph import build_graph

    print("validation -- Sa at four colours, where rho = 7:", flush=True)
    S, ok = shrink(build_graph(degrey.build_Sa()), 4)
    print(f"  {len(S)} vertices, forcing={ok}\n", flush=True)

    print("de Grey's G at five colours:", flush=True)
    g = degrey.build_G()
    S, ok = shrink(g, 5)
    print(f"  rho(G,5) <= {len(S)} of {g.n}" if ok else f"  {len(S)} not forcing")
    if ok:
        print(f"  against the 63 a core of three needs: "
              f"{'WITHIN' if len(S) <= 63 else 'still above'}")
