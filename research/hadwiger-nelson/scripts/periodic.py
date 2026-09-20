"""A strictly stronger gate than blocking.

A coset colouring is the special case n = 5 of something much larger.  For ANY
finite abelian group G and homomorphism phi from the edge module M to G with
0 not in phi(D), colouring v by the colour of phi(v - v0) in a proper
colouring of the Cayley graph Cay(G, phi(D)) is proper on the graph, because
adjacent vertices differ by an element of D and phi sends it to a nonzero
connection element.  So

    chi(Gamma) <= chi(Cay(G, phi(D)))   for every such phi,

and therefore a 6-chromatic unit-distance graph must have EVERY finite abelian
quotient Cayley graph of chromatic number at least 6.

Blocking is the case G = Z/5: there phi(D) misses 0, so the Cayley graph is
the complete graph K_5 and its chromatic number is exactly 5.  Which is why
"no phi to Z/5" was the right condition -- but only the first of a family.
At G = Z/n for n > 5 the Cayley graph is usually NOT complete, and its
chromatic number can still be 5 or less, so the screen keeps biting after
blocking stops.

Run here at every modulus up to 40, with phi sampled at random: a single hit
proves the graph 5-colourable and no amount of blocking can save it.
"""
import sys, random, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver


def cayley_chromatic(n, S, cap=5):
    """chi(Cay(Z/n, S)) if it is at most cap, else None."""
    S = {s % n for s in S} | {(-s) % n for s in S}
    S.discard(0)
    E = [(v, (v + s) % n) for v in range(n) for s in S if v < (v + s) % n]
    for k in range(2, cap + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return None


def periodic_screen(vecs, moduli=range(5, 41), trials=4000, seed=0,
                    cap=5, report=None):
    """Look for a periodic 5-colouring: a modulus and a phi that work.

    Returns (n, phi, chi) on the first hit, or None if every sample failed --
    which is evidence, not proof, since phi is sampled rather than enumerated.
    """
    rng = random.Random(seed)
    r = len(vecs[0])
    for n in moduli:
        best = None
        for _ in range(trials):
            c = [rng.randrange(n) for _ in range(r)]
            S = set()
            ok = True
            for d in vecs:
                t = sum(x * y for x, y in zip(c, d)) % n
                if t == 0:
                    ok = False
                    break
                S.add(t)
            if not ok:
                continue
            k = cayley_chromatic(n, S, cap)
            if k is not None and (best is None or k < best):
                best = k
            if k is not None and k <= cap:
                return n, c, k
        if report:
            report(n, best)
    return None


if __name__ == "__main__":
    from hn.degrey import build_Sa
    from hn.graph import build_graph
    from hn.homcol import (denominator_29_directions, edge_vectors,
                           has_homomorphism)

    t0 = time.time()
    print("control -- Sa, which is 4-chromatic and has a coset colouring:",
          flush=True)
    v = edge_vectors(build_graph(build_Sa()))
    hit = periodic_screen(v, moduli=[5], trials=200)
    print(f"  screen at n = 5: {'5-colourable' if hit else 'no hit'}"
          f"  (blocking says {'no' if has_homomorphism(v, 5)[0] else 'yes'} "
          f"coset colouring)  [{time.time()-t0:.0f}s]", flush=True)

    print("\nthe denominator-29 directions, which BLOCK at n = 5:", flush=True)
    d29 = denominator_29_directions()
    assert has_homomorphism(d29, 5)[0] is None
    hit = periodic_screen(d29, moduli=range(6, 41), trials=600,
                          report=lambda n, b: print(
                              f"    n = {n:2d}: best Cayley chi seen "
                              f"{b}", flush=True))
    print(f"  stronger screen: "
          + (f"5-COLOURABLE at n = {hit[0]}, Cayley chi = {hit[2]}"
             if hit else "no periodic 5-colouring found")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
