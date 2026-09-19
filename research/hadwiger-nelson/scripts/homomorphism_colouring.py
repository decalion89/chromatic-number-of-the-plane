"""Homomorphism colourings, by SAT, in any dimension.

A map phi from the edge-vector module to Z/n that is nonzero on every edge
vector IS a proper n-colouring: c(p) = phi(p), and adjacent points differ by
an edge vector, so their colours differ. It is the most structural colouring
there is -- constant on cosets of ker phi -- and its existence means no subset
of the module can ever need more than n colours.

So a necessary condition for a unit-distance graph to be 6-chromatic is that
its edge-vector module admit NO homomorphism to Z/5 avoiding the edge vectors.
That is a finite check and, encoded as SAT over the coordinates, one that
works in dimension 32 where brute force cannot.

Each coordinate takes a value in Z/n; each edge vector d contributes the
constraint sum_i phi_i d_i != 0 mod n, encoded by carrying a running partial
sum through the coordinates.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from math import gcd
from pysat.formula import IDPool
from pysat.solvers import Solver
from hn.degrey import build_G, build_Sa
from hn.graph import build_graph
from hn.mixed import three_hexagon_gadget


def edge_vectors(g):
    out = set()
    for a, b in g.edges():
        d = (g.vertices[b].x - g.vertices[a].x,
             g.vertices[b].y - g.vertices[a].y)
        out.add(tuple(d[0].c) + tuple(d[1].c))
    den = 1
    for v in out:
        for q in v:
            den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
    return sorted({tuple(int(Fr(q) * den) for q in v) for v in out})


def has_homomorphism(vecs, n):
    dim = len(vecs[0])
    pool = IDPool()
    cls = []

    def phi(i, val):
        return pool.id(("phi", i, val))

    for i in range(dim):
        cls.append([phi(i, v) for v in range(n)])
        for a in range(n):
            for b in range(a + 1, n):
                cls.append([-phi(i, a), -phi(i, b)])

    for t, d in enumerate(vecs):
        idx = [i for i in range(dim) if d[i] % n]
        if not idx:
            return None, "an edge vector is divisible by n: impossible"

        def s(j, val):
            return pool.id(("s", t, j, val))

        # running partial sum after the first j coordinates of idx
        for val in range(n):
            cls.append([-phi(idx[0], val), s(0, (val * d[idx[0]]) % n)])
        for j in range(1, len(idx)):
            for prev in range(n):
                for val in range(n):
                    cls.append([-s(j - 1, prev), -phi(idx[j], val),
                                s(j, (prev + val * d[idx[j]]) % n)])
        cls.append([-s(len(idx) - 1, 0)])       # the total must not be zero

    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            return None, "no homomorphism"
        m = set(sv.get_model())
        return tuple(next(v for v in range(n) if phi(i, v) in m)
                     for i in range(dim)), "found"


t0 = time.time()
for name, g in (("three-hexagon gadget", build_graph(three_hexagon_gadget()[2])),
                ("Sa", build_graph(build_Sa())),
                ("de Grey's G", build_G())):
    vecs = edge_vectors(g)
    print(f"{name}: {g}, {len(vecs)} edge vectors, dim {len(vecs[0])}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for n in (4, 5, 6, 7):
        phi, why = has_homomorphism(vecs, n)
        print(f"    Z/{n}: {why}" + (f"  phi = {phi}" if phi else "")
              + f"  [{time.time()-t0:.0f}s]", flush=True)
