"""Blocking and folding pull against each other: measure both on real graphs.

A step set can block every coset colouring, or it can fold a ball into a
high-chromatic graph, and the two want opposite things.  Blocking asks the
edge directions to cover PG(r-1,5), which takes at least six of them spread
through a module of rank r; folding asks for RELATIONS among those steps,
because Z-independent steps build a tree and a tree is bipartite.  The
quantity that governs folding is therefore the corank |S| - r of the relation
lattice -- and blocking pushes r up while folding pushes it down.

Each graph gets four numbers: distinct edge directions, the rank of the module
they generate, the corank, and the chromatic number.
"""
import sys
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from math import gcd
from hn.homcol import edge_vectors, has_homomorphism
from hn.graph import build_graph
from hn.geometry import Point, Rotation
from hn.field import Field
from hn import degrey
from pysat.solvers import Solver

FLD = Field((3, 11))


def moser_spindle():
    """Two unit rhombi sharing the origin, spun until the far vertices meet."""
    h = FLD.rational(Fraction(1, 2))
    u = Point(h, FLD.sqrt(3) * h)
    one = Point(FLD.one(), FLD.zero())
    rhombus = [Point(FLD.zero(), FLD.zero()), one, u, one + u]
    rho = Rotation(FLD.rational(Fraction(5, 6)),
                   FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
    return rhombus + [rho(p) for p in rhombus[1:]]


def tri_lattice(r=5):
    """A patch of the Eisenstein lattice: the classic 3-colourable example."""
    h = FLD.rational(Fraction(1, 2))
    u = Point(h, FLD.sqrt(3) * h)
    one = Point(FLD.one(), FLD.zero())
    return [one.scaled(a) + u.scaled(b)
            for a in range(-r, r + 1) for b in range(-r, r + 1)]


def rank(vs):
    """Rank over Q by fraction-free elimination."""
    M = [list(map(Fraction, v)) for v in vs]
    if not M:
        return 0
    r = 0
    for c in range(len(M[0])):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r


def directions(vs):
    """Distinct edge directions, identifying v with -v."""
    out = set()
    for v in vs:
        g = 0
        for t in v:
            g = gcd(g, abs(t))
        w = tuple(t // g for t in v) if g > 1 else tuple(v)
        out.add(min(w, tuple(-t for t in w)))
    return sorted(out)


def chi(g, hi=6):
    for k in range(2, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
        for a, b in g.edges():
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return f">{hi}"


cases = [("Moser spindle", build_graph(moser_spindle())),
         ("triangular patch", build_graph(tri_lattice())),
         ("de Grey S", build_graph(degrey.build_S())),
         ("de Grey Sa", build_graph(degrey.build_Sa())),
         ("de Grey Sb", build_graph(degrey.build_Sb())),
         ("de Grey Y", build_graph(degrey.build_Y())),
         ("de Grey G", degrey.build_G())]

print(f"{'graph':18} {'n':>5} {'m':>6} {'|S|':>5} {'rank':>5} "
      f"{'corank':>7} {'chi':>4}  blocked")
print("-" * 68)
for name, g in cases:
    S = directions(edge_vectors(g))
    r = rank(S)
    phi, _ = has_homomorphism(S, 5)
    print(f"{name:18} {g.n:5} {g.m:6} {len(S):5} {r:5} {len(S)-r:7} "
          f"{str(chi(g)):>4}  {'NO' if phi else 'YES'}", flush=True)
