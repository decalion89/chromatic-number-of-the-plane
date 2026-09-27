"""Unit-distance graphs over Z[zeta_n], built by walking the unit steps.

`hn.cyclotomic` gives the arithmetic; this gives the graphs.  A walk of s steps
from the origin along the 2n roots of unity lands in Z[zeta_n], every step is
a unit step by construction, and the exact test for an edge is
``norm2(a - b) == 1`` with no floating point anywhere in the decision.

For n = 3 this reproduces the Eisenstein walk the rest of the package uses.
For n = 15 it is a different world: 30 steps instead of 6, and the reachable
set is dense in the plane rather than a lattice, so a walk keeps finding new
points at every radius instead of filling a fixed grid.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from typing import Dict, List, Sequence, Set, Tuple

from .cyclotomic import CycloRing
from .graph import UnitDistanceGraph

_TOL = 1e-7


class CycloPoint:
    """A point of Z[zeta_n], carrying its ring so equality stays exact."""

    __slots__ = ("ring", "c", "_z")

    def __init__(self, ring: CycloRing, c: Tuple[int, ...]):
        self.ring = ring
        self.c = c
        self._z = None

    @property
    def z(self) -> complex:
        if self._z is None:
            self._z = self.ring.to_complex(self.c)
        return self._z

    def __add__(self, other) -> "CycloPoint":
        return CycloPoint(self.ring, self.ring.add(self.c, other.c))

    def __sub__(self, other) -> "CycloPoint":
        return CycloPoint(self.ring, self.ring.sub(self.c, other.c))

    def dist2(self, other) -> Tuple[int, ...]:
        return self.ring.norm2(self.ring.sub(self.c, other.c))

    def is_unit_apart(self, other) -> bool:
        return self.ring.is_unit_apart(self.c, other.c)

    def __eq__(self, other) -> bool:
        return isinstance(other, CycloPoint) and self.c == other.c

    def __hash__(self) -> int:
        return hash(self.c)

    def __repr__(self) -> str:
        return f"CycloPoint({self.c})"


def cyclo_origin(ring: CycloRing) -> CycloPoint:
    return CycloPoint(ring, ring.zero())


def cyclo_units(ring: CycloRing) -> List[CycloPoint]:
    return [CycloPoint(ring, u) for u in ring.unit_steps()]


def cyclo_walk(ring: CycloRing, steps: int, cap: int = 60000) -> List[CycloPoint]:
    """Every point reachable from the origin in at most `steps` unit steps."""
    units = cyclo_units(ring)
    seen: Dict[Tuple[int, ...], None] = {ring.zero(): None}
    frontier = [cyclo_origin(ring)]
    for _ in range(steps):
        nxt = []
        for p in frontier:
            for u in units:
                q = p + u
                if q.c not in seen:
                    seen[q.c] = None
                    nxt.append(q)
                    if len(seen) >= cap:
                        return [CycloPoint(ring, c) for c in seen]
        frontier = nxt
        if not frontier:
            break
    return [CycloPoint(ring, c) for c in seen]


def build_cyclo_graph(points: Sequence[CycloPoint]) -> UnitDistanceGraph:
    """Float hash proposes the pairs, exact ring arithmetic decides them."""
    pts = list(points)
    n = len(pts)
    buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, p in enumerate(pts):
        z = p.z
        buckets[(int(z.real // 1.0), int(z.imag // 1.0))].append(i)
    adj: List[Set[int]] = [set() for _ in range(n)]
    for (bx, by), idxs in buckets.items():
        near: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                near.extend(buckets.get((bx + dx, by + dy), ()))
        for i in idxs:
            zi = pts[i].z
            for j in near:
                if j <= i:
                    continue
                d = zi - pts[j].z
                if abs(abs(d) - 1.0) > _TOL:
                    continue
                if pts[i].is_unit_apart(pts[j]):     # exact, decides the edge
                    adj[i].add(j)
                    adj[j].add(i)
    return UnitDistanceGraph(pts, adj)


def cyclo_graph(n: int, steps: int, cap: int = 60000) -> UnitDistanceGraph:
    ring = CycloRing(n)
    return build_cyclo_graph(cyclo_walk(ring, steps, cap))


# -- constructions over Q(zeta_n) --------------------------------------------

def cyclo_point(field, c) -> CycloPoint:
    return CycloPoint(field, c)


def moser_spindle_cyclotomic(field):
    """The Moser spindle as seven points of Q(zeta_n), 11 | n.

    Two unit rhombi sharing the origin.  A rhombus 0, 1, u, 1+u with u a
    primitive sixth root has its far vertex at distance sqrt(3), and the
    spindle rotation rho = (5 + sqrt(-11))/6 satisfies |1 - rho|^2 = 1/3, so
    the two far vertices land exactly one apart.  Nothing here needs a real
    quadratic field: the rotation is a multiplication.
    """
    from .cyclotomic import moser_rotation

    one = field.one()
    u = field.neg(field.mul(field.zeta(field.n // 3), field.zeta(field.n // 3)))
    rhombus = [field.zero(), one, u, field.add(one, u)]
    rho = moser_rotation(field)
    pts = rhombus + [field.mul(rho, q) for q in rhombus[1:]]
    return [CycloPoint(field, c) for c in pts]


def cyclo_generated(field, seeds, gens, rounds: int, cap: int = 40000,
                    radius: float = 0.0):
    """Close a seed set under a list of maps, each given as a field element.

    A generator ``(kind, value)`` is either ``("add", v)`` for a translation by
    v or ``("mul", v)`` for a rotation about the origin by v.  Rotations of
    order 11 and 33 are what Q(zeta_33) brings that no multiquadratic field
    has, so they belong in the generator list beside the unit steps.
    """
    # Z[zeta_n] is dense for phi(n) > 2, so growth spreads thin unless it is
    # held in: 20000 points let loose covered a wide disc at 4.4 edges each and
    # forced nothing.  Capping the radius spends the same budget inside the
    # ball that holds the targets, where the constraint has to come from.
    #
    # The complex value rides along rather than being recomputed.  The radius
    # test needs it for every candidate, and `to_complex` is a degree-long sum
    # -- twenty complex multiply-adds at n = 33 -- while a translation shifts
    # it by a constant and a rotation multiplies it by one, both O(1).  The
    # embedding is a ring homomorphism, so carrying it costs nothing in
    # exactness: it decides only which points to *try*, and every point kept is
    # exact.
    zg = [(kind, v, field.to_complex(v)) for kind, v in gens]
    seen = {p.c: p.z for p in seeds}
    frontier = list(seen.items())
    for _ in range(rounds):
        nxt = []
        for c, zc in frontier:
            for kind, v, zv in zg:
                if kind == "add":
                    zq = zc + zv
                    if radius and abs(zq) > radius:
                        continue
                    q = field.add(c, v)
                else:
                    zq = zc * zv
                    if radius and abs(zq) > radius:
                        continue
                    q = field.mul(c, v)
                if q not in seen:
                    seen[q] = zq
                    nxt.append((q, zq))
                    if len(seen) >= cap:
                        return _points(field, seen)
        frontier = nxt
        if not frontier:
            break
    return _points(field, seen)


def _points(field, seen) -> List[CycloPoint]:
    out = []
    for c, z in seen.items():
        p = CycloPoint(field, c)
        p._z = z                       # already known; skip the degree-long sum
        out.append(p)
    return out


# -- the integer fast path ----------------------------------------------------

def common_denominator(steps) -> int:
    """Least d with every step's coefficients in (1/d) Z."""
    from math import lcm

    d = 1
    for s in steps:
        for x in s:
            d = lcm(d, Fraction(x).denominator)
    return d


def scaled_graph_walk(field, steps, rounds: int, radius: float,
                      cap: int = 400000, seeds=None, centres=None):
    """Grow by translation with integer coordinates instead of fractions.

    Hashing a long tuple of Fractions once per candidate is what generation
    actually spends its time on -- not the arithmetic.  Every step has a
    bounded denominator, so one common denominator clears them all and a point
    becomes a tuple of ints: addition is int addition and the hash is an int
    hash.  Exactly the lesson `hn.fast.IntBasis` already learned for the
    multiquadratic side.

    `seeds` and `centres` are what make an *asymmetric* ball possible.  A ball
    a rotation preserves never forces a proper subset of a target orbit, which
    is what every core equal to its orbit size was saying; growing from several
    centres breaks that, and at order 15 it took the core from 15 down to 5.

    Returns the scaled integer coordinates, their complex values, and the
    denominator, so a caller can wire the edges without leaving integers.
    """
    den = common_denominator(steps)
    ints = [tuple(int(Fraction(x) * den) for x in s) for s in steps]
    zs = [field.to_complex(s) for s in steps]
    if seeds is None:
        seeds = [tuple(Fraction(0) for _ in range(field.degree))]
    iseeds = [tuple(int(Fraction(x) * den) for x in c) for c in seeds]
    if centres is None:
        centres = [0j]
    seen = {c: field.to_complex(s) for c, s in zip(iseeds, seeds)}
    frontier = list(seen.items())
    for _ in range(rounds):
        nxt = []
        stop = False
        for c, zc in frontier:
            for iv, zv in zip(ints, zs):
                zq = zc + zv
                if radius and all(abs(zq - m) > radius for m in centres):
                    continue
                q = tuple(a + b for a, b in zip(c, iv))
                if q not in seen:
                    seen[q] = zq
                    nxt.append((q, zq))
                    if len(seen) >= cap:
                        stop = True
                        break
            if stop:
                break
        if stop or not nxt:
            break
        frontier = nxt
    return seen, ints, den


def step_edges(field, points, steps) -> UnitDistanceGraph:
    """Edges by table lookup: p ~ p + u for u a unit step in the set.

    A float spatial hash asks every nearby pair, which is fine for a lattice
    and hopeless for a dense ring: 109000 points inside radius 2 put ~8600 in
    every unit cell, so each point proposes tens of thousands of candidates.
    A set closed under a known step list does not need the question asked --
    its edges *are* the steps, found by one dictionary lookup each.

    Sound but not complete.  Q(zeta_n) has modulus-one elements outside any
    finite step list, so two points can be one apart with their difference
    absent from `steps`, and that edge is missed.  The result is a subgraph of
    the true unit-distance graph on these points, which is still a
    unit-distance graph: a colouring bound proved on it holds for the plane,
    and only the sharpness is lost, never the soundness.
    """
    index = {p.c: i for i, p in enumerate(points)}
    adj: List[Set[int]] = [set() for _ in points]
    for u in steps:
        uc = tuple(Fraction(x) for x in u)
        for i, p in enumerate(points):
            j = index.get(field.add(p.c, uc))
            if j is not None and j != i:
                adj[i].add(j)
                adj[j].add(i)
    return UnitDistanceGraph(list(points), adj)


def scaled_graph(field, steps, rounds: int, radius: float, cap: int = 400000,
                 seeds=None, centres=None):
    """Walk and wire in one pass, entirely in integer coordinates.

    Splitting the two costs more than the walk itself: 7921 points and 198
    steps is 1.57 million lookups, and rebuilding a long Fraction tuple for
    each one took 54 seconds against 0 for the walk.
    """
    seen, ints, den = scaled_graph_walk(field, steps, rounds, radius, cap,
                                        seeds, centres)
    order = list(seen)
    index = {c: i for i, c in enumerate(order)}
    adj: List[Set[int]] = [set() for _ in order]
    for iv in ints:
        for i, c in enumerate(order):
            j = index.get(tuple(a + b for a, b in zip(c, iv)))
            if j is not None and j != i:
                adj[i].add(j)
                adj[j].add(i)
    inv = Fraction(1, den)
    pts = []
    for c in order:
        q = CycloPoint(field, tuple(x * inv for x in c))
        q._z = seen[c]
        pts.append(q)
    return UnitDistanceGraph(pts, adj)
