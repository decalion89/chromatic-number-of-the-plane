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


def cyclo_generated(field, seeds, gens, rounds: int, cap: int = 40000):
    """Close a seed set under a list of maps, each given as a field element.

    A generator ``(kind, value)`` is either ``("add", v)`` for a translation by
    v or ``("mul", v)`` for a rotation about the origin by v.  Rotations of
    order 11 and 33 are what Q(zeta_33) brings that no multiquadratic field
    has, so they belong in the generator list beside the unit steps.
    """
    seen = {p.c: None for p in seeds}
    frontier = [p.c for p in seeds]
    for _ in range(rounds):
        nxt = []
        for c in frontier:
            for kind, v in gens:
                q = field.add(c, v) if kind == "add" else field.mul(c, v)
                if q not in seen:
                    seen[q] = None
                    nxt.append(q)
                    if len(seen) >= cap:
                        return [CycloPoint(field, x) for x in seen]
        frontier = nxt
        if not frontier:
            break
    return [CycloPoint(field, c) for c in seen]
