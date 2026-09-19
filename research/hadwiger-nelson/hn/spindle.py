"""The spindling argument, automated.

This is de Grey's method, stated so a machine can run it.

    Let G be a unit-distance graph containing points p and q with
    |p - q|^2 = d, and suppose *every* proper k-colouring of G gives p and q
    the same colour.  Let rho be the rotation about p by arccos(1 - 1/(2d)),
    which moves q to a point q' with |q - q'| = 1.  Then in G union rho(G):

        colour(p) = colour(q)    (G is a copy of G)
        colour(p) = colour(q')   (rho(G) is a copy of G)
        => colour(q) = colour(q'),  yet q and q' are adjacent.

    So G union rho(G) has no proper k-colouring at all.

The whole search therefore reduces to a question a SAT solver answers directly:
is there a pair (p, q) at a spindle-able distance that every k-colouring is
forced to paint the same colour?  A graph can be k-colourable and still contain
such a pair -- which is exactly why searching balls for a non-k-colourable
subgraph finds nothing, while spindling one of those same balls succeeds.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pysat.formula import CNF
from pysat.solvers import Solver

from .coloring import DEFAULT_SOLVER
from .field import Field
from .geometry import Point, Rotation, rotation_joining
from .graph import UnitDistanceGraph, build_graph

__all__ = ["ForcedPairFinder", "spindle_union", "candidate_pairs_from",
           "ForcedDisjunctionFinder", "targets_at_third", "triple_spindle_union",
           "spindle_union_auto",
           "SeparationTest"]


class ForcedPairFinder:
    """One CNF, many 'must p and q share a colour?' queries.

    For each candidate pair a selector s_pq is added with the clauses
        (~s_pq | ~x[p][c] | ~x[q][c])   for every colour c,
    which under the assumption s_pq demands that the colour sets of p and q be
    disjoint.  UNSAT under that assumption means no proper colouring separates
    them, i.e. the pair is forced monochromatic.  (Any proper colouring gives a
    satisfying assignment with singleton colour sets, and those are disjoint
    exactly when the colours differ, so the two statements coincide.)
    """

    def __init__(self, graph: UnitDistanceGraph, k: int, pairs: Sequence[Tuple[int, int]]):
        self.graph = graph
        self.k = k
        self.pairs = list(pairs)
        n = graph.n
        self.n = n
        self.cnf = CNF()
        for v in range(n):
            self.cnf.append([self.x(v, c) for c in range(k)])
        for u, v in graph.edges():
            for c in range(k):
                self.cnf.append([-self.x(u, c), -self.x(v, c)])
        clique = graph.find_clique(min(3, k))
        if clique:
            for i, v in enumerate(clique):
                self.cnf.append([self.x(v, i)])
                for c in range(k):
                    if c != i:
                        self.cnf.append([-self.x(v, c)])
        self.sel: Dict[Tuple[int, int], int] = {}
        for idx, (p, q) in enumerate(self.pairs):
            s = 1 + n * k + idx
            self.sel[(p, q)] = s
            for c in range(k):
                self.cnf.append([-s, -self.x(p, c), -self.x(q, c)])

    def x(self, v: int, c: int) -> int:
        return 1 + v * self.k + c

    def run(
        self, solver: str = DEFAULT_SOLVER, verbose: bool = True, stop_at_first: bool = False
    ) -> List[Tuple[int, int]]:
        """Return the pairs that are forced monochromatic."""
        s = Solver(name=solver, bootstrap_with=self.cnf)
        forced: List[Tuple[int, int]] = []
        try:
            if not s.solve():
                if verbose:
                    print("    graph is already not k-colourable", flush=True)
                return []
            for (p, q) in self.pairs:
                if s.solve(assumptions=[self.sel[(p, q)]]) is False:
                    forced.append((p, q))
                    if verbose:
                        print(f"    FORCED monochromatic: {p} ~ {q}", flush=True)
                    if stop_at_first:
                        break
        finally:
            s.delete()
        return forced


def candidate_pairs_from(
    graph: UnitDistanceGraph,
    pivot: int,
    field: Field,
    max_per_distance: int = 2,
    max_d2: float = 60.0,
) -> List[Tuple[int, int]]:
    """Pairs (pivot, q) that can actually be spindled.

    A pair is usable when |pivot - q|^2 is rational *and* the resulting
    rotation arccos(1 - 1/(2 d2)) has its sine inside the field -- otherwise
    the rotated copy would leave the ring and the arithmetic would stop being
    exact.  Distances need not be integers: 4/3 and 7/3 occur in this ring and
    spindle perfectly well.

    The ball is symmetric under the 6-fold lattice rotation, so equidistant
    targets are largely interchangeable; `max_per_distance` keeps only a few
    representatives of each distance instead of all six or twelve.
    """
    from fractions import Fraction

    p = graph.vertices[pivot]
    seen_d2: Dict[Fraction, int] = {}
    out: List[Tuple[int, int]] = []
    for j, q in enumerate(graph.vertices):
        if j == pivot:
            continue
        d2 = p.dist2(q)
        if not d2.is_rational():
            continue
        val = d2.c[0]
        if val <= 0 or val > max_d2:
            continue
        if seen_d2.get(val, 0) >= max_per_distance:
            continue
        try:
            rotation_joining(val, field)
        except ValueError:
            continue  # the spindle rotation would need a larger field
        seen_d2[val] = seen_d2.get(val, 0) + 1
        out.append((pivot, j))
    return out


def spindle_union(
    graph: UnitDistanceGraph, pivot: int, target: int, field: Field
) -> UnitDistanceGraph:
    """G union rho(G), rho the rotation about `pivot` that sends `target` to
    distance exactly 1 from itself.

    The result is not k-colourable whenever (pivot, target) is a forced
    monochromatic pair for k colours in G.
    """
    p = graph.vertices[pivot]
    q = graph.vertices[target]
    d2 = p.dist2(q)
    if not d2.is_rational():
        raise ValueError(f"|p-q|^2 = {d2} is not rational; no spindle rotation")
    rot = rotation_joining(d2.c[0], field)
    turn = rot.about(p)
    image = [turn(v) for v in graph.vertices]
    if not q.is_unit_apart(turn(q)):
        raise AssertionError("spindle rotation did not place the target at distance 1")
    return build_graph(list(graph.vertices) + image)


# --- the three-copy pigeonhole spindle -----------------------------------
#
# The two-copy argument needs a pair that is forced monochromatic outright.
# One extra copy weakens what has to be forced, and the geometry allows it at
# exactly one distance.
#
# Rotated images of a point at distance r from the pivot must be pairwise at
# distance 1 for the contradiction to land.  On a circle of radius r you can
# inscribe two such points whenever r >= 1/2, but *three* only when the circle
# is the circumcircle of a unit equilateral triangle -- that is, r = 1/sqrt(3),
# or d^2 = 1/3, and then the rotations are by 120 and 240 degrees.
#
# With three copies, pigeonhole does the rest.  Suppose every k-colouring of G
# makes p monochromatic with q1 *or* with q2, both at distance 1/sqrt(3) from p.
# In G u rho(G) u rho^2(G) each copy forces one of the two, so two of the three
# copies force the same q; their images of that q are two vertices of the unit
# triangle, hence adjacent, yet both share p's colour.  Contradiction.
#
# The hypothesis is a disjunction rather than a single forced pair, which is a
# far weaker thing to ask of G -- and weaker hypotheses are what a search can
# actually find.

class ForcedDisjunctionFinder:
    """Finds {q1, q2} such that every k-colouring ties p to one of them.

    Encoded with one selector per candidate pair: asserting it demands a
    colouring in which p differs from *both*, so UNSAT under that assumption
    is exactly the disjunction being forced.
    """

    def __init__(self, graph: UnitDistanceGraph, k: int, pivot: int, targets: Sequence[int]):
        self.graph = graph
        self.k = k
        self.pivot = pivot
        self.targets = list(targets)
        n = graph.n
        self.n = n
        self.cnf = CNF()
        for v in range(n):
            self.cnf.append([self.x(v, c) for c in range(k)])
        for u, v in graph.edges():
            for c in range(k):
                self.cnf.append([-self.x(u, c), -self.x(v, c)])
        clique = graph.find_clique(min(3, k))
        if clique:
            for i, v in enumerate(clique):
                self.cnf.append([self.x(v, i)])
                for c in range(k):
                    if c != i:
                        self.cnf.append([-self.x(v, c)])
        self.sel: Dict[Tuple[int, int], int] = {}
        nxt = 1 + n * k
        for a in range(len(self.targets)):
            for b in range(a + 1, len(self.targets)):
                q1, q2 = self.targets[a], self.targets[b]
                s = nxt
                nxt += 1
                self.sel[(q1, q2)] = s
                for c in range(k):
                    self.cnf.append([-s, -self.x(pivot, c), -self.x(q1, c)])
                    self.cnf.append([-s, -self.x(pivot, c), -self.x(q2, c)])

    def x(self, v: int, c: int) -> int:
        return 1 + v * self.k + c

    def run(
        self, solver: str = DEFAULT_SOLVER, verbose: bool = True, stop_at_first: bool = True
    ) -> List[Tuple[int, int]]:
        s = Solver(name=solver, bootstrap_with=self.cnf)
        found: List[Tuple[int, int]] = []
        try:
            if not s.solve():
                if verbose:
                    print("    graph is already not k-colourable", flush=True)
                return []
            for pair, lit in self.sel.items():
                if s.solve(assumptions=[lit]) is False:
                    found.append(pair)
                    if verbose:
                        print(f"    FORCED disjunction: p ~ {pair[0]} or p ~ {pair[1]}", flush=True)
                    if stop_at_first:
                        break
        finally:
            s.delete()
        return found


def targets_at_third(graph: UnitDistanceGraph, pivot: int) -> List[int]:
    """Vertices at squared distance 1/3 from `pivot` -- the only radius where
    three rotated copies are pairwise adjacent."""
    p = graph.vertices[pivot]
    third = Fraction(1, 3)
    out = []
    for j in range(graph.n):
        if j == pivot:
            continue
        d2 = p.dist2(graph.vertices[j])
        if d2.is_rational() and d2.c[0] == third:
            out.append(j)
    return out


def triple_spindle_union(graph: UnitDistanceGraph, pivot: int) -> UnitDistanceGraph:
    """G union rho(G) union rho^2(G), rho the 120-degree rotation about `pivot`.

    Not k-colourable whenever some pair of vertices at distance 1/sqrt(3) from
    the pivot carries a forced disjunction for k colours.
    """
    from .geometry import ROT60

    p = graph.vertices[pivot]
    rho = ROT60 ** 2                      # 120 degrees
    turn = rho.about(p)
    verts = list(graph.vertices)
    image1 = [turn(v) for v in verts]
    image2 = [turn(v) for v in image1]
    return build_graph(verts + image1 + image2)


class SeparationTest:
    """One query decides a whole family of disjunctions.

    Give each target q a selector s_q asserting that p and q take disjoint
    colour sets.  Assume all of them at once:

      SAT   -- some colouring separates p from every target, so *no*
               disjunction over any subset of them is forced.  Dead end, and
               one query establishes it instead of O(|Q|^2).
      UNSAT -- the disjunction over the whole target set is forced, and the
               solver's UNSAT core hands back a subset that already forces it.
               A core of size <= 2 is what the three-copy argument needs.
    """

    def __init__(self, graph: UnitDistanceGraph, k: int, pivot: int, targets: Sequence[int]):
        self.graph = graph
        self.k = k
        self.pivot = pivot
        self.targets = list(targets)
        n = graph.n
        self.n = n
        self.cnf = CNF()
        for v in range(n):
            self.cnf.append([self.x(v, c) for c in range(k)])
        for u, v in graph.edges():
            for c in range(k):
                self.cnf.append([-self.x(u, c), -self.x(v, c)])
        clique = graph.find_clique(min(3, k))
        if clique:
            for i, v in enumerate(clique):
                if v == pivot or v in self.targets:
                    continue  # never pin a vertex the query is about
                self.cnf.append([self.x(v, i)])
        self.sel: Dict[int, int] = {}
        for idx, q in enumerate(self.targets):
            s = 1 + n * k + idx
            self.sel[q] = s
            for c in range(k):
                self.cnf.append([-s, -self.x(pivot, c), -self.x(q, c)])

    def x(self, v: int, c: int) -> int:
        return 1 + v * self.k + c

    def run(self, solver: str = DEFAULT_SOLVER, shrink_rounds: int = 40, subset=None):
        """Return (separable, core).

        `separable` True means a colouring separates the pivot from all
        targets.  Otherwise `core` is a target subset whose disjunction is
        forced, shrunk by re-solving on the previous core and then by trying
        each single deletion.
        """
        s = self._solver(solver)
        try:
            want = self.targets if subset is None else list(subset)
            assumptions = [self.sel[q] for q in want]
            if s.solve(assumptions=assumptions):
                return True, []
            core_lits = s.get_core() or assumptions
            back = {v: q for q, v in self.sel.items()}
            core = [back[l] for l in core_lits]
            for _ in range(shrink_rounds):
                res = s.solve(assumptions=[self.sel[q] for q in core])
                if res:
                    break
                new = [back[l] for l in (s.get_core() or [])]
                if not new or len(new) >= len(core):
                    break
                core = new
            for q in list(core):
                if len(core) <= 1:
                    break
                trial = [t for t in core if t != q]
                if s.solve(assumptions=[self.sel[t] for t in trial]) is False:
                    core = trial
            return False, core
        finally:
            pass

    def _solver(self, solver: str):
        if getattr(self, "_s", None) is None:
            self._s = Solver(name=solver, bootstrap_with=self.cnf)
        return self._s

    def close(self):
        if getattr(self, "_s", None) is not None:
            self._s.delete()
            self._s = None


def _squarefree_factors(n: int) -> List[int]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def spindle_union_auto(graph: UnitDistanceGraph, pivot: int, target: int):
    """`spindle_union`, enlarging the field on demand.

    The earlier search only considered targets whose spindle rotation already
    lived in the graph's field.  That restriction was never needed and cost a
    large share of the candidates -- in de Grey's graph the missing radicals
    are sqrt17, sqrt13 and sqrt19, together accounting for more candidate pairs
    than sqrt2 does by two orders of magnitude.  Whether a pair is forced is
    decided by SAT on the graph alone; the field only has to be wide enough to
    write the rotated copy down, and widening it is free.

    Returns (graph, field).
    """
    from .field import Field, embed
    from .geometry import Point, required_radical

    p, q = graph.vertices[pivot], graph.vertices[target]
    d2 = p.dist2(q)
    if not d2.is_rational():
        raise ValueError(f"|p-q|^2 = {d2} is not rational")
    val = d2.c[0]
    if val < Fraction(1, 4):
        raise ValueError(f"d^2 = {val} < 1/4: no rotation separates such a pair by 1")
    needed = _squarefree_factors(required_radical(val))
    field = graph.vertices[0].field
    missing = [g for g in needed if g not in field.gens and g > 1]
    big = Field(tuple(sorted(set(field.gens) | set(missing)))) if missing else field
    verts = [Point(embed(v.x, big), embed(v.y, big)) for v in graph.vertices]
    rot = rotation_joining(val, big)
    turn = rot.about(verts[pivot])
    image = [turn(v) for v in verts]
    if not verts[target].is_unit_apart(image[target]):
        raise AssertionError("spindle rotation did not land the target at distance 1")
    return build_graph(verts + image), big
