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

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pysat.formula import CNF
from pysat.solvers import Solver

from .coloring import DEFAULT_SOLVER
from .field import Field
from .geometry import Point, Rotation, rotation_joining
from .graph import UnitDistanceGraph, build_graph

__all__ = ["ForcedPairFinder", "spindle_union", "candidate_pairs_from"]


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
