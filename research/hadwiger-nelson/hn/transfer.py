"""Gluing copies along an interface instead of at a point.

The spindle argument joins its copies at the pivot: every copy is an isometry
fixing one point, so every chosen image carries that point's colour, and the
union fails when no choice survives.  Everything in this package amplifies
that way, and `hn.transversal` says how far it goes -- capacity 2 over any
multiquadratic field, 5 at best anywhere.

Periodic structures are usually attacked the other way, by a transfer matrix.
Take a small set W of points and ask which colourings of it a k-colouring of
the configuration can realise.  Take a second set W' congruent to W, and ask
which pairs (colouring of W, colouring of W') are *jointly* realisable.  That
relation is a transfer matrix, and chaining copies -- gluing each copy's W' to
the next copy's W, which are literally the same points -- iterates it.  If the
reachable set empties after m steps, a chain of m copies has no k-colouring.

Two things make it sound.

The patterns are actual colour tuples, not partitions.  Glued points are the
same points, so their colours agree exactly; only one global colour
permutation may be quotiented out, and that is done once at the start rather
than per step.

And extra edges only help.  Copies placed by isometries may meet elsewhere and
create unit-distance edges the relation never accounted for, but an edge only
removes colourings.  The computed reachable set is therefore a superset of the
true one: if it empties, the true one has emptied already.
"""
from __future__ import annotations

from itertools import product
from typing import Dict, List, Optional, Sequence, Set, Tuple

from pysat.formula import CNF
from pysat.solvers import Solver

from .graph import UnitDistanceGraph

DEFAULT_SOLVER = "cd19"

Pattern = Tuple[int, ...]


def _base_cnf(graph: UnitDistanceGraph, k: int) -> CNF:
    cnf = CNF()
    for v in range(graph.n):
        cnf.append([1 + v * k + c for c in range(k)])
    for u, v in graph.edges():
        for c in range(k):
            cnf.append([-(1 + u * k + c), -(1 + v * k + c)])
    return cnf


def realisable_patterns(graph: UnitDistanceGraph, k: int, W: Sequence[int],
                        solver: str = DEFAULT_SOLVER) -> Set[Pattern]:
    """Colourings of W that some k-colouring of the graph induces."""
    W = list(W)
    out: Set[Pattern] = set()
    with Solver(name=solver, bootstrap_with=_base_cnf(graph, k)) as s:
        for pat in product(range(k), repeat=len(W)):
            if s.solve(assumptions=[1 + v * k + c for v, c in zip(W, pat)]):
                out.add(pat)
    return out


def transfer_relation(graph: UnitDistanceGraph, k: int, W: Sequence[int],
                      Wp: Sequence[int], solver: str = DEFAULT_SOLVER
                      ) -> Dict[Pattern, Set[Pattern]]:
    """Which colourings of W' can accompany each colouring of W.

    One solver, one assumption list per pair, so the whole matrix costs
    k^(2|W|) incremental calls rather than that many formulas.
    """
    W, Wp = list(W), list(Wp)
    rel: Dict[Pattern, Set[Pattern]] = {}
    with Solver(name=solver, bootstrap_with=_base_cnf(graph, k)) as s:
        for p in product(range(k), repeat=len(W)):
            ass_p = [1 + v * k + c for v, c in zip(W, p)]
            if not s.solve(assumptions=ass_p):
                continue
            reach: Set[Pattern] = set()
            for q in product(range(k), repeat=len(Wp)):
                ass = ass_p + [1 + v * k + c for v, c in zip(Wp, q)]
                if s.solve(assumptions=ass):
                    reach.add(q)
            if reach:
                rel[p] = reach
    return rel


def chain_length(rel: Dict[Pattern, Set[Pattern]],
                 start: Optional[Set[Pattern]] = None,
                 limit: int = 64) -> Optional[int]:
    """Copies needed before no pattern survives, or None if it never empties.

    Forward reachability on the transfer relation. The set can only shrink or
    cycle, so `limit` steps past a repeat is already conclusive.
    """
    cur = set(rel) if start is None else set(start) & set(rel)
    seen = []
    for step in range(1, limit + 1):
        nxt: Set[Pattern] = set()
        for p in cur:
            nxt |= rel.get(p, set())
        nxt &= set(rel)
        if not nxt:
            return step
        if nxt in seen:
            return None
        seen.append(nxt)
        cur = nxt
    return None


def congruent_pairs(graph: UnitDistanceGraph, W: Sequence[int],
                    limit: int = 200) -> List[Tuple[int, ...]]:
    """Vertex tuples with the same pairwise distances as W, in order.

    A chain needs W' congruent to W so that one isometry carries the copy's
    interface onto the next copy's, and the glued points really are the same
    points.
    """
    W = list(W)
    want = [[graph.vertices[a].dist2(graph.vertices[b]) for b in W] for a in W]
    out: List[Tuple[int, ...]] = []
    n = len(W)

    def extend(pick: List[int]) -> None:
        if len(out) >= limit:
            return
        i = len(pick)
        if i == n:
            if tuple(pick) != tuple(W):
                out.append(tuple(pick))
            return
        for v in range(graph.n):
            if v in pick:
                continue
            ok = True
            for j, u in enumerate(pick):
                if graph.vertices[v].dist2(graph.vertices[u]) != want[i][j]:
                    ok = False
                    break
            if ok:
                extend(pick + [v])
                if len(out) >= limit:
                    return

    extend([])
    return out
