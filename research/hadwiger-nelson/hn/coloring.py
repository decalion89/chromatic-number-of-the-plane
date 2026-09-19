"""Deciding k-colourability with SAT, and shrinking the witness.

Encoding (the standard direct encoding):
  x[v,c]  -- vertex v may take colour c
  (1) every active vertex takes at least one colour
  (2) for every edge and every colour, not both endpoints take it

At-most-one-colour clauses are deliberately omitted.  They are unnecessary:
if the colour sets of adjacent vertices are disjoint, picking any colour from
each vertex's set yields a proper colouring.  Leaving them out shrinks the
formula by k(k-1)/2 clauses per vertex.

Vertex selectors a[v] make the formula *reusable*: solving under assumptions
{a[v] : v in S} decides the induced subgraph on S, so one formula answers
every subgraph question and the solver's UNSAT core hands back a set of
vertices that is already non-k-colourable.
"""

from __future__ import annotations

import threading
import time
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from pysat.formula import CNF
from pysat.solvers import Solver

from .graph import UnitDistanceGraph

__all__ = ["ColoringInstance", "is_k_colorable", "find_uncolorable_core", "chromatic_number"]

DEFAULT_SOLVER = "cd19"  # CaDiCaL 1.9.5


class ColoringInstance:
    """A reusable CNF for 'is the induced subgraph on S k-colourable?'."""

    def __init__(self, graph: UnitDistanceGraph, k: int, symmetry_break: bool = True):
        self.graph = graph
        self.k = k
        n = graph.n
        self.nv = n
        # x[v,c] -> 1 + v*k + c ;  a[v] -> 1 + n*k + v
        self.cnf = CNF()
        self._clique: List[int] = []
        for v in range(n):
            self.cnf.append([-self.a(v)] + [self.x(v, c) for c in range(self.k)])
        for u, v in graph.edges():
            for c in range(self.k):
                self.cnf.append([-self.x(u, c), -self.x(v, c)])
        if symmetry_break:
            # A unit-distance graph in the plane has no K4 (four mutually
            # unit-distant points do not exist), so the largest clique we can
            # pin is a triangle.  Pinning it removes k!/(k-3)! of the colour
            # permutations for free.
            clique = graph.find_clique(min(3, self.k))
            if clique:
                self._clique = clique
                for i, v in enumerate(clique):
                    self.cnf.append([self.x(v, i)])
                    for c in range(self.k):
                        if c != i:
                            self.cnf.append([-self.x(v, c)])

    def x(self, v: int, c: int) -> int:
        return 1 + v * self.k + c

    def a(self, v: int) -> int:
        return 1 + self.nv * self.k + v

    def decode(self, model: Sequence[int]) -> List[int]:
        pos = set(l for l in model if l > 0)
        colors = []
        for v in range(self.nv):
            chosen = -1
            for c in range(self.k):
                if self.x(v, c) in pos:
                    chosen = c
                    break
            colors.append(chosen)
        return colors

    # -- solving ---------------------------------------------------------
    def solve(
        self,
        subset: Optional[Iterable[int]] = None,
        timeout: Optional[float] = None,
        solver: str = DEFAULT_SOLVER,
        with_proof: bool = False,
    ) -> Tuple[Optional[bool], Optional[List[int]], Optional[List[int]], Optional[List[str]]]:
        """Return (sat, colouring, core, proof).

        sat is None on timeout.  `core` is a vertex subset that is already
        non-k-colourable.  `proof` is DRAT, only when with_proof and UNSAT.
        """
        verts = list(range(self.nv)) if subset is None else sorted(set(subset))
        assumptions = [self.a(v) for v in verts]
        # forbid the excluded vertices so their variables cannot carry weight
        s = Solver(name=solver, bootstrap_with=self.cnf, with_proof=with_proof)
        timer = None
        if timeout:
            timer = threading.Timer(timeout, s.interrupt)
            timer.daemon = True
            timer.start()
        try:
            if timeout:
                res = s.solve_limited(assumptions=assumptions, expect_interrupt=True)
            else:
                res = s.solve(assumptions=assumptions)
            if res is None:
                return None, None, None, None
            if res:
                return True, self.decode(s.get_model()), None, None
            core_lits = s.get_core() or assumptions
            core = sorted(l - 1 - self.nv * self.k for l in core_lits)
            proof = s.get_proof() if with_proof else None
            return False, None, core, proof
        finally:
            if timer:
                timer.cancel()
            s.delete()


def is_k_colorable(
    graph: UnitDistanceGraph, k: int, timeout: Optional[float] = None
) -> Tuple[Optional[bool], Optional[List[int]]]:
    """Reduce by k-core, then decide.  Returns (answer, colouring)."""
    core = graph.k_core(k)
    if core.n == 0:
        return True, graph.greedy_coloring()
    inst = ColoringInstance(core, k)
    sat, coloring, _, _ = inst.solve(timeout=timeout)
    if sat is not True:
        return sat, None
    # extend the core colouring back to the whole graph (degree < k vertices
    # are coloured last, in reverse removal order)
    return True, _extend_to_full(graph, core, coloring, k)


def _extend_to_full(graph, core, core_coloring, k) -> List[int]:
    idx = {}
    for i, p in enumerate(core.vertices):
        idx[p] = i
    full = [-1] * graph.n
    pending = []
    for v, p in enumerate(graph.vertices):
        if p in idx:
            full[v] = core_coloring[idx[p]]
        else:
            pending.append(v)
    changed = True
    while changed:
        changed = False
        for v in list(pending):
            used = {full[w] for w in graph.adj[v] if full[w] >= 0}
            free = [c for c in range(k) if c not in used]
            if free:
                full[v] = free[0]
                pending.remove(v)
                changed = True
    return full


def chromatic_number(
    graph: UnitDistanceGraph, lo: int = 1, hi: int = 8, timeout: Optional[float] = None
) -> Tuple[Optional[int], Optional[List[int]]]:
    """Smallest k for which the graph is k-colourable."""
    best_coloring = None
    for k in range(lo, hi + 1):
        sat, coloring = is_k_colorable(graph, k, timeout=timeout)
        if sat is None:
            return None, None
        if sat:
            return k, coloring
        best_coloring = None
    return None, best_coloring


def find_uncolorable_core(
    graph: UnitDistanceGraph,
    k: int,
    timeout: Optional[float] = None,
    rounds: int = 60,
    greedy: bool = True,
    verbose: bool = True,
) -> Optional[UnitDistanceGraph]:
    """Shrink `graph` to a small subgraph that is still not k-colourable.

    Two phases.  First, iterated UNSAT cores: re-solving restricted to the
    previous core usually yields a strictly smaller one, and this converges
    fast.  Second, greedy vertex deletion: try dropping each remaining vertex
    and keep the deletion whenever the rest is still uncolourable.  The result
    is vertex-minimal with respect to single deletions.
    """
    work = graph.k_core(k)
    if work.n == 0:
        return None
    inst = ColoringInstance(work, k, symmetry_break=False)
    sat, _, core, _ = inst.solve(timeout=timeout)
    if sat is not False:
        return None
    subset = set(core)
    for r in range(rounds):
        sat, _, core, _ = inst.solve(subset=subset, timeout=timeout)
        if sat is not False:
            break
        if len(core) >= len(subset):
            break
        subset = set(core)
        if verbose:
            print(f"    core round {r+1}: {len(subset)} vertices")

    if greedy:
        order = sorted(subset, key=lambda v: len(work.adj[v]))
        for v in order:
            if v not in subset:
                continue
            trial = subset - {v}
            sat, _, core, _ = inst.solve(subset=trial, timeout=timeout)
            if sat is False:
                subset = set(core) if len(core) < len(trial) else trial
        if verbose:
            print(f"    after greedy deletion: {len(subset)} vertices")

    sub = work.induced(subset)
    # k-core again: greedy deletion can strand low-degree vertices
    sub = sub.k_core(k).largest_component()
    return sub
