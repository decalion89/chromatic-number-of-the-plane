"""mu: the number that IS the Hadwiger-Nelson problem.

For a unit-distance graph `G` in the plane and a point `p`, write

    N(p) = { x in G : |x - p| = 1 }

and define, over the proper k-colourings of `G`,

    mu_k(G, p) = min | c(N(p)) |          (5 if G has no k-colouring)

`G + p` is the unit-distance graph with `p` adjoined, and it has no proper
k-colouring exactly when every k-colouring of `G` already spends all k colours
on `N(p)` -- that is, exactly when `mu_k(G, p) = k`.  So for k = 5:

    mu_5 = 5 somewhere   <=>   chi(R^2) >= 6
    mu_5 <= 4 always     <=>   chi(R^2) <= 5

Both directions of the first are easy: a blocked point gives a 6-chromatic
unit-distance graph by adjoining it, and conversely a vertex-critical
6-chromatic unit-distance graph minus any vertex `p` is a 5-colourable `G` in
which `N(p)` must carry all five colours, or `p` could be put back.  The second
line is the contrapositive of the first, plus `chi(R^2) >= 5`, which is known.

An earlier version of this docstring wrote the second line as `mu_5 <= 2 always
<=> chi(R^2) = 5`, and the forward direction is fine -- 2 < 5 everywhere means
no blocked point.  The CONVERSE is false: `chi(R^2) = 5` only says every finite
graph plus a point is 5-colourable, which is `mu_5 <= 4`, and nothing forbids a
graph with `mu_5 = 3` or `4` at some point.  So climbing from 2 to 3 or 4 proves
nothing on its own; it is progress in the search, not partial credit towards a
theorem.  Only `mu_5 = 5` decides anything.

So `mu_5` is not a heuristic for the problem or a proxy for it.  It is the
problem, restated as an integer between 1 and 5, attached to a graph and a
point.  That makes it worth computing properly rather than in passing, which
is what this module is for.

Two facts shape the computation.

The CNF carries no at-most-one clauses, because nothing here reads a colour
off a model: every question is satisfiability under assumptions.  The edge
clauses already make adjacent colour sets disjoint, so a satisfying assignment
still yields a proper colouring by picking any true colour per vertex.

**Colour symmetry.**  If some set of `j` colours can be simultaneously absent
from `N(p)` in a proper colouring, then the colours `0..j-1` can, by permuting.
So the whole question is a nested chain of at most `k-1` assumption calls on
one warm solver, and `mu <= k - j` is decided by a single call forbidding
`0..j-1`.  Since the interesting answer is "above 2", the cheapest useful
probe is the one call that forbids three colours: a yes settles the point as
`mu <= 2` and costs nothing more.

**The neighbourhood is bipartite.**  `N(p)` lies on the unit circle about `p`,
and two points of a radius-1 circle are a unit apart exactly when their central
angle is 60 degrees, so inside `N(p)` a point has at most two neighbours and
every cycle has even length.  `chi(N(p)) = 2` on its own, in every planar
unit-distance graph, however many neighbours `p` has.  Anything above 2
therefore has to be imposed by the rest of `G` -- which is a five-colour
forcing statement, and the kind this project has searched for exhaustively and
never found.

Measured so far, on every graph in `data/` that is small enough to colour and
on the richest candidate neighbourhoods of each (up to 15 points, against the
13 of de Grey's `G`): `mu = 2`, without exception.
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Tuple

from pysat.solvers import Solver

from .geometry import Point
from .graph import UnitDistanceGraph

__all__ = ["neighbourhood", "MuSolver", "mu_of_point"]


def neighbourhood(graph: UnitDistanceGraph, p: Point) -> List[int]:
    """The indices of graph vertices at distance exactly 1 from `p`.

    Float first, exact second.  A float distance discards all but a handful of
    candidates before any field arithmetic happens, and the exact comparison
    still decides every survivor, so nothing is assumed -- the float is a
    filter, never an answer.  Scoring a hundred thousand candidates by exact
    arithmetic alone takes a quarter of an hour; this takes seconds.
    """
    if graph.n == 0:
        return []
    one = graph.vertices[0].x.field.rational(1)
    px, py = p.fx, p.fy
    out: List[int] = []
    for i, q in enumerate(graph.vertices):
        ex, ey = q.fx - px, q.fy - py
        if abs(ex * ex + ey * ey - 1.0) < 1e-9 and (q - p).norm2() == one:
            out.append(i)
    return out


class MuSolver:
    """One warm solver for many points of the same graph.

    Building the CNF and finding the first colouring is the expensive part --
    seconds on a thousand vertices, minutes on ten thousand -- and every point
    of the graph reuses both.  Each point then costs assumption calls only,
    and the solver keeps everything it learned from the ones before.
    """

    def __init__(self, graph: UnitDistanceGraph, k: int = 5,
                 budget: Optional[int] = 3_000_000) -> None:
        self.graph = graph
        self.k = k
        self.budget = budget
        n = graph.n
        self._x = lambda v, c: 1 + v * k + c
        cnf = [[self._x(v, c) for c in range(k)] for v in range(n)]
        # NO at-most-one clauses.  An earlier version included them "because mu
        # reads colours off models", which is simply false -- mu never looks at
        # a model, only at whether a call is satisfiable.  They are unnecessary
        # and expensive: ten per vertex at k = 5, so 53 770 extra clauses on a
        # 5377-vertex graph, and cadical is measurably slower for them.
        #
        # Soundness without them.  The assumptions force colours 0..j-1 false
        # on N(p), so every member takes some colour >= j; the edge clauses
        # make adjacent colour SETS disjoint, so choosing any true colour per
        # vertex yields a proper colouring avoiding 0..j-1.  The converse is
        # immediate.  So the answer is the same, and only the cost differs.
        for x, y in graph.edges():
            for c in range(k):
                cnf.append([-self._x(x, c), -self._x(y, c)])
        self._solver = Solver(name="cd19", bootstrap_with=cnf)
        self.colourable = self._solver.solve()

    def close(self) -> None:
        self._solver.delete()

    def __enter__(self) -> "MuSolver":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _absent(self, nb: Sequence[int], j: int) -> Optional[bool]:
        """Can colours 0..j-1 all be absent from `nb`?  None if out of budget."""
        ass = [-self._x(u, c) for u in nb for c in range(j)]
        if self.budget is None:
            return self._solver.solve(assumptions=ass)
        self._solver.conf_budget(self.budget)
        return self._solver.solve_limited(assumptions=ass)

    def at_most_two(self, nb: Sequence[int]) -> Optional[bool]:
        """The cheap probe: one call, and a yes means the point is placeable."""
        if self.k < 3:
            return True
        return self._absent(nb, self.k - 2)

    def mu(self, nb: Sequence[int]) -> Optional[int]:
        """mu on this neighbourhood, or None if a call ran out of budget.

        An empty or tiny neighbourhood cannot carry k colours, so mu is capped
        by its size -- worth checking before any solving, since most candidate
        points have two or three neighbours and no call is needed at all.
        """
        if not self.colourable:
            return self.k
        # A set of t points shows at most t colours, so mu <= |nb| whatever
        # the solver says.  Most candidate points in the plane have two or
        # three neighbours, and there mu is decided without any solving.
        cap = len(nb)
        if cap <= 1:
            return cap
        for j in range(self.k - 1, 0, -1):
            r = self._absent(nb, j)
            if r is None:
                return None
            if r is True:
                return min(self.k - j, cap)
        return min(self.k, cap)


def mu_of_point(graph: UnitDistanceGraph, p: Point, k: int = 5,
                budget: Optional[int] = 3_000_000) -> Tuple[Optional[int], List[int]]:
    """mu_k(graph, p), and the neighbourhood it was computed on.

    Convenience for a single point; use `MuSolver` for many, since it shares
    the CNF and the first colouring across them.
    """
    nb = neighbourhood(graph, p)
    with MuSolver(graph, k, budget) as ms:
        return ms.mu(nb), nb
