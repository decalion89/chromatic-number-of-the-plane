"""Independence ratio: a finite certificate for a measurable lower bound.

Every search in this package looks for one binary fact -- a graph with no
proper 5-colouring -- and until it appears there is nothing to show. There is
a continuous quantity that gives the same conclusion and that any finite graph
reports a value for.

Let S be a measurable set avoiding distance 1, of upper density d, and let G
be a finite unit-distance graph on n vertices. Average over rigid motions
sigma of the plane: the expected size of sigma(V) ∩ S is n d, and that
intersection is an independent set of a copy of G, so it never exceeds
alpha(G). Hence

    m_1(R^2) <= alpha(G) / n

for every finite unit-distance graph G, where m_1 is the largest density a
1-avoiding set can have. And if five measurable classes cover the plane one of
them has density at least 1/5, so

    alpha(G) / n < 1/5   for a single finite G   =>   chi_m(R^2) >= 6.

That is a certificate, not a search: one graph, one independent-set
computation, and the number either clears 0.2 or does not. The Moser spindle
gives 2/7 = 0.2857. The best known bound on m_1 is about 0.2470, from Fourier
methods rather than from a graph, so the graph route has room in it.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from .graph import UnitDistanceGraph

DEFAULT_SOLVER = "cd19"


def _independent_at_least(graph: UnitDistanceGraph, size: int, pool: IDPool):
    """CNF asserting an independent set of at least `size` vertices."""
    cnf = CNF()
    xs = [pool.id(("v", v)) for v in range(graph.n)]
    for u, v in graph.edges():
        cnf.append([-xs[u], -xs[v]])
    card = CardEnc.atleast(lits=xs, bound=size, vpool=pool,
                           encoding=EncType.seqcounter)
    cnf.extend(card.clauses)
    return cnf, xs


def has_independent_set(graph: UnitDistanceGraph, size: int,
                        solver: str = DEFAULT_SOLVER):
    """Is there an independent set this large?  Returns it, or None."""
    pool = IDPool()
    cnf, xs = _independent_at_least(graph, size, pool)
    with Solver(name=solver, bootstrap_with=cnf) as s:
        if not s.solve():
            return None
        model = set(l for l in s.get_model() if l > 0)
        return [v for v in range(graph.n) if xs[v] in model]


def independence_number(graph: UnitDistanceGraph, lo: int = 0,
                        hi: Optional[int] = None,
                        solver: str = DEFAULT_SOLVER) -> Tuple[int, List[int]]:
    """alpha(G) exactly, by binary search on the cardinality bound."""
    if hi is None:
        hi = graph.n
    best: List[int] = []
    lo = max(lo, 0)
    while lo <= hi:
        mid = (lo + hi) // 2
        found = has_independent_set(graph, mid, solver)
        if found is None:
            hi = mid - 1
        else:
            best = found
            lo = mid + 1
    return len(best), best


def independence_ratio(graph: UnitDistanceGraph, **kw) -> float:
    a, _ = independence_number(graph, **kw)
    return a / graph.n


def measurable_bound(graph: UnitDistanceGraph, **kw) -> float:
    """The lower bound on the measurable chromatic number this graph gives.

    n / alpha(G), since a class of density d needs alpha/n >= d and k classes
    covering the plane force some d >= 1/k.
    """
    a, _ = independence_number(graph, **kw)
    return graph.n / a if a else float("inf")
