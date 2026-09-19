"""Which pairs of points must differ, and which must agree, in every colouring.

The two-orbit block is k-independent. Its six copies rule out the choice
pattern by geometry alone -- the a-orbit of a 1/3 circle is a triangle, so at
most one copy per orbit takes the first leg, and the second rotation pairs the
b-images off, so at most three take the second; four are needed and six copies
close it whatever k is. So chi(R^2) >= 6 needs exactly one thing this package
does not have: a pivot whose core at five colours is a genuinely joint pair,
with one leg at squared distance 1/3.

Searching for that inside existing graphs found nothing, across 880 pivots.
The one that worked at three colours was built, not found, out of a triangle:
three points that take three different colours in every 3-colouring, with the
pivot adjacent to one and the other two as the legs. The same recipe at five
colours needs five points taking five different colours in every 5-colouring.

That is the object this module looks for. Call a set of vertices a *rainbow*
if every colouring gives its members pairwise different colours; a k-rainbow
in a k-colourable graph is what plays the triangle's role. Since the maximum
clique of a plane unit-distance graph is three, a rainbow of five cannot be
made of edges: at least some of its pairs must be forced apart at a distance
other than one, which is a real question about the graph and not about the
plane.

The queries are cheap. The plain colouring formula is symmetric under
permuting colours, so

    u and v can share a colour   <=>  x[u][0] and x[v][0] is satisfiable
    u and v can differ           <=>  x[u][0] and x[v][1] is satisfiable

and each is one incremental solve under two assumptions, on a formula built
once. Forced-same pairs matter as much as forced-different ones: a pair forced
to agree, followed by one unit step, is a pair forced to differ -- which is how
a rainbow gets built at a distance the plane does not hand you.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .graph import UnitDistanceGraph

__all__ = ["ColourRelations", "forced_different_pairs", "extend_rainbow",
           "find_rainbow", "min_colours_on", "pressure", "max_core_bound",
           "centre_pressure"]


class ColourRelations:
    """Forced-same and forced-different queries on one persistent solver.

    Vacuous unless the graph is k-colourable: with no colouring at all every
    query is UNSAT and every pair comes back forced both ways. `colourable`
    records the answer so callers can refuse to report nonsense.
    """

    def __init__(self, graph: UnitDistanceGraph, k: int):
        from pysat.solvers import Solver

        self.graph, self.k = graph, k
        cls = [[self._x(v, c) for c in range(k)] for v in range(graph.n)]
        for u, v in graph.edges():
            for c in range(k):
                cls.append([-self._x(u, c), -self._x(v, c)])
        self._s = Solver(name="cd15", bootstrap_with=cls)
        self.colourable = self._s.solve()
        self.calls = 0

    def _x(self, v: int, c: int) -> int:
        return 1 + v * self.k + c

    def _ask(self, u: int, cu: int, v: int, cv: int) -> bool:
        self.calls += 1
        return self._s.solve(assumptions=[self._x(u, cu), self._x(v, cv)])

    def can_share(self, u: int, v: int) -> bool:
        return self._ask(u, 0, v, 0)

    def can_differ(self, u: int, v: int) -> bool:
        return self._ask(u, 0, v, 1)

    def different(self, u: int, v: int) -> bool:
        """True when no colouring gives u and v the same colour."""
        return not self.can_share(u, v)

    def same(self, u: int, v: int) -> bool:
        """True when no colouring gives u and v different colours."""
        return not self.can_differ(u, v)

    def close(self) -> None:
        self._s.delete()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def forced_different_pairs(rel: ColourRelations, u: int,
                           candidates: Iterable[int]) -> List[int]:
    """Every candidate that cannot share u's colour. Edges are free."""
    adj = rel.graph.adj[u]
    return [v for v in candidates
            if v != u and (v in adj or rel.different(u, v))]


def extend_rainbow(rel: ColourRelations, seed: Sequence[int],
                   candidates: Iterable[int]) -> List[int]:
    """Candidates forced different from every member of the seed."""
    out = list(candidates)
    for s in seed:
        out = forced_different_pairs(rel, s, out)
        if not out:
            break
    return [v for v in out if v not in set(seed)]


def find_rainbow(graph: UnitDistanceGraph, k: int, size: int,
                 seeds: Optional[Sequence[Sequence[int]]] = None,
                 candidates: Optional[Sequence[int]] = None,
                 report=None) -> Optional[List[int]]:
    """A set of `size` vertices pairwise forced to differ, or None.

    Triangles are the natural seeds: three vertices already pairwise forced
    apart, for free, by adjacency. Growing from them turns an n-choose-5 search
    into one pass over the triangles and a small clique problem on what each
    one leaves.
    """
    rel = ColourRelations(graph, k)
    try:
        if not rel.colourable:
            raise ValueError(f"graph is not {k}-colourable; every query is "
                             f"vacuously UNSAT")
        pool = list(range(graph.n) if candidates is None else candidates)
        if seeds is None:
            seeds = [[a, b, c] for a in range(graph.n)
                     for b in graph.adj[a] if b > a
                     for c in graph.adj[a] & graph.adj[b] if c > b]
        for i, seed in enumerate(seeds):
            if len(seed) >= size:
                return list(seed)[:size]
            rest = extend_rainbow(rel, seed, pool)
            if report:
                report(i, seed, rest, rel.calls)
            need = size - len(seed)
            if len(rest) < need:
                continue
            found = _clique(rel, rest, need)
            if found is not None:
                return list(seed) + found
        return None
    finally:
        rel.close()


def _clique(rel: ColourRelations, pool: Sequence[int],
            need: int) -> Optional[List[int]]:
    """A forced-different clique of the given size inside the pool."""
    if need == 0:
        return []
    if need == 1:
        return [pool[0]] if pool else None
    for i, u in enumerate(pool):
        rest = forced_different_pairs(rel, u, pool[i + 1:])
        if len(rest) >= need - 1:
            sub = _clique(rel, rest, need - 1)
            if sub is not None:
                return [u] + sub
    return None


# -- colour pressure, and what it forbids ---------------------------------
#
# Let W be a k-colourable unit-distance graph and p one of its vertices. The
# only thing that ever restricts c(p) is its neighbourhood, so the colours p
# can take are exactly the complement of c(N(p)), and the useful invariant is
#
#     pressure(p) = min |c(N(p))|   over all k-colourings of W,
#
# which `min_colours_on` computes in at most k incremental solves.
#
# THEOREM. If p has a core of size r -- a set T of r vertices with
# c(p) in c(T) in every k-colouring -- then pressure(p) >= k - r.
#
# Proof. Take a colouring with |c(N(p))| = pressure(p). At least
# k - pressure(p) colours are free for p, while c(T) offers at most r values.
# If k - pressure(p) > r, some free colour lies outside c(T); recolouring p to
# it keeps the colouring proper, because properness at p asks only that its
# colour avoid c(N(p)), and it leaves every other vertex, T included, alone.
# That contradicts T being a core. []
#
# The contrapositive is what makes it useful: a pivot whose neighbourhood can
# be squeezed into k - r - 1 colours has no core of size r, and no search will
# ever find one there. Two corollaries decide this package's two regimes.
#
#   A core of one -- the classical spindle -- needs pressure k-1, so c(p) is
#   determined outright. A core of two -- what the two-orbit block closes --
#   needs pressure k-2, which at five colours means three.
#
# Measured on de Grey's G at k = 5, pressure is exactly 2 at all 1581
# vertices, the two hubs of degree 60 included. So no vertex of G has a core
# of size two at five colours, and none has a core of size one either. That is
# a proof of what 880 pivots of searching reported one negative at a time.
#
# Where pressure can come from is equally sharp. Two points of a unit circle
# are one apart exactly when they subtend 60 degrees, and the 60-degree orbit
# closes after six steps, so every component of a circle's own graph is a path
# or a 6-cycle: bipartite. A unit circle therefore never needs a third colour
# on its own account, and any pressure above two is entirely the ambient
# graph's doing -- it has to rule out every 2-colouring of the circle, which
# is a joint condition on the whole circle and not a statement about any one
# pair. Radius 1 is what is special here, not circles in general: at radius
# 1/sqrt(3) the chord-1 angle is 120 degrees and the orbit is a triangle,
# odd, which is the leg the block needs and the one such distance Niven's
# theorem leaves on a rational radius.
#
# SECOND THEOREM (criticality blocks forcing). If chi(G - u) < chi(G) = k,
# then u belongs to no forced pair at k colours.
#
# Proof, both directions at once. Take a (k-1)-colouring of G - u on colours
# 1..k-1 and give u the colour k, which is proper and makes u the only vertex
# carrying it; so c(u) differs from every c(v), and no pair through u is
# forced same. For the other direction let v be non-adjacent to u and recolour
# v to k as well: every neighbour of v lies in G - u, since u is not one of
# them, so all of them still carry colours 1..k-1 and the colouring stays
# proper. Now c(u) = c(v), so no non-adjacent pair through u is forced
# different either. []
#
# The same colouring proves more, and this is the sharpest form. In it p is
# the ONLY vertex carrying the kth colour, so c(p) lies outside c(T) for every
# set T at once:
#
#   COROLLARY. A k-vertex-critical graph has no core of any size, at any
#   vertex. Not one, not two, not thirty-four.
#
# de Grey's G is 5-vertex-critical, so it is completely inert at five colours:
# every pivot is separable from every set of targets simultaneously, and no
# search over it can ever return a forced anything. That accounts for all 880
# pivots, in one line, before any solver runs.
#
# It also says exactly what a candidate graph must not be. W = G union f(G)
# escapes: removing one vertex still leaves a whole 5-chromatic copy, so W is
# not vertex-critical, no pivot can be given a colour of its own, and the full
# target set IS a core. How far that core shrinks is then the real question,
# and pressure answers it from below -- measured at 2 on these unions, so
# their smallest possible core is three, never two. Between them the two theorems
# account for every negative this package recorded at k = 5, and they say what
# a construction has to look like: the forced vertices must be ones whose
# removal leaves the graph still k-chromatic, which for a graph built on top
# of G means points added outside it. That is what the Moser spindle does at
# three colours and what de Grey's construction does at four.
#
# Note what the first theorem does NOT give. Pressure k-2 says c(p) is
# confined to two colours; it does not say which two, and the two-orbit block
# needs them named by fixed vertices y and z. Pressure >= k-2 is necessary for
# a joint core of two, not sufficient, and `hn.spindle.SeparationTest` is
# still what decides the naming.

def min_colours_on(rel: "ColourRelations", vertices: Sequence[int]) -> int:
    """The least number of colours c(S) can use, over all k-colourings.

    Returns 0 for the empty set. By colour symmetry only the initial segment
    of each size needs testing, so this is at most k solves.
    """
    S = list(dict.fromkeys(vertices))
    if not S:
        return 0
    for t in range(1, rel.k + 1):
        rel.calls += 1
        drop = [-rel._x(s, c) for s in S for c in range(t, rel.k)]
        if rel._s.solve(assumptions=drop):
            return t
    return rel.k + 1        # unreachable while the graph is k-colourable


def pressure(rel: "ColourRelations", p: int) -> int:
    """min |c(N(p))| over all k-colourings: what p's neighbourhood forces."""
    return min_colours_on(rel, sorted(rel.graph.adj[p]))


def max_core_bound(rel: "ColourRelations", p: int) -> int:
    """The largest core size the theorem still permits at p.

    A core of size r needs pressure(p) >= k - r, so no core smaller than
    k - pressure(p) can exist. Returned as that lower bound: 5 - 2 = 3 on
    every vertex of de Grey's G, which is why no pair there is ever forced.
    """
    return rel.k - pressure(rel, p)


def centre_pressure(rel: "ColourRelations", vertices: Sequence[int]) -> dict:
    """What a centre with this unit-circle set rules in and out.

    `smallest_possible_core` is the theorem's bound; `verdict` reads it. Note
    that reaching k-2 only makes a joint core of two possible -- the pair
    still has to be named, and that is a separation query, not this one.
    """
    t = min_colours_on(rel, vertices)
    return {
        "circle": len(set(vertices)),
        "pressure": t,
        "free_for_centre": rel.k - t,
        "smallest_possible_core": rel.k - t,
        "verdict": ("uncolourable: chi >= k+1" if t >= rel.k else
                    "centre forced; a core of one is possible" if t == rel.k - 1
                    else "a joint core of two is possible" if t == rel.k - 2
                    else f"no core smaller than {rel.k - t}"),
    }
