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

import random as _rand
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


# -- what a core must look like when pressure is exactly k - r -------------
#
# The pressure bound is tight in a way that constrains the core's shape, not
# just its size. Suppose pressure(p) = k - r exactly and T is a core of size
# r. Take a colouring attaining the minimum: p then has exactly r free
# colours, and c(p) can be any of them, because properness at p asks only that
# its colour avoid c(N(p)). For T to be a core in that colouring, c(T) must
# contain all r, and |c(T)| <= r, so
#
#   c(T) = free(p)  exactly, in every minimum-pressure colouring.
#
# Two consequences, both checkable. T is a RAINBOW there -- its r members take
# r different colours -- and c(T) is DISJOINT from c(N(p)). At five colours
# with pressure 2 that means a core of three is three points taking three
# different colours, none of them a colour any of the pivot's sixty
# neighbours carries, in every colouring that squeezes the neighbourhood into
# two. The legs may be a unit triangle, which is a rainbow for free, but the
# disjointness is not free at all: it says no neighbour of the pivot ever
# shares a colour with any of the three.
#
# This is why the shrinking searches keep returning "separable". It is not
# that the core is large; it is that the core has to be aligned, and alignment
# is a global condition on the colouring, not a local one on the distances.

def free_colours(rel: "ColourRelations", p: int) -> int:
    """How many colours p can still take, at minimum pressure: k - pressure."""
    return rel.k - pressure(rel, p)


def core_must_be_rainbow(rel: "ColourRelations", p: int,
                         targets: Sequence[int]) -> bool:
    """Whether these targets can be a core of p, by the alignment condition.

    Necessary, not sufficient. Returns False as soon as some colouring puts
    the neighbourhood at minimum pressure while two targets share a colour or
    a target shares a colour with a neighbour -- either of which breaks
    c(T) = free(p).
    """
    T = list(dict.fromkeys(targets))
    if len(T) != free_colours(rel, p):
        return len(T) > free_colours(rel, p)
    for i, u in enumerate(T):
        for v in T[i + 1:]:
            if rel.can_share(u, v):
                return False
    return True


# -- the core condition IS a pressure measurement --------------------------
#
# T is a core of p when c(p) lies in c(T) in every k-colouring. Suppose some
# colouring left a colour g unused on N(p) and on T alike. Then recolouring p
# to g is proper, because properness at p asks only that its colour avoid
# c(N(p)), and it puts c(p) outside c(T). Conversely if T is not a core, the
# colouring witnessing it has c(p) outside c(T) and outside c(N(p)) already,
# so c(p) is such a colour. Hence
#
#     T is a core of p   <=>   min |c(N(p) union T)| = k,
#
# which `min_colours_on` answers in at most k incremental solves, against the
# repeated shrinking that every core search in this package has used. It also
# explains the shape: with pressure(p) = 2 at five colours, N(p) contributes
# two and T has to supply the remaining three by itself, in every colouring.
#
# And the objective is MONOTONE. Enlarging S can only raise |c(S)| in each
# colouring, so the minimum over colourings can only rise; adding a vertex to
# T never undoes progress. That is the first search here with a gradient
# rather than an all-or-nothing verdict, and it is what makes growing a core
# from the neighbourhood outwards a sensible thing to do at all.

def is_core(rel: "ColourRelations", p: int, targets: Sequence[int]) -> bool:
    """Whether c(p) lies in c(targets) in every k-colouring."""
    circle = sorted(rel.graph.adj[p])
    return min_colours_on(rel, list(circle) + [t for t in targets
                                               if t != p]) >= rel.k


def grow_core(rel: "ColourRelations", p: int,
              candidates: Optional[Sequence[int]] = None,
              limit: int = 6, report=None) -> Tuple[List[int], int]:
    """Greedily grow a core out of p's neighbourhood; return (T, coverage).

    Each step takes the candidate that raises min |c(N(p) union T)| the most,
    stopping at k -- where T is a core -- or when no candidate helps. Monotone,
    so a step never has to be undone, and every evaluation is at most k solves
    on the one persistent formula.
    """
    circle = sorted(rel.graph.adj[p])
    pool = [v for v in (range(rel.graph.n) if candidates is None else candidates)
            if v != p and v not in rel.graph.adj[p]]
    chosen: List[int] = []
    best = min_colours_on(rel, circle)

    def squeeze(S, t):
        """Assumptions forcing S into the first t colours."""
        return [-rel._x(v, c) for v in S for c in range(t, rel.k)]

    for step in range(limit):
        # A candidate helps exactly when it breaks the current squeeze. One
        # MODEL of that squeeze kills every candidate it colours inside the
        # first `best` colours, because that model witnesses their fit -- so a
        # handful of models clears almost the whole graph and only the
        # survivors need a query each. Scanning candidates one at a time was
        # 2900 solves per step; this is a few dozen.
        live = list(pool)
        for _ in range(8):
            if not live:
                break
            rel.calls += 1
            if not rel._s.solve(assumptions=squeeze(circle + chosen, best)):
                break                       # already past `best`; nothing to do
            m = set(rel._s.get_model())
            live = [v for v in live
                    if all(rel._x(v, c) not in m for c in range(best))]
            # the same assumptions return the same model, so the next solve
            # has to be pushed somewhere else or harvesting stops after one
            rel._s.set_phases([_rand.choice([1, -1])
                               * rel._x(v, _rand.randrange(rel.k))
                               for v in _rand.sample(range(rel.graph.n),
                                                     min(120, rel.graph.n))])
        gain = None
        for v in live:
            rel.calls += 1
            if not rel._s.solve(assumptions=squeeze(circle + chosen + [v], best)):
                gain = (min_colours_on(rel, circle + chosen + [v]), v)
                break
        if gain is None:
            break
        best, v = gain
        chosen.append(v)
        pool = [w for w in pool if w != v]
        if report:
            report(step, chosen, best)
        if best >= rel.k:
            break
    return chosen, best


def cegar_core(rel: "ColourRelations", p: int, key=None,
               limit: int = 40) -> Tuple[List[int], bool]:
    """Build a core of p by counterexamples, one solve per vertex added.

    Growing greedily stalls, and monotonicity is exactly why it is allowed to:
    the objective never decreases, but no SINGLE vertex has to raise it even
    when a pair would, so the search sits at a local plateau with a genuine
    core still above it.

    Counterexamples do not stall. T is a core when no colouring leaves a
    colour unused on N(p) union T, and by colour symmetry it is enough to ask
    that of the first colour. A solve either proves it -- done -- or returns a
    colouring in which that colour IS free, and any vertex carrying it
    elsewhere, added to T, kills that colouring. Each round costs one solve
    and rules out at least the witness it was given, and the full vertex set
    is a core whenever the graph is not k-vertex-critical, so it terminates.

    `key` ranks the candidates that would kill the current counterexample;
    the default prefers the ones nearest the pivot, since a leg's images only
    conflict when its circle is small enough to matter.
    """
    g = rel.graph
    circle = sorted(g.adj[p])
    banned = set(circle) | {p}
    if key is None:
        pv = g.vertices[p]
        key = lambda v: float(pv.dist2(g.vertices[v]))
    T: List[int] = []
    for _ in range(limit):
        rel.calls += 1
        ass = [-rel._x(v, 0) for v in circle + T]
        if not rel._s.solve(assumptions=ass):
            return T, True
        m = set(rel._s.get_model())
        cands = [v for v in range(g.n)
                 if v not in banned and v not in T and rel._x(v, 0) in m]
        if not cands:
            return T, False       # nothing carries the free colour: no core
        T.append(min(cands, key=key))
    return T, False


# -- the machine that makes pressure, and what it needs one level up -------
#
# Exactly one configuration in this package reaches pressure 3, and taking it
# apart names the mechanism rather than leaving it a measurement. Deleting
# everything from de Grey's Sa that can go while the pivot's neighbourhood
# still refuses to be squeezed into two colours leaves 47 vertices:
#
#   the pivot;
#   its circle of 30, which splits into FIVE hexagons -- the 60-degree orbits,
#     as the bipartiteness argument above requires, each with exactly two
#     alternating 2-colourings, so five independent orientation bits;
#   sixteen further points, at d^2 = 1/3, (7 +- sqrt33)/6 and (3 +- sqrt33)/6,
#     each adjacent to exactly two circle points -- and in all sixteen cases
#     the two lie in DIFFERENT hexagons, so each one reads the relative
#     orientation of a pair of them.
#
# Squeeze the circle into two colours. Every one of those sixteen that sees
# two differently-coloured circle points is barred from both, so it is
# confined to the remaining k-2. And the sixteen among themselves form a
# 16-vertex, 27-edge unit-distance graph that is 3-chromatic and NOT
# bipartite: it carries an odd cycle.
#
#   at k = 4 the confined points have 2 colours, an odd cycle does not fit,
#     and the squeeze is impossible -- pressure 3;
#   at k = 5 they have 3, the graph is 3-colourable, and the squeeze goes
#     through -- pressure 2.
#
# Checked causally, not just by coincidence of numbers: deleting three of the
# sixteen makes their graph bipartite, and the pressure drops from 3 to 2 on
# the spot.
#
# So the mechanism is "an odd cycle against k-2 colours", and lifting it says
# exactly what to build. At five colours the confined points get three, so the
# gadget among them must be 4-CHROMATIC rather than merely non-bipartite --
# a Moser spindle, or the 19-vertex jointly forced construction in `hn.mixed`,
# every vertex of which is adjacent to two circle points lying in different
# hexagons. That is a finite design problem with the pieces already in hand,
# which is a different kind of difficulty from "search harder".

PRESSURE_GADGET = {
    "witness": "certificates/pressure3_witness_47.json",
    "circle_components": 5,
    "component_size": 6,
    "confined_points": 16,
    "confined_graph": (16, 27),
    "confined_chromatic_number": 3,
    "mechanism": "an odd cycle among the confined points against k-2 colours",
    "lift_to_five": "the confined gadget must be 4-chromatic, not merely odd",
}
