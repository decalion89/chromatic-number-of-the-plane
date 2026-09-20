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

import collections
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
# RETRACTED, and this one mattered. This package asserted throughout that de
# Grey's G is 5-vertex-critical, and used the corollary to explain why nothing
# is ever forced on it at five colours. The assertion is FALSE.
#
#     G - 1420 is still 5-chromatic.
#
# Vertex 1420 has degree 4. Removing it and asking for a 4-colouring with a
# triangle pinned to 0, 1, 2 comes back UNSAT after 1581 seconds. So a proper
# subgraph of G on 1580 vertices is already 5-chromatic, G is not
# vertex-critical, and the corollary does not apply to it.
#
# What was offered as evidence before was separability -- which is a
# CONSEQUENCE of criticality, not a proof of it -- and that was the error: a
# consequence was read backwards. The warning sign was there to be noticed, in
# that others have published 5-chromatic unit-distance graphs of several
# hundred vertices.
#
# The theorem above is unaffected: a k-vertex-critical graph does have no core
# of any size. It simply says nothing about G.
#
# What survives, because it was measured rather than inferred: pressure is
# exactly 2 at all 1581 vertices; no forced pair turned up in 29930 queries;
# no forced-different non-edge in 40539 pairs; and 1201 vertices -- the union
# of the 133 highest-degree closed neighbourhoods -- are not rainbow-forcing.
# Those results stand, and they are now UNEXPLAINED: criticality was the
# explanation and it is gone.
#
# What changes: rho(G,5) <= 1580, since G - 1420 is a 5-chromatic subgraph and
# every 5-chromatic subgraph is rainbow-forcing. The claim that rho = n on G
# rested on criticality and goes with it.
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
# The mechanism is LIST COLOURING, and that is the whole of it.
#
# Squeezing the circle into two colours hands each auxiliary point a list --
# the colours its circle neighbours leave it. Seeing two differently-coloured
# circle points gives a list of k-2; seeing two of the same colour gives k-1.
# So the auxiliary graph does not face a single palette, it faces a list
# assignment, and the squeeze is refused exactly when that assignment admits
# no proper colouring.
#
# Measured over all 32 orientations of the five hexagons:
#
#   k = 4: lists of size 2 and 3, and the auxiliary graph is NOT
#          list-colourable in 32 of 32 orientations -- pressure 3;
#   k = 5: the same lists become size 3 and 4, and it IS list-colourable in
#          all 32 -- pressure 2.
#
# The odd cycle among the confined points is the special case where every
# list is the same pair, and it is load-bearing where it applies: deleting
# three of the sixteen makes their graph bipartite and the pressure drops from
# 3 to 2 at once. But it covers only half the orientations. In the other
# sixteen the confined set is bipartite, and the minimal refusal there is 22
# vertices using ELEVEN auxiliaries of which only six are confined -- the
# other five carry lists of size 3 and sit at d^2 = 1/3. A frame with only
# "confined" and "free" in it cannot express that; lists can.
#
# The lift is now an exact and standard question. At five colours the
# auxiliary graph must fail to be list-colourable with lists of sizes 3 and 4
# -- a choosability question, not a chromatic one. A 4-chromatic subgraph all
# of whose vertices carry the SAME list of three is one sufficient way, which
# is why the Moser spindle and the 19-vertex jointly forced construction are
# the pieces to try; they are not the only way.

PRESSURE_GADGET = {
    "witness": "certificates/pressure3_witness_47.json",
    "circle_components": 5,
    "component_size": 6,
    "confined_points": 16,
    "confined_graph": (16, 27),
    "confined_chromatic_number": 3,
    "mechanism": ("the auxiliary graph is not list-colourable with the lists "
                  "the squeezed circle hands it: k-2 where it sees two "
                  "differently-coloured circle points, k-1 where it sees two "
                  "of the same"),
    "orientations_refused_at_k4": 32,
    "orientations_refused_at_k5": 0,
    "lift_to_five": ("an auxiliary graph that is not list-colourable with "
                     "lists of sizes 3 and 4 -- choosability, not chromatic "
                     "number"),
}


def circle_hexagons(graph, pivot: int):
    """The pivot's circle split into its 60-degree orbits, with parities.

    Returns (components, parity, component_of). Every component is a path or
    a 6-cycle -- adjacency on a unit circle means exactly 60 degrees -- so a
    2-colouring of the circle is a choice of colour pair plus one orientation
    bit per component, and `parity` fixes the reference alternation.
    """
    circle = set(graph.adj[pivot])
    seen, comps = set(), []
    for s in sorted(circle):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in graph.adj[x] & circle:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(set(comp))
    par, comp_of = {}, {}
    for ci, comp in enumerate(comps):
        s = min(comp)
        par[s], front = 0, [s]
        for v in comp:
            comp_of[v] = ci
        while front:
            x = front.pop()
            for y in graph.adj[x] & comp:
                if y not in par:
                    par[y] = 1 - par[x]
                    front.append(y)
    return comps, par, comp_of


def induced_lists(graph, pivot: int, k: int, bits: int):
    """The list each non-circle vertex gets when the circle is squeezed.

    `bits` picks an orientation per hexagon; the circle then uses colours 0
    and 1, and every other vertex keeps the colours none of its circle
    neighbours took -- k-2 of them where it sees both, k-1 where it sees one.
    """
    comps, par, comp_of = circle_hexagons(graph, pivot)
    circle = set(graph.adj[pivot])
    fixed = {v: par[v] ^ ((bits >> comp_of[v]) & 1) for v in circle}
    out = {}
    for v in range(graph.n):
        if v == pivot or v in circle:
            continue
        used = {fixed[u] for u in graph.adj[v] & circle}
        out[v] = [c for c in range(k) if c not in used]
    return out, len(comps)


def list_colourable(graph, vertices: Sequence[int], lists: dict) -> bool:
    """Whether the induced subgraph admits a colouring respecting the lists."""
    from pysat.solvers import Solver

    vs = list(vertices)
    pos = {v: i for i, v in enumerate(vs)}
    var, n = {}, 0
    for v in vs:
        if not lists.get(v):
            return False
        for c in lists[v]:
            n += 1
            var[(v, c)] = n
    cls = [[var[(v, c)] for c in lists[v]] for v in vs]
    for i, u in enumerate(vs):
        for w in vs[i + 1:]:
            if w in graph.adj[u]:
                for c in set(lists[u]) & set(lists[w]):
                    cls.append([-var[(u, c)], -var[(w, c)]])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


def minimise_core(rel: "ColourRelations", p: int,
                  targets: Sequence[int]) -> List[int]:
    """Drop what is not needed: counterexample construction is greedy.

    `cegar_core` adds whichever vertex kills the current counterexample, so
    what it returns is a core but rarely a small one. Each removal is one
    `is_core` query, which is at most k solves, so shrinking is cheap next to
    building -- and shrinking to size matters, because the blockable sizes are
    exactly one, two and three.
    """
    cur = list(targets)
    changed = True
    while changed and len(cur) > 1:
        changed = False
        for t in list(cur):
            trial = [x for x in cur if x != t]
            if trial and is_core(rel, p, trial):
                cur, changed = trial, True
                break
    return cur


# -- and adding one point to a critical graph is provably useless ----------
#
# There is a real tension in what a core needs. Rigid graphs have few
# colourings and so small cores -- but the most rigid are the vertex-critical
# ones, and criticality kills cores outright. Loose graphs escape criticality,
# but their cores are huge: a counterexample can park the free colour
# anywhere, and hundreds of rounds of counterexample-killing close nothing on
# the unions G union f(G).
#
# The obvious resolution is to add ONE point to a k-critical graph. In
# W = G + v every old vertex u still has W - u containing G - u, which is
# (k-1)-colourable, so u stays critical; only v has W - v = G, which is
# k-chromatic. So v is the unique non-critical vertex, the unique possible
# pivot, and G keeps every bit of its rigidity.
#
# It does not work, and the reason is exact.
#
# THEOREM. Let G be k-vertex-critical and W = G + v. Then every core of v
# contains all of V(G) \ N(v) -- so the minimal core has size n - deg(v),
# which is the whole graph but for the pivot's neighbours.
#
# Proof. Fix any u in V(G) \ N(v). G - u is (k-1)-colourable, and v has at
# most deg(v) neighbours, all in G - u, so W - u is (k-1)-colourable too;
# colour it with 1..k-1 and give u the colour k. Now u is the ONLY vertex
# carrying k, and none of v's neighbours does, so recolouring v to k is
# proper. In that colouring c(v) = k = c(u) and nothing else has it, so a
# target set avoiding u fails to contain c(v) and is not a core. Hence every
# core contains u. []
#
# Verified on the Moser spindle, which is 4-critical: adding any of the exact
# points one away from two of its vertices gives a pivot of degree 2 whose
# minimal core is 5 -- exactly the five vertices outside its neighbourhood,
# every one of them forced in.
#
# So the tension does not resolve that way. A usable core needs a host that is
# k-chromatic, has the pivot non-critical, AND has no vertex that can be left
# alone in a colour -- which means no k-critical subgraph avoiding the pivot
# can be (k-1)-coloured alongside it. That is a third condition, and nothing
# in this package satisfies all three.

def forced_into_every_core(rel: "ColourRelations", p: int,
                           candidates: Optional[Sequence[int]] = None
                           ) -> List[int]:
    """Vertices no core of p can omit.

    u is forced in exactly when dropping it breaks the core condition on
    everything else -- which happens whenever some colouring leaves u alone in
    a colour that p may take.
    """
    others = [v for v in range(rel.graph.n) if v != p]
    pool = others if candidates is None else list(candidates)
    return [u for u in pool
            if not is_core(rel, p, [x for x in others if x != u])]


# -- what the three conditions really ask for: unique colourability --------
#
# The three conditions a usable core needs -- k-chromatic, pivot non-critical,
# no vertex leavable alone in a colour -- have a classical name between them.
#
# A core of size r asks that c(p) lie in c(T) always, with r as small as
# possible; and what makes that hard is that colour CLASSES move. T has to
# meet every class that could be free at p, and in a loose graph the classes
# are large and mobile, so T is large. Pin the classes and T shrinks to one
# representative each:
#
#   in a UNIQUELY k-colourable graph -- one whose k-colouring is unique up to
#   permuting colours -- a pivot of pressure q has a core of exactly k - q,
#   namely one vertex from each class its neighbourhood misses.
#
# At k = 5 with the pressure 2 that every graph here measures, that is a core
# of exactly THREE, which is exactly the size the blocking bounds allow. So
# the remaining object has a name:
#
#   A UNIQUELY 5-COLOURABLE UNIT-DISTANCE GRAPH, with a pivot whose three
#   free-class representatives sit at radii and angles matching one of the
#   27344 blocking patterns, gives chi(R^2) >= 6.
#
# The mechanism is visible one level down. The triangular lattice IS uniquely
# 3-colourable -- its colouring is the Eisenstein residue modulo (1 - omega)
# -- and measuring a pivot there gives pressure 2 and a core of size ONE, at
# squared distance 3. That is the classical rhombus forcing, recovered as a
# special case: sqrt(3) forces two points to agree at three colours precisely
# because the lattice's classes cannot move.
#
# It also closes the loop with criticality. A graph is uniquely k-colourable
# exactly when the forced-same relation has k classes, and a k-vertex-critical
# graph has no forced-same pair at all -- so a k-critical graph is never
# uniquely k-colourable beyond the trivial case. The two theorems are the two
# ends of the same axis, and the object wanted sits at the far end from
# de Grey's G.

UNIQUE_COLOURING_TARGET = {
    "statement": ("a uniquely 5-colourable unit-distance graph gives a core "
                  "of size 5 - pressure at every pivot, which is 3 wherever "
                  "pressure is 2 -- exactly the blockable size"),
    "witness_one_level_down": ("the triangular lattice is uniquely "
                               "3-colourable; a pivot there has pressure 2 "
                               "and a core of ONE, at d^2 = 3, which is the "
                               "classical rhombus forcing"),
    "why_criticality_is_the_opposite_end": ("uniquely k-colourable means the "
                                            "forced-same relation has k "
                                            "classes; a k-critical graph has "
                                            "no forced-same pair at all"),
}


# -- why the difficulty is exactly where it is ----------------------------
#
# Two facts already proved above combine into a ladder that says where this
# whole method works and where it stops, with no search involved.
#
# The unit circle around any point is bipartite -- adjacency there means
# exactly 60 degrees and the 60-degree orbit is a 6-cycle -- so a
# neighbourhood containing an edge has pressure at least 2, for free, and
# pressure above 2 is entirely the ambient graph's doing. Call 2 the FREE
# pressure. The pressure theorem then makes the smallest possible core
#
#     r >= k - pressure(p)  =  k - 2   when nothing beyond the circle helps,
#
# and blocking, by the counting and misalignment results, reaches r <= 3 and
# no further:
#
#   k = 3: free core 1. A core of one is a forced pair, blocked by the
#          classical spindle. Rigidity is FREE here, since 2 = k-1 means the
#          circle alone determines the centre's colour -- which is exactly why
#          the triangular lattice is uniquely 3-colourable and why sqrt(3)
#          forces agreement.
#   k = 4: free core 2, blocked by counting. de Grey's construction.
#   k = 5: free core 3, blocked ONLY by misalignment, and only with
#          reflections. The free pressure exactly saturates the blocking
#          bound -- there is no slack anywhere.
#   k = 6: free core 4, and nothing blocks four. The method is not hard here,
#          it is impossible.
#
# So the reason five is the frontier is not that nobody has searched hard
# enough. It is the last value of k at which the free pressure of a unit
# circle and the largest blockable core still meet, and they meet exactly.

FREE_PRESSURE = 2
LADDER = {3: "core 1, classical spindle; rigidity free since 2 = k-1",
          4: "core 2, blocked by counting; de Grey",
          5: "core 3, blocked only by misalignment with reflections -- exact",
          6: "core 4, nothing blocks four: the method is impossible"}


# -- how far a graph is from unique colourability, in one number ----------
#
# THEOREM. A uniquely k-colourable graph has pressure exactly k-1 at EVERY
# vertex.
#
# Proof. Pressure is at most k-1 always, since v's own colour never appears in
# N(v). If some v had pressure at most k-2, a minimising colouring would leave
# it two free colours; switching v between them moves v to a different class
# and so gives a genuinely different partition, not a permutation of the same
# one. []
#
# That turns "how close is this graph to uniquely k-colourable" into a
# measurement, vertex by vertex, with a gradient rather than a verdict. And it
# says exactly how far the problem is from its target:
#
#   k = 3: pressure 2 is the FREE pressure of a unit circle, so unique
#          3-colourability costs nothing extra -- the triangular lattice.
#   k = 4: pressure 3 is one above free. de Grey's Sa reaches it at some
#          vertices, which is where its forcing comes from.
#   k = 5: pressure 4 is two above free, and nothing measured here reaches
#          even 3 -- every graph tried comes back at a flat 2.
#
# The gap is 2 against 4, stated in the same units as everything else.

def unique_colouring_defect(rel: "ColourRelations") -> dict:
    """Vertices at pressure k-1, and the shortfall at the rest.

    Zero defect everywhere is necessary for unique k-colourability, never
    sufficient; but the count is monotone in added points, so it is something
    a search can climb.
    """
    hist = collections.Counter()
    for v in range(rel.graph.n):
        if rel.graph.adj[v]:
            hist[rel.k - 1 - pressure(rel, v)] += 1
    total = sum(hist.values())
    return {
        "defect_histogram": dict(sorted(hist.items())),
        "at_k_minus_one": hist.get(0, 0),
        "vertices": total,
        "uniquely_colourable_possible": hist.get(0, 0) == total,
    }


# -- the invariant everything reduces to ----------------------------------
#
# Strip the pivot out of the core condition and one graph invariant is left:
#
#     rho(W, k) = the least size of a set that uses all k colours in EVERY
#                 k-colouring -- a rainbow-forcing set.
#
# Every bound above is a statement about it.
#
#   In a k-vertex-critical graph rho = n. For any u, colour W - u with k-1 and
#   give u the kth; a set omitting u then misses that colour. So every vertex
#   is needed, which IS the inertness of de Grey's G.
#   In a uniquely k-colourable graph rho = k, one representative per class.
#   Measured: 3 on a triangular-lattice patch at three colours, and 6 on de
#   Grey's Sa at four.
#
# And it converts into cores directly, with the pivot doing the work.
#
# THEOREM. If S is a rainbow-forcing set and p is any point, then S \ N(p) is
# a core of p.
#
# Proof. N(p) union (S \ N(p)) contains S, so it uses every colour in every
# colouring, which is the core condition. []
#
# So the core's size is |S| - |S intersect N(p)| and the requirement is just
# that p be adjacent to at least |S| - 3 elements of S:
#
#   |S| = 4: one element -- trivially placed.
#   |S| = 5: two -- the two circle intersections of any pair less than two
#            apart, which always exist.
#   |S| = 6: three, which needs them concyclic at radius EXACTLY one, and
#            then p is their circumcentre.
#
# One trap, met on the first attempt and worth stating: p must not itself lie
# in S. A minimal forcing set loses the property when any element is dropped,
# so a circumcentre that happens to be a forcing vertex removes its own target
# and the core evaporates -- measured, on Sa, as a perfectly good circumcentre
# of degree 30 whose core was not one.

def forcing_set(rel: "ColourRelations", order: Optional[Sequence[int]] = None,
                limit: int = 400) -> Tuple[List[int], bool]:
    """A set using all k colours in every colouring, built by counterexamples.

    Returns (S, closed). `closed` is False when the limit was reached, in
    which case S is a lower bound witness rather than a forcing set.
    """
    pool = list(range(rel.graph.n) if order is None else order)
    S: List[int] = []
    for _ in range(limit):
        rel.calls += 1
        if not rel._s.solve(assumptions=[-rel._x(v, 0) for v in S]):
            return S, True
        m = set(rel._s.get_model())
        cand = [v for v in pool if v not in S and rel._x(v, 0) in m]
        if not cand:
            return S, False
        S.append(cand[0])
    return S, False


def shrink_forcing_set(rel: "ColourRelations",
                       S: Sequence[int]) -> List[int]:
    cur, changed = list(S), True
    while changed and len(cur) > 1:
        changed = False
        for t in list(cur):
            trial = [x for x in cur if x != t]
            if trial and min_colours_on(rel, trial) >= rel.k:
                cur, changed = trial, True
                break
    return cur


def core_from_forcing_set(rel: "ColourRelations", p: int,
                          S: Sequence[int]) -> List[int]:
    """S minus N(p), which is a core of p whenever S is forcing and p not in S."""
    nb = rel.graph.adj[p]
    return [v for v in S if v != p and v not in nb]


# -- why stacking copies can never work, as a theorem ---------------------
#
# S is rainbow-forcing exactly when it meets every colour class of every
# colouring -- equivalently, every independent set I with W - I still
# (k-1)-colourable. That reading kills the whole family of constructions this
# package kept returning to.
#
# THEOREM. Let W split as A, B and a shared part, with W - a - b
# (k-1)-colourable for every non-adjacent a in A, b in B. Then every
# rainbow-forcing set contains all of A or all of B, so rho >= min(|A|, |B|).
#
# Proof. For such a pair, colour W - a - b with k-1 and give a and b the kth
# colour; they are non-adjacent, so it is proper, and {a, b} is then a colour
# class. A forcing set must therefore contain a or b. Ranging over all pairs,
# it is a vertex cover of the complete bipartite graph between A and B, and
# the only vertex covers of that are A and B. []
#
# The hypothesis is exactly what a union of two k-critical graphs supplies:
# deleting one vertex from each copy leaves both (k-1)-colourable. Verified on
# two Moser spindles glued at two vertices -- all 22 non-adjacent cross pairs
# are killable, and the forcing set comes back as all five of one side plus
# both shared vertices, rho = 7 against the bound's 5.
#
# WITHDRAWN AS APPLIED TO G. The hypothesis asks that W - a - b be
# (k-1)-colourable for every non-adjacent cross pair, which for copies of de
# Grey's G means G - a must be 4-colourable -- exactly the vertex-criticality
# that is retracted above, since G - 1420 is still 5-chromatic. So the bound
# rho >= 1357 for G union (G + t) does not follow, and the conclusion that
# "every union in this package was dead before it was built" is withdrawn with
# it. The theorem stands for genuinely critical A and B; G is not one.
#
# Monotonicity points the other way in any case: rho(W) <= rho(G) whenever G
# sits inside W, so unions can only lower it.

def cross_pair_bound(graph, k: int, A: Sequence[int], B: Sequence[int],
                     sample: int = 8) -> dict:
    """Check the hypothesis on a sample and report the bound it gives."""
    import random

    from .coloring import is_k_colorable
    from .graph import build_graph

    rng = random.Random(0)
    pairs, killable = 0, 0
    for a in rng.sample(list(A), min(sample, len(A))):
        for b in rng.sample(list(B), min(sample, len(B))):
            if b in graph.adj[a]:
                continue
            pairs += 1
            sub = build_graph([q for m, q in enumerate(graph.vertices)
                               if m != a and m != b])
            if is_k_colorable(sub, k - 1)[0]:
                killable += 1
    return {
        "pairs_tested": pairs,
        "killable": killable,
        "hypothesis_holds_on_sample": pairs > 0 and killable == pairs,
        "bound": min(len(A), len(B)),
    }


# -- when rho is n, and therefore when it can be small --------------------
#
# Two disjoint reasons make rho equal n, and between them they account for
# every graph in this package at five colours.
#
#   k > chi(W). Every vertex is then removable: chi(W - u) <= chi(W) <= k-1,
#     so u can be left alone in the kth colour and every forcing set needs it.
#     The three-hexagon gadget is 4-chromatic, so at k = 5 it measures
#     rho > 123 with all twelve sampled vertices critical -- useless, for a
#     reason that has nothing to do with its geometry.
#   W k-vertex-critical. Same conclusion by the same colouring. de Grey's G.
#
# So rho can only be small where chi(W) = k EXACTLY and W is not
# vertex-critical. The three-hexagon gadget has both at four colours -- rho 9,
# zero of twelve sampled vertices critical. At five colours the reading has
# changed: G is NOT critical after all, so it meets both conditions too, and
# the unions no longer fall to the cross-pair theorem.
#
# Adding copies makes it worse, not better. With three pairwise-overlapping
# copies a forcing set must hit every cross TRIPLE, so it is a vertex cover of
# a complete tripartite 3-uniform hypergraph and has to swallow two whole
# parts rather than one. The family degrades with every copy added.
#
# What is left is the same object from the other side. Few realisable colour
# classes is what makes rho small, and a graph with exactly one 5-colouring up
# to permutation has exactly five of them, giving rho = 5. The target is a
# uniquely 5-colourable unit-distance graph, which by the pressure theorem
# needs pressure 4 -- every neighbourhood using four colours -- at every
# single vertex.

RHO_IS_N_WHEN = ("k > chi(W), since every vertex is then removable; or W is "
                 "k-vertex-critical, by the same colouring")


# -- the confined set is indexed by a cut, and that bounds it -------------
#
# Squeeze the pivot's circle into colours {0,1} and an auxiliary seeing both
# gets the list {2,3,4}, size 3, the SAME list for every such auxiliary, while
# one seeing a single colour gets a list of size 4. (Two circles meet in at
# most two points, so no auxiliary ever sees three: the list sizes are exactly
# 3 and 4.) The confined set is where the squeeze is resisted, in every
# orientation, since a 2-colouring the design does not control is an escape.
#
# STATED CAREFULLY, because an earlier version of this note overreached. The
# confined set failing to be list-colourable is the LOCAL obstruction -- the
# one the k = 4 construction actually runs on -- not an equivalent of pressure
# 3. A 3-colourable confined set does not by itself give a colouring of the
# whole graph, and pressure 3 could in principle come from deeper structure
# instead. What follows explains the measurements and says what the local
# mechanism can and cannot do; it is not a proof that the pressure must be 2.
#
#   WHAT IT IS NOT. "The auxiliary graph is 4-chromatic" is not the condition.
#   Measured on de Grey's Sa at its best pivot: chi(A) = 4 already, degeneracy
#   4, and the pressure at five colours is still 2. The lists are what differ.
#
# Writing bits[a] for the parity of hexagon a, the point at position i in it
# takes colour (i + bits[a]) mod 2, so the auxiliary u_i + v_j from hexagons a
# and b is confined exactly when i + j + bits[a] + bits[b] is odd. Two
# consequences, both measured rather than assumed.
#
#   SAME HEXAGON: the bits cancel, so those auxiliaries are confined in EVERY
#   orientation. They are the points at sqrt 3 from the pivot, and their graph
#   has maximum degree ONE -- a matching, 9 disjoint edges on 18 points. They
#   can contribute 2 to a chromatic number and never more.
#
#   ACROSS HEXAGONS: only eps_ab = bits[a] XOR bits[b] matters, and eps is a
#   CUT of K_t, so eps_ab + eps_bc + eps_ac = 0 for every triple. The 2^t
#   orientations give only 2^(t-1) distinct confined sets, and no design can
#   choose the pair-parities independently: of the 8 patterns at t = 3 exactly
#   4 arise, and every missing one breaks the triangle identity.
#
# Measured maxima of chi(confined set), over all orientations:
#
#     2 hexagons   48 auxiliaries    chi in [1, 1]
#     3 hexagons  126 auxiliaries    chi in [2, 3]
#     de Grey Sa, 5 hexagons          chi in [2, 3]   (4 orientations at 2)
#
# Four is what is needed, in every orientation, and nothing here reaches it.

CONFINED_SET_CONDITION = (
    "the local mechanism behind pressure 3 is the confined set's list "
    "problem, not the auxiliary graph's chromatic number -- chi(A) is already "
    "4 on Sa and buys nothing, because only the auxiliaries seeing both "
    "colours share a list; whether pressure 3 could come from deeper "
    "structure instead is open here"
)

CONFINED_SET_IS_A_CUT = (
    "confinement of u_i + v_j depends on i + j + bits[a] + bits[b], so "
    "same-hexagon auxiliaries are always confined and form a matching, while "
    "cross-hexagon ones see only the cut eps = bits[a] XOR bits[b]: 2^(t-1) "
    "distinct confined sets, not 2^C(t,2)"
)


# -- why k = 4 works and k = 5 does not, in one number -------------------
#
# The confined set is 2-DEGENERATE. Measured over de Grey's own angle family,
# hexagons at k*theta/2, taking the worst orientation at each size:
#
#     hexagons   auxiliaries   confined   max degree   degeneracy
#         3          126          72          4            2
#         4          240         132          4            2
#         5          390         210          4            2
#         6          576         306          4            2
#         7          798         420          4            2
#
# The confined set grows nearly sixfold and the degeneracy does not move. And
# degeneracy is exactly the statistic the list-colouring argument turns on:
#
#     d-degenerate  =>  (d+1)-CHOOSABLE   (greedy, peeling low-degree first)
#
# so a 2-degenerate confined set is 3-choosable. That settles both colour
# counts at once, in opposite directions:
#
#   k = 4. The lists have size 2. A 2-degenerate graph need not be
#   2-choosable -- an odd cycle is 2-degenerate and has list chromatic number
#   3 -- so the lists can fail, and de Grey's Sa reaches pressure 3.
#
#   k = 5. The lists have size 3. Every 2-degenerate graph IS 3-choosable, so
#   the lists always complete, and the pressure is 2. No number of hexagons
#   changes this, because none of them changes the degeneracy.
#
# So the frontier at five colours is not a search that has not yet succeeded.
# The mechanism the method runs on lives in the gap between "2-degenerate is
# not always 2-choosable" and "2-degenerate is always 3-choosable", and at
# five colours it falls on the wrong side of it by exactly one. Reaching
# pressure 3 through this mechanism needs a confined set of degeneracy at
# least 3 -- a statement about the geometry of Minkowski sums of hexagons, not
# about how hard the search is run. Reaching it some OTHER way is not excluded
# by any of this; nothing measured here does.

CONFINED_SET_IS_TWO_DEGENERATE = {
    "measured_sizes": [72, 132, 210, 306, 420],
    "degeneracy": 2,
    "max_degree": 4,
    "consequence_k4": "lists of size 2; 2-degenerate is not 2-choosable "
                      "(an odd cycle is not), so pressure 3 is reachable",
    "consequence_k5": "lists of size 3; 2-degenerate IS 3-choosable by "
                      "Erdos-Rubin-Taylor, so the lists always complete and "
                      "the pressure is 2, at every size",
    "what_would_break_it": "a confined set of degeneracy >= 3",
}


# -- how much of the design space that measurement actually covers --------
#
# Two auxiliaries q = u + v and q' = u' + v' from the SAME pair of hexagons a
# and b differ by D1 + D2 with D1 = r_a(w^i - w^i'), D2 = r_b(w^j - w^j'), and
#
#     |q - q'|^2 = |D1|^2 + |D2|^2 + 2 Re(t z) = 1,   t = r_a conj(r_b),
#
# a LINE in t meeting the unit circle twice: for two hexagons the good angles
# are a finite set, and enumerating all 31 x 31 pairs of hexagon differences
# gives 48 solutions, EIGHT distinct modulo 60 degrees. Measured at each:
#
#     degeneracy 0: 1 angle      degeneracy 2: 3 angles
#     degeneracy 1: 2 angles     degeneracy 3 or more: NONE
#
# One case is free and useless: D2 = 0, same v and adjacent u, holds at every
# angle -- but confinement of u_i + v_j turns on i + j, adjacent i differ by
# one, so those two auxiliaries never share a parity. Every always-available
# edge joins a confined auxiliary to a non-confined one, which is why
# Pythagorean rotations leave the confined set with no edges at all.
#
# WHAT THIS DOES NOT SHOW, stated because the obvious reading is wrong. The
# enumeration is complete for TWO hexagons only. With three, an edge can join
# auxiliaries from DIFFERENT pairs -- q from (a,b) and q' from (a,c) -- which
# is two free angles against one equation, a CURVE rather than a discrete set.
# de Grey's theta/2 is exactly such a case: it is absent from the eight, and
# chi(confined) is 1 at two hexagons and jumps to 3 at three. So the flat
# degeneracy measured as hexagons are added describes ONE CURVE through a
# continuous space, not the space. It explains why de Grey's family stops; it
# does not close the route.

TWO_HEXAGON_ANGLES_ARE_FINITE = {
    "raw_solutions": 48,
    "distinct_mod_60": 8,
    "best_degeneracy": 2,
    "free_case": "D2 = 0 holds at every angle but always joins a confined "
                 "auxiliary to a non-confined one, so it never helps",
    "caveat": "complete for two hexagons only; with three the edges can cross "
              "hexagon pairs, giving two angles against one equation -- a "
              "curve, and de Grey's theta/2 lives on it",
}


# -- scanning the two-parameter space de Grey's family is a curve in ------
#
# With three hexagons the confined set's edges come from equations relating
# TWO free angles, so the design space is a surface rather than a list. It can
# still be scanned, because a configuration is rich exactly where many
# equations hold at once. Fixing alpha1, an edge between q = u_i + v_j from
# the pair (0,1) and q' = u_i' + w_l from (0,2) needs
#
#     |(w^i - w^i') + e^(i alpha1) w^j - e^(i alpha2) w^l| = 1,
#
# which for A = the first two terms and C = A conj(w^l) reads
# Re(C e^(-i alpha2)) = |A|^2 / 2, so alpha2 = arg C -+ arccos(|A|/2) whenever
# |A| <= 2. Each (i, i', j, l) gives up to two values, 2160 in all, and the
# rich configurations are the HISTOGRAM PEAKS of those values.
#
# The method reproduces de Grey's own configuration as a check: at
# alpha1 = theta/2 the value alpha2 = theta carries multiplicity 144, and the
# configuration measures 126 auxiliaries, 306 edges, confined degeneracy 2 --
# agreeing exactly with the field arithmetic.
#
# Scanned over 900 values of alpha1 and the six richest alpha2 at each, 5400
# configurations: nothing exceeds degeneracy 2.
#
# TWO MEASUREMENT BUGS worth recording, because both made the scan blind
# rather than wrong-looking. Rounding auxiliary coordinates to 1e-6 while
# testing distances against 1e-7 found 50 of de Grey's 306 edges; the keys are
# now rounded to 1e-9, where the arithmetic really carries 1e-15. And the
# largest histogram peak is always alpha2 = alpha1, which is the same six
# points, so coincident offsets are dropped before ranking rather than after.

THREE_HEXAGON_SCAN = {
    "free_parameters": 2,
    "candidate_equations": 2160,
    "validated_on": "de Grey's alpha1 = theta/2, alpha2 = theta: multiplicity "
                    "144, 126 auxiliaries, 306 edges, degeneracy 2",
    "configurations_evaluated": 5400,
    "best_degeneracy": 2,
    "needed": 3,
}


# -- what the confined set looks like, and what it never does -------------
#
# On de Grey's own three-hexagon configuration every confined auxiliary has
# degree 0, 2 or 4 -- never odd -- with 42 of 72 isolated, 36 edges, and
# components of cycle rank 1 or 4. That looked like a law, and it is not: over
# 1200 configurations sampled across the three-hexagon surface, 1276
# orientations carry an odd confined degree and degree 5 occurs. The even
# degrees come from the reflective symmetry of that particular offset family,
# not from the construction.
#
# What does hold, across 1200 configurations and 9600 orientations:
#
#     largest confined degree seen:  5
#     largest DEGENERACY seen:       2
#
# Never 3, which is what lists of size 3 would need to fail. Stated as a
# measurement, since nothing here proves it: the surface is continuous and
# only its histogram peaks were evaluated.

# Extending the peak method one hexagon further: given the offsets already
# placed, an edge joining an auxiliary that uses the NEW hexagon to one that
# does not reads |r_new w^m - X| = 1 with X the old auxiliary minus the old
# circle point it pairs with -- the same line-meets-circle problem. Scanned
# over three free angles, 360 four-hexagon configurations: still 2. So the
# cap holds over a complete enumeration at two hexagons, 9600 orientations at
# three, 360 configurations at four, and de Grey's own family out to seven.

CONFINED_DEGENERACY_NEVER_THREE = {
    "configurations": 1200,
    "four_hexagon_configurations": 360,
    "orientations": 9600,
    "max_degree_seen": 5,
    "max_degeneracy_seen": 2,
    "needed": 3,
    "status": "measurement, not proof -- peaks of a continuous surface",
    "de_grey_is_special": "its confined degrees are all even, a symmetry of "
                          "that offset family; odd degrees occur elsewhere",
}


# -- and the obvious generalisation is vacuous ----------------------------
#
# The confined set is one subgraph; the honest object is everything outside
# the circle, where squeezing into two colours gives
#
#     L(v) = k - |colours v sees on the circle|,  so 3, 4 or 5 at k = 5,
#
# and greedy completes the list colouring unless some subgraph S has
# deg_S(v) >= L(v) for every v in S. That maximal surviving subgraph, the
# L-CORE, would be the exact greedy obstruction. Measured, over every
# orientation, taking the SMALLEST core each graph achieves:
#
#     de Grey Sa   k=4: 366    k=5: 360     (397 vertices)
#     de Grey Sb   k=4: 366    k=5: 360
#     de Grey Y    k=4: 730    k=5: 718     (791 vertices)
#     gadget       k=4: 108    k=5: 108     (127 vertices)
#
# Non-empty everywhere and at both colour counts, because these graphs have
# average degree around ten and almost nothing peels. So the greedy criterion
# separates nothing: it does not explain why the pressure is 3 at four colours
# and 2 at five, and it yields no new sufficient condition.
#
# That is why the degeneracy result above is stated about the confined set's
# MECHANISM and not about the pressure. The mechanism is genuinely sufficient
# when it fails -- an odd cycle on one shared 2-list is infeasible outright,
# which is what k = 4 runs on -- and at five colours it never fails. Whether
# pressure 3 could come from somewhere else stays open, and nothing measured
# here decides it in either direction.

L_CORE_IS_NOT_THE_OBSTRUCTION = {
    "smallest_l_core": {"Sa": {4: 366, 5: 360}, "Sb": {4: 366, 5: 360},
                        "Y": {4: 730, 5: 718}, "gadget": {4: 108, 5: 108}},
    "verdict": "non-empty at both colour counts, so greedy peeling separates "
               "nothing and gives no sufficient condition; the confined-set "
               "mechanism remains sufficient-when-it-fails, and only that",
}


# -- the missing object, stated precisely ---------------------------------
#
# Reading the theorems together narrows the gap more than any one of them
# does, and it is worth writing down what is actually left.
#
# The pressure theorem says a core of size r at p needs pressure(p) >= k - r.
# At five colours a core of THREE therefore needs pressure >= 2 -- and the
# measured pressure is exactly 2, everywhere, on every graph here. So pressure
# does not exclude a core of three. It never did; the flat 2 was read as a
# wall when it is only a floor.
#
# Blocking a core of three is solved: 27344 misalignment configurations do it,
# and the counting bound (r <= 2) together with the pressure floor (r >= 3)
# pins the size at exactly three from both sides.
#
# So neither pressure nor blocking is what is missing. What is missing is:
#
#   A 5-chromatic graph W with a pivot p that is NOT essential -- chi(W - p)
#   still 5, so the criticality corollary cannot give p a colour of its own --
#   carrying an actual core of three: min |c(N(p) u T)| = 5 for some |T| = 3.
#
# The first half is easy and known: W = G u (G + t) works, because removing
# any vertex leaves a whole 5-chromatic copy. The second half has never been
# found on anything tested. That single measurement is the gap.
#
# Note what it is NOT. It is not a bigger search over de Grey's G, which the
# criticality corollary settles in one line. It is not more pressure, which is
# already sufficient at this size. It is not a better blocking configuration,
# of which there are 27344. It is one core, of three, at five colours.

WHAT_IS_MISSING = {
    "not_pressure": "a core of three needs pressure >= k - r = 2, and the "
                    "measured pressure is exactly 2 everywhere -- a floor "
                    "that is met, not a wall",
    "not_blocking": "27344 misalignment configurations block three targets; "
                    "counting caps r at 2 and pressure floors it at 3, so the "
                    "size is pinned at exactly three",
    "missing": "a 5-chromatic graph with a NON-essential pivot carrying an "
               "actual core of three at five colours",
    "where_to_look": "unions G u (G + t): removing any vertex leaves a whole "
                     "5-chromatic copy, so no pivot gets a colour of its own",
}


# -- and all of it is one number after all --------------------------------
#
# A core of size r at p means min |c(N(p) u T)| = k: every colouring uses all
# k colours on that set. That is exactly what rho, the rainbow-forcing number,
# measures -- the least size of a set using all k colours in every colouring.
# So a core of size r at p forces
#
#     rho <= deg(p) + r,
#
# and on de Grey's G, whose maximum degree is 60, a core of three needs
# rho <= 63. Measured there: rho = n = 1581, because a k-vertex-critical graph
# has rho = n. And the cross-pair theorem says a union of critical graphs has
# rho >= min(|A|,|B|), so G u f(G) has rho >= 1581 too, however heavily it
# folds. The folded union was therefore ruled out before it was built:
# 2373 vertices, 789 of them folded, pressure 2 at every hub, and no core in
# 60 steps at any of four pivots. Confirmation, not discovery.
#
# So the ladder, the pressure, the cores, the criticality corollary and the
# cross-pair bound are all the same statement about rho, and the target is:
#
#     A 5-chromatic unit-distance graph with rho <= 63.
#
# For comparison rho = k exactly when the graph is uniquely k-colourable, so
# what is wanted is a unit-distance graph that is NEARLY uniquely 5-colourable
# -- and every graph measured in this package has rho = n or rho >= n/2.
# That is the whole distance to the goal, in one number.

RHO_IS_THE_WHOLE_GAP = {
    "identity": "a core of size r at p means N(p) u T uses all k colours in "
                "every colouring, so rho <= deg(p) + r",
    "degrey_needs": 63,
    "degrey_has": "at most 1580, since G - 1420 is a 5-chromatic subgraph; "
                  "the earlier 1581 rested on a criticality claim now "
                  "retracted, and 1201 vertices are measured not forcing",
    "why": "vertex-critical graphs have rho = n -- but G is NOT one, so that "
           "route to rho = n is gone; the cross-pair theorem still gives "
           "unions of critical graphs rho >= min(|A|,|B|)",
    "folded_union_checked": {"vertices": 2373, "folded": 789, "pressure": 2,
                             "core_in_60_steps": False},
    "target": "a 5-chromatic unit-distance graph with rho <= 63, which is to "
              "say one that is nearly uniquely 5-colourable",
    "rho_of_uniquely_colourable": 5,
}


# -- the frontier in one measurement --------------------------------------
#
# The rho reading is only worth anything if it comes out right where the
# construction actually succeeds. At four colours de Grey's Sa reaches
# pressure 3 with cores down to 7 at a pivot of degree 30, which predicts
# rho(Sa, 4) <= 37. Measured on the same 397 points, one colour apart:
#
#     Moser spindle, k = 4:   rho = 7 = n     (4-vertex-critical, so rho = n)
#     de Grey Sa,    k = 4:   rho = 7         on 397 vertices
#     de Grey Sa,    k = 5:   the same construction runs past 358 without
#                             closing
#
# Seven vertices out of 397 use all four colours in every colouring. That is
# the rigidity the whole spindle method runs on, and the prediction is met
# with room to spare. Add ONE colour to the same graph and the construction
# that closed at seven does not close at 358.
#
# Stated carefully: the 358 bounds the greedy path, not rho itself, since a
# different set of that size might still be forcing. What it does show is that
# the counterexample construction -- the same one, on the same points --
# behaves in completely different regimes at four and five colours.
#
# That is the frontier, in one measurement, on one graph.

RHO_JUMPS_AT_FIVE = {
    "graph": "de Grey Sa, 397 vertices, maximum degree 30",
    "rho_at_4": 7,
    "greedy_at_5_did_not_close_by": 358,
    "predicted_bound_at_4": 37,
    "moser_spindle_at_4": 7,
    "caveat": "358 bounds the greedy path, not rho: a different set of that "
              "size might still be forcing",
}


# -- small rho does NOT need a small critical subgraph --------------------
#
# One direction is immediate: a k-chromatic subgraph uses all k colours in
# every k-colouring of the whole graph, so
#
#     rho(G, k) <= |H|  for every k-chromatic subgraph H of G.
#
# At four colours that would explain everything -- Sa contains Moser spindles,
# seven vertices, 4-critical -- and it is the obvious reading of rho(Sa,4) = 7.
#
# It is the WRONG reading. The seven vertices the construction returns are
#
#     S = [0, 3, 4, 6, 7, 8, 9],
#
# which induce FIVE edges, have chromatic number 3, and include three isolated
# vertices. Not a spindle, not critical, not 4-chromatic. Re-derived from a
# CNF built from scratch rather than through ColourRelations: for each of the
# four colours, no proper 4-colouring of Sa leaves that colour off S. And
# minimal -- all seven single deletions break it.
#
# So the forcing is AMBIENT. It is carried by the rest of the 397 vertices,
# not by anything inside the set, and the converse of the theorem is false.
#
# That matters for what is left. A core of three at a degree-60 pivot needs
# rho <= 63. Had small rho required a small k-chromatic subgraph, this would
# be asking for a 5-chromatic unit-distance graph on 63 vertices, against a
# published record of around five hundred -- hopeless. It does not. It is
# asking for ambient rigidity of exactly the kind Sa already exhibits at four
# colours, with a set that is nearly edgeless.

SMALL_RHO_IS_AMBIENT = {
    "theorem": "rho(G,k) <= |H| for every k-chromatic subgraph H",
    "converse": False,
    "witness": {"graph": "de Grey Sa", "k": 4, "set": [0, 3, 4, 6, 7, 8, 9],
                "induced_edges": 5, "induced_chi": 3, "isolated": 3,
                "minimal": True,
                "verified": "independent CNF, all four colours checked"},
    "consequence": "rho <= 63 at five colours does not ask for a 63-vertex "
                   "5-chromatic unit-distance graph; it asks for ambient "
                   "rigidity, which Sa already has at four",
}


# -- and being at the boundary is the hardest case, not a free one --------
#
# "Pressure does not exclude a core of three" is true and was worth saying,
# but it is easy to read as more room than it gives. The inequality
# pressure >= k - r is necessary; sitting exactly ON it is the worst place to
# be.
#
# A core of size r at p means every colouring uses all k colours on
# N(p) u T. The circle supplies at least pressure(p) of them. So in the
# colourings that squeeze the circle down to its minimum, T has to supply
#
#     k - pressure(p)   colours, all of them new,
#
# which at k = 5 with pressure 2 and r = 3 means THREE targets carrying three
# colours, all different, none of them the circle's two -- in every such
# colouring. Three vertices forced to be rainbow and colour-disjoint from the
# circle is close to the original problem restated.
#
# Room appears only where pressure EXCEEDS k - r. That is exactly what happens
# at four colours: de Grey's Sa has pressure 3 against k - r = 4 - 7 < 0, so
# the targets have slack everywhere, and the measured cores go down to 7 while
# rho comes out at 7 on 397 vertices. The spindle method has always run on
# that slack.
#
# So the honest form is: pressure 2 leaves a core of three possible and gives
# it no margin at all. The place to look is a pivot with pressure 3 at five
# colours, which is the thing the degeneracy measurements say this angle
# family cannot produce.

BOUNDARY_PRESSURE_HAS_NO_MARGIN = (
    "pressure >= k - r is necessary, and sitting exactly on it means the r "
    "targets must carry all k - pressure remaining colours, rainbow and "
    "disjoint from the circle's, in every squeezing colouring; margin needs "
    "pressure strictly above k - r, which is what four colours has and five "
    "does not"
)


# -- and how far that is from happening on G ------------------------------
#
# rho(Sa,4) = 7 works because the pressure at the pivot is 3: N(p) always
# carries three colours, so p together with a small piece of its circle
# already forces all four. ONE pivot suffices at four colours.
#
# At five it cannot. Pressure is 2 everywhere, so {p} u N(p) forces three
# colours at most -- the pivot's own plus the circle's two -- and three is not
# five. Measured on de Grey's G, where every test is a single SAT call and
# comes back in under two seconds:
#
#     one hub's closed neighbourhood      61 vertices   not forcing
#     adjacent hub pairs                  38 vertices   not forcing
#     the 133 highest-degree closed neighbourhoods, united:
#                                       1201 vertices   NOT FORCING
#
# Seventy-six per cent of the graph, and a 5-colouring still exists that
# leaves one colour off all of it. Against seven vertices out of 397 at four
# colours.
#
# What that does and does not prove. A superset of a forcing set is forcing,
# so a non-forcing set contains NO 5-chromatic subgraph: any 5-chromatic
# subgraph of G must use a vertex outside those 1201. It does not bound rho
# itself, since some other set of that size might force.
#
# It also bears on the criticality question from an unexpected side. If G had
# a small 5-chromatic subgraph it would sit in the dense part, and this says
# it does not -- which is consistent with G being vertex-critical after all,
# and with the published smaller 5-chromatic graphs being separate
# constructions rather than subgraphs of this one.

ONE_PIVOT_CANNOT_REACH_FIVE = {
    "mechanism_at_4": "pressure 3 makes {p} u N(p) force all four colours, "
                      "which is why rho(Sa,4) = 7",
    "at_5": "pressure 2 makes {p} u N(p) force three, and three is not five",
    "measured_not_forcing": {"one_hub": 61, "adjacent_pair": 38,
                             "133_neighbourhoods": 1201},
    "graph_size": 1581,
    "consequence": "a superset of a forcing set is forcing, so no 5-chromatic "
                   "subgraph of G fits inside those 1201 vertices",
}


# -- why rho is cheap to bound below and expensive to bound above ---------
#
# The two directions cost completely different amounts, and it is worth saying
# so, because the asymmetry decides which measurements this package can make.
#
# Showing a set is NOT forcing is SAT: exhibit a proper k-colouring that
# leaves some colour off it. Small sets have many such colourings and the
# solver finds one in under two seconds even at 1201 vertices.
#
# Showing a set IS forcing is UNSAT: prove that no such colouring exists. For
# the whole vertex set that is exactly "G is not 4-colourable", the instance
# kissat needed 522 seconds for. For V - {v} it is "chi(G - v) < k", which is
# the vertex-criticality question. There is no cheap direction.
#
# So peeling from the top, which would give an upper bound on rho, costs a
# hard UNSAT per vertex; growing from the bottom, which gives lower-bound
# evidence, is nearly free. This package therefore has:
#
#     1201 vertices measured NOT forcing      (cheap, certain)
#     rho(G,5) <= 1581                        (trivial, since chi(G) = 5)
#
# and closing the gap between them is the same computation as deciding
# vertex-criticality. Both are running; neither is cheap.

RHO_BOUNDS_ARE_ASYMMETRIC = (
    "not-forcing is SAT and costs seconds; forcing is UNSAT and costs the "
    "same as the 4-colourability proof, so lower-bound evidence is nearly "
    "free and any upper bound below n is as hard as vertex-criticality"
)


# -- the retraction, as a fact rather than a note -------------------------

DEGREY_G_IS_NOT_VERTEX_CRITICAL = {
    "claim_retracted": "de Grey's G is 5-vertex-critical",
    "witness": 1420,
    "witness_degree": 4,
    "result": "G - 1420 admits no proper 4-colouring (UNSAT, 1581 s, with a "
              "triangle pinned to 0,1,2 as in the main certificate)",
    "so": "a proper subgraph on 1580 vertices is already 5-chromatic",
    "why_i_believed_it": "separability was measured and read backwards -- it "
                         "is a consequence of criticality, not a proof",
    "unaffected": "the theorem that a k-vertex-critical graph has no core; it "
                  "simply does not apply to G",
    "now_unexplained": ["pressure exactly 2 at all 1581 vertices",
                        "no forced pair in 29930 queries",
                        "no forced-different non-edge in 40539 pairs",
                        "1201 vertices not rainbow-forcing"],
    "verification": "two solvers, formulas rebuilt independently; drat-trim "
                    "is not installed here so the proof is unchecked",
}


# -- rho only goes DOWN when structure is added --------------------------
#
# If G sits inside W, every proper k-colouring of W restricts to one of G, so
# a set forcing in G forces in W:
#
#     rho(W) <= rho(G)   whenever G is contained in W.
#
# Which reverses the strategy this package had been following. The way to a
# small rho is a BIGGER graph, not a smaller one.
#
# That was hidden while the cross-pair theorem was thought to apply to unions
# of copies of de Grey's G. The theorem is correct, but its hypothesis asks
# that W - a - b be (k-1)-colourable for every non-adjacent cross pair, which
# for copies of G means G - a must be 4-colourable -- exactly the
# vertex-criticality that is now retracted. G - 1420 is 5-chromatic, so the
# hypothesis fails and the bound rho >= 1357 for G u (G + t) is withdrawn.
# Unions are back on the table, and monotonicity says they can only help.
#
# Verified where rho can be computed exactly by brute force: the Moser spindle
# has rho = 7 at four colours, and embedding it in a 13-vertex graph built by
# rotating it keeps rho at 7 with the same witness -- never rising.

RHO_IS_MONOTONE = (
    "rho(W) <= rho(G) whenever G sits inside W, since colourings of W "
    "restrict; so the route to a small rho is a larger graph, and the "
    "cross-pair bound that seemed to forbid unions needed the criticality "
    "that is now retracted"
)


# -- unioning copies drives rho down to k + 1 -----------------------------
#
# Monotonicity says rho(W) <= rho(G) when G sits inside W, but <= is not <,
# and the whole strategy rests on the difference. Measured at four colours,
# where rho is cheap: Sa unioned with rotated copies of itself through the
# half-Moser angle, each graph strictly containing the last, minimal forcing
# set computed the same way at every size.
#
#     copies      n       m    minimal rho
#          1    397    1974          7
#          2    619    3324          5
#          3    829    4638          5
#          7   1645    9558          5
#
# It DROPS, from 7 to 5, and then holds. Five is one above the floor, since
# rho >= k always. So enlarging works, and it works immediately -- one extra
# copy captures whatever rigidity there is, and further copies add nothing.
#
# The five vertices are the interesting part. On Sa u rot(Sa) they are
# [107, 208, 502, 568, 618]: five points spanning seven edges, three
# overlapping triangles, and 3-CHROMATIC -- 208 can take 107's colour and 618
# can take 502's. A 3-chromatic set of five points forces all FOUR colours,
# because of the 614 vertices around it. The same ambient mechanism as the
# seven on Sa, now sharper: rho = k + 1 with a set that is locally
# unremarkable.
#
# If five colours behave the same way, rho would fall to 6 -- an order of
# magnitude below the 63 a core of three needs. That is the experiment the
# whole strategy now turns on, and at five colours it costs one hard UNSAT:
# assuming every selector is exactly the 4-colourability instance.

UNION_DRIVES_RHO_TO_K_PLUS_ONE = {
    "k": 4,
    "trend": {1: 7, 2: 5, 3: 5, 7: 5},
    "floor": "rho >= k, so 5 is one above it",
    "witness": [107, 208, 502, 568, 618],
    "witness_edges": 7,
    "witness_chi": 3,
    "reading": "a 3-chromatic set of five points forces all four colours; the "
               "rigidity is entirely ambient",
}


# -- but the drop is not a general pattern, and one guard is needed -------
#
# rho -> k + 1 under unioning is a fact about Sa, not about unioning. Measured
# on two other bases, with the same method:
#
#   triangular patch, k = 3:  rho = 3 at one copy -- already the FLOOR, since
#     rho = k exactly when the graph is uniquely k-colourable, and the lattice
#     is: its 3-colouring is the coset colouring and nothing else. Unions have
#     nothing left to give.
#
#   Moser spindle, k = 4:  rho = 7 = n at every size. It is 4-vertex-critical,
#     so rho = n, and the copies add six to fifteen vertices -- too little
#     ambient structure to change anything.
#
# So the drop needs a graph in BETWEEN: not uniquely colourable, where rho is
# already minimal, and not critical, where rho is n. Sa at four colours is
# exactly that, and de Grey's G at five is not critical either -- which is
# what makes the measurement worth the hard solve.
#
# THE GUARD. rho means nothing unless proper k-colourings exist. Without that
# check a union that stops being k-colourable reports every set as forcing for
# want of a counterexample: the triangular patch returned rho = 1 at three
# colours, below the floor rho >= k, because two copies at the Moser angle are
# already 4-chromatic. That is the spindle construction, arriving as a bug.
# Checked on the result it threatened: every Sa union up to seven copies and
# 1645 vertices IS 4-colourable, so the 7 -> 5 drop is real.

RHO_DROP_NEEDS_THE_MIDDLE = {
    "uniquely_colourable": "rho = k already, nothing to gain "
                           "(triangular patch at k=3, rho=3)",
    "vertex_critical": "rho = n, nothing helps "
                       "(Moser spindle at k=4, rho=7=n at every size)",
    "in_between": "Sa at k=4, rho 7 -> 5 under unioning",
    "vacuity_guard": "a union that stops being k-colourable reports every set "
                     "as forcing; the triangular patch gave rho = 1, below "
                     "the floor, because two copies at the Moser angle are "
                     "4-chromatic",
    "checked": "every Sa union to 1645 vertices is 4-colourable, so the drop "
               "is real",
}


# -- what a small rho would actually hand over ---------------------------
#
# A core of size r at p means N(p) u T is forcing. So a forcing set S of ANY
# shape can be padded into one: pick a pivot p, set T = S minus N(p) minus p,
# and N(p) u T contains S, hence forces. The core at p is then |T| <= |S|.
#
#     rho = 6  =>  a core of at most 6 at EVERY pivot,
#                  and of 3 at any pivot adjacent to three of the six.
#
# Which is the whole argument: cores of three are blockable, 27344
# configurations do it, and the counting bound caps r at 3 from the other
# side. So rho = 6 on a 5-chromatic unit-distance graph would close it.
#
# That is exactly what unioning delivered at four colours. rho = 5 = k + 1 on
# Sa u rot(Sa), verified independently of the code that found it: the CNF
# rebuilt from scratch, all four colours checked, all five single deletions
# breaking it, the union confirmed 4-colourable so the question is not vacuous
# -- and the five points inducing seven edges with chromatic number 3.
#
# Whether five colours behave the same way is the open measurement. It costs
# one hard UNSAT, and unavoidably so: any proof that a set forces in a
# 5-chromatic graph contains a proof that the graph is not 4-colourable, since
# a 4-colouring would leave the fifth colour off every set at once.

SMALL_RHO_CLOSES_IT = {
    "padding": "N(p) u T contains S for T = S minus N(p), so a forcing set of "
               "any shape gives a core of at most |S| at every pivot",
    "rho_6_gives": "a core of at most 6 everywhere, and of 3 at any pivot "
                   "adjacent to three of the six -- which is blockable",
    "achieved_at_k4": {"graph": "Sa u rot(Sa)", "n": 619, "rho": 5,
                       "witness": [107, 208, 502, 568, 618],
                       "verified": "independent CNF, four colours, all five "
                                   "deletions, non-vacuity confirmed"},
    "cost_at_k5": "one hard UNSAT, unavoidably: forcing in a 5-chromatic "
                  "graph implies not 4-colourable",
}


# -- and it is not a fact about Sa: all three pieces do it ---------------
#
# Each of de Grey's 4-chromatic pieces, alone and with ONE rotated copy
# through the half-Moser angle, rho measured identically and guarded against
# vacuity:
#
#     Sa alone   397 vertices   rho = 7      Sa + copy    619   rho = 5
#     Sb alone   397            rho = 7      Sb + copy    619   rho = 5
#     Y  alone   791            rho = 7      Y  + copy   1233   rho = 5
#
# Three independent bases, all exactly 7 alone and 5 with a copy -- k + 1 in
# every case. The drop is a property of unioning a non-critical, not-uniquely-
# colourable graph with a rotated copy of itself, not a property of Sa.
#
# A PREDICTION WITH A CONSEQUENCE. If de Grey's G behaves the same way at five
# colours, rho(G,5) is single figures. Padding then gives a core of at most
# that size at EVERY pivot -- and the forward core construction found no core
# in 60 steps at four pivots of the folded union. Both cannot be right.
#
# The greedy construction is the likelier to be wrong: its own docstring notes
# that no single vertex has to raise the objective even when a pair would, so
# it sits on plateaus with genuine cores above it. But the alternative is that
# rho(G,5) is large and the four-colour pattern simply does not transfer. The
# measurement settles it, and costs one hard UNSAT.

RHO_DROP_ON_ALL_THREE = {
    "alone": {"Sa": 7, "Sb": 7, "Y": 7},
    "with_one_rotated_copy": {"Sa": 5, "Sb": 5, "Y": 5},
    "k": 4,
    "reading": "k + 1 in every case, so the drop is about unioning, not Sa",
    "tension": "if it transfers to five colours, padding gives tiny cores "
               "everywhere -- but cegar_core found none in 60 steps on the "
               "folded union, so one of the two is wrong",
}


# -- the tension resolves, and against the optimistic reading ------------
#
# Padding is a proof, not a search: on Sa u rot(Sa) at four colours, with
# rho = 5 and S = [107, 208, 502, 568, 618], the set T = S minus N(p) minus p
# is a core of p for any pivot p NOT IN S. Verified:
#
#     pivot   0 (deg 36)   padding gives a core of 4, confirmed;  cegar finds 7
#     pivot  25 (deg 22)   a core of 5, confirmed;                cegar finds 10
#     pivot  26 (deg 22)   a core of 5, confirmed;                cegar finds 9
#     pivot  27 (deg 22)   a core of 5, confirmed;                cegar finds 10
#     pivot 107 (in S)     padding gives 1, NOT a core -- correctly
#     pivot 208 (in S)     padding gives 2, NOT a core -- correctly
#
# Two things. The padding argument needs p outside S, since otherwise
# N(p) u T no longer contains S; the last two rows are that condition failing
# and being caught rather than believed.
#
# And the greedy core builder WORKS. It finds cores of 7 to 10 where 4 and 5
# exist -- suboptimal, never failing. So its verdict on the folded union at
# five colours, no core at all in 60 steps at four pivots, is meaningful
# rather than an artefact of the method. Which resolves the tension the wrong
# way: the four-colour pattern probably does NOT transfer, and rho(G,5) is
# probably large.
#
# Stated as probable, not settled: the direct measurement is what decides it.

CEGAR_IS_SUBOPTIMAL_NOT_BLIND = {
    "padding_needs": "p outside S, else N(p) u T no longer contains S",
    "verified_cores": {0: 4, 25: 5, 26: 5, 27: 5},
    "cegar_found": {0: 7, 25: 10, 26: 9, 27: 10},
    "reading": "it overshoots but never fails, so 'no core in 60 steps' on "
               "the folded union at five colours is evidence, not an artefact",
    "consequence": "the four-colour pattern probably does not transfer and "
                   "rho(G,5) is probably large -- probable, not settled",
}


# -- two routes to a small rho, and only one can reach 63 ----------------
#
# STRUCTURAL. A k-chromatic subgraph uses all k colours in every k-colouring,
# so rho <= its size. At four colours the smallest unit-distance one is the
# Moser spindle: seven vertices. rho(Sa,4) = 7 sits exactly on that bound.
#
# AMBIENT. Unioning with a rotated copy takes rho to 5, BELOW the structural
# bound, with a set containing no 4-chromatic subgraph at all -- five points,
# seven edges, chromatic number 3.
#
#     Sa            397 vertices   4-chromatic subgraph 30   rho = 7
#     Sa u rot(Sa)  619 vertices   4-chromatic subgraph 11   rho = 5
#
# (Both subgraph figures come from an UNSAT core and are not minimal -- the
# true minimum is the spindle's 7 -- but rho is below even the core's answer
# in the second row, and below the spindle's 7 as well.)
#
# Only the ambient route matters at five colours. There the structural bound
# is the smallest 5-chromatic unit-distance graph, around five hundred
# vertices in the published record, which is eight times the 63 a core of
# three needs. Ambient rigidity is the only mechanism that could close that
# gap, and at four colours it does beat the structural bound -- by 7/5 on the
# union, and by 30/7 against what the same core method finds.

TWO_ROUTES_TO_SMALL_RHO = {
    "structural": "rho <= |smallest k-chromatic subgraph|; at k=4 that is the "
                  "Moser spindle's 7, and rho(Sa,4) = 7 sits on it",
    "ambient": "unioning takes rho to 5, below the structural bound, with a "
               "set that is 3-chromatic and nearly edgeless",
    "at_k5": "the structural bound is around 500 vertices, eight times the "
             "63 needed, so only the ambient route could ever close it",
}


# -- the pivot's circle cannot be enriched: G is saturated ---------------
#
# A core of three at p needs N(p) u T forcing, and supersets of forcing sets
# force -- so a BIGGER circle makes the condition easier. de Grey's G has
# maximum degree 60 while its edge module carries 134 unit steps, so the
# circle looks like it has room to double.
#
# It has none. Of the 74 module steps from the pivot that do not already land
# on a vertex of G, every single one lands on a point adjacent to NOTHING ELSE
# -- contact count exactly 1, all 74. They are pendants: free, unconstrained,
# and worth nothing to the pressure. Measured:
#
#     +  0 circle points   degree  60   pressure at k=5 = 2
#     + 10                 degree  70   pressure = 2
#     + 30                 degree  90   pressure = 2
#     + 74                 degree 134   pressure = 2
#
# Degree more than doubles and the pressure does not move. So de Grey's G
# already contains every constrained point on its pivot's circle: within its
# own module the construction is saturated, and enriching the circle is not a
# lever that exists there. It would need points outside the module, which
# means a larger field.

CIRCLE_IS_SATURATED = {
    "pivot_degree": 60,
    "module_unit_steps": 134,
    "new_points_available": 74,
    "contacts_each": 1,
    "pressure_at_degree": {60: 2, 70: 2, 90: 2, 134: 2},
    "reading": "every new circle point is a pendant, so G already holds every "
               "constrained point on its pivot's circle",
}
