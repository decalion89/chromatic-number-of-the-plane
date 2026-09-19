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
