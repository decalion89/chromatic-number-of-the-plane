"""Forced colour relations, and the ladder chi >= 6 reduces to."""
import pytest

from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.forced import (ColourRelations, centre_pressure, find_rainbow,
                       min_colours_on)
from hn.graph import build_graph
from hn.mixed import joint_core_union


def test_edges_are_forced_different():
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        assert rel.colourable
        for u, v in list(g.edges())[:20]:
            assert rel.different(u, v)
            assert not rel.can_share(u, v)
    finally:
        rel.close()


def test_min_colours_recovers_the_chromatic_number():
    """On the whole vertex set the measurement is just chi."""
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        assert min_colours_on(rel, range(g.n)) == 4
        assert min_colours_on(rel, []) == 0
        assert min_colours_on(rel, [0]) == 1
    finally:
        rel.close()


def test_vacuous_queries_are_refused():
    """A graph with no k-colouring makes every query UNSAT; say so."""
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 3)
    try:
        assert not rel.colourable
    finally:
        rel.close()
    with pytest.raises(ValueError):
        find_rainbow(g, 3, 3)


def test_unit_circle_graphs_are_bipartite():
    """Two points of a unit circle are adjacent only at 60 degrees, and the
    60-degree orbit closes after six steps, so every component of a circle
    graph is a path or a 6-cycle and no circle ever needs three colours on its
    own account. Any third colour has to be forced by the ambient graph --
    which is why every neighbourhood of de Grey's G comes back at exactly two,
    however large: the two hubs of degree 60 included.
    """
    for g in (build_graph(joint_core_union()), build_G()):
        for v in range(g.n):
            circle = set(g.adj[v])
            colour = {}
            for a in sorted(circle):
                if a in colour:
                    continue
                colour[a], stack = 0, [a]
                while stack:
                    x = stack.pop()
                    for y in g.adj[x] & circle:
                        if y not in colour:
                            colour[y] = 1 - colour[x]
                            stack.append(y)
                        else:
                            assert colour[y] != colour[x], \
                                f"circle of {v} is not bipartite at {x},{y}"


@pytest.mark.slow
def test_degrey_G_has_no_local_pressure_at_five():
    """Every one of the 1581 neighbourhoods is 2-colourable at k=5.

    The ladder needs three: at three the pivot's colour is confined to two,
    which is the joint core the two-orbit block closes. This pins the measured
    fact that G supplies none of it.
    """
    g = build_G()
    rel = ColourRelations(g, 5)
    try:
        assert rel.colourable
        seen = {min_colours_on(rel, sorted(g.adj[v])) for v in range(g.n)}
        assert seen == {2}
    finally:
        rel.close()


def test_centre_pressure_reads_the_ladder():
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        hub = max(range(g.n), key=lambda v: len(g.adj[v]))
        r = centre_pressure(rel, sorted(g.adj[hub]))
        assert r["pressure"] + r["free_for_centre"] == 4
        assert r["smallest_possible_core"] == 4 - r["pressure"]
    finally:
        rel.close()


# -- the pressure theorem -------------------------------------------------

def test_pressure_bound_matches_the_construction():
    """A core of r needs pressure >= k - r; the built core of two saturates it.

    On the four points at three colours the pivot sees one neighbour, so
    pressure is 1 = k - 2, the theorem permits a core of two, and a separation
    query confirms one. The bound is tight exactly where the construction is.
    """
    from hn.forced import max_core_bound, pressure
    from hn.mixed import joint_core_configuration
    from hn.spindle import SeparationTest

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    ia, ib = g.vertices.index(pts[2]), g.vertices.index(pts[3])
    rel = ColourRelations(g, 3)
    try:
        assert pressure(rel, bp) == 1
        assert max_core_bound(rel, bp) == 2
    finally:
        rel.close()
    st = SeparationTest(g, 3, bp, [ia, ib])
    try:
        assert not st.run(subset=[ia, ib])[0]
    finally:
        st.close()


def test_pressure_forbids_small_cores_where_it_is_low():
    """The contrapositive, checked directly: squeeze the neighbourhood into
    k - r - 1 colours and every set of r vertices fails to be a core."""
    from itertools import combinations

    from hn.forced import pressure
    from hn.spindle import SeparationTest

    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        hub = max(range(g.n), key=lambda v: len(g.adj[v]))
        assert pressure(rel, hub) == 1          # so no core of size 1 or 2
    finally:
        rel.close()
    others = [v for v in range(g.n) if v != hub and v not in g.adj[hub]]
    st = SeparationTest(g, 4, hub, others)
    try:
        for pair in list(combinations(others, 2))[:40]:
            assert st.run(subset=list(pair))[0], f"{pair} should not be a core"
    finally:
        st.close()


# -- criticality makes a graph inert ---------------------------------------

def test_critical_graphs_have_no_core_at_all():
    """Colour G - p with k-1 and give p the kth: p is then alone in its
    colour, so c(p) lies outside c(T) for EVERY target set at once.

    The Moser spindle is 4-critical, so at four colours every one of its
    vertices is separable from all its non-neighbours simultaneously. This is
    the corollary that accounts for all 880 pivots of fruitless searching on
    de Grey's G, which is 5-critical, before any solver runs.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.spindle import SeparationTest

    pts, _doc = load_certificate("certificates/moser_spindle_no3coloring.json")
    g = build_graph(pts)
    assert not is_k_colorable(g, 3)[0] and is_k_colorable(g, 4)[0]
    for p in range(g.n):
        rest = build_graph([q for j, q in enumerate(pts) if j != p])
        assert is_k_colorable(rest, 3)[0], f"{p} not critical"
        others = [v for v in range(g.n) if v != p and v not in g.adj[p]]
        if not others:
            continue
        st = SeparationTest(g, 4, p, others)
        try:
            assert st.run(subset=others)[0], f"pivot {p} should be separable"
        finally:
            st.close()


def test_union_of_two_copies_is_not_critical_and_has_a_core():
    """The escape the corollary points at: in W = H union f(H), removing a
    vertex leaves a whole chi-chromatic copy, so no vertex is critical, no
    pivot can be given a colour of its own, and the full target set is a core.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.geometry import Point
    from hn.spindle import SeparationTest

    pts, _doc = load_certificate("certificates/moser_spindle_no3coloring.json")
    dx = pts[1].x - pts[0].x
    dy = pts[1].y - pts[0].y
    both = list(dict.fromkeys(pts + [Point(q.x + dx, q.y + dy) for q in pts]))
    w = build_graph(both)
    assert is_k_colorable(w, 4)[0]
    for v in range(w.n):
        rest = build_graph([q for j, q in enumerate(both) if j != v])
        if is_k_colorable(rest, 3)[0]:
            break
    else:
        # no vertex is critical, so some pivot must fail to be separable
        hub = max(range(w.n), key=lambda v: len(w.adj[v]))
        others = [v for v in range(w.n) if v != hub and v not in w.adj[hub]]
        st = SeparationTest(w, 4, hub, others)
        try:
            assert not st.run(subset=others)[0], "the full set must be a core"
        finally:
            st.close()


def test_general_block_agrees_with_the_two_sat_one():
    from hn.mixed import (blocks_targets, blocks_two_targets,
                          joint_core_configuration, joint_core_copies)

    E, pts, rho, sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    pair = [g.vertices.index(pts[2]), g.vertices.index(pts[3])]
    copies = joint_core_copies(E, rho, sigma)
    assert blocks_targets(g, bp, pair, copies) is True
    assert blocks_two_targets(g, bp, pair, copies) is True
    # one orbit alone cannot close it: three copies leave an escape
    assert not blocks_targets(g, bp, pair, copies[:3])


# -- the alignment condition ----------------------------------------------

def test_core_must_be_a_rainbow_disjoint_from_the_neighbourhood():
    """When pressure(p) = k - r exactly, a core of size r has a shape.

    In a colouring attaining the minimum the pivot has exactly r free colours
    and can take any of them, so c(T) must contain all r while holding at most
    r -- hence c(T) = free(p): the targets are a rainbow AND none of them ever
    shares a colour with a neighbour of the pivot. The built pair satisfies
    both, and its legs are the two non-pivot corners of a unit triangle, so
    they are forced apart by adjacency.
    """
    from hn.forced import core_must_be_rainbow, free_colours, pressure
    from hn.mixed import joint_core_configuration

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    legs = [g.vertices.index(pts[2]), g.vertices.index(pts[3])]
    rel = ColourRelations(g, 3)
    try:
        assert pressure(rel, bp) == 1
        assert free_colours(rel, bp) == 2 == len(legs)
        assert not rel.can_share(*legs)
        assert core_must_be_rainbow(rel, bp, legs)
        # a repeated target is one distinct point, fewer than the free
        # colours, so it cannot be a core and the test says so
        assert not core_must_be_rainbow(rel, bp, [legs[0], legs[0]])
    finally:
        rel.close()


# -- cores built forwards --------------------------------------------------

def test_core_condition_is_a_pressure_measurement():
    """T is a core of p exactly when N(p) union T uses every colour.

    If some colouring left a colour free on both, recolouring p to it stays
    proper -- properness at p asks only that its colour avoid c(N(p)) -- and
    puts c(p) outside c(T). The two characterisations must agree with the
    separation test on the built construction, where the answer is known.
    """
    from hn.forced import is_core, min_colours_on
    from hn.mixed import joint_core_configuration
    from hn.spindle import SeparationTest

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    y, z = g.vertices.index(pts[2]), g.vertices.index(pts[3])
    rel = ColourRelations(g, 3)
    try:
        assert is_core(rel, bp, [y, z])
        assert not is_core(rel, bp, [y])
        assert not is_core(rel, bp, [z])
        circle = sorted(g.adj[bp])
        assert min_colours_on(rel, circle + [y, z]) == 3
        assert min_colours_on(rel, circle + [y]) < 3
    finally:
        rel.close()
    st = SeparationTest(g, 3, bp, [y, z])
    try:
        assert not st.run(subset=[y, z])[0]      # forced, i.e. a core
        assert st.run(subset=[y])[0]             # not forced alone
        assert st.run(subset=[z])[0]
    finally:
        st.close()


def test_cegar_core_finds_the_built_pair():
    from hn.forced import cegar_core, is_core
    from hn.mixed import joint_core_configuration

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    rel = ColourRelations(g, 3)
    try:
        T, ok = cegar_core(rel, bp)
        assert ok and is_core(rel, bp, T)
        assert set(T) == {g.vertices.index(pts[2]), g.vertices.index(pts[3])}
    finally:
        rel.close()


@pytest.mark.slow
def test_cegar_core_fails_on_a_critical_graph():
    """A k-vertex-critical graph has no core, so the construction must run out.

    This is the criticality corollary as an executable statement: de Grey's G
    is 5-vertex-critical, so no amount of counterexample-killing can ever
    close a core there, and the builder has to come back empty.
    """
    from hn.forced import cegar_core

    g = build_G()
    rel = ColourRelations(g, 5)
    try:
        _T, ok = cegar_core(rel, 0, limit=40)
        assert not ok, "a critical graph must yield no core"
    finally:
        rel.close()


# -- the mechanism behind the only pressure 3 ------------------------------

def test_pressure_gadget_is_an_odd_cycle_against_k_minus_two():
    """Take the 47-vertex witness apart and check the machine causally.

    The circle splits into five hexagons, every one of the sixteen further
    points is adjacent to exactly two circle points lying in DIFFERENT
    hexagons, and those sixteen form a graph that is 3-chromatic but not
    bipartite. Squeezing the circle into two colours confines them to k-2, so
    the odd cycle bites at four colours and not at five -- and deleting three
    of them to kill the odd cycle drops the pressure from 3 to 2, which is the
    causal half rather than a coincidence of numbers.
    """
    import itertools

    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.forced import PRESSURE_GADGET, min_colours_on

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    circle = set(g.adj[piv])
    extras = [v for v in range(g.n) if v != piv and v not in circle]
    assert len(circle) == 30 and len(extras) == PRESSURE_GADGET["confined_points"]

    comps, seen = [], set()
    for s in sorted(circle):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in g.adj[x] & circle:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(set(comp))
    assert len(comps) == 5 and all(len(c) == 6 for c in comps)

    for v in extras:
        nb = [c for c in circle if c in g.adj[v]]
        assert len(nb) == 2
        assert sum(1 for c in comps if nb[0] in c) == 1
        assert not any(nb[0] in c and nb[1] in c for c in comps), \
            "both circle neighbours in one hexagon would read nothing"

    sub = build_graph([g.vertices[v] for v in extras])
    assert (sub.n, sub.m) == PRESSURE_GADGET["confined_graph"]
    assert is_k_colorable(sub, 3)[0] and not is_k_colorable(sub, 2)[0]

    def pressure_of(keep, k):
        h = build_graph([g.vertices[v] for v in sorted(keep)])
        j = h.vertices.index(g.vertices[piv])
        rel = ColourRelations(h, k)
        try:
            return min_colours_on(rel, sorted(h.adj[j])) if rel.colourable else None
        finally:
            rel.close()

    assert pressure_of(range(g.n), 4) == 3
    assert pressure_of(range(g.n), 5) == 2
    for drop in itertools.combinations(extras, 3):
        rest = [v for v in extras if v not in drop]
        if is_k_colorable(build_graph([g.vertices[v] for v in rest]), 2)[0]:
            assert pressure_of(set(range(g.n)) - set(drop), 4) == 2, \
                "killing the odd cycle must kill the pressure"
            break
    else:
        raise AssertionError("no three extras make the gadget bipartite")


def test_the_odd_cycle_covers_only_half_the_orientations():
    """The correction to the tempting one-line story.

    "An odd cycle against k-2 colours" is load-bearing but incomplete: over
    the 32 orientations of the five hexagons the confined set is bipartite in
    sixteen of them, so those are killed by something longer-range, through
    the confined points' other neighbours rather than among themselves.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    circle = set(g.adj[piv])
    extras = [v for v in range(g.n) if v != piv and v not in circle]
    seen, comps = set(), []
    for s in sorted(circle):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in g.adj[x] & circle:
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
            for y in g.adj[x] & comp:
                if y not in par:
                    par[y] = 1 - par[x]
                    front.append(y)

    bipartite = 0
    for bits in range(1 << len(comps)):
        def col(v):
            return par[v] ^ ((bits >> comp_of[v]) & 1)
        barred = [v for v in extras
                  if len({col(c) for c in g.adj[v] & circle}) == 2]
        sub = build_graph([g.vertices[v] for v in barred])
        if sub.n == 0 or is_k_colorable(sub, 2)[0]:
            bipartite += 1
    assert bipartite == 16, f"expected half, got {bipartite}/32"


def test_the_mechanism_is_list_colouring():
    """Named completely: the auxiliary graph versus the lists the circle gives.

    Squeezing the pivot's circle into two colours hands every other vertex a
    list -- k-2 colours where it sees two differently-coloured circle points,
    k-1 where it sees two of the same. Over all 32 orientations of the five
    hexagons the auxiliary graph fails to be list-colourable in ALL of them at
    four colours, and in NONE at five, where the same lists grow to sizes 3
    and 4. That is exactly the measured pressure, 3 and then 2.
    """
    from hn.certify import load_certificate
    from hn.forced import circle_hexagons, induced_lists, list_colourable

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    comps, _par, _co = circle_hexagons(g, piv)
    assert len(comps) == 5 and all(len(c) == 6 for c in comps)
    extras = [v for v in range(g.n) if v != piv and v not in g.adj[piv]]

    for k, expected in ((4, 32), (5, 0)):
        refused = 0
        sizes = set()
        for bits in range(1 << len(comps)):
            lists, _ = induced_lists(g, piv, k, bits)
            sizes |= {len(lists[v]) for v in extras}
            if not list_colourable(g, extras, lists):
                refused += 1
        assert refused == expected, f"k={k}: {refused} refusals"
        assert sizes == {k - 2, k - 1}
