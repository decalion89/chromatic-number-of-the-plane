"""The 5-chromatic unit-distance graphs in Q(sqrt3, sqrt11, sqrt247).

Everything here is recomputed from the JSON files alone -- the exact
coordinates -- so the files are the objects under test, not the scripts that
made them.  A transcription error shows up as a wrong edge count or a graph
that turns out to be 4-colourable.

The construction is a ladder of two glues:

  * Sa, de Grey's 39 points under the 12-element group <rot60, reflect>, has
    NO pair of vertices forced to share a colour in every 4-colouring.
  * Glue Sa to a rotated image of itself about a VERTEX.  The glue circle is
    that vertex's own unit circle -- not capped, and the rotation is free.
    The union has forced-equal pairs, all at squared distance 64/9.  The glue
    manufactured them out of a carrier that had none.
  * Spindle one.  The rotation about the pivot by 2*arcsin(3/16) carries the
    far endpoint to distance exactly 1 from itself; the pivot is fixed, so both
    copies agree there, forced equality carries that to the endpoint and its
    image, and those two are adjacent.

The angle needs sqrt(2223 * 16384) = 384 sqrt(247), and 247 = 13 * 19 -- a
radical de Grey's field does not contain.

Three sizes are recorded, from the same recipe with more care taken each time:
1139 from Sa itself, 951 from Sa peeled to its 327 highest-degree vertices, and
803 after cutting the union down to the part that still forces.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as Fr

import pytest

from hn.coloring import is_k_colorable
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

HERE = os.path.dirname(__file__)
FIELD = Field((3, 11, 247))

# file, vertices, edges
GRAPHS = [
    ("five_247.json", 1139, 6475),
    ("five_247_b.json", 951, 5171),
    ("five_247_c.json", 803, 4065),
]


def _load(name):
    with open(os.path.join(HERE, os.pardir, "data", name)) as fh:
        d = json.load(fh)
    pts = [Point(FIELD.element([Fr(a, b) for a, b in x]),
                 FIELD.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


@pytest.fixture(scope="module")
def loaded():
    return {name: _load(name) for name, _, _ in GRAPHS}


@pytest.mark.parametrize("name,n,m", GRAPHS)
def test_vertices_and_edges_recomputed_exactly(loaded, name, n, m):
    d, pts = loaded[name]
    assert tuple(d["field_generators"]) == (3, 11, 247)
    assert d["n"] == n
    assert len(pts) == n and len(set(pts)) == n
    g = build_graph(pts)
    assert g.n == n
    assert sum(len(a) for a in g.adj) // 2 == m
    for u, v in g.edges():
        assert pts[u].dist2(pts[v]) == 1


@pytest.mark.parametrize("name,n,m", GRAPHS)
def test_chromatic_number_is_five(loaded, name, n, m):
    _, pts = loaded[name]
    g = build_graph(pts)
    assert is_k_colorable(g, 4, timeout=1800)[0] is False
    ok, colouring = is_k_colorable(g, 5, timeout=1800)
    assert ok is True
    for u, v in g.edges():
        assert colouring[u] != colouring[v]


@pytest.mark.parametrize("name,n,m", GRAPHS)
def test_does_not_embed_in_de_greys_field(loaded, name, n, m):
    """Squared distances are invariant under every isometry of the plane, so
    the field they generate is an invariant of the graph and not of this
    drawing of it.  Each graph has a squared distance with a
    sqrt(741) = sqrt3 sqrt247 component, and Q(sqrt3,sqrt5,sqrt7,sqrt11)
    contains no square root of 13*19.  So no congruent copy of any of them
    lives in the field de Grey's G needs."""
    _, pts = loaded[name]
    carries = [i for i in range(FIELD.dim)
               if FIELD._prod[i] % 13 == 0 or FIELD._prod[i] % 19 == 0]
    assert carries
    for i in range(len(pts)):
        for j in range(i + 1, min(i + 40, len(pts))):
            if any(pts[i].dist2(pts[j]).c[k] != 0 for k in carries):
                return
    pytest.fail("no squared distance carries sqrt(247)")


def test_the_spindle_closes(loaded):
    """The three vertices the first construction turns on: the pivot is at
    distance 8/3 from the forced-equal endpoint, and that endpoint is one unit
    from its own image.  Those two are the contradiction."""
    d, pts = loaded["five_247.json"]
    sp = d["spindle"]
    pivot, u, v = sp["pivot"], sp["forced_pair_endpoint"], sp["its_image"]
    assert pts[pivot].dist2(pts[u]) == Fr(64, 9)
    assert pts[u].dist2(pts[v]) == 1
    assert v in build_graph(pts).adj[u]


def test_the_spindle_angle_needs_247():
    """cos t = 1 - 1/(2*(64/9)) = 119/128 and sin t = 384*sqrt(247)/16384.
    Nothing smaller than sqrt(13*19) will do."""
    cos = Fr(1) - Fr(1, 2) / Fr(64, 9)
    assert cos == Fr(119, 128)
    r = (1 - cos * cos)
    n = r.numerator * r.denominator
    d = 2
    while d * d <= n:
        while n % (d * d) == 0:
            n //= d * d
        d += 1
    assert n == 247


# ---------------------------------------------------------------------------
# The one with a group of its own.  Every graph above has isometry group of
# order 1, because the glue rotation is about one vertex and the spindle about
# one pivot.  Doing both over whole orbits instead -- six glues at once, then
# six spindles at once -- keeps the union invariant, since g rot_w g^-1 =
# rot_{g(w)} and Sa is already fixed by the group.  The carrier that results
# has mean degree 13.34 and 153 forced pairs where the sequential chain
# reached 52, and its spindle is the first 5-chromatic unit-distance graph
# here that carries a symmetry.

SYM = os.path.join(HERE, os.pardir, "data", "five_symmetric.json")


@pytest.fixture(scope="module")
def symmetric():
    with open(SYM) as fh:
        d = json.load(fh)
    pts = [Point(FIELD.element([Fr(a, b) for a, b in x]),
                 FIELD.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


def test_symmetric_graph_size(symmetric):
    d, pts = symmetric
    assert d["n"] == 7141 and d["m"] == 47682
    assert len(pts) == 7141 and len(set(pts)) == 7141
    g = build_graph(pts)
    assert g.n == 7141
    assert sum(len(a) for a in g.adj) // 2 == 47682


def test_symmetric_graph_is_c6_invariant(symmetric):
    """The group is about the origin by construction, so checking the
    60-degree rotation maps the set onto itself is O(n), where searching for an
    unknown centre would be O(n^2).  The reflection is not in it: the group is
    C6, of order 6, not D6."""
    from hn.geometry import _rot60
    _, pts = symmetric
    S = set(pts)
    r = _rot60(FIELD)
    assert all(r(p) in S for p in pts)
    assert all(r(r(p)) in S for p in pts)
    assert not all(Point(p.x, -p.y) in S for p in pts)


def test_symmetric_graph_refuses_four_colours(symmetric):
    _, pts = symmetric
    g = build_graph(pts)
    for u, v in list(g.edges())[:200]:
        assert pts[u].dist2(pts[v]) == 1
    assert is_k_colorable(g, 4, timeout=1800)[0] is False


# ---------------------------------------------------------------------------
# The other forced orbit, spent alone.  The symmetric carrier's 153 forced
# pairs fall in two orbits; the second, at squared distance 64/3, has
# cos 125/128 and sin sqrt(759)/128 with 759 = 3*11*23, so its spindle lives in
# Q(sqrt3, sqrt11, sqrt23).  Same vertex count as the first and six fewer
# edges: a different graph, not a redrawing.

F23 = Field((3, 11, 23))
SRC23 = os.path.join(HERE, os.pardir, "data", "five_23.json")


@pytest.fixture(scope="module")
def graph23():
    with open(SRC23) as fh:
        d = json.load(fh)
    pts = [Point(F23.element([Fr(a, b) for a, b in x]),
                 F23.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


def test_third_field_size_and_edges(graph23):
    d, pts = graph23
    assert tuple(d["field_generators"]) == (3, 11, 23)
    assert d["n"] == 7141 and d["m"] == 47676
    g = build_graph(pts)
    assert g.n == 7141
    assert sum(len(a) for a in g.adj) // 2 == 47676
    for u, v in list(g.edges())[:200]:
        assert pts[u].dist2(pts[v]) == 1


def test_third_field_is_c6_invariant(graph23):
    from hn.geometry import _rot60
    _, pts = graph23
    S = set(pts)
    r = _rot60(F23)
    assert all(r(p) in S for p in pts)


def test_third_field_refuses_four_colours(graph23):
    _, pts = graph23
    assert is_k_colorable(build_graph(pts), 4, timeout=1800)[0] is False


def test_the_second_spindle_angle_needs_23():
    """cos t = 1 - (1/2)/(64/3) = 125/128, sin t = sqrt(759)/128, and
    759 = 3*11*23, so this rotation is outside Q(sqrt3,sqrt11,sqrt247)."""
    cos = Fr(1) - Fr(1, 2) / Fr(64, 3)
    assert cos == Fr(125, 128)
    r = 1 - cos * cos
    n = r.numerator * r.denominator
    d = 2
    while d * d <= n:
        while n % (d * d) == 0:
            n //= d * d
        d += 1
    assert n == 759 == 3 * 11 * 23
