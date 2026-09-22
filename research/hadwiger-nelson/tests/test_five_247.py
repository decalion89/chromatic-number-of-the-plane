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
807 after cutting the union down to the part that still forces.
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
    ("five_247_c.json", 807, 4091),
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
