"""The 1139-vertex 5-chromatic unit-distance graph in Q(sqrt3, sqrt11, sqrt247).

Everything here is recomputed from `data/five_247.json` alone -- the exact
coordinates -- so the file is the object under test, not the script that made
it.  A transcription error shows up as a wrong edge count or a graph that turns
out to be 4-colourable.

The construction is a ladder of two glues:

    Sa   de Grey's 39 points under the 12-element group <rot60, reflect>,
         397 points, and NO pair of them is forced to share a colour in every
         4-colouring.
    H    Sa u rho1(Sa), rho1 the 60-degree rotation about the VERTEX Sa[25].
         The glue circle is the unit circle about that vertex: 20 points, and
         it is not capped.  H has eight forced-equal pairs, all at distance
         8/3.  The glue manufactured them.
    Z    H u rho2(H), rho2 the rotation about H[157] by 2*arcsin(3/16), which
         carries H[327] to distance exactly 1 from itself.  rho2 fixes H[157],
         so both copies agree there; forced equality carries that agreement to
         H[327] and to its image; and those two are adjacent.

The angle needs sqrt(2223*16384) = 384*sqrt(247), and 247 = 13*19.
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

DATA = os.path.join(os.path.dirname(__file__), os.pardir, "data", "five_247.json")
FIELD = Field((3, 11, 247))


@pytest.fixture(scope="module")
def record():
    with open(DATA) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def points(record):
    out = []
    for x, y in record["points"]:
        out.append(Point(FIELD.element([Fr(a, b) for a, b in x]),
                         FIELD.element([Fr(a, b) for a, b in y])))
    return out


@pytest.fixture(scope="module")
def graph(points):
    return build_graph(points)


def test_field_is_the_one_recorded(record):
    assert tuple(record["field_generators"]) == (3, 11, 247)
    assert 247 == 13 * 19


def test_vertex_count(record, points):
    assert record["n"] == 1139
    assert len(points) == 1139
    assert len(set(points)) == 1139


def test_edges_recomputed_exactly(graph):
    assert graph.n == 1139
    assert sum(len(a) for a in graph.adj) // 2 == 6475


def test_every_edge_is_exactly_unit(graph, points):
    for u, v in graph.edges():
        assert points[u].dist2(points[v]) == 1


def test_not_four_colourable(graph):
    ok, _ = is_k_colorable(graph, 4, timeout=1800)
    assert ok is False


def test_five_colourable(graph):
    ok, colouring = is_k_colorable(graph, 5, timeout=1800)
    assert ok is True
    for u, v in graph.edges():
        assert colouring[u] != colouring[v]


def test_does_not_embed_in_de_greys_field(points):
    """Squared distances are invariant under every isometry of the plane, so
    the field they generate is an invariant of the graph and not of this
    drawing of it.  Some squared distance of Z has a sqrt(741) = sqrt3*sqrt247
    component, and Q(sqrt3,sqrt5,sqrt7,sqrt11) contains no square root of
    13*19.  So no congruent copy of Z lives in the field de Grey's G needs."""
    carries_247 = [m for m in range(FIELD.dim)
                   if FIELD._prod[m] % 13 == 0 or FIELD._prod[m] % 19 == 0]
    assert carries_247, "the field must have basis elements involving 247"
    for i in range(len(points)):
        for j in range(i + 1, min(i + 40, len(points))):
            d2 = points[i].dist2(points[j])
            if any(d2.c[m] != 0 for m in carries_247):
                return
    pytest.fail("no squared distance carries sqrt(247): Z would embed after all")


def test_the_spindle_closes(record, points):
    """The three vertices the construction turns on.  Z interleaves H[i] with
    its image, so the indices are recorded rather than assumed: the pivot is
    at distance 8/3 from the forced-equal endpoint, and that endpoint is one
    unit from its own image.  Those two are the contradiction."""
    sp = record["spindle"]
    pivot, u, v = sp["pivot"], sp["forced_pair_endpoint"], sp["its_image"]
    assert points[pivot].dist2(points[u]) == Fr(64, 9)
    assert points[u].dist2(points[v]) == 1
    assert v in build_graph(points).adj[u]


def test_the_spindle_angle_needs_247():
    """cos t = 1 - 1/(2*(64/9)) = 119/128 and sin t = 384*sqrt(247)/16384.
    Nothing smaller than sqrt(13*19) will do."""
    d2 = Fr(64, 9)
    cos = Fr(1) - Fr(1, 2) / d2
    assert cos == Fr(119, 128)
    s2 = 1 - cos * cos
    n = s2.numerator * s2.denominator
    r = n
    d = 2
    while d * d <= r:
        while r % (d * d) == 0:
            r //= d * d
        d += 1
    assert r == 247
