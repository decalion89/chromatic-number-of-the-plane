"""The dense 5-chromatic graphs, recomputed from their stored coordinates.

Until the tuned chain, the two halves of this project could not be combined.
Stacking orbits of glue centres takes a carrier to mean degree 18.47 and every
such carrier is 4-colourable in under a second -- density never bought a
chromatic number.  Spending a spindle buys one, but the spindle adds a sparse
second copy, so every 5-chromatic graph here sat between 10.1 and 13.8.

Tuning the composite forced pair to distance exactly 1 removes the spindle: the
pair becomes an EDGE, so the union of the carrier with one rotated copy already
refuses four, and that union is two copies of the carrier and keeps its degree.

    cos phi = 1/(2 D^2) - 1 = -119/128        sin phi = 3 sqrt(247)/128

with D^2 = 64/9 the forced distance -- and sqrt(247) is the radical the
ordinary spindle at 64/9 needs anyway, so the whole construction stays inside
the carrier's own field Q(sqrt3, sqrt11, sqrt247).

These two files are the result at 16.81 and 18.48.  They are the densest
5-chromatic unit-distance graphs in the project, and the reason they matter is
negative: free@5 measures 16.40 % and 13.80 % on them, against 15.73 % at
degree 10.12, so an 83 % rise in mean degree moves it not at all.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as Fr

import pytest

from hn.coloring import is_k_colorable
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

HERE = os.path.dirname(__file__)
FIELD = Field((3, 11, 247))
D2 = Fr(64, 9)

# file, n, m, mean degree to two places, measured free@5
GRAPHS = [
    ("five_dense_2.json", 6925, 58213, 16.81, 0.16404),
    ("five_dense_10.json", 12469, 115189, 18.48, 0.13802),
]


def _load(name):
    with open(os.path.join(HERE, os.pardir, "data", name)) as fh:
        d = json.load(fh)
    pts = [Point(FIELD.element([Fr(a, b) for a, b in x]),
                 FIELD.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


def test_the_tuning_angle_is_exact_and_lands_on_one():
    """cos phi is rational, sin phi carries sqrt247, and the composite pair
    lands at distance exactly 1 -- which is why no spindle is needed."""
    cos = FIELD.rational(Fr(1) / (2 * D2) - 1)
    assert cos == FIELD.rational(Fr(-119, 128))
    sin = FIELD.sqrt(247) * FIELD.rational(Fr(3, 128))
    assert cos * cos + sin * sin == FIELD.rational(Fr(1))
    u = Point(FIELD.rational(Fr(8, 3)), FIELD.zero())     # |u|^2 = 64/9
    assert u.norm2() == FIELD.rational(D2)
    r = Rotation(cos, sin)(u)
    assert Point(u.x + r.x, u.y + r.y).norm2() == FIELD.rational(Fr(1))


@pytest.mark.parametrize("name,n,m,deg,free5", GRAPHS)
def test_counts_match_the_stored_coordinates(name, n, m, deg, free5):
    d, pts = _load(name)
    assert len(pts) == len(set(pts)), "duplicate points in the file"
    g = build_graph(pts)
    assert g.n == n == d["n"]
    edges = sum(len(a) for a in g.adj) // 2
    assert edges == m == d["m"]
    assert round(2.0 * edges / g.n, 2) == deg


@pytest.mark.parametrize("name,n,m,deg,free5", GRAPHS)
def test_they_are_two_copies_of_one_carrier(name, n, m, deg, free5):
    """The chain doubles the vertex count exactly once, minus the shared pair.

    tau carries v to q, so the two copies already share that point; everything
    else is new.  An odd count would mean the copies had merged somewhere and
    the construction was not what it claims.
    """
    d, pts = _load(name)
    assert d["n"] % 2 == 1, "two copies sharing one point is odd"


@pytest.mark.parametrize("name,n,m,deg,free5", GRAPHS)
def test_they_refuse_four(name, n, m, deg, free5):
    """The whole claim.  The composite pair is forced equal and is an edge.

    is_k_colorable returns a PAIR, (answer, colouring).  `not` on a two-element
    tuple is always False, so `assert not is_k_colorable(g, 4)` fails whatever
    the solver says -- which is how this test first read, and it reported the
    graphs as broken when the bug was here.
    """
    _, pts = _load(name)
    g = build_graph(pts)
    answer, colouring = is_k_colorable(g, 4)
    assert answer is False, f"expected no 4-colouring, solver said {answer}"
    assert colouring is None


@pytest.mark.parametrize("name,n,m,deg,free5", GRAPHS)
def test_free_at_five_did_not_move(name, n, m, deg, free5):
    """Recorded, not recomputed -- a 5-colouring of these costs 2906 s and
    8921 s.  The point is the flatness: 15.73 % at degree 10.12, and still
    between 12.5 % and 16.4 % at 16.81 and 18.48."""
    assert 0.12 < free5 < 0.17
