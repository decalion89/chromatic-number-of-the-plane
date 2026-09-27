"""The tuned family: 5-chromatic graphs whose field is a parameter.

Every other construction here takes the field it is given -- the radii pick the
radicals.  This one chooses.

Let H force c(v) = c(q) in every 4-colouring, with |v - q| = D, and let tau be
the isometry carrying v to q whose rotational part is the rotation by phi:

        tau(p) = q + R_phi (p - v)

tau(H) is a congruent copy, so it forces c(tau v) = c(tau q), that is
c(q) = c(tau q).  Forcing is transitive, so H u tau(H) forces c(v) = c(tau q),
a NEW forced pair, and with u = v - q,

        v - tau(q) = u + R_phi u,      |u + R_phi u|^2 = 2 D^2 (1 + cos phi)

The composite distance therefore sweeps [0, 2D] as phi turns, and the angle
that lands on a chosen d is

        cos phi = d^2 / (2 D^2) - 1

which is RATIONAL whenever d^2 and D^2 are.  Only sin phi carries a radical,
and that radical is the entire arithmetic cost.  With D^2 = 64/9, writing
d^2 = p/q,

        cos phi = (9p - 128q) / (128q)
        sin^2 phi = [(128q)^2 - (9p - 128q)^2] / (128q)^2
                  = 9 p (256q - 9p) / (128q)^2

so the cost is sqrt(p (256q - 9p)) and choosing p and q chooses the field.

Two consequences are recorded in data/:

  * d^2 = 1 makes the composite pair ADJACENT, so the chained union already
    refuses four and no spindle is needed -- 2041 points instead of 4061.
    Its radical is sqrt(247), which is also what the ordinary spindle at 64/9
    needs, so tuning to 1 is that spindle in disguise.  A sanity check, not a
    new field.
  * d^2 = 1/3 needs sqrt(759) = sqrt3 sqrt11 sqrt23 and d^2 = 2 needs
    sqrt(476) = 2 sqrt7 sqrt17, both outside Q(sqrt3,sqrt5,sqrt7,sqrt11),
    while d^2 = 4 needs sqrt(880) = 4 sqrt5 sqrt11 and lands back inside it by
    a route de Grey did not take.

Everything below is recomputed from the stored coordinates.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as Fr

import pytest

from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

HERE = os.path.dirname(__file__)
D2 = Fr(64, 9)                       # the forced distance every case starts from

# file, p, q (target d^2 = p/q), radical, new generator, spindled?
CASES = [
    ("five_tuned_1_1.json", 1, 1, 247, 247, False),
    ("five_tuned_1_3.json", 1, 3, 759, 23, True),
    ("five_tuned_2_1.json", 2, 1, 476, 17, True),
    ("five_tuned_4_1.json", 4, 1, 880, 11, True),   # back inside de Grey's field
    # d^2 = 16: radical 16 sqrt7, spindled by de Grey's 4e rotation -- the whole graph lives in
    # Q(sqrt3, sqrt7, sqrt11), i.e. in the CM field Q(sqrt-3, sqrt-7, sqrt-11) (no sqrt5, no sqrt247)
    ("five_tuned_16_1_3_7_11.json", 16, 1, 1792, 7, True),
]


def _load(name):
    with open(os.path.join(HERE, os.pardir, "data", name)) as fh:
        d = json.load(fh)
    field = Field(tuple(d["field_generators"]))
    pts = [Point(field.element([Fr(a, b) for a, b in x]),
                 field.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, field, pts


@pytest.mark.parametrize("name,p,q,rad,gen,spindled", CASES)
def test_the_radical_is_the_formula(name, p, q, rad, gen, spindled):
    assert rad == p * (256 * q - 9 * p)
    d = json.load(open(os.path.join(HERE, os.pardir, "data", name)))
    assert d["radical"] == rad
    assert d["target_d2"] == [p, q]
    assert gen in d["field_generators"]


@pytest.mark.parametrize("name,p,q,rad,gen,spindled", CASES)
def test_the_rotation_is_exact_and_lands_on_the_target(name, p, q, rad, gen,
                                                       spindled):
    """cos phi is rational, sin phi carries the radical, and the composite
    distance is EXACTLY the target -- checked on a vector of the right length
    rather than on the graph, since the identity is about phi alone."""
    field = Field(tuple(json.load(
        open(os.path.join(HERE, os.pardir, "data", name)))["field_generators"]))
    target = Fr(p, q)
    cos = field.rational(Fr(9 * p - 128 * q, 128 * q))
    sin = field.sqrt(rad) * field.rational(Fr(3, 128 * q))
    assert cos * cos + sin * sin == field.rational(Fr(1))
    assert cos == field.rational(target / (2 * D2) - 1)
    # u is any vector of length D; take it along the axis, |u|^2 = 64/9
    u = Point(field.rational(Fr(8, 3)), field.zero())
    assert u.norm2() == field.rational(D2)
    r = Rotation(cos, sin)(u)
    composite = Point(u.x + r.x, u.y + r.y)
    assert composite.norm2() == field.rational(target)


@pytest.mark.parametrize("name,p,q,rad,gen,spindled", CASES)
def test_the_stored_graph_matches_its_own_counts(name, p, q, rad, gen, spindled):
    d, field, pts = _load(name)
    assert len(pts) == len(set(pts)), "duplicate points in the file"
    g = build_graph(pts)
    assert g.n == d["n"]
    assert sum(len(a) for a in g.adj) // 2 == d["m"]


@pytest.mark.parametrize("name,p,q,rad,gen,spindled",
                         [c for c in CASES if c[4] in (23, 17)])
def test_does_not_embed_in_de_greys_field(name, p, q, rad, gen, spindled):
    """Squared distances are invariant under every isometry of the plane, so
    the field they generate is an invariant of the GRAPH, not of the drawing.
    A squared distance with a nonzero coefficient on a basis element carrying
    23 or 17 cannot occur in Q(sqrt3,sqrt5,sqrt7,sqrt11)."""
    d, field, pts = _load(name)
    carries = [i for i in range(field.dim) if field._prod[i] % gen == 0]
    assert carries
    # The file lists the carrier first, then its tau-image, then the spindled
    # copy.  Distances INSIDE any one of those blocks are the carrier's own,
    # since tau and the spindle are isometries -- the new radical only shows up
    # ACROSS the blocks.  So sample across the whole list, not from the front.
    step = max(1, len(pts) // 120)
    probe = list(range(0, len(pts), step))
    seen = False
    for a, i in enumerate(probe):
        for j in probe[a + 1:]:
            e = (pts[i] - pts[j]).norm2()
            if any(e.c[k] != 0 for k in carries):
                seen = True
                break
        if seen:
            break
    assert seen, f"no squared distance carries sqrt({gen})"


def test_tuning_to_one_needs_no_spindle():
    """d^2 = 1 puts the composite pair at distance exactly 1, so it is an EDGE
    and the chained union refuses four on its own -- half the vertices."""
    d, field, pts = _load("five_tuned_1_1.json")
    assert d["mechanism"] == "composed forcing tuned to 1"
    assert d["n"] == 2041
    spindled = json.load(open(os.path.join(HERE, os.pardir, "data",
                                           "five_tuned_1_3.json")))
    assert spindled["n"] > d["n"]


def test_the_composite_distance_sweeps_the_whole_interval():
    """2 D^2 (1 + cos phi) runs from 0 to 4 D^2 as cos phi runs over [-1, 1],
    so every rational target below 2D is reachable -- and the ones that are
    not integers are reachable too, which is what makes 1/3 available."""
    field = Field((3, 5, 7, 11))
    for p, q in ((1, 1), (2, 1), (4, 1), (1, 3), (2, 3), (16, 1)):
        target = Fr(p, q)
        assert target < 4 * D2
        cos = Fr(9 * p - 128 * q, 128 * q)
        assert -1 <= cos <= 1
        assert 2 * D2 * (1 + cos) == target
        rad = p * (256 * q - 9 * p)
        assert Fr(9 * rad, (128 * q) ** 2) == 1 - cos * cos
    # and 4 D^2 itself is the far end: cos phi = 1, tau is a translation
    assert Fr(9 * 4 * D2.numerator, 128 * D2.denominator) != 0
    far = 4 * D2
    assert 2 * D2 * (1 + 1) == far


def test_the_d2_16_graph_lives_in_the_cm_field_q_sqrt_minus_3_minus_7_minus_11():
    """Every point x + iy has x in span(1, sqrt21, sqrt33, sqrt77) and y in
    span(sqrt3, sqrt7, sqrt11, sqrt231): the complex coordinate lies in Q(sqrt-3, sqrt-7, sqrt-11),
    whose non-split places lie above 17, 41, 83, ... (scripts/fieldscreen.py)."""
    d, field, pts = _load("five_tuned_16_1_3_7_11.json")
    assert tuple(field.gens) == (3, 7, 11)
    re_, im_ = {1, 21, 33, 77}, {3, 7, 11, 231}
    for p in pts:
        assert all(c == 0 or field._prod[i] in re_ for i, c in enumerate(p.x.c))
        assert all(c == 0 or field._prod[i] in im_ for i, c in enumerate(p.y.c))
