"""The spindling arguments, checked on cases whose answers are known."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fractions import Fraction

import pytest

from hn.coloring import is_k_colorable
from hn.field import QSQRT3_11 as F
from hn.geometry import ROT60, Point, eisenstein, origin
from hn.graph import build_graph
from hn.spindle import (ForcedPairFinder, SeparationTest, candidate_pairs_from,
                        spindle_union, targets_at_third, triple_spindle_union)


def unit_triangle_about_origin():
    """The unit equilateral triangle inscribed in the circle of radius 1/sqrt(3),
    plus its centre.  The centre is the pivot; the three corners are the targets."""
    r = F.sqrt(3).inverse()
    q = Point(r, F.zero())
    rho = ROT60 ** 2
    return [origin(), q, rho(q), rho(rho(q))]


def test_the_circle_of_radius_one_over_sqrt3_carries_a_unit_triangle():
    """The geometric fact the three-copy argument rests on."""
    pts = unit_triangle_about_origin()
    centre, corners = pts[0], pts[1:]
    for c in corners:
        assert centre.dist2(c) == Fraction(1, 3)
    for i in range(3):
        for j in range(i + 1, 3):
            assert corners[i].dist2(corners[j]) == 1


def test_three_colours_tie_the_centre_to_a_corner():
    """The triangle uses all three colours, so the centre must repeat one."""
    g = build_graph(unit_triangle_about_origin())
    assert g.m == 3                      # the centre is adjacent to nothing
    targets = targets_at_third(g, 0)
    assert len(targets) == 3
    sep, core = SeparationTest(g, 3, 0, targets).run()
    assert sep is False                  # no 3-colouring separates the centre from all
    assert len(core) == 3                # and all three corners are needed


def test_four_colours_leave_the_centre_free():
    """With a spare colour the same configuration forces nothing -- which is
    why the plane is hard at k=4 and easy at k=3."""
    g = build_graph(unit_triangle_about_origin())
    sep, core = SeparationTest(g, 4, 0, targets_at_third(g, 0)).run()
    assert sep is True
    assert core == []


def test_core_of_size_three_is_rejected_as_unusable():
    """Only cores of size 1 (two copies) or 2 (three copies, at d^2 = 1/3) can
    be spindled; four points pairwise at distance 1 do not fit on a circle."""
    g = build_graph(unit_triangle_about_origin())
    _, core = SeparationTest(g, 3, 0, targets_at_third(g, 0)).run()
    assert len(core) == 3
    assert len(core) > 2                 # the search must skip this, not spindle it


def test_triple_spindle_images_are_pairwise_adjacent():
    g = build_graph(unit_triangle_about_origin())
    spun = triple_spindle_union(g, 0)
    target = g.vertices[1]
    rho = ROT60 ** 2
    images = [target, rho(target), rho(rho(target))]
    for i in range(3):
        for j in range(i + 1, 3):
            assert images[i].dist2(images[j]) == 1
            assert images[i] in spun.vertices and images[j] in spun.vertices


def test_separation_test_agrees_with_pairwise_forced_search():
    """A core of size 1 is exactly a forced pair, by either route."""
    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    forced = ForcedPairFinder(rhombus, 3, candidate_pairs_from(rhombus, 0, F)).run(verbose=False)
    assert len(forced) == 1
    _, q = forced[0]
    sep, core = SeparationTest(rhombus, 3, 0, [q]).run()
    assert sep is False and core == [q]


def test_two_copy_spindle_still_works_after_the_refactor():
    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    forced = ForcedPairFinder(rhombus, 3, candidate_pairs_from(rhombus, 0, F)).run(verbose=False)
    spun = spindle_union(rhombus, forced[0][0], forced[0][1], F)
    assert is_k_colorable(spun, 3)[0] is False
