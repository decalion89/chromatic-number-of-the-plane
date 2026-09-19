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
                        spindle_union, spindle_union_auto, targets_at_third,
                        triple_spindle_union)


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


# -- the field is a choice, not a constraint -------------------------------

def test_embedding_preserves_arithmetic():
    from hn.field import Field, embed

    big = Field((2, 3, 7, 11, 13))
    x = F.sqrt(33) * 2 + F.rational(5)
    y = F.sqrt(3) - F.rational(1)
    bx, by = embed(x, big), embed(y, big)
    assert (bx * by).c[0] == (x * y).c[0]
    assert embed(x + y, big) == bx + by
    assert bx * bx == embed(x * x, big)


def test_embedding_refuses_a_field_that_lacks_a_radical():
    from hn.field import Field, embed

    with pytest.raises(ValueError):
        embed(F.sqrt(11), Field((2, 3)))


def test_spindle_union_auto_matches_the_fixed_field_version():
    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    fixed = spindle_union(rhombus, 0, 3, F)
    auto, field = spindle_union_auto(rhombus, 0, 3)
    assert (auto.n, auto.m) == (fixed.n, fixed.m) == (7, 11)
    assert is_k_colorable(auto, 3)[0] is False


def test_spindle_union_auto_widens_the_field_when_needed():
    """A pair whose spindle needs a radical the graph's field lacks is not a
    reason to skip it -- forcing is decided by SAT on the graph, and the field
    only has to be wide enough to write the rotated copy down."""
    from fractions import Fraction

    from hn.field import Field
    from hn.geometry import Point, required_radical

    # d^2 = 25/9 needs sqrt(91) = sqrt7 * sqrt13, absent from Q(sqrt3, sqrt11)
    assert required_radical(Fraction(25, 9)) == 91
    base = Field((3, 11))
    p = Point(base.zero(), base.zero())
    q = Point(base.rational(Fraction(5, 3)), base.zero())
    g = build_graph([p, q])
    spun, field = spindle_union_auto(g, 0, 1)
    assert 7 in field.gens and 13 in field.gens
    assert 3 in field.gens and 11 in field.gens        # the originals survive
    assert spun.n == 3                                 # pivot fixed, target imaged


def test_spindle_rejects_targets_closer_than_half():
    from fractions import Fraction

    from hn.field import Field
    from hn.geometry import Point

    base = Field((3,))
    g = build_graph([Point(base.zero(), base.zero()),
                     Point(base.rational(Fraction(1, 4)), base.zero())])
    with pytest.raises(ValueError):
        spindle_union_auto(g, 0, 1)                    # no rotation separates them by 1


# -- a gradient where there was none ---------------------------------------

def test_separation_difficulty_reports_effort_per_query():
    """Asking "is this pair forced?" gives yes or no and no sense of distance.
    The solver's conflict count supplies the missing gradient: zero conflicts
    means wide open, many means nearly forced."""
    from fractions import Fraction

    from hn.spindle import SeparationDifficulty

    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    targets = [3]
    td = SeparationDifficulty(rhombus, 3, 0, targets)
    try:
        rows = td.effort_ranking({Fraction(3): targets})
    finally:
        td.close()
    assert len(rows) == 1
    d2, separable, conflicts, decisions, core = rows[0]
    assert d2 == 3
    assert separable is False          # the rhombus forces its far tips
    assert core == [3]
    assert conflicts >= 0 and decisions >= 0


def test_effort_ranking_puts_the_hardest_first():
    from fractions import Fraction

    from hn.spindle import SeparationDifficulty

    g = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1),
                     eisenstein(2, 0), eisenstein(2, 1), eisenstein(1, 2)])
    p = g.vertices[0]
    groups = {}
    for j in range(1, g.n):
        d2 = p.dist2(g.vertices[j])
        if d2.is_rational() and d2.c[0] >= Fraction(1, 4):
            groups.setdefault(d2.c[0], []).append(j)
    td = SeparationDifficulty(g, 4, 0, sorted({j for js in groups.values() for j in js}))
    try:
        rows = td.effort_ranking(groups)
    finally:
        td.close()
    conflicts = [r[2] for r in rows]
    assert conflicts == sorted(conflicts, reverse=True)
