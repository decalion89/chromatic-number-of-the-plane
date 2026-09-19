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


def test_effort_measures_constraint_not_size():
    """The control the gradient needs to be worth anything.

    A bigger formula can cost more conflicts for no reason but its size, which
    would make the hill climb an artefact.  Adding a far-away translate doubles
    the vertex count while adding no constraint between the copies, so whatever
    it buys is the size effect alone.  Measured on de Grey's graph at k=5:

        G alone                        1581 vertices   score   88
        G + a translate 100 units away 3162 vertices   score  129   (x1.5)
        G u rho(G), genuinely tighter  3008 vertices   score 4216   (x48)

    Size buys half again; tightening buys forty-eight times.  The score tracks
    constraint.

    Checked here structurally on a small graph, since at four vertices the
    conflict counts are too small for a ratio to mean anything: a far copy must
    leave the original's forcing exactly as it was, and must itself force
    nothing.
    """
    from fractions import Fraction

    from hn.geometry import Point
    from hn.spindle import SeparationDifficulty

    base = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    shift = Point(F.rational(100), F.zero())
    translated = build_graph(list(base.vertices) + [v + shift for v in base.vertices])

    assert translated.n == 2 * base.n
    assert translated.m == 2 * base.m          # no edge crosses between copies

    def rows(g, k):
        p = g.vertices[0]
        groups = {}
        for j in range(1, g.n):
            d2 = p.dist2(g.vertices[j])
            if d2.is_rational() and d2.c[0] >= Fraction(1, 4):
                groups.setdefault(d2.c[0], []).append(j)
        td = SeparationDifficulty(g, k, 0, sorted({j for js in groups.values() for j in js}))
        try:
            return {r[0]: r[1] for r in td.effort_ranking(groups)}
        finally:
            td.close()

    before, after = rows(base, 3), rows(translated, 3)
    # the rhombus still forces its far tips, and the near distances are unchanged
    assert before[Fraction(3)] is False
    assert after[Fraction(3)] is False
    for d2, separable in before.items():
        assert after[d2] == separable
    # and every distance reaching the far copy forces nothing
    for d2, separable in after.items():
        if d2 > 100:
            assert separable is True


# -- keeping the search local, so the machine lasts -------------------------

def test_local_ball_keeps_the_pivot_and_its_neighbourhood():
    from hn.generate import hex_ball
    from hn.spindle import local_ball

    g = build_graph(hex_ball(4))
    centre = next(i for i, v in enumerate(g.vertices)
                  if v.x.is_zero() and v.y.is_zero())
    ball, bp = local_ball(g, centre, 1.5)
    assert ball.n < g.n
    assert ball.vertices[bp] == g.vertices[centre]
    # every kept vertex really is within the radius
    p = ball.vertices[bp]
    for v in ball.vertices:
        assert (v.fx - p.fx) ** 2 + (v.fy - p.fy) ** 2 <= 1.5 ** 2 + 1e-9
    # and every unit neighbour of the centre survived
    assert len(ball.adj[bp]) == len(g.adj[centre])


def test_core_preserves_forcing_not_just_colourability():
    """A vertex of degree below k is always colourable last, so the colourings
    of G restricted to G - v are exactly those of G - v.  Forcing among the
    survivors is therefore unchanged -- which is what licenses peeling the
    graph between rounds of the hill climb."""
    from hn.generate import hex_ball
    from hn.spindle import SeparationTest, core_preserving_forcing

    rhombus = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    g = build_graph(rhombus + hex_ball(3))
    tip = next(j for j in range(g.n) if g.vertices[j] == eisenstein(1, 1))
    before = SeparationTest(g, 3, 0, [tip]).run()
    res = core_preserving_forcing(g, 3, 0, [tip])
    assert res is not None
    core, pivot, targets = res
    after = SeparationTest(core, 3, pivot, targets).run()
    assert before[0] == after[0] is False       # forced before, forced after
    assert core.n <= g.n


def test_core_protects_the_pivot_and_targets_from_peeling():
    from hn.spindle import core_preserving_forcing

    g = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    res = core_preserving_forcing(g, 5, 0, [3])   # k above every degree here
    assert res is not None
    core, pivot, targets = res
    assert core.vertices[pivot] == g.vertices[0]
    assert targets and core.vertices[targets[0]] == g.vertices[3]


def test_effort_ranking_honours_a_conflict_budget():
    """A query approaching forcing costs unboundedly much, so a round without a
    cap can simply never return -- and a round that never returns reports
    nothing at all.  Exhausting the budget comes back as separable None, which
    is the strongest signal the measure gives short of UNSAT."""
    from fractions import Fraction

    from hn.spindle import SeparationDifficulty

    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    td = SeparationDifficulty(rhombus, 3, 0, [3])
    try:
        rows = td.effort_ranking({Fraction(3): [3]}, conflict_budget=100000)
    finally:
        td.close()
    # this one is genuinely forced and cheap, so the budget does not bite
    assert rows[0][1] is False
    assert rows[0][4] == [3]


def test_triple_spindle_works_over_a_larger_field():
    """It is called exactly when a narrowed core reaches size 2, on graphs
    living in Q(sqrt3, sqrt5, sqrt7, sqrt11) -- so building its 120-degree
    rotation from the module-level ROT60, which is over Q(sqrt3, sqrt11),
    would fail at the one moment it matters."""
    from fractions import Fraction

    from hn.geometry import DEGREY_FIELD, Point, _rot60
    from hn.spindle import triple_spindle_union

    F4 = DEGREY_FIELD
    r = F4.sqrt(3).inverse()
    q = Point(r, F4.zero())
    rho = _rot60(F4) ** 2
    g = build_graph([Point(F4.zero(), F4.zero()), q, rho(q), rho(rho(q))])
    spun = triple_spindle_union(g, 0)
    assert spun.n >= g.n
    assert spun.vertices[0].field == F4
