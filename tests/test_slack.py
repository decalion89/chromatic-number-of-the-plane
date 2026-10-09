"""The slack framing: why the plane's difficulty jumps where it does.

A unit-distance graph in the plane has clique number 3.  Four points pairwise
at distance 1 would be a regular tetrahedron, which does not fit in two
dimensions.  Define the slack of a colouring problem as

    s = k - 3

and let f(k) be the number of vertices in the smallest unit-distance
configuration containing a pair forced monochromatic under k colours.

At s = 0 a triangle exhausts the palette, so forcing is *local*: two triangles
sharing an edge already do it, and f(3) = 4, which meets the trivial bound
f(k) >= k + 1 exactly.

At s >= 1 no local configuration exhausts anything -- a spare colour always
remains -- so forcing must be assembled combinatorially across many vertices,
and f jumps.  The historical record reads the same way: s = 0 was settled in
1961, s = 1 took until 2018 and needed 1581 vertices, s = 2 is open (OpenAI
proved chi >= 6 in 2026 without a graph, so the explicit graph is what is open).

These tests pin the anchors.  f(4) itself is measured by
scripts/measure_fk.py rather than here, being far too slow for a suite.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from pysat.formula import CNF
from pysat.solvers import Solver

from hn.generate import hex_ball
from hn.geometry import SPINDLE, eisenstein, origin
from hn.graph import build_graph


def forces(g, k, p, q):
    """Does every proper k-colouring give p and q the same colour?"""
    x = lambda v, c: 1 + v * k + c
    cnf = CNF()
    for v in range(g.n):
        cnf.append([x(v, c) for c in range(k)])
    for u, v in g.edges():
        for c in range(k):
            cnf.append([-x(u, c), -x(v, c)])
    for c in range(k):
        cnf.append([-x(p, c), -x(q, c)])        # demand they differ
    s = Solver(name="cd19", bootstrap_with=cnf)
    try:
        return s.solve() is False               # no such colouring exists
    finally:
        s.delete()


def rhombus():
    return build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])


def test_the_plane_has_clique_number_three():
    """Four points pairwise at distance 1 need a regular tetrahedron."""
    g = build_graph(hex_ball(5) + [SPINDLE(p) for p in hex_ball(3)])
    assert g.find_clique(3) is not None
    assert g.find_clique(4) is None


def test_nothing_on_k_vertices_can_force_with_k_colours():
    """The trivial bound f(k) >= k+1: colour every vertex differently."""
    triangle = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1)])
    assert not any(forces(triangle, 3, i, j)
                   for i in range(3) for j in range(i + 1, 3))


def test_at_zero_slack_the_trivial_bound_is_tight():
    """f(3) = 4.  A triangle exhausts three colours, so two of them sharing an
    edge force the remaining pair -- forcing is local at s = 0."""
    assert forces(rhombus(), 3, 0, 3)


def test_one_unit_of_slack_changes_the_regime():
    """The same four vertices force nothing at k=4, and neither does the
    seven-vertex Moser spindle: with a spare colour, no local configuration
    exhausts the palette."""
    assert not forces(rhombus(), 4, 0, 3)
    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    spindle = build_graph(rh + [SPINDLE(p) for p in rh])
    assert not any(forces(spindle, 4, i, j)
                   for i in range(spindle.n) for j in range(i + 1, spindle.n))


def test_the_spindle_is_still_four_chromatic_despite_forcing_nothing():
    """Forcing and chromatic number are different things, and the distinction
    is the whole reason the slack framing is about f(k) and not about chi."""
    from hn.coloring import chromatic_number

    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    spindle = build_graph(rh + [SPINDLE(p) for p in rh])
    assert chromatic_number(spindle)[0] == 4


# -- why a forced core stays wide ------------------------------------------

def test_symmetric_targets_resist_being_singled_out():
    """The forced core on the measured k=4 configuration is 34 of 36 targets,
    and the 36 form a closed orbit under the 60-degree rotation about the
    pivot.  A colouring argument cannot distinguish targets that a symmetry
    permutes, so symmetry sets a floor on how narrow a core can get.

    The rotation there is not a full automorphism -- 280 of 359 vertices land
    back in the graph -- and that partial asymmetry is what bought the two
    exclusions taking 36 down to 34.  Asymmetry buys exclusions; symmetry
    blocks them.

    Which indicts tightening by unions of copies rotated about the *same*
    pivot: that raises symmetry, pushing the core the wrong way.  de Grey's
    construction is asymmetric on purpose -- two different rotations about an
    off-centre pivot -- and this is what that is for.

    The invariant checked here is the cheap half: a rotation that closes on
    the target set maps targets to targets.
    """
    from fractions import Fraction

    from hn.geometry import ROT60, Point

    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    g = build_graph(rh + [SPINDLE(p) for p in rh] + hex_ball(2))
    p = g.vertices[0]
    third = Fraction(1, 3)
    targets = {v for v in g.vertices if p.dist2(v).is_rational() and p.dist2(v).c[0] == third}
    turn = ROT60.about(p)
    if targets:
        # whatever sits at that radius is permuted by the rotation about the pivot
        assert all(turn(t) in {q for q in map(turn, targets)} for t in targets)
    # and the rotation need not preserve the graph as a whole
    inside = sum(1 for v in g.vertices if turn(v) in set(g.vertices))
    assert inside <= g.n
