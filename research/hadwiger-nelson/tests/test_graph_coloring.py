"""Graph construction, reductions, and colouring, checked against known values."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from hn.coloring import ColoringInstance, chromatic_number, find_uncolorable_core, is_k_colorable
from hn.field import QSQRT3_11 as F
from hn.generate import hex_ball, unit_vectors
from hn.geometry import SPINDLE, eisenstein, origin
from hn.graph import build_graph
from hn.spindle import ForcedPairFinder, candidate_pairs_from, spindle_union


def moser_spindle():
    rhombus = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    return build_graph(rhombus + [SPINDLE(p) for p in rhombus])


def test_moser_spindle_shape():
    g = moser_spindle()
    assert (g.n, g.m) == (7, 11)
    assert sorted(g.degrees()) == [3, 3, 3, 3, 3, 3, 4]


def test_moser_spindle_is_four_chromatic():
    g = moser_spindle()
    assert is_k_colorable(g, 3)[0] is False
    assert is_k_colorable(g, 4)[0] is True
    assert chromatic_number(g)[0] == 4


def test_triangular_lattice_is_three_chromatic():
    g = build_graph(hex_ball(4))
    assert chromatic_number(g)[0] == 3


def test_reported_colouring_is_actually_proper():
    g = moser_spindle()
    ok, coloring = is_k_colorable(g, 4)
    assert ok
    for u, v in g.edges():
        assert coloring[u] != coloring[v]


def test_k_core_preserves_colourability():
    """A vertex of degree < k can always be coloured last, so the k-core
    decides k-colourability for the whole graph."""
    g = build_graph(hex_ball(3) + [SPINDLE(p) for p in hex_ball(2)])
    for k in (3, 4):
        assert is_k_colorable(g, k)[0] == is_k_colorable(g.k_core(k), k)[0]


def test_planar_unit_distance_graphs_have_no_k4():
    """Four mutually unit-distant points do not exist in the plane, so colour
    symmetry breaking can only ever pin a triangle."""
    g = build_graph(hex_ball(3))
    assert g.find_clique(4) is None
    assert g.find_clique(3) is not None


def test_uncolorable_core_is_found_and_minimal_ish():
    g = moser_spindle()
    core = find_uncolorable_core(g, 3, verbose=False)
    assert core is not None
    assert is_k_colorable(core, 3)[0] is False
    assert core.n <= g.n


def test_forced_pair_discovers_the_rhombus_lemma():
    """The classic step, found rather than assumed: a unit rhombus forces its
    two far tips to share a colour whenever only 3 colours are available."""
    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    pairs = candidate_pairs_from(rhombus, 0, F)
    forced = ForcedPairFinder(rhombus, 3, pairs).run(verbose=False)
    assert forced, "no forced pair found in the unit rhombus"
    p, q = forced[0]
    assert rhombus.vertices[p].dist2(rhombus.vertices[q]) == 3


def test_spindling_a_forced_pair_removes_all_colourings():
    """The whole de Grey move, end to end, at k=3: rhombus -> Moser spindle."""
    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    assert is_k_colorable(rhombus, 3)[0] is True          # G alone is fine
    forced = ForcedPairFinder(rhombus, 3, candidate_pairs_from(rhombus, 0, F)).run(verbose=False)
    spun = spindle_union(rhombus, forced[0][0], forced[0][1], F)
    assert (spun.n, spun.m) == (7, 11)
    assert is_k_colorable(spun, 3)[0] is False            # the union is not
