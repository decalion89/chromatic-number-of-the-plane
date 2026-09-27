"""The vectorised integer path must agree with exact Fraction arithmetic.

Everything at scale runs through `hn.fast`; it is only trustworthy insofar as
it reproduces the slow, obviously-correct path exactly.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from hn.certify import exact_edges
from hn.fast import IntBasis, fast_edges, fast_edges_complete, fast_graph_lazy, fast_walk
from hn.generate import unit_vectors, walk_ball
from hn.graph import build_graph


def test_int_rows_round_trip_through_exact_points():
    U = unit_vectors(m_max=2)
    b = IntBasis.covering(U)
    rows = b.rows(U)
    assert b.points(rows) == U


def test_every_generated_unit_vector_verifies_exactly():
    for m in (1, 2, 3):
        U = unit_vectors(m_max=m)
        b = IntBasis.covering(U)
        assert all(b.is_unit_vector(r) for r in b.rows(U))
        assert all(p.norm2() == 1 for p in U)


def test_fast_walk_matches_the_exact_walk():
    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    fast = set(b.points(fast_walk(b, b.rows(U), steps=3, radius=2.0)))
    slow = set(p for p in walk_ball(3, m_max=1, radius=2.0))
    assert fast == slow


def test_complete_edges_match_brute_force():
    """The grid-based finder must return exactly the O(n^2) answer."""
    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    rows = fast_walk(b, b.rows(U), steps=2, radius=1.6)
    pts = b.points(rows)
    assert sorted(fast_edges_complete(b, rows)) == sorted(exact_edges(pts))


def test_partial_edges_are_a_subset_of_the_truth():
    """Lookup-by-unit-vector may miss edges; it must never invent one."""
    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    rows = fast_walk(b, b.rows(U), steps=3, radius=2.0)
    partial = set(fast_edges(b, rows, b.rows(U)))
    truth = set(exact_edges(b.points(rows)))
    assert partial <= truth


def test_lazy_vertices_survive_the_k_core():
    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    g = fast_graph_lazy(b, fast_walk(b, b.rows(U), steps=4, radius=2.0))
    core = g.k_core(4)
    assert hasattr(core.vertices, "rows")
    assert len(core.vertices.rows) == core.n
    assert core.vertices[0].field == b.field


def test_overflow_guard_is_wired_up():
    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    rows = fast_walk(b, b.rows(U), steps=3, radius=2.0)
    assert 0.0 < b.overflow_headroom(rows) < 1e-6
