"""The 24-point graph with circular chromatic number exactly four.

Read from data/tight_four.json and re-derived here, so a reader who doubts
the claim can run this rather than take it.  The file stores only exact
rationals in the field basis; the edges, the chromatic number and the
refusal of every circular clique below 4 are all recomputed.

Why the object is interesting: the Moser spindle is 4-chromatic with
chi_c = 7/2, half a colour loose, while this one has no slack at all.  The
measurements in hn.homcol say the palette caps that drive de Grey's
construction appear exactly where a carrier is tight in that sense, so a
small tight graph is the template for what is missing one level up.
"""
import json
import os
from fractions import Fraction

import pytest

from hn.fast import IntBasis, fast_edges_complete
from hn.field import Field
from hn.geometry import Point

DATA = os.path.join(os.path.dirname(__file__), "..", "data",
                    "tight_four.json")


def _load():
    with open(DATA) as fh:
        doc = json.load(fh)
    K = Field(tuple(doc["field"]["generators"]))
    pts = []
    for row in doc["points"]:
        half = len(row) // 2
        x = K.element([Fraction(s) for s in row[:half]])
        y = K.element([Fraction(s) for s in row[half:]])
        pts.append(Point(x, y))
    return doc, pts


def _edges(pts):
    basis = IntBasis.covering(pts)
    rows = basis.rows(pts)
    assert basis.overflow_headroom(rows) < 1.0
    return sorted(set((min(a, b), max(a, b))
                      for a, b in fast_edges_complete(basis, rows)))


def _maps_to_circular_clique(n, edges, p, q):
    """Is there a homomorphism to K(p/q)?  The centre is pinned to kill the
    p rotations, which is sound because the clique is vertex-transitive."""
    from pysat.solvers import Solver
    cnf = [[1 + v * p + j for j in range(p)] for v in range(n)]
    for a, b in edges:
        for j in range(p):
            for d in range(-(q - 1), q):
                cnf.append([-(1 + a * p + j), -(1 + b * p + (j + d) % p)])
    cnf += [[-(1 + j)] for j in range(1, p)] + [[1]]
    solver = Solver(name="cd15", bootstrap_with=cnf)
    ok = solver.solve()
    solver.delete()
    return ok


def _colourable(n, edges, k):
    from pysat.solvers import Solver
    cnf = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in edges:
        for c in range(k):
            cnf.append([-(1 + a * k + c), -(1 + b * k + c)])
    solver = Solver(name="cd15", bootstrap_with=cnf)
    ok = solver.solve()
    solver.delete()
    return ok


def test_file_describes_the_graph_it_claims():
    doc, pts = _load()
    assert len(pts) == doc["vertices"] == 24
    edges = _edges(pts)
    assert len(edges) == doc["edges"] == 54
    assert edges == [tuple(e) for e in doc["unit_edges"]]


def test_chromatic_number_is_four():
    doc, pts = _load()
    edges = _edges(pts)
    assert not _colourable(len(pts), edges, 3)
    assert _colourable(len(pts), edges, 4)


@pytest.mark.parametrize("p,q", [(35, 9), (31, 8), (27, 7), (23, 6),
                                 (19, 5), (15, 4), (11, 3), (7, 2)])
def test_refuses_the_ratios_just_below_four(p, q):
    """chi_c = 4 needs every ratio below 4 refused; these are the largest
    ones with small denominators, and the widest net is in the scripts."""
    doc, pts = _load()
    edges = _edges(pts)
    assert Fraction(p, q) < 4
    assert not _maps_to_circular_clique(len(pts), edges, p, q)


def test_it_does_map_at_four_so_the_value_is_exactly_four():
    doc, pts = _load()
    edges = _edges(pts)
    assert _maps_to_circular_clique(len(pts), edges, 4, 1)
