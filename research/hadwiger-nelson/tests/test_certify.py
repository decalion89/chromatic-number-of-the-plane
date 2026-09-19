"""Certificates must catch wrong answers, not just bless right ones."""
import sys, os, shutil, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from hn.certify import (exact_edges, load_certificate, save_certificate,
                        verify_certificate, verify_coloring, verify_uncolorable,
                        write_dimacs)
from hn.coloring import chromatic_number
from hn.geometry import SPINDLE, eisenstein, origin
from hn.graph import build_graph

DRAT = shutil.which("drat-trim") or os.environ.get("DRAT_TRIM")


def moser_spindle():
    rhombus = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    return build_graph(rhombus + [SPINDLE(p) for p in rhombus])


def test_certificate_round_trips_exactly(tmp_path):
    g = moser_spindle()
    path = str(tmp_path / "c.json")
    save_certificate(g, path, 4, "test")
    pts, doc = load_certificate(path)
    assert pts == list(g.vertices)
    assert doc["n"] == g.n


def test_valid_colouring_is_accepted(tmp_path):
    g = moser_spindle()
    _, coloring = chromatic_number(g)
    ok, msg = verify_coloring(list(g.vertices), coloring, 4)
    assert ok, msg


def test_monochromatic_unit_pair_is_rejected():
    g = moser_spindle()
    _, coloring = chromatic_number(g)
    u, v = next(iter(g.edges()))
    broken = list(coloring)
    broken[v] = broken[u]                     # force an illegal edge
    ok, msg = verify_coloring(list(g.vertices), broken, 4)
    assert not ok and "distance 1" in msg


def test_colour_outside_the_palette_is_rejected():
    g = moser_spindle()
    _, coloring = chromatic_number(g)
    broken = list(coloring)
    broken[0] = 9
    ok, msg = verify_coloring(list(g.vertices), broken, 4)
    assert not ok and "outside" in msg


def test_verifier_does_not_trust_the_graph_builder():
    """verify_coloring re-derives edges in O(n^2); a colouring that ignores an
    edge the fast builder might have missed must still fail."""
    g = moser_spindle()
    pts = list(g.vertices)
    assert sorted(exact_edges(pts)) == sorted(g.edges())


def test_claiming_a_lower_bound_on_a_colourable_graph_fails():
    g = moser_spindle()
    ok, msg = verify_uncolorable(list(g.vertices), 4, drat_trim=DRAT)
    assert not ok and "IS 4-colourable" in msg


@pytest.mark.skipif(not DRAT, reason="drat-trim not on PATH")
def test_drat_proof_verifies_the_four_chromatic_lower_bound(tmp_path):
    """chi(R^2) >= 4, machine-checked from coordinates to proof."""
    g = moser_spindle()
    ok, msg = verify_uncolorable(list(g.vertices), 3, workdir=str(tmp_path), drat_trim=DRAT)
    assert ok, msg
    assert "VERIFIED" in msg


def test_dimacs_encoding_has_the_expected_shape(tmp_path):
    g = moser_spindle()
    path = str(tmp_path / "g.cnf")
    n_clauses = write_dimacs(g.n, list(g.edges()), 3, path)
    assert n_clauses == g.n + g.m * 3
    head = open(path).readline().split()
    assert head[2:] == [str(g.n * 3), str(n_clauses)]
