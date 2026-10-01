"""Planes over real quadratic fields that need four colours (notes/quadratic_planes.md): exact checks, no solver."""
import copy
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import verify_quadratic_planes as vq  # noqa: E402

D = os.path.join(ROOT, "data", "quadratic_planes")
FIELDS = sorted(int(f[1:-5]) for f in os.listdir(D) if f.startswith("q") and f.endswith(".json"))


def load(name):
    with open(os.path.join(D, name)) as fh:
        return json.load(fh)


def test_fields():
    assert {11, 23, 35, 47, 59, 71, 119, 131, 155, 179, 191, 239, 359, 431} <= set(FIELDS)
    for d in FIELDS:
        assert d % 4 == 3 and d % 3 == 2            # d = 11 mod 12: the only real quadratic fields that can need 4


@pytest.mark.parametrize("d", FIELDS)
def test_graph_is_a_unit_distance_graph(d):
    g = load(f"q{d}.json")
    assert g["d"] == d
    ok, msg = vq.check_graph(g)
    assert ok, msg


@pytest.mark.parametrize("d", FIELDS)
def test_colourings(d):
    g = load(f"q{d}.json")
    assert vq.check_colouring(g, g["four_colouring"], 4)
    ok, msg = vq.check_critical(g)
    assert ok, msg


@pytest.mark.parametrize("d", FIELDS)
def test_stored_cnf_and_logs(d):
    g = load(f"q{d}.json")
    ok, msg = vq.check_cnf(g, os.path.join(D, f"q{d}.cnf"))
    assert ok, msg
    for enc in (0, 1):
        with open(os.path.join(D, f"q{d}.logs", f"encoding{enc}.kissat.log")) as fh:
            assert "s UNSATISFIABLE" in fh.read()
        with open(os.path.join(D, f"q{d}.logs", f"encoding{enc}.drat-trim.log")) as fh:
            assert "s VERIFIED" in fh.read()


def test_upper_bounds():
    planes = load("finite_planes.json")["planes"]
    assert vq.plane_colouring_ok(7, planes["7"]["colouring"]) and planes["7"]["colours"] == 4
    assert vq.plane_colouring_ok(11, planes["11"]["colouring"]) and planes["11"]["colours"] == 5
    assert vq.upper_bound(11, planes)[0] == 4 and vq.upper_bound(23, planes)[0] == 4
    assert vq.upper_bound(59, planes)[0] == 4 and vq.upper_bound(71, planes)[0] == 4
    assert vq.upper_bound(119, planes)[0] == 4 and vq.upper_bound(191, planes)[0] == 4    # 119 = 0, 191 = 3^2 mod 7
    assert vq.upper_bound(239, planes)[0] == 4 and vq.upper_bound(35, planes)[0] == 4
    assert vq.upper_bound(359, planes)[0] == 4 and vq.upper_bound(431, planes)[0] == 4
    assert vq.upper_bound(179, planes)[0] == 4 and vq.upper_bound(155, planes)[0] == 4     # 155 = 1 mod 7
    assert vq.upper_bound(131, planes)[0] == 4                                               # 131 = 3 mod 8: Fischer
    assert vq.upper_bound(47, planes)[0] == 5


def test_unit_length_is_exact():
    # (sqrt11/6, 5/6) has length 1 in Q(sqrt11)^2: 11/36 + 25/36 = 1; scaled by D = 6
    assert vq.unit(11, 6, (0, 0, 0, 0), (0, 1, 5, 0))
    assert not vq.unit(11, 6, (0, 0, 0, 0), (0, 1, 4, 0))
    assert not vq.unit(11, 6, (0, 0, 0, 0), (1, 1, 5, 0))     # cross term 2 a b sqrt11 does not vanish


def test_tampering_is_caught():
    g = copy.deepcopy(load("q11.json"))
    g["points"][3][0] += 1
    assert not vq.check_graph(g)[0]
    g = load("q11.json")
    bad = list(g["four_colouring"])
    a, b = g["edges"][0]
    bad[b] = bad[a]
    assert not vq.check_colouring(g, "".join(bad), 4)
    g = copy.deepcopy(load("q11.json"))
    v = "5"
    col = list(g["critical_3_colourings"][v])
    a, b = next(e for e in g["edges"] if 5 not in e)
    col[b] = col[a]
    g["critical_3_colourings"][v] = "".join(col)
    assert not vq.check_critical(g)[0]
