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
    assert {11, 23, 35, 47, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 455, 599, 611, 791, 911, 935,
            959} <= set(FIELDS)
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
    assert vq.upper_bound(95, planes)[0] == 4                                                # 95 = 2^2 mod 7
    assert vq.upper_bound(131, planes)[0] == 4                                               # 131 = 3 mod 8: Fischer
    assert vq.upper_bound(251, planes)[0] == 4                                               # 251 = 3 mod 8: Fischer
    assert vq.upper_bound(455, planes)[0] == 4                                               # 455 = 5 7 13 = 0 mod 7
    assert vq.upper_bound(935, planes)[0] == 4                                               # 935 = 2^2 mod 7
    assert vq.upper_bound(263, planes)[0] == 4 and vq.upper_bound(599, planes)[0] == 4       # both 2^2 mod 7
    assert vq.upper_bound(959, planes)[0] == 4                                               # 959 = 7 137
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


def test_cake_lpr_log_covers_every_field():
    """data/quadratic_planes/cake_lpr_checks.txt: the checker run with kissat, drat-trim and cake_lpr has a
    VERIFIED UNSAT line for every field and ends with exit status 0"""
    log = open(os.path.join(D, "cake_lpr_checks.txt"), encoding="utf-8").read()
    for d in FIELDS:
        assert f"Q(sqrt{d}): ok: kissat UNSATISFIABLE, drat-trim VERIFIED and cake_lpr VERIFIED UNSAT" in log, d
    assert "\n# exit status 0" in log and "\nCONFIRMED: " in log


sys.path.insert(0, os.path.join(ROOT, "lean", "tools"))
import field_lean  # noqa: E402

LEAN = os.path.join(ROOT, "lean")


def test_lean_files_match_the_data():
    """lean/Sqrt{d}.lean (the Lean proof of chi(Q(sqrt d)^2) = 4) is the file lean/tools/field_lean.py writes from
    data/quadratic_planes/q{d}.json, for every field with a Lean proof, so each formal proof is about the published
    graph; and the 4-colouring of F_7^2 in lean/QuadraticPlanes.lean is that of finite_planes.json."""
    import subprocess
    r = subprocess.run([sys.executable, os.path.join(LEAN, "tools", "field_lean.py"), "--check"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_lean_fields_are_built_and_checked():
    """Every field of field_lean.FIELDS has its file, its library in lakefile.toml (built by default), its line in
    PrintAxioms.lean and axioms.expected, and its kernel replay in the Lean workflow; and no other Sqrt{d}.lean
    exists."""
    lake = open(os.path.join(LEAN, "lakefile.toml"), encoding="utf-8").read()
    default = lake.split("defaultTargets = [", 1)[1].split("]", 1)[0]
    pa = open(os.path.join(LEAN, "PrintAxioms.lean"), encoding="utf-8").read()
    ax = open(os.path.join(LEAN, "axioms.expected"), encoding="utf-8").read()
    wf = open(os.path.join(ROOT, ".github", "workflows", "lean.yml"), encoding="utf-8").read()
    replayed = wf.split("for m in ", 1)[1].split(";", 1)[0].split()
    for d in field_lean.FIELDS:
        m = f"Sqrt{d}"
        assert os.path.exists(os.path.join(LEAN, m + ".lean")), m
        assert f'name = "{m}"' in lake and f'"{m}"' in default, m
        assert f"import {m}\n" in pa and f"#print axioms {m}.chromaticNumber_eq_four\n" in pa, m
        assert (f"'{m}.chromaticNumber_eq_four' depends on axioms: [propext, Classical.choice, Quot.sound]"
                in ax.splitlines()), m
        assert m in replayed, m
    assert "QuadraticPlanes" in replayed and '"QuadraticPlanes"' in default
    on_disk = {int(f[4:-5]) for f in os.listdir(LEAN) if f.startswith("Sqrt") and f.endswith(".lean")}
    assert on_disk == set(field_lean.FIELDS)


@pytest.mark.parametrize("d", field_lean.FIELDS)
def test_lrat_proof_is_for_the_stored_formula(d):
    """data/quadratic_planes/q{d}.lrat, which Sqrt{d}.lean checks, refers only to the clauses of q{d}.cnf and earlier
    lemmas, and ends with the empty clause."""
    with open(os.path.join(D, f"q{d}.cnf")) as fh:
        header = fh.readline().split()
    nclauses = int(header[3])
    with open(os.path.join(D, f"q{d}.lrat")) as fh:
        lines = [l.split() for l in fh if l.strip()]
    added = set()
    empty = False
    for t in lines:
        if t[1] == "d":
            continue
        k = int(t[0])
        assert k > nclauses and k not in added
        z = t.index("0", 1)
        hints = [int(x) for x in t[z + 1:-1]]
        assert all(0 < abs(h) <= nclauses or abs(h) in added for h in hints)
        added.add(k)
        empty = empty or z == 1
    assert empty
