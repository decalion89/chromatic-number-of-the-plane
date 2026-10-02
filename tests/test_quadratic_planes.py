"""Planes over real quadratic fields that need four colours (notes/quadratic_planes.md): exact checks, no solver."""
import copy
import hashlib
import json
import os
import re
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
    assert {11, 23, 35, 47, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851,
            911, 935, 959} <= set(FIELDS)
    for d in FIELDS:
        assert d % 4 == 3 and d % 3 == 2            # d = 11 mod 12: the only real quadratic fields that can need 4


@pytest.mark.parametrize("d", FIELDS)
def test_graph_is_a_unit_distance_graph(d):
    g = load(f"q{d}.json")
    assert g["d"] == d
    ok, msg = vq.check_graph(g)
    assert ok, msg
    assert msg.endswith("no triangle, girth 4"), msg       # the paper says every graph has girth 4


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



def _circle_eigenvalues(p, level):
    """the eigenvalues of Cay((Z/p^level)^2, unit circle) for the characters of exact level `level`, as fractions
    of the degree: one character per norm class (rotations act transitively on primitive vectors of a given norm)"""
    import math
    m = p ** level
    T = [(x, y) for x in range(m) for y in range(m) if (x * x + y * y) % m == 1]
    reps = {}
    for x in range(m):
        for y in range(m):
            N = (x * x + y * y) % m
            if N % p and N not in reps:
                reps[N] = (x, y)
    return [sum(math.cos(2 * math.pi * ((a * u + b * v) % m) / m) for u, v in T) / len(T) for a, b in reps.values()]


def test_padic_planes_hoffman():
    """notes §6, Q(sqrt47): Hoffman's ratio -mu/(1 - mu) of the level-1 plane F_p^2 is below 1/4 for p = 23, 31, 43,
    47 (no 4-colouring), not for p = 11, 19; at level 2 the new eigenvalues are at most 2/(p + 1) of the degree,
    so the least eigenvalue, and with it Hoffman's ratio, is that of level 1"""
    ratio = {}
    for p in (11, 19, 23, 31, 43, 47):
        mu = min(_circle_eigenvalues(p, 1))
        ratio[p] = -mu / (1 - mu)
    assert ratio[11] > 0.28 and ratio[19] > 0.27
    assert all(ratio[p] < 0.25 for p in (23, 31, 43, 47))
    assert 0.2499 < ratio[31] < 0.25
    for p in (11, 19):
        new = _circle_eigenvalues(p, 2)
        assert max(abs(x) for x in new) <= 2 / (p + 1) + 1e-12
        assert min(new) > min(_circle_eigenvalues(p, 1))


def _no_four_colouring(m, pts):
    """True if the subgraph of Cay((Z/m)^2, unit circle) induced on pts has no proper 4-colouring (CaDiCaL, with the
    colours of one unit triangle fixed, which loses nothing)"""
    from pysat.solvers import Solver
    idx = {tuple(q): i for i, q in enumerate(pts)}
    assert len(idx) == len(pts)
    T = [(a, b) for a in range(m) for b in range(m) if (a * a + b * b) % m == 1]
    E = {(min(i, j), max(i, j)) for i, (x, y) in enumerate(pts) for a, b in T
         for j in [idx.get(((x + a) % m, (y + b) % m))] if j is not None}
    adj = {}
    for u, v in E:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    tri = next((u, v, min(adj[u] & adj[v])) for u, v in sorted(E) if adj[u] & adj[v])
    s = Solver(name="cadical153")
    for v in range(len(pts)):
        s.add_clause([4 * v + c + 1 for c in range(4)])
    for u, v in E:
        for c in range(4):
            s.add_clause([-(4 * u + c + 1), -(4 * v + c + 1)])
    for c, v in enumerate(tri):
        s.add_clause([4 * v + c + 1])
    return not s.solve()


def test_eleven_adic_levels_need_five_colours():
    """notes §6, Q(sqrt47): data/quadratic_planes/padic11.json. The 69 points of level 1 and the 244 points of level 2
    have no proper 4-colouring; every level-2 point lies above a level-1 point; and the logs of the level-3 check
    (the preimage of the level-2 points, 29 524 points) say UNSATISFIABLE and VERIFIED"""
    g = load("padic11.json")
    assert g["p"] == 11 and len(g["level1"]) == 69 and len(g["level2"]) == 244
    l1 = {tuple(q) for q in g["level1"]}
    assert all((x % 11, y % 11) in l1 for x, y in g["level2"])
    assert _no_four_colouring(11, g["level1"])
    assert _no_four_colouring(121, g["level2"])
    logs = os.path.join(D, "padic11.logs")
    check = open(os.path.join(logs, "level3.check.txt"), encoding="utf-8").read()
    assert "level 3: 29524 points, 1689039 edges, 1452 unit vectors" in check
    assert "s UNSATISFIABLE" in open(os.path.join(logs, "level3.kissat.log"), encoding="utf-8").read()
    assert "s VERIFIED" in open(os.path.join(logs, "level3.drat-trim.log"), encoding="utf-8").read()

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
    VERIFIED UNSAT line for every field, for its own encoding and for the stored formula, and ends with exit
    status 0; and it checked the stored formulas as they are now (the sha256 it printed for each q{d}.cnf)"""
    log = open(os.path.join(D, "cake_lpr_checks.txt"), encoding="utf-8").read()
    for d in FIELDS:
        digest = hashlib.sha256(open(os.path.join(D, f"q{d}.cnf"), "rb").read()).hexdigest()
        assert re.search(rf"Q\(sqrt{d}\): ok: q{d}\.cnf is the 3-colouring formula of this graph \(\d+ clauses, "
                         rf"sha256 {digest}\)", log), (d, "the log does not check this q{d}.cnf")
        assert f"Q(sqrt{d}): ok: kissat UNSATISFIABLE, drat-trim VERIFIED and cake_lpr VERIFIED UNSAT" in log, d
        # the stored formula itself: the hypothesis of the Lean files of field_lean.COND_FIELDS
        assert (f"Q(sqrt{d}): ok: kissat UNSATISFIABLE, drat-trim VERIFIED and cake_lpr VERIFIED UNSAT (LRAT) on "
                f"q{d}.cnf itself") in log, d
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
    """Every field of field_lean.FIELDS and field_lean.COND_FIELDS has its file, its library in lakefile.toml (built
    by default), its line in PrintAxioms.lean and axioms.expected, and its kernel replay in the Lean workflow; every
    published graph has its file, and no other Sqrt{d}.lean exists."""
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
    for d in field_lean.COND_FIELDS:
        m = f"Sqrt{d}"
        thm = f"{m}.not_colorable_three_of_unsatisfiable" if d == 47 else f"{m}.chromaticNumber_eq_four_of_unsatisfiable"
        assert os.path.exists(os.path.join(LEAN, m + ".lean")), m
        assert f'name = "{m}"' in lake and f'"{m}"' in default, m
        assert f"import {m}\n" in pa and f"#print axioms {thm}\n" in pa, m
        assert f"'{thm}' depends on axioms: [propext, Classical.choice, Quot.sound]" in ax.splitlines(), m
        assert m in replayed, m
    for m in ("QuadraticPlanes", "ColouringFormula", "PadicPlanes"):
        assert m in replayed and f'"{m}"' in default, m
    assert "import PadicPlanes\n" in pa
    for t in ("padicSeven_chromaticNumber", "padicThree_chromaticNumber", "padicTwo_chromaticNumber"):
        assert f"#print axioms PadicPlanes.{t}\n" in pa
        assert f"'PadicPlanes.{t}' depends on axioms: [propext, Classical.choice, Quot.sound]" in ax.splitlines()
    on_disk = {int(f[4:-5]) for f in os.listdir(LEAN) if f.startswith("Sqrt") and f.endswith(".lean")}
    assert on_disk == set(field_lean.FIELDS) | set(field_lean.COND_FIELDS) == set(FIELDS)
    assert not set(field_lean.FIELDS) & set(field_lean.COND_FIELDS)


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


def test_padic_planes_table():
    """note §5, planes over Q_p: chi(Q_2^2) = 2, chi(Q_3^2) = 3, chi(Q_7^2) = 4, and the bounds up to p = 83"""
    import subprocess, sys
    out = subprocess.run([sys.executable, os.path.join(ROOT, "data", "quadratic_planes", "scripts", "padic_planes.py")],
                         capture_output=True, text=True, check=True).stdout
    assert "exact: chi(Q_2^2) = 2, chi(Q_3^2) = 3, chi(Q_7^2) = 4" in out
    rows = {int(l.split("|")[0]): [c.strip() for c in l.split("|")] for l in out.splitlines() if l[:1].isdigit()}
    assert rows[7][2] == "4" and rows[7][3] == "4" and "Q(sqrt11)" in rows[7][4]
    assert rows[11][2] == "5" and rows[11][3] == "4" and rows[19][2] == "5" and rows[19][3] == "4"
    assert rows[83][3] == "5" and all(rows[p][3] == "4" for p in (23, 31, 43, 47, 59, 67, 71, 79))


def test_padic_reach():
    """note §5: the 27 certified fields reach every prime p = 3 (mod 4) from 7 on below 10^5 (padic_reach.c finds the
    first one they miss, 2 129 503 819, which is prime, 3 mod 4, and has no d as a nonzero square)"""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "padic_reach", os.path.join(ROOT, "data", "quadratic_planes", "scripts", "padic_reach.py"))
    pr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pr)
    ds, count = pr.check(10 ** 5)
    assert len(ds) == 27 and count == 4808
    assert 11 in pr.reached(7, ds) and pr.reached(pr.P0, ds) == []


def test_padic_measurable_bounds():
    """note §5, measurable colourings of Q_p^2: the interval-arithmetic bounds 5, 5, 6, 7 at p = 23, 31, 59, 71"""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "padic_measurable", os.path.join(ROOT, "data", "quadratic_planes", "scripts", "padic_measurable.py"))
    pm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pm)
    assert {p: pm.bound(p)[2] for p in (7, 11, 23, 31, 59, 71)} == {7: 3, 11: 4, 23: 5, 31: 5, 59: 6, 71: 7}
    assert float(pm.bound(31)[1].a) > 4.00004

