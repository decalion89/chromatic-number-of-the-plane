"""Finite witnesses at four colours (notes/circular_planes.md §6.9, data/number_fields/circular/at_four): the
characters that bound kappa of the unit vectors of Q(sqrt 59)^2 with denominators 210 and 1050, and the referee's
check of Lemma F13 (a tight square or a character) on small finite abelian groups."""
import gzip
import hashlib
import json
import lzma
import os
import shutil
import subprocess
import sys
from fractions import Fraction as Fr

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIRC = os.path.join(ROOT, "data", "number_fields", "circular")
AT4 = os.path.join(CIRC, "at_four")
sys.path.insert(0, AT4)
from check_theta import least_margin, units  # noqa: E402


def _theta(name):
    with open(os.path.join(AT4, name)) as fh:
        return json.load(fh)["theta"]


def test_theta_210_keeps_every_unit_vector_above_a_quarter():
    """kappa(U_210) >= 15/59 > 1/4: no graph built from these vectors has chi_c = 4."""
    n, m = least_margin(59, 210, _theta("theta59_210.json"))
    assert n == 108 and m == Fr(15, 59)


@pytest.mark.slow
def test_theta_1050_margin_is_a_quarter():
    n, m = least_margin(59, 1050, _theta("theta59_1050.json"))
    assert n == 300 and m == Fr(1, 4)


def test_brute_force_units_match_the_enumerator():
    sys.path.insert(0, os.path.join(ROOT, "data", "quadratic_planes", "winding"))
    from kappaD import units as kunits
    assert sorted(units(59, 210)) == sorted(kunits(59, 210))


def test_tight_four_cycle_that_is_not_a_square():
    out = subprocess.run([sys.executable, "example_z8.py"], cwd=os.path.join(AT4, "indep"),
                         capture_output=True, text=True, check=True).stdout
    assert "tight 4-cycle 0-2-4-6: True" in out and "tight squares: []" in out
    assert "A=B everywhere: True" in out and "kappa(S) = 1/4" in out


@pytest.mark.parametrize("group", [["2", "4"], ["3", "3"], ["8"]])
def test_lemma_F13_on_small_groups(group):
    out = subprocess.run([sys.executable, "at4check_py.py", *group], cwd=os.path.join(AT4, "indep"),
                         capture_output=True, text=True, check=True, timeout=900).stdout
    last = out.strip().splitlines()[-1]
    assert last.startswith("# summary") and last.endswith("BAD=0")


def test_q3_11_vectors_and_character_at_a_quarter():
    """the 27 vectors of Q(sqrt3, sqrt11)^2 are those of the certificate; theta311 has least distance exactly 1/4."""
    out = subprocess.run([sys.executable, "q3_11.py", "cert311_open_4.json.gz"], cwd=AT4,
                         capture_output=True, text=True, check=True, timeout=600).stdout
    assert "pairwise distinct up to sign" in out
    assert "the certificate's units are these vectors up to sign: True" in out
    assert "least distance to Z over the 27 vectors: 1/4" in out


@pytest.mark.parametrize("checker", ["check_open.py", "check_open_indep.py"])
def test_q3_11_open_certificate_both_checkers(checker):
    """kappa <= 1/4 for the 27 vectors: no character maps them into (1/4, 3/4)."""
    out = subprocess.run([sys.executable, checker, os.path.join("at_four", "cert311_open_4.json.gz")], cwd=CIRC,
                         capture_output=True, text=True, check=True, timeout=1800).stdout
    assert out.startswith(("VERIFIED", "ACCEPTED")) and "3,11" in out and "27 units" in out
    assert "12073 nodes" in out


def test_checkers_reject_a_wrong_unit_over_q3_11(tmp_path):
    import gzip
    with gzip.open(os.path.join(AT4, "cert311_open_4.json.gz"), "rt") as fh:
        C = json.load(fh)
    C["units"][0] = [x + (1 if k == 0 else 0) for k, x in enumerate(C["units"][0])]
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(C))
    for cmd in (["check_open.py", str(bad)], ["check_open_indep.py", str(bad)]):
        res = subprocess.run([sys.executable] + cmd, cwd=CIRC, capture_output=True, text=True, timeout=600)
        assert res.returncode != 0 or not res.stdout.startswith(("VERIFIED", "ACCEPTED")), cmd


FW4 = os.path.join(CIRC, "finite_witness")
DRAT = shutil.which("drat-trim") or os.environ.get("DRAT_TRIM")
W4 = ["check_witness4.py", "witness_q3_11.json.gz", "witness_q3_11.cnf.gz"]


def test_witness_q3_11_checker(tmp_path):
    """The explicit graph with chi_c = 4 over Q(sqrt3, sqrt11) (finite_witness/README.md): 1874 points, all 8085 unit
    pairs as edges, a proper 4-colouring, 4992 cycles of lengths divisible by 4, and the stored formula is the one
    check_witness4.py builds (drat-trim refutes it in the slow test below)."""
    for n in W4:
        shutil.copy(os.path.join(FW4, n), tmp_path / n)
    out = subprocess.run([sys.executable, "check_witness4.py", W4[1], W4[2]], cwd=tmp_path, capture_output=True,
                         text=True, check=True, timeout=900).stdout
    assert "(1) Q(sqrt3, sqrt11): 1874 points, 8085 edges, every edge at distance exactly 1" in out
    assert "not listed as edges: 0 (induced)" in out
    assert "(2) the colouring is a proper 4-colouring" in out and "(3) 4992 cycles" in out
    assert "sha256 febd6b3d96bb8bd43962882ab3bf09bea13a9f8cb3ff563f2691dc340d3dd26f" in out


def test_witness_q3_11_referee_encoding(tmp_path):
    """The referee's program (indep_W4/), which shares no code with check_witness4.py, recomputes every unit pair and
    writes the second encoding, the formula that kissat, drat-trim and cake_lpr refuted (indep_W4/results/)."""
    r = subprocess.run([sys.executable, os.path.join(FW4, "indep_W4", "ref_check4.py"),
                        os.path.join(FW4, "witness_q3_11.json.gz"), str(tmp_path), "--no-sat"],
                       capture_output=True, text=True, timeout=900)
    assert "8085 unit-distance pairs" in r.stdout and "declared edge set == recomputed edge set" in r.stdout
    assert "REFEREE: INCOMPLETE (tasks 1-4 passed" in r.stdout and r.returncode == 2
    assert hashlib.sha256((tmp_path / "ref4_main.cnf").read_bytes()).hexdigest() == \
        "1f28389f2aec69a069146c82b7c1a7870a5955fa72a380e21a8d3335efcdcde6"


def test_witness_q3_11_checker_rejects_a_moved_point(tmp_path):
    """A point moved by 1/84 in one coordinate breaks an edge: check_witness4.py rejects the file."""
    shutil.copy(os.path.join(FW4, W4[0]), tmp_path / W4[0])
    W = json.load(gzip.open(os.path.join(FW4, W4[1]), "rt"))
    v = W["edges"][0][1]
    W["points"][v][0] += 1
    with gzip.open(tmp_path / "bad.json.gz", "wt") as f:
        json.dump(W, f)
    r = subprocess.run([sys.executable, W4[0], "bad.json.gz"], cwd=tmp_path, capture_output=True, text=True,
                       timeout=900)
    assert r.returncode != 0 and "not a unit-distance edge" in r.stderr


W4B = ["check_witness4.py", "witness_q2_3.json.gz", "witness_q2_3.cnf.gz"]


def test_witness_q2_3_checker(tmp_path):
    """The second explicit graph with chi_c = 4, over Q(sqrt2, sqrt3) (finite_witness/README.md): 1657 points over
    the basis (1, sqrt2, sqrt3, sqrt6) with denominator 36, all 6238 unit pairs as edges, a proper 4-colouring, 6062
    cycles of lengths 4 and 8, and the stored formula is the one check_witness4.py builds."""
    for n in W4B:
        shutil.copy(os.path.join(FW4, n), tmp_path / n)
    out = subprocess.run([sys.executable, "check_witness4.py", W4B[1], W4B[2]], cwd=tmp_path, capture_output=True,
                         text=True, check=True, timeout=900).stdout
    assert "(1) Q(sqrt2, sqrt3): 1657 points, 6238 edges, every edge at distance exactly 1" in out
    assert "not listed as edges: 0 (induced)" in out
    assert "(2) the colouring is a proper 4-colouring" in out and "(3) 6062 cycles" in out and "(lengths 4..8" in out
    assert "sha256 cd6382b7729db29c82b046dfd047b74a98b42bba41048007906aebf3ed5beff2" in out


def test_witness_q2_3_referee_encoding(tmp_path):
    """The referee's program for the second witness (indep_W4b/), written from the file format alone, recomputes all
    1 371 996 pairs and writes its own encoding, the formula that kissat, drat-trim and cake_lpr refuted."""
    r = subprocess.run([sys.executable, os.path.join(FW4, "indep_W4b", "ref_check_z.py"), "--no-sat",
                        "--cnf", str(tmp_path / "ref_z.cnf"), "--results", str(tmp_path / "res")],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout[-2000:]
    assert "RESULT: PASS (tasks 1-4, solver runs skipped)" in r.stdout
    assert hashlib.sha256((tmp_path / "ref_z.cnf").read_bytes()).hexdigest() == \
        "a8e9d4bdee31892b03d2f6bbe4e865f38414434cd72be38768e3a65bf2483553"


def test_witness_q2_3_checker_rejects_a_moved_point(tmp_path):
    """A point moved by 1/36 in one coordinate breaks an edge: check_witness4.py rejects the file."""
    shutil.copy(os.path.join(FW4, W4B[0]), tmp_path / W4B[0])
    W = json.load(gzip.open(os.path.join(FW4, W4B[1]), "rt"))
    v = W["edges"][0][1]
    W["points"][v][0] += 1
    with gzip.open(tmp_path / "bad.json.gz", "wt") as f:
        json.dump(W, f)
    r = subprocess.run([sys.executable, W4B[0], "bad.json.gz"], cwd=tmp_path, capture_output=True, text=True,
                       timeout=900)
    assert r.returncode != 0 and "not a unit-distance edge" in r.stderr


@pytest.mark.slow
@pytest.mark.skipif(not DRAT, reason="drat-trim not found (PATH or DRAT_TRIM)")
def test_witness_q2_3_proof(tmp_path):
    """drat-trim verifies the stored DRAT proof for the second witness: chi_c = chi = 4 over Q(sqrt2, sqrt3)."""
    for n in W4B:
        shutil.copy(os.path.join(FW4, n), tmp_path / n)
    with lzma.open(os.path.join(FW4, "witness_q2_3.drat.xz"), "rb") as src, \
            open(tmp_path / "witness_q2_3.drat", "wb") as dst:
        shutil.copyfileobj(src, dst)
    out = subprocess.run([sys.executable, W4B[0], W4B[1], W4B[2], "witness_q2_3.drat", DRAT], cwd=tmp_path,
                         capture_output=True, text=True, check=True, timeout=3600).stdout
    assert "(5) drat-trim: s VERIFIED" in out and "chi_c(H) = chi(H) = 4" in out


KISSAT = shutil.which("kissat") or os.environ.get("KISSAT")


@pytest.mark.skipif(not KISSAT, reason="kissat not found (PATH or KISSAT)")
def test_effective_construction_on_a_distance_graph():
    """The remark after Lemma F17 (note, end of 6.9), run by the referee's program (indep_W4/remark_check.py) for
    U = D = {2, 3, 5, 6} in Z, where kappa(D) <= 1/4: the relations rho_p from the rational cones leave the period
    system without a solution, and every proper 4-colouring of the graph G on their chains has a tight cycle, by a SAT
    check that does not use the lemma; the chains of the basis alone give a 4-colouring without tight cycles."""
    r = subprocess.run([sys.executable, os.path.join(FW4, "indep_W4", "remark_check.py"), "2", "3", "5", "6"],
                       capture_output=True, text=True, timeout=600, env=dict(os.environ, KISSAT=KISSAT))
    assert r.returncode == 0, r.stderr
    assert "|Pi| = 4" in r.stdout and "(iii) p in a box of 448 candidates satisfying every range on S: 0" in r.stdout
    assert "G_S: |H| = 26, |E| = 86, CNF 926 vars 6310 clauses -> ['s UNSATISFIABLE']" in r.stdout
    assert "G_basis_only: |H| = 9, |E| = 15, CNF 138 vars 576 clauses -> ['s SATISFIABLE']" in r.stdout
    assert "RESULT: every proper 4-colouring of G has a tight cycle" in r.stdout
    assert "(23,6)-colourability of G" in r.stdout and "['s UNSATISFIABLE'] ; 4-colourable: ['s SATISFIABLE']" in r.stdout


@pytest.mark.slow
@pytest.mark.skipif(not DRAT, reason="drat-trim not found (PATH or DRAT_TRIM)")
def test_witness_q3_11_proof(tmp_path):
    """drat-trim verifies the stored DRAT proof that every proper 4-colouring of the witness has a tight listed
    cycle, so chi_c = chi = 4."""
    for n in W4:
        shutil.copy(os.path.join(FW4, n), tmp_path / n)
    with lzma.open(os.path.join(FW4, "witness_q3_11.drat.xz"), "rb") as src, \
            open(tmp_path / "witness_q3_11.drat", "wb") as dst:
        shutil.copyfileobj(src, dst)
    out = subprocess.run([sys.executable, W4[0], W4[1], W4[2], "witness_q3_11.drat", DRAT], cwd=tmp_path,
                         capture_output=True, text=True, check=True, timeout=3600).stdout
    assert "(5) drat-trim: s VERIFIED" in out and "chi_c(H) = chi(H) = 4" in out


# The subgraphs of H4 and H4' in the Cayley graphs of their construction's unit vectors (finite_witness/README.md):
# file stem, field, points, edges, generators, cycles, unit pairs among the points, sha256 of the formula
CAY = {"q3_11": ("witness_q3_11_cayley", "Q(sqrt3, sqrt11)", 1874, 7887, 27, 3389, 8085,
                 "b394ee9342158c22354cb00fa0d7b373aa98b48d33a0795f62993ecf6b531230"),
       "q2_3": ("witness_q2_3_cayley", "Q(sqrt2, sqrt3)", 1657, 6199, 60, 5264, 6238,
                "802921acb1e509aa3ab79f3ba9efaed50b036cb3d9d7277027ee5bb918b0b85b")}


@pytest.mark.parametrize("which", ["q3_11", "q2_3"])
def test_cayley_witness_checker(which, tmp_path):
    """The same points as H4 (H4'), with only the pairs that differ by one of the construction's unit vectors or its
    negative as edges: check_witness4.py checks the generators, that the edges are exactly these pairs, that every point
    is joined to the origin, the colouring, the cycles and the stored formula."""
    name, field, n, e, g, c, allpairs, sha = CAY[which]
    for f in ("check_witness4.py", name + ".json.gz", name + ".cnf.gz"):
        shutil.copy(os.path.join(FW4, f), tmp_path / f)
    out = subprocess.run([sys.executable, "check_witness4.py", name + ".json.gz", name + ".cnf.gz"], cwd=tmp_path,
                         capture_output=True, text=True, check=True, timeout=900).stdout
    assert f"(1) {field}: {n} points, {e} edges, every edge at distance exactly 1" in out
    assert f"unit-distance pairs among the points: {allpairs}; not listed as edges: {allpairs - e} (not induced)" in out
    assert f"(1b) {g} generators U" in out and f"{n} of the {n} points are joined to the fixed vertex (the origin)" in out
    assert "(2) the colouring is a proper 4-colouring" in out and f"(3) {c} cycles" in out and f"sha256 {sha}" in out


def _units_q2_3(D):
    """the 120 vectors zeta_24^j w^l (j mod 24, |l| <= 2), w = (1 + 2 sqrt(-2))/3, over (1, sqrt2, sqrt3, sqrt6)/D"""
    def mul(p, q):
        return (p[0] * q[0] + 2 * p[1] * q[1] + 3 * p[2] * q[2] + 6 * p[3] * q[3],
                p[0] * q[1] + p[1] * q[0] + 3 * (p[2] * q[3] + p[3] * q[2]),
                p[0] * q[2] + p[2] * q[0] + 2 * (p[1] * q[3] + p[3] * q[1]),
                p[0] * q[3] + p[3] * q[0] + p[1] * q[2] + p[2] * q[1])

    def cmul(z, w):
        re = tuple(s - t for s, t in zip(mul(z[0], w[0]), mul(z[1], w[1])))
        im = tuple(s + t for s, t in zip(mul(z[0], w[1]), mul(z[1], w[0])))
        return re, im
    q = Fr(1, 4)
    zeta, one = ((0, q, 0, q), (0, -q, 0, q)), ((1, 0, 0, 0), (0, 0, 0, 0))
    w, wb = ((Fr(1, 3), 0, 0, 0), (0, Fr(2, 3), 0, 0)), ((Fr(1, 3), 0, 0, 0), (0, Fr(-2, 3), 0, 0))
    pw = [cmul(wb, wb), wb, one, w, cmul(w, w)]
    out, z = set(), one
    for _ in range(24):
        for v in pw:
            u = [Fr(t) * D for t in cmul(z, v)[0] + cmul(z, v)[1]]
            assert all(t.denominator == 1 for t in u)
            out.add(tuple(int(t) for t in u))
        z = cmul(z, zeta)
    return out


def test_cayley_witness_generators():
    """The stored generators are the construction's vectors up to sign: the 27 of at_four/q3_11.py, and the 60
    vectors zeta_24^j w^l (one of each pair +-u), rebuilt here."""
    from q3_11 import units_311
    sign = lambda S: {min(tuple(u), tuple(-t for t in u)) for u in S}
    W = json.load(gzip.open(os.path.join(FW4, "witness_q3_11_cayley.json.gz"), "rt"))
    assert len(W["generators"]) == 27 and sign(W["generators"]) == sign(units_311())
    W = json.load(gzip.open(os.path.join(FW4, "witness_q2_3_cayley.json.gz"), "rt"))
    U = _units_q2_3(W["denominator"])
    assert len(U) == 120 and len(W["generators"]) == 60 and sign(W["generators"]) == sign(U)


def test_cayley_witness_checker_rejects_a_foreign_edge(tmp_path):
    """Adding to the Cayley subgraph one of the 198 unit pairs of H4 in another direction breaks the generator check."""
    shutil.copy(os.path.join(FW4, "check_witness4.py"), tmp_path / "check_witness4.py")
    W = json.load(gzip.open(os.path.join(FW4, "witness_q3_11_cayley.json.gz"), "rt"))
    H = json.load(gzip.open(os.path.join(FW4, "witness_q3_11.json.gz"), "rt"))
    assert W["points"] == H["points"]
    have = {tuple(e) for e in W["edges"]}
    extra = next(e for e in H["edges"] if (min(e), max(e)) not in have)
    W["edges"].append(extra)
    with gzip.open(tmp_path / "bad.json.gz", "wt") as f:
        json.dump(W, f)
    r = subprocess.run([sys.executable, "check_witness4.py", "bad.json.gz"], cwd=tmp_path, capture_output=True,
                       text=True, timeout=900)
    assert r.returncode != 0 and "the edges are not the pairs of points that differ by an element" in r.stderr


@pytest.mark.slow
@pytest.mark.skipif(not DRAT, reason="drat-trim not found (PATH or DRAT_TRIM)")
@pytest.mark.parametrize("which", ["q3_11", "q2_3"])
def test_cayley_witness_proof(which, tmp_path):
    """drat-trim verifies the stored DRAT proof: every proper 4-colouring of the Cayley subgraph has a tight listed
    cycle, so chi_c = chi = 4."""
    name = CAY[which][0]
    for f in ("check_witness4.py", name + ".json.gz", name + ".cnf.gz"):
        shutil.copy(os.path.join(FW4, f), tmp_path / f)
    with lzma.open(os.path.join(FW4, name + ".drat.xz"), "rb") as src, open(tmp_path / (name + ".drat"), "wb") as dst:
        shutil.copyfileobj(src, dst)
    out = subprocess.run([sys.executable, "check_witness4.py", name + ".json.gz", name + ".cnf.gz", name + ".drat",
                          DRAT], cwd=tmp_path, capture_output=True, text=True, check=True, timeout=3600).stdout
    assert "(5) drat-trim: s VERIFIED" in out and "chi_c(H) = chi(H) = 4" in out


def test_quarter_character_residues_and_maps():
    """Proposition F15': in Z_2[zeta_12]/4 exactly 24 residues have r rbar = 1 (zeta^k and (1 + 2i) zeta^k), they are
    the residues of the norm-one elements, and exactly twelve additive maps to Z/4 vanish on none of them."""
    import quarter_character as qc
    n1, nunits, same, chars, rot, kern = qc.part1()
    assert (n1, nunits, same, len(chars), rot, kern) == (24, 192, True, 12, True, True)
    assert qc.L0 in chars


def test_quarter_character_on_unit_vectors():
    """c = lambda(x + iy) is 1, 2 or 3 on every unit vector with the given denominators over Q(sqrt59) and Q(sqrt83)
    (and agrees there with the closed formula x1 + 3 y1 + s (3 y2 - x2) mod 4)."""
    import quarter_character as qc
    out = qc.part2_quadratic(59, [210, 1050])
    assert out[210] == {0: 0, 1: 34, 2: 40, 3: 34} and out[1050] == {0: 0, 1: 96, 2: 108, 3: 96}
    out = qc.part2_quadratic(83, [210, 630])
    assert all(cnt[0] == 0 for cnt in out.values()) and sum(out[630].values()) == 60


def test_quarter_character_q3_11():
    """the same over Q(sqrt3, sqrt11): the 27 vectors and the edges of witness_q3_11 and of its Cayley subgraph."""
    import quarter_character as qc
    res = qc.part2_q3_11()
    assert set(res) == {"27", "witness_q3_11", "witness_q3_11_cayley"}
    assert all(cnt[0] == 0 for cnt in res.values())
    assert sum(res["witness_q3_11"].values()) == 8085 and sum(res["witness_q3_11_cayley"].values()) == 7887


def test_quarter_character_has_no_analogue_over_q2_3():
    """over Q_2(sqrt2, sqrt3): 384 norm-one residues modulo 4, and no additive map to Z/4 is nonzero on all of them."""
    import quarter_character as qc
    assert qc.part3_q2_3() == (384, 0)
