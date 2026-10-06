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
