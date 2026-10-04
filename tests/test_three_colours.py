"""Theorem B (notes/three_colours_number_fields.md): the limit test of Lemma B1 on known and new fields."""
import os, subprocess, sys, shutil
import pytest

D = os.path.join(os.path.dirname(__file__), "..", "data", "number_fields", "three_colours")
pytestmark = pytest.mark.skipif(shutil.which("gp") is None, reason="PARI/GP not installed")


def run(*args, timeout=1200):
    return subprocess.run([sys.executable, *args], cwd=D, capture_output=True, text=True, check=True,
                          timeout=timeout).stdout


def test_known_quadratic_cases():
    out = run("known_cases.py")
    for d in (11, 23, 35, 47, 59):
        assert f"{d} " in out
    lines = {int(l.split()[0]): l for l in out.splitlines() if l.strip()}
    for d in (11, 23, 35, 47, 59):
        assert "(False, None)" in lines[d]
    for d in (3, 7):
        assert "(True," in lines[d]


def test_sqrt2_sqrt7_needs_four_colours():
    out = run("check_saved.py", "V_sqrt2_sqrt7.json")
    assert "INFEASIBLE" in out


def test_sqrt2_sqrt7_in_gp():
    gp = shutil.which("gp")
    v = os.path.join(D, "_v27.gp"); r = os.path.join(D, "_run27.gp")
    try:
        run("to_gp.py", "V_sqrt2_sqrt7.json", "_v27.gp")
        with open(r, "w") as fh:
            fh.write(open(v).read() + open(os.path.join(D, "check.gp")).read())
        out = subprocess.run([gp, "-q", "-s", "2000000000", r], capture_output=True, text=True, check=True,
                             timeout=1200).stdout
        assert "rank of L = 8 (expected 8)" in out and "INFEASIBLE" in out
    finally:
        for f in (v, r):
            if os.path.exists(f):
                os.remove(f)


def test_two_roots_certificate():
    out = run("check_mq.py", "cert_two_roots_2_7_N25.json")
    assert out.startswith("VERIFIED") and "50 unit vectors" in out


def test_criterion_examples():
    gp = shutil.which("gp")
    out = subprocess.run([gp, "-q", "classify_examples.gp"], cwd=D, capture_output=True, text=True, check=True,
                         timeout=600).stdout
    want = {"Q(sqrt2,sqrt7)": "[0, 0]", "Q(sqrt3,sqrt5)": "[0, 0]", "Q(sqrt2,sqrt3)": "[0, 0]",
            "Q(c7,sqrt7)": "[0, 0]", "Q(sqrt11)": "[0, 0]", "Q(sqrt3,sqrt7)": "[0, 1]", "Q(sqrt7)": "[0, 1]",
            "Q(sqrt3)": "[0, 1]", "Q(c7,sqrt2)": "[1, 0]", "Q(sqrt2)": "[1, 0]", "Q(cbrt2)": "[1, 1]"}
    lines = dict(l.split(": [a, b] = ") for l in out.splitlines() if ": [a, b] = " in l)
    for k, v in want.items():
        assert lines[k] == v, (k, lines.get(k))


CIRC = os.path.join(os.path.dirname(__file__), "..", "data", "number_fields", "circular")


def test_circular_seven_halves_certificate():
    out = subprocess.run([sys.executable, "check_open.py", "cert_sqrt11_open_7_2_N25.json.gz"], cwd=CIRC,
                         capture_output=True, text=True, check=True, timeout=1200).stdout
    assert out.startswith("VERIFIED") and "(2/7, 5/7)" in out and "70 units" in out


def test_seven_adic_residue_colouring():
    # a + bi -> 2a + 3b takes only the values 2..5 on the norm-one elements of F_49 (Proposition C1, p = 7)
    p = 7
    vals = set()
    for a in range(p):
        for b in range(p):
            if (a * a + b * b) % p == 1:
                vals.add((2 * a + 3 * b) % p)
    assert vals == {2, 3, 4, 5}
