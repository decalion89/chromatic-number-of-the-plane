"""Theorem D (notes/circular_planes.md §5, papers/three-colours Section 9): the exact computations behind the gap
above 3, rerun from the stored programs (data/number_fields/circular/probe/ and a referee's probe/indep2/)."""
import os, re, subprocess, sys
from collections import Counter
from fractions import Fraction

import pytest

P = os.path.join(os.path.dirname(__file__), "..", "data", "number_fields", "circular", "probe")


def run(cwd, *args, timeout=900):
    return subprocess.run([sys.executable, "-B", *args], cwd=cwd, capture_output=True, text=True, check=True,
                          timeout=timeout).stdout


def test_lemma_values_polygons_and_padic_points():
    """Lemma 8 congruences, the base-case identities, the polygons P and X(eps), the points c_k (k <= 15)."""
    out = run(os.path.join(P, "indep2"), "check1_basic.py")
    assert out.rstrip().endswith("ALL OK")
    assert "c_k checked exactly for k <= 15" in out


def test_padic_points_second_program():
    out = run(P, "far_points.py")
    lines = [l for l in out.splitlines() if l.startswith("k=")]
    assert len(lines) == 10 and all("kappa = 3/10" in l for l in lines)


def test_induction_thresholds():
    out = run(P, "lemma_lp.py", "6")
    pat = r"k=(\d+): .*?r1_C = 1/2 - s_C = (\d+/\d+).*?\n.*?r1_Q = 1/3 - eps_Q = (\d+/\d+)"
    got = {int(m.group(1)): (Fraction(m.group(2)), Fraction(m.group(3))) for m in re.finditer(pat, out)}
    assert got[2] == (Fraction(17, 56), Fraction(17, 56))
    assert got[3] == (Fraction(333, 1106), Fraction(166, 553))
    for k in (4, 5, 6):
        assert got[k] == (Fraction(3303, 10981), Fraction(13183, 43924))


def test_probe_at_theta_star():
    """At theta* = 0.3001611 only the five components of Proposition 6 remain from level 17 on."""
    out = run(P, "levels3.py", "3001611/10000000", "19")
    counts = [int(c) for c in re.findall(r"^level \d+ N=\d+: (\d+) comps", out, re.M)]
    assert counts == [13, 33, 57, 81, 97, 109, 125, 133, 105, 89, 73, 49, 33, 25, 17, 9, 5, 5, 5]
    for k in (17, 18, 19):
        line = next(l for l in out.splitlines() if l.startswith(f"level {k} "))
        assert "'X'" not in line


@pytest.mark.slow
def test_base_case_classes():
    """The 25 classes of the base case: least s with the shifts of the proof (about a minute)."""
    out = run(os.path.join(P, "indep2"), "check2_basecase.py")
    vals = Counter(re.findall(r"least s \(centred shifts\) = (\S+)", out))
    assert vals == Counter({"0": 1, "11/56": 8, "1/4": 4, "2/7": 4, "3/8": 4, "1/6": 4})
