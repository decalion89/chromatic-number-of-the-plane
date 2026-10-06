"""Finite witnesses at four colours (notes/circular_planes.md §6.9, data/number_fields/circular/at_four): the
characters that bound kappa of the unit vectors of Q(sqrt 59)^2 with denominators 210 and 1050, and the referee's
check of Lemma F13 (a tight square or a character) on small finite abelian groups."""
import json
import os
import subprocess
import sys
from fractions import Fraction as Fr

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AT4 = os.path.join(ROOT, "data", "number_fields", "circular", "at_four")
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
