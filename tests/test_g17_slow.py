"""Part A of scripts/g17_alpha.py solved again, by CaDiCaL instead of kissat (marked slow): no independent set
of 58 points of G_17 contains a point together with its whole circle. The DRAT-checked kissat runs are in
certificates/g17_part_a_checks.txt."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pysat.solvers import Solver

import g17_alpha as g


@pytest.mark.slow
def test_part_a_is_unsatisfiable():
    """Each formula A_c has no model."""
    for c in g.INDEPENDENT_CIRCLES:
        text, _ = g.formula_a(c)
        clauses = [list(map(int, l.split()[:-1])) for l in text.split("\n")[1:] if l]
        with Solver(name="cd19", bootstrap_with=clauses) as s:
            assert s.solve() is False, c
