"""The decompositions of 2 and 3 recomputed with PARI/GP (scripts/decompositions.gp).

notes/local_colourings.md uses them in sections 8, 10 and 11: the places over 2 of Q(sqrt3, sqrt11) and
Q(sqrt2, sqrt3) and their inertness in L(i), and, for L = Q(sqrt3, sqrt q), a place over 2 with residue
field F_2 exactly when q = 2 or q = 1, 3 (mod 8), and one over 3 with residue field F_3 exactly when
q = 1 (mod 3). The test runs the script and checks every line it prints. It is skipped when gp is not
installed; the continuous integration installs it.
"""
import ast
import os
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

pytestmark = pytest.mark.skipif(shutil.which("gp") is None, reason="PARI/GP (gp) is not installed")


def _lines():
    out = subprocess.run(["gp", "-q", os.path.join("scripts", "decompositions.gp")], cwd=ROOT,
                         capture_output=True, text=True, timeout=600, check=True).stdout
    rows = []
    for line in out.splitlines():
        gens, rest = line.split("] ", 1)
        p, dec = rest.split(" ", 1)
        rows.append((tuple(ast.literal_eval(gens + "]")), int(p), [tuple(d) for d in ast.literal_eval(dec)]))
    return rows


def test_the_two_theorems():
    got = {(v, p): d for v, p, d in _lines()}
    assert got[(3, 11), 2] == [(2, 1), (2, 1)]      # L = Q(sqrt3, sqrt11): 2 = P1^2 P2^2, residue fields F_2
    assert got[(-1, 3, 11), 2] == [(2, 2), (2, 2)]  # both inert in L(i), residue fields F_4
    assert got[(2, 3), 2] == [(4, 1)]               # L = Q(sqrt2, sqrt3): 2 = P^4, residue field F_2
    assert got[(-1, 2, 3), 2] == [(4, 2)]           # inert in L(i) = Q(zeta24), residue field F_4


def test_two_square_roots():
    rows = [(v[1], p, d) for v, p, d in _lines() if len(v) == 2 and v[0] == 3]
    qs = sorted({q for q, _, _ in rows})
    assert qs == [q for q in range(2, 75) if q != 3 and all(q % r for r in range(2, q))]
    for q, p, d in rows:
        assert sum(e * f for e, f in d) == 4
        if p == 2:
            assert (min(f for _, f in d) == 1) == (q == 2 or q % 8 in (1, 3)), q
        if p == 3:
            assert (min(f for _, f in d) == 1) == (q % 3 == 1), q


NINE = [(3, 7), (3, 7, 13), (3, 7, 13, 19), (3, 5), (3, 11), (3, 5, 7), (3, 7, 11), (3, 13, 17),
        (3, 7, 13, 23)]


def test_more_square_roots():
    """Residue degree 1 over 3 exactly when every generator other than 3 is 1 (mod 3)."""
    got = {v: d for v, p, d in _lines() if p == 3}
    for v in NINE:
        d = got[v]
        assert sum(e * f for e, f in d) == 2 ** len(v)
        assert {f for _, f in d} == {1 if all(g % 3 == 1 for g in v[1:]) else 2}, v
