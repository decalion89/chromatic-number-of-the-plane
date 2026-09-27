"""Independent sets of the anisotropic plane G_17 (scripts/g17_alpha.py): the graph, its independent circles,
the 57-point rosette, the totalizer encoding, and the formulas of part A, whose DRAT check is recorded in
certificates/g17_part_a_checks.txt. tests/test_g17_slow.py solves them again with CaDiCaL."""
import hashlib
import itertools
import os
import re
import sys

HN = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(HN, "scripts"))

from pysat.solvers import Solver

import g17_alpha as g


def test_the_graph():
    assert len(g.POINTS) == 289 and len(g.UNITS) == 18 and len(g.EDGES) == 289 * 18 // 2
    # z -> l z + a and z -> l conj(z) + a, N(l) = 1, map edges to edges
    mul = lambda u, v: ((u[0] * v[0] + g.N0 * u[1] * v[1]) % g.Q, (u[0] * v[1] + u[1] * v[0]) % g.Q)
    edges = set(g.EDGES)
    for l in g.UNITS[:4]:
        for conj in (0, 1):
            for a in [(0, 0), (1, 5), (7, 3)]:
                img = {z: g.INDEX[g.add(mul(l, (z[0], -z[1] % g.Q) if conj else z), a)] for z in g.POINTS}
                assert all(tuple(sorted((img[g.POINTS[u]], img[g.POINTS[v]]))) in edges for u, v in g.EDGES)


def test_the_independent_circles():
    assert g.INDEPENDENT_CIRCLES == [4, 5, 9, 11, 12, 14, 15]
    assert all(len(g.circle(c)) == 18 for c in range(1, 17))
    assert {c: len(g.region(c)) for c in g.INDEPENDENT_CIRCLES} == {4: 90, 5: 108, 9: 90, 11: 90, 12: 108,
                                                                  14: 90, 15: 90}


def test_the_rosette():
    r = g.ROSETTE
    assert len(set(r)) == 57 and g.independent(r)
    assert (0, 0) in r and set(g.circle(12)) <= set(r)
    counts = {c: sum(1 for z in r if z != (0, 0) and g.norm(z) == c) for c in range(1, 17)}
    assert {c: k for c, k in counts.items() if k} == {3: 6, 4: 12, 5: 8, 6: 4, 9: 6, 12: 18, 14: 2}
    # the other 38 points lie in the region of the circle N = 12
    assert set(r) - set(g.circle(12)) - {(0, 0)} <= set(g.region(12))


def test_the_totalizer_is_exact():
    """With the literals fixed, the clauses are satisfiable exactly when at least t of them are true."""
    for n in range(1, 9):
        for t in range(1, n + 1):
            clauses, _ = g.totalizer(list(range(1, n + 1)), t, n)
            with Solver(name="cd19", bootstrap_with=clauses) as s:
                for x in itertools.product((0, 1), repeat=n):
                    ok = s.solve(assumptions=[i + 1 if b else -(i + 1) for i, b in enumerate(x)])
                    assert ok == (sum(x) >= t), (n, t, x)


def test_part_a_formulas_are_those_of_the_certificate_log():
    logged = {}
    for line in open(os.path.join(HN, "certificates", "g17_part_a_checks.txt")):
        m = re.match(r"A_(\d+): .* sha256 ([0-9a-f]{64}); kissat UNSAT .* drat-trim VERIFIED", line)
        if m:
            logged[int(m.group(1))] = m.group(2)
    assert sorted(logged) == g.INDEPENDENT_CIRCLES
    for c in g.INDEPENDENT_CIRCLES:
        text, _ = g.formula_a(c)
        assert hashlib.sha256(text.encode()).hexdigest() == logged[c]


def test_the_rosette_breaks_one_clause_of_formula_b():
    """Of the clauses of formula B on the vertex variables (edges, circles, vertex 0), the rosette breaks
    exactly one: its centre with its whole circle N = 12. Part A is what rules such sets out."""
    s = {g.INDEX[z] + 1 for z in g.ROSETTE}
    text = g.formula_b()
    body = [list(map(int, l.split()[:-1])) for l in text.split("\n")[1:] if l]
    top = int(text.split()[2])
    original = [c for c in body if all(abs(l) <= 289 for l in c)]
    broken = [c for c in original if not any((l > 0) == (abs(l) in s) for l in c)]
    # the only violated clause among those on the vertex variables: centre 0 with its whole circle N = 12
    assert len(broken) == 1 and sorted(-l - 1 for l in broken[0]) == sorted([0] + [g.INDEX[z] for z in g.circle(12)])
    assert top > 289

