"""alpha(G_13) = 36 (notes/g13.md): the graph, the encodings of the formulas, and the certificate, without running
kissat or drat-trim.

- The graph: 169 vertices, 14 unit vectors, 1 183 edges, the independent circles N = 1, 6, 7, 9, 11; the 4 732 maps
  z -> l z + a, l conj(z) + a (N(l) = 1) are automorphisms; the 15 known sets of 36 points are independent.
- The encodings of scripts/g13/g13cnf.py, checked by brute force against their definitions on small inputs: the
  totalizer (at least t, at most t), the counter whose outputs are forced both ways, the lex-leader chains (also
  with the automorphisms of G_13 themselves), value precedence.
- The certificate: the formulas of part A, written again by scripts/g13/enum_cert.py, have the SHA-256 of the log;
  the leaves of part B cover every assignment, and formula E37_B under each leaf has the SHA-256 of a drat-trim
  VERIFIED line of the log.
"""
import hashlib
import itertools
import os
import random
import sys

import pytest
from pysat.solvers import Solver

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "g13"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import g13                                    # noqa: E402
import g13cnf as gc                           # noqa: E402
import verify_g13 as vg                       # noqa: E402

SET36 = vg.KNOWN_G13[0]


def test_the_graph():
    assert g13.Q == 13 and len(g13.POINTS) == 169 and len(g13.UNITS) == 14 and len(g13.EDGES) == 1183
    assert all(len(set(a)) == 14 for a in g13.ADJ)
    assert all(v in g13.ADJ[w] for v in range(169) for w in g13.ADJ[v])
    maps = g13.automorphisms()
    assert len({tuple(m) for m in maps}) == len(maps) == 4732
    edges = set(g13.EDGES)
    for m in maps:
        assert sorted(m) == list(range(169))
        assert all(tuple(sorted((m[a], m[b]))) in edges for a, b in g13.EDGES)
    circles = {c: g13.circle(c) for c in range(1, 13)}
    assert all(len(p) == 14 for p in circles.values())
    assert [c for c, p in circles.items() if not any(b in g13.ADJ[g13.INDEX[a]] for a in p for b in
                                                        [g13.INDEX[z] for z in p])] == [1, 6, 7, 9, 11]
    for S in vg.KNOWN_G13:
        assert len(set(S)) == 36 and not any(b in g13.ADJ[a] for a, b in itertools.combinations(S, 2))


def check(cnf, cases):
    """cases: (assumptions, whether the formula has a model with them)"""
    if [] in cnf.clauses:                     # totalizer_atleast with t > n adds the empty clause
        assert not any(want for _, want in cases)
        return
    with Solver(name="cd19", bootstrap_with=cnf.clauses) as s:
        for assumptions, want in cases:
            assert s.solve(assumptions=assumptions) == want, assumptions


def fix(lits, bits):
    return [l if b else -l for l, b in zip(lits, bits)]


@pytest.mark.parametrize("n", range(1, 8))
def test_the_totalizer(n):
    lits = list(range(1, n + 1))
    for t in range(0, n + 2):
        for enc, holds in ((gc.totalizer_atleast, lambda k: k >= t), (gc.totalizer_atmost, lambda k: k <= t)):
            cnf = gc.CNF()                    # auxiliary variables start above the 845 colour variables
            enc(cnf, lits, t)
            check(cnf, [(fix(lits, bits), holds(sum(bits))) for bits in itertools.product((0, 1), repeat=n)])


@pytest.mark.parametrize("n", range(1, 8))
def test_the_counter_outputs_are_forced_both_ways(n):
    lits = list(range(1, n + 1))
    for t in range(1, n + 1):
        cnf = gc.CNF()
        out = gc.unary_count(cnf, lits, t)
        assert len(out) == t
        cases = []
        for bits in itertools.product((0, 1), repeat=n):
            right = [o if sum(bits) >= j + 1 else -o for j, o in enumerate(out)]
            cases.append((fix(lits, bits) + right, True))
            cases += [(fix(lits, bits) + [-r], False) for r in right]
        check(cnf, cases)


def test_the_lex_leader_chains():
    rng = random.Random(20260927)
    for _ in range(400):
        n = rng.randint(2, 7)
        perms = [rng.sample(range(n), n) for _ in range(rng.randint(1, 3))]
        order = rng.sample(range(n), n)
        L = rng.randint(1, n)
        lits = list(range(1, n + 1))
        cnf = gc.CNF()
        gc.lex_leader_sets(cnf, lits, perms, order, L)
        check(cnf, [(fix(lits, bits),
                     all([bits[v] for v in order[:L]] >= [bits[p[v]] for v in order[:L]] for p in perms))
                    for bits in itertools.product((0, 1), repeat=n)])


def test_the_lex_leader_chains_with_the_group_of_g13():
    """the 4 732 automorphisms on the first 25 positions of lex_order(): the images of a set satisfy the chains
    exactly when their 25-prefix is the largest in the orbit (every image with that prefix, and 150 others at
    random, for three sets)"""
    maps, order = g13.automorphisms(), gc.lex_order()
    cnf = gc.CNF()
    lits = [1 + v for v in range(169)]
    gc.lex_leader_sets(cnf, lits, maps, order, 25)          # a group: the inverses are the same maps
    prefix = lambda T: [int(v in T) for v in order[:25]]
    rng = random.Random(13)
    cases = []
    for S in (rng.sample(range(169), 3), rng.sample(range(169), 8), SET36):
        images = sorted({frozenset(m[v] for v in S) for m in maps}, key=sorted)
        best = max(prefix(T) for T in images)
        top = [T for T in images if prefix(T) == best]
        rest = [T for T in images if prefix(T) != best]
        for T in top + rng.sample(rest, min(150, len(rest))):
            cases.append((fix(lits, [int(v in T) for v in range(169)]), T in top))
    check(cnf, cases)


@pytest.mark.parametrize("colours", [[1, 2, 3, 4], [1, 2, 3]])
def test_value_precedence(colours):
    for m in range(1, 6):
        cnf = gc.CNF()
        gc.value_precedence(cnf, colours, list(range(m)))
        cases = []
        for col in itertools.product(range(5), repeat=m):
            first = {c: col.index(c) if c in col else None for c in range(5)}
            want = all(first[b] is None or (first[a] is not None and first[a] < first[b])
                       for a, b in zip(colours, colours[1:]))
            cases.append(([gc.X(v, c) if col[v] == c else -gc.X(v, c) for v in range(m) for c in range(5)], want))
        check(cnf, cases)


def test_part_a_formulas_are_the_certified_ones():
    logged = vg.logged(vg.LOG_A)
    assert sorted(logged) == sorted(f"E37_A{c}.cnf" for c in vg.CIRCLES)
    for c in vg.CIRCLES:
        assert hashlib.sha256(vg.formula_a(c).encode()).hexdigest() == logged[f"E37_A{c}.cnf"], c


def test_every_leaf_of_part_b_is_certified():
    cubes = vg.read_cubes(vg.CUBES)
    assert len(cubes) == 4822 and vg.covers(cubes)
    logged = vg.logged(vg.LOG_B, vg.gzip.open)
    assert sorted(logged) == sorted(f"E37_B_leaf{i}" for i in range(len(cubes)))
    digest = vg.LeafDigests(vg.formula_b().split("\n", 1))
    assert all(digest(cube) == logged[f"E37_B_leaf{i}"] for i, cube in enumerate(cubes))


def test_the_cover_proof_resolves_to_the_empty_clause():
    """the LRAT proof of the cover formula that verify_g13.py hands to cake_lpr: each step is the resolvent of its
    two hints on the variable the tree branches on, and the last step is the empty clause"""
    cubes = vg.read_cubes(vg.CUBES)
    cnf, lrat, _ = vg.cover_proof(cubes)
    clauses = {k + 1: frozenset(map(int, line.split()[:-1])) for k, line in enumerate(cnf.splitlines()[1:])}
    last = None
    for line in lrat.splitlines():
        head, hints = line.split(" 0 ", 1)
        step, *lits = map(int, head.split())
        a, b = map(int, hints.split()[:2])
        pivot = [l for l in clauses[a] if -l in clauses[b]]
        assert len(pivot) == 1
        assert (clauses[a] | clauses[b]) - {pivot[0], -pivot[0]} == frozenset(lits)
        clauses[step], last = frozenset(lits), step
    assert clauses[last] == frozenset()
