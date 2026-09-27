"""The colourings behind the upper bounds for the small finite planes (notes/local_colourings.md, sections 4, 8
and 12), stored in data/small_plane_colourings.json and checked here on every edge of the rebuilt graphs:
- G_q for q = 3, 5, 7, 13: the anisotropic plane over F_q, x^2 - n y^2 with n the least non-residue;
- G_q for q = 4, 8: Cay(F_{4^f}, mu_{2^f + 1}), the plane at an inert place over 2 with residue field F_{2^f};
- H_q for q <= 16, q != 13: the hyperbola graphs of tests/test_split_places.py (H_13 has a linear colouring there).
The colourings were found by a SAT solver; the lower bounds that need no solver are checked too: G_3 has an odd
cycle, G_4 is the Clebsch graph with independence number 5 < 16/3, G_5 has alpha <= 7 < 25/3 by Hoffman's bound
in interval arithmetic, and G_8 contains K_4 (the subfield F_4, since mu_3 is in mu_9)."""
import itertools
import json
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from test_split_places import hyperbola_graph  # noqa: E402

DATA = json.load(open(os.path.join(HERE, "..", "data", "small_plane_colourings.json")))


def gq_prime(q):
    n = next(a for a in range(2, q) if pow(a, (q - 1) // 2, q) == q - 1)
    U = [(x, y) for x in range(q) for y in range(q) if (x * x - n * y * y) % q == 1]
    assert len(U) == q + 1
    E = {tuple(sorted((q * x + y, q * ((x + a) % q) + (y + b) % q))) for x in range(q) for y in range(q)
         for a, b in U}
    return q * q, sorted(E)


def gq_char2(f):
    m = 2 * f
    poly = {4: 0b10011, 6: 0b1000011}[m]            # x^4 + x + 1, x^6 + x + 1

    def mul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> m:
                a ^= poly
        return r

    def power(a, e):
        r = 1
        for _ in range(e):
            r = mul(r, a)
        return r

    n = 1 << m
    U = [a for a in range(1, n) if power(a, 2 ** f + 1) == 1]
    assert len(U) == 2 ** f + 1
    return n, sorted({tuple(sorted((a, a ^ u))) for a in range(n) for u in U}), power


def graph(kind, q):
    if kind == "H_q":
        V, _, E = hyperbola_graph(q)
        return len(V), E
    if q in (4, 8):
        n, E, _ = gq_char2(q.bit_length() - 1)
        return n, E
    return gq_prime(q)


# (graph, q): the number of colours of the stored colouring, which is the chromatic number
CASES = [("G_q", 3, 3), ("G_q", 4, 4), ("G_q", 5, 4), ("G_q", 7, 4), ("G_q", 8, 4), ("G_q", 13, 6)] + \
        [("H_q", q, k) for q, k in [(2, 2), (3, 3), (4, 4), (5, 3), (7, 4), (8, 4), (9, 3), (11, 4), (16, 4)]]


@pytest.mark.parametrize("kind, q, k", CASES)
def test_the_stored_colouring_is_proper(kind, q, k):
    n, E = graph(kind, q)
    col = DATA[kind][str(q)]
    assert len(col) == n and set(col) == set(range(k))
    assert all(col[a] != col[b] for a, b in E)


def test_the_stored_colourings_are_all_checked():
    assert {(kind, int(q)) for kind in ("G_q", "H_q") for q in DATA[kind]} == {(k, q) for k, q, _ in CASES}


def test_lower_bounds_without_a_solver():
    # G_3 is not bipartite: an odd cycle
    n, E = gq_prime(3)
    side, adj = {0: 0}, {v: set() for v in range(n)}
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    stack, odd = [0], False
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in side:
                side[w] = 1 - side[v]
                stack.append(w)
            elif side[w] == side[v]:
                odd = True
    assert odd
    # G_4 is the Clebsch graph: 5-regular on 16 vertices, triangle-free, independence number 5 < 16/3
    n, E, _ = gq_char2(2)
    adj = {v: set() for v in range(n)}
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    assert all(len(adj[v]) == 5 for v in adj) and not any(adj[a] & adj[b] for a, b in E)
    alpha = max(r for r in range(1, 7) for S in itertools.combinations(range(n), r)
                if all(b not in adj[a] for a, b in itertools.combinations(S, 2)))
    assert alpha == 5 and 3 * alpha < 16
    # G_5: Hoffman's bound in interval arithmetic
    from finite_hoffman import hoffman
    d, lmin, amax = hoffman(5, "inert")
    assert d == 6 and amax == 7 and 3 * amax < 25
    # G_8 contains K_4: the subfield F_4 = {a : a^4 = a} of F_64, whose differences lie in F_4^* = mu_3, in mu_9
    n, E, power = gq_char2(3)
    F4 = [a for a in range(n) if power(a, 4) == a]
    assert len(F4) == 4
    assert all(tuple(sorted(pair)) in set(E) for pair in itertools.combinations(F4, 2))
