"""Homomorphism colourings: the structural screen, and what it says."""
from hn.degrey import build_G
from hn.graph import build_graph
from hn.homcol import edge_vectors, has_homomorphism, screen
from hn.mixed import joint_core_union, three_hexagon_gadget


def _verify(phi, vecs, n):
    """Check phi directly rather than trusting the solver."""
    return all(sum(p * x for p, x in zip(phi, d)) % n for d in vecs)


def test_edge_vectors_include_a_real_extension_s_radical_part():
    """`.c` on a RealExtElement gives only its base component.

    Using it alone silently truncated the edge vectors of anything living over
    sqrt(v) and ran the homomorphism search on the wrong module. The joint-core
    union lives over Q(sqrt3, sqrt11)(sqrt v), so its coordinates need both
    halves: dimension 16, not 8.
    """
    g = build_graph(joint_core_union())
    vecs = edge_vectors(g)
    assert len(vecs[0]) == 16


def test_common_content_is_divided_out():
    """Clearing denominators multiplies every vector by one integer, and if
    that integer shares a factor with n the whole set looks divisible by n --
    a scaling artefact reported as an obstruction until the content was
    removed. After removing it, the vectors have no common factor."""
    from math import gcd

    for pts in (three_hexagon_gadget()[2], joint_core_union()):
        vecs = edge_vectors(build_graph(pts))
        content = 0
        for v in vecs:
            for x in v:
                content = gcd(content, abs(x))
        assert content == 1


def test_everything_here_is_five_colourable_by_cosets():
    """A homomorphism from the edge-vector module to Z/5 that avoids every edge
    vector IS a proper 5-colouring, c(p) = phi(p), constant on cosets of its
    kernel. It exists for all three graphs, which is why every rigidity
    measurement in this package came back flat -- a coset colouring composes
    with the automorphisms of Z/5 and with any translation, so colour classes
    move freely and forcing sets are enormous.

    It is also a screen: such a graph can never be 6-chromatic.
    """
    for name, g in (("gadget", build_graph(three_hexagon_gadget()[2])),
                    ("joint core", build_graph(joint_core_union())),
                    ("de Grey G", build_G())):
        vecs = edge_vectors(g)
        phi, why = has_homomorphism(vecs, 5)
        assert phi is not None, f"{name}: {why}"
        assert _verify(phi, vecs, 5), f"{name}: phi does not avoid every vector"
        assert "never 6-chromatic" in screen(g, 5)["verdict"]


def test_rank_two_lattices_can_never_block_a_coset_colouring():
    """A plane lattice's unit vectors occupy exactly three of the six
    projective points mod 5, and which three is decided by the square class of
    their norm. Three hyperplanes cover 13 of the 25 points of (Z/5)^2, so a
    homomorphism always survives -- blocking needs rank at least 4.
    """
    def proj(v, p=5):
        v = (v[0] % p, v[1] % p)
        if v == (0, 0):
            return None
        for lam in range(1, p):
            w = ((v[0] * lam) % p, (v[1] * lam) % p)
            if w[0] == 1 or (w[0] == 0 and w[1] == 1):
                return w
        return v

    def norm(a, b):
        return a * a - a * b + b * b

    squares, nonsquares = set(), set()
    for N in (1, 3, 7, 13, 21, 49, 91):
        vecs = [(a, b) for a in range(-12, 13) for b in range(-12, 13)
                if norm(a, b) == N]
        assert vecs
        cls = {proj(v) for v in vecs} - {None}
        assert len(cls) == 3, f"norm {N} gave {len(cls)} classes"
        (squares if N % 5 in (1, 4) else nonsquares).add(frozenset(cls))
    assert len(squares) == 1 and len(nonsquares) == 1
    assert squares != nonsquares
    assert set(next(iter(squares))) == {(0, 1), (1, 0), (1, 1)}
    assert set(next(iter(nonsquares))) == {(1, 2), (1, 3), (1, 4)}


def test_unit_directions_can_never_contain_a_projective_line():
    """The cheapest blocking set in PG(r-1,5) is a line, six points, and unit
    vectors cannot supply one.

    They all satisfy Q(d) = N for a single N, so their directions carry a form
    value in one square class; and no projective line over F_5 has all six of
    its points in one square class. Checked exhaustively over every binary
    form and every target: the maximum is five, reached only when the form
    degenerates to rank one, the sixth point being its radical.
    """
    p = 5
    squares = {(t * t) % p for t in range(1, p)}

    def direction(v):
        for lam in range(1, p):
            w = tuple((t * lam) % p for t in v)
            for i in range(len(w)):
                if w[i] == 1 and all(w[j] == 0 for j in range(i)):
                    return w
        return None

    best = 0
    for a in range(p):
        for b in range(p):
            for c in range(p):
                pts = {}
                for x in range(p):
                    for y in range(p):
                        if (x, y) == (0, 0):
                            continue
                        d = direction((x, y))
                        if d is None:
                            continue
                        pts.setdefault(d, set()).add(
                            (a * x * x + b * x * y + c * y * y) % p)
                for N in range(1, p):
                    cls = {(N * s) % p for s in squares}
                    best = max(best, sum(1 for vals in pts.values()
                                         if vals & cls))
    assert best == 5, f"a line should never be fully covered; got {best}"


def test_cyclotomic_step_set_admits_no_coset_five_colouring():
    """Q(zeta_7) blocks, where every multiquadratic family fails.

    Multiquadratic fields give at most 3^t directions against PG(2t-1,5)'s
    (5^2t - 1)/4 points, because 5 splits there and the norm-one group mod 5
    is a product of factors of size 2 or 3. A cyclic Galois group with 5
    inert avoids that: in Q(zeta_7), 5 is a primitive root mod 7, O/5 is
    F_5^6, and the norm-one subgroup has 126 elements.

    Modulus-one elements come free from Hilbert 90 -- u = alpha/conj(alpha) --
    and 63 of their directions admit no homomorphism to Z/5 at all. Checked by
    brute force over all 15625 maps rather than on the solver's word.
    """
    import itertools
    import random

    from hn.homcol import CYCLOTOMIC_BLOCKING

    D, N = 6, 7

    def mul(x, y):
        raw = [0] * (2 * D + 1)
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    if b:
                        raw[i + j] += a * b
        for k in range(2 * D, D - 1, -1):
            c = raw[k]
            if c:
                raw[k] = 0
                for j in range(D):
                    raw[k - D + j] -= c
        return tuple(raw[:D])

    def conj(x):
        out = [0] * D
        for i, a in enumerate(x):
            if not a:
                continue
            k = (-i) % N
            if k < D:
                out[k] += a
            else:
                for j in range(D):
                    out[j] -= a
        return tuple(out)

    def inv_mod5(x):
        e, r, b = 5 ** D - 2, tuple([1] + [0] * (D - 1)), tuple(t % 5 for t in x)
        while e:
            if e & 1:
                r = tuple(t % 5 for t in mul(r, b))
            b = tuple(t % 5 for t in mul(b, b))
            e >>= 1
        return r

    def direction(v):
        v = tuple(t % 5 for t in v)
        if not any(v):
            return None
        for lam in range(1, 5):
            w = tuple((t * lam) % 5 for t in v)
            for i in range(D):
                if w[i] == 1 and all(w[j] == 0 for j in range(i)):
                    return w
        return None

    random.seed(7)
    seen, vecs = set(), []
    for _ in range(4000):
        a = tuple(random.randint(-3, 3) for _ in range(D))
        if not any(a):
            continue
        ca = conj(a)
        if all(t % 5 == 0 for t in ca):
            continue
        u = tuple(t % 5 for t in mul(a, inv_mod5(ca)))
        d = direction(u)
        if d and d not in seen:
            seen.add(d)
            vecs.append(u)

    assert len(vecs) == CYCLOTOMIC_BLOCKING["directions"] == 63
    assert not any(all(t % 5 == 0 for t in v) for v in vecs)
    survivors = sum(1 for f in itertools.product(range(5), repeat=D)
                    if all(sum(p * t for p, t in zip(f, v)) % 5 for v in vecs))
    assert survivors == 0, f"{survivors} coset colourings survive"
