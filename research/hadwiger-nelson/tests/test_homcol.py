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
