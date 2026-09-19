"""Gluing copies along an interface: the transfer-matrix construction."""
from hn.geometry import SPINDLE, eisenstein, origin
from hn.graph import build_graph
from hn.transfer import (chain_length, congruent_pairs, realisable_patterns,
                         transfer_relation)


def moser():
    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    return build_graph(rh + [SPINDLE(p) for p in rh])


def test_realisable_patterns_respect_the_edges():
    """No realisable pattern gives two adjacent interface points one colour."""
    g = moser()
    W = [0, 1, 2]
    pats = realisable_patterns(g, 4, W)
    assert 0 < len(pats) < 4 ** 3
    for pat in pats:
        for i, u in enumerate(W):
            for j, v in enumerate(W):
                if j > i and v in g.adj[u]:
                    assert pat[i] != pat[j]


def test_no_pattern_is_realisable_without_a_colouring():
    """The spindle is 4-chromatic, so at k=3 the interface has nothing."""
    assert realisable_patterns(moser(), 3, [0, 1, 2]) == set()


def test_congruent_sets_really_are_congruent():
    g = moser()
    W = [0, 1, 2]
    for Wp in congruent_pairs(g, W, limit=8):
        for i, a in enumerate(W):
            for j, b in enumerate(W):
                assert (g.vertices[a].dist2(g.vertices[b])
                        == g.vertices[Wp[i]].dist2(g.vertices[Wp[j]]))


def test_a_colourable_graph_never_empties_its_chain():
    """The spindle 4-colours, so no chain of copies of it can fail that way."""
    g = moser()
    W = [0, 1, 2]
    for Wp in congruent_pairs(g, W, limit=4):
        rel = transfer_relation(g, 4, W, Wp)
        assert rel
        assert chain_length(rel) is None


def test_an_empty_relation_empties_at_once():
    assert chain_length({}) == 1
