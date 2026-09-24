"""The unit-distance graph on Q(sqrt-3, sqrt-11, sqrt-247) is 5-colourable by reduction at 11.

The place of K+ above 11 does not split in K, every unit vector reduces into the norm-one group N1
(12 elements) of F_121, and a 5-colouring of the finite plane Cay(F_121, N1) pulls back to the
whole field (hn/adelic.py, notes/rigidity.md section 9)."""
import json
from fractions import Fraction as Fr

from hn.adelic import finite_plane_11_colouring, reduce11
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"


def test_the_finite_plane_over_F11_is_5_but_not_4_colourable():
    from pysat.solvers import Solver
    col, N1 = finite_plane_11_colouring()
    assert len(N1) == 12 and len(col) == 121
    for (a, b), c in col.items():
        for (x, y) in N1:
            assert col[((a + x) % 11, (b + y) % 11)] != c
    V = sorted(col)
    idx = {v: i for i, v in enumerate(V)}
    with Solver(name="cadical195") as s:
        var = lambda v, c: 4 * v + c + 1
        for i in range(121):
            s.add_clause([var(i, c) for c in range(4)])
        for (a, b) in V:
            for (x, y) in N1:
                i, j = idx[(a, b)], idx[((a + x) % 11, (b + y) % 11)]
                if i < j:
                    for c in range(4):
                        s.add_clause([-var(i, c), -var(j, c)])
        assert not s.solve()


def test_the_803_graph_is_coloured_by_reduction_at_both_places_above_11():
    d = json.load(open(f"{ROOT}/data/five_247_c.json"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P)
    col, N1 = finite_plane_11_colouring()
    N1 = set(N1)
    for place in (1, -1):
        r = [reduce11(p, place) for p in P]
        for i, j in g.edges():
            assert ((r[j][0] - r[i][0]) % 11, (r[j][1] - r[i][1]) % 11) in N1
            assert col[r[i]] != col[r[j]]
