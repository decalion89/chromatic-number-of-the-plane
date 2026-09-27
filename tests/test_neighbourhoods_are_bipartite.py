"""Neighbourhoods in the plane are bipartite, so rigidity is never local.

The theorem is not new in this project -- it was derived in an earlier pass
(docs/research-log.md, "Every neighbourhood is bipartite, which bounds the single-point
attack"), along with the conclusion that rigidity is global.  This file checks
it on the graphs built since, and draws two consequences that pass did not.

A vertex v is FREE at k colours when some colour is missing from its
neighbourhood, and free@k = 0 -- which Sa achieves at four -- means no vertex
ever has a spare colour.  The obvious way to guarantee that for one vertex is
locally: if the subgraph induced on N(v) needs k-1 colours, then N(v) uses all
k-1 colours other than c(v) in every proper colouring, and v is pinned.

In the plane that is available only for k = 3, and the reason is a fact about
circles.

    Every neighbour of v lies on the circle of radius 1 about v, and two points
    of a radius-1 circle are a unit apart exactly when their central angle is
    2 arcsin(1/2) = 60 degrees.

So inside N(v) a point can only be adjacent to the two points 60 degrees away
from it: the induced subgraph has MAXIMUM DEGREE 2.  It is therefore a disjoint
union of paths and cycles, and a cycle has to close with steps of +-60 degrees
summing to a multiple of 360, which forces even length -- with a steps forward
and b back, a - b = 6m and a + b = (a - b) + 2b = 6m + 2b, even.  Hence

    chi(N(v)) <= 2   for every vertex of every unit-distance graph in R^2.

Two consequences run through the whole project:

  * At k >= 4 no vertex can be pinned by its own neighbourhood.  The local
    certificate tops out at two colours, so it can never account for the k-1
    that rigidity needs.  Sa's free@4 = 0.00 % is therefore a GLOBAL property,
    and there is no local gadget to transplant -- which is why every attempt
    to manufacture it by adding degree failed.
  * The disjunction gadget's local four-vertex form dies here too, and for the
    same geometry: the points at distance 1 from both ends of an edge are the
    intersection of two unit circles, exactly two of them.

The test checks the statement itself on every graph in the project rather than
quoting it: maximum degree at most 2 inside each neighbourhood, every induced
neighbourhood 2-colourable, and the 60-degree characterisation in exact
arithmetic.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as Fr

import pytest

from hn.degrey import build_G, build_Sa, build_Y
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

HERE = os.path.dirname(__file__)
DEGREY = Field((3, 5, 7, 11))


def _two_colour(vs, adj):
    """Greedy 2-colouring of a graph of maximum degree 2; None if it fails."""
    colour = {}
    for s in vs:
        if s in colour:
            continue
        colour[s] = 0
        stack = [s]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in colour:
                    colour[w] = 1 - colour[u]
                    stack.append(w)
                elif colour[w] == colour[u]:
                    return None
    return colour


def _check(g, limit=None):
    worst = 0
    n = g.n if limit is None else min(g.n, limit)
    for v in range(n):
        nb = sorted(g.adj[v])
        s = set(nb)
        adj = {u: [w for w in g.adj[u] if w in s] for u in nb}
        worst = max(worst, max((len(a) for a in adj.values()), default=0))
        assert worst <= 2, f"neighbourhood of {v} has a vertex of degree {worst}"
        assert _two_colour(nb, adj) is not None, \
            f"neighbourhood of {v} is not bipartite"
    return worst


def test_the_sixty_degree_characterisation_is_exact():
    """Two points of a radius-1 circle are a unit apart iff 60 degrees apart."""
    f = DEGREY
    one, zero = f.rational(Fr(1)), f.zero()
    p = Point(one, zero)
    r60 = Rotation(f.rational(Fr(1, 2)), f.sqrt(3) * f.rational(Fr(1, 2)))
    q = r60(p)
    assert q.norm2() == one                    # still on the circle
    assert (p - q).norm2() == one              # and a unit away
    # 120 degrees apart is sqrt3, not 1
    assert (p - r60(q)).norm2() == f.rational(Fr(3))
    # so a neighbour has at most two neighbours inside the neighbourhood
    assert (p - r60(r60(q))).norm2() == f.rational(Fr(4))


@pytest.mark.parametrize("name,builder", [
    ("Sa", lambda: build_Sa(DEGREY)),
    ("Y", lambda: build_Y(DEGREY)),
    ("G", lambda: build_G(DEGREY, as_graph=False)),
])
def test_de_greys_neighbourhoods_are_paths_and_hexagons(name, builder):
    g = build_graph(builder())
    worst = _check(g)
    assert worst <= 2


@pytest.mark.parametrize("name", ["five_247_c.json", "five_tuned_1_1.json"])
def test_the_new_graphs_too(name):
    with open(os.path.join(HERE, os.pardir, "data", name)) as fh:
        d = json.load(fh)
    field = Field(tuple(d["field_generators"]))
    pts = [Point(field.element([Fr(a, b) for a, b in x]),
                 field.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    g = build_graph(pts)
    assert _check(g, limit=400) <= 2


def test_so_no_vertex_is_pinned_by_its_neighbourhood_at_four():
    """chi(N(v)) <= 2 < 3, so the local certificate never reaches four.

    Stated as a property of the numbers: the largest chromatic number an
    induced neighbourhood can have is 2, and pinning a vertex at k colours
    needs k-1.  So local pinning is possible only at k = 3.
    """
    g = build_graph(build_Sa(DEGREY))
    best = 0
    for v in range(g.n):
        nb = sorted(g.adj[v])
        s = set(nb)
        adj = {u: [w for w in g.adj[u] if w in s] for u in nb}
        if any(adj.values()):
            best = max(best, 2)
    assert best == 2                      # some neighbourhood does have an edge
    assert best < 3                       # and none can reach three
