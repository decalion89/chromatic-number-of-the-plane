#!/usr/bin/env python3
"""The two-orbit block on a pair that is forced only jointly.

Every forced pair this package found by search turned out to be forced one leg
at a time: the 409-vertex certificate blocks a pair whose sqrt(3) leg is
already forced alone, so a core of 1, and two copies of the ordinary spindle
would close it too. The block was valid and unnecessary. What was missing was
an instance where neither leg is forced by itself and the pair is forced only
together -- the case the classical argument cannot reach, because it has no
single target to spindle.

Here is one, built rather than found:

    x = (0, 0)      y = (1, 0)      z = (1/2, sqrt(3)/2)      unit triangle
    p = (5/6, -sqrt(11)/6)

    |p - x|^2 = 1          p is adjacent to x
    |p - y|^2 = 1/3
    |p - z|^2 = (7 + sqrt(33)) / 6

In any 3-colouring x, y, z take three different colours and p differs from x,
so p repeats y's colour or z's. Neither alone: both completions exist. The
core is exactly two, and it is mixed -- 1/3 on one leg, (7 + sqrt(33))/6 on
the other -- which is what lets it past the capacity ceiling that caps every
same-distance core at 2 and every same-distance block at 3 colours.

The block needs a rotation closing each circle. On the 1/3 circle that is 120
degrees, Niven's only odd-order radius. On the other,

    cos t = (-5 + 3 sqrt(33)) / 16,  sin^2 t = (-66 + 30 sqrt(33)) / 256,

and sin t lives in no multiquadratic field: that sin^2 has a negative
conjugate, and multiquadratic fields are totally real. Adjoining it once,
over Q(sqrt 3, sqrt 11), is what `hn.realext` is for.

Six copies -- the three rotations by 120 degrees and the same three composed
with t -- and the union has no 3-colouring. Nineteen points, thirty-three
edges. Not a record for chi >= 4; the Moser spindle does it in seven. It is
the mechanism that is new, and this is the smallest k where it can be shown to
work at all.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import exact_edges, save_certificate, verify_certificate
from hn.coloring import is_k_colorable
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from hn.mixed import (blocks_two_targets, compose_rotations,
                      joint_core_configuration, joint_core_copies,
                      joint_core_union)
from hn.spindle import SeparationTest

DRAT = os.environ.get("HN_DRAT") or None
OUT = os.environ.get(
    "HN_OUT", os.path.join(HN_DIR, "certificates",
                           "genuine_pair_19_no3coloring.json"))


def main() -> None:
    E, (p, x, y, z), rho, sigma = joint_core_configuration()
    one, third = E.rational(1), E.rational(Fraction(1, 3))
    print("the four points")
    print(f"  |p-x|^2 = {p.dist2(x)}   adjacent: {p.dist2(x) == one}")
    print(f"  |p-y|^2 = {p.dist2(y)}   = 1/3: {p.dist2(y) == third}")
    print(f"  |p-z|^2 = {p.dist2(z)}")
    print(f"  x,y,z unit triangle: "
          f"{x.dist2(y) == one and y.dist2(z) == one and z.dist2(x) == one}")

    base4 = build_graph([p, x, y, z])
    bp = base4.vertices.index(p)
    ia, ib = base4.vertices.index(y), base4.vertices.index(z)
    st = SeparationTest(base4, 3, bp, [ia, ib])
    try:
        print("\nthe core, by SAT on the four points")
        print(f"  y alone separable : {st.run(subset=[ia])[0]}   "
              f"(True = not forced on its own)")
        print(f"  z alone separable : {st.run(subset=[ib])[0]}")
        print(f"  pair forced       : {not st.run(subset=[ia, ib])[0]}")
    finally:
        st.close()

    ident = Rotation(E.rational(1), E.zero())
    copies = joint_core_copies(E, rho, sigma)
    print(f"\nthe six copies")
    print(f"  rho^3 = 1            : "
          f"{compose_rotations(compose_rotations(rho, rho), rho) == ident}")
    print(f"  sigma is a rotation  : "
          f"{sigma.cos * sigma.cos + sigma.sin * sigma.sin == one}")
    print(f"  sigma closes the big circle: "
          f"{z.dist2(sigma.about(p)(z)) == one}")
    print(f"  rho closes the 1/3 circle  : "
          f"{y.dist2(rho.about(p)(y)) == one}")
    print(f"  blocks (2-SAT)       : "
          f"{blocks_two_targets(base4, bp, [ia, ib], copies)}")

    pts = joint_core_union()
    g = build_graph(pts)
    edges = exact_edges(pts)
    ok3, _ = is_k_colorable(g, 3)
    ok4, col4 = is_k_colorable(g, 4)
    print(f"\nthe union: {g}")
    print(f"  edges re-derived exactly : {len(edges)}  "
          f"agrees with the hash: {set(edges) == set(g.edges())}")
    print(f"  3-colourable: {ok3}    4-colourable: {ok4}")

    crit = list(range(len(pts)))
    for i in sorted(crit, key=lambda i: len(g.adj[i])):
        trial = [j for j in crit if j != i]
        if not is_k_colorable(build_graph([pts[j] for j in trial]), 3)[0]:
            crit = trial
    print(f"  a 3-critical subgraph: {len(crit)} vertices")

    notes = {
        "mechanism": "two-orbit block on a genuinely joint core of two",
        "pivot": "p = (5/6, -sqrt(11)/6)",
        "legs": ["|p-y|^2 = 1/3", "|p-z|^2 = (7 + sqrt(33))/6"],
        "why_it_is_new": (
            "neither leg is forced on its own -- both singleton separations "
            "are satisfiable -- so no single spindle closes this pair, and "
            "the classical argument has no target. The pair is forced only "
            "jointly, and six copies of it around the pivot exhaust the "
            "colours. The earlier 409-vertex certificate blocked a pair "
            "whose 1/3 leg was already forced alone, a core of one, so it "
            "did not need the block; this one does."),
        "field": ("Q(sqrt 3, sqrt 11)(sqrt v), v = (-66 + 30 sqrt 33)/256. "
                  "The extension is unavoidable: v has a negative conjugate, "
                  "so sqrt(v) lies in no totally real field, and every "
                  "multiquadratic field is totally real."),
        "rotations": ["cos = -1/2, sin = sqrt(3)/2 (order 3)",
                      "cos = (-5 + 3 sqrt 33)/16, sin = sqrt(v)"],
        "critical_subgraph_order": len(crit),
        "not_a_record": ("the Moser spindle reaches chi >= 4 on seven "
                         "vertices; what is new here is the mechanism, not "
                         "the size"),
    }
    path = os.path.abspath(OUT)
    save_certificate(g, path, 3, "not 3-colourable, by the two-orbit block on "
                     "a pair forced only jointly", notes=notes)
    print(f"\nwrote {path}")
    ok, msg = verify_certificate(path, drat_trim=DRAT)
    print(f"  {'VERIFIED' if ok else 'FAILED'}: {msg}")


if __name__ == "__main__":
    main()
