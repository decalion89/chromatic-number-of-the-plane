#!/usr/bin/env python3
"""Construct the forced pair instead of searching for it.

A pair {a,b} is forced exactly when the graph plus the edges (p,a),(p,b) has
no k-colouring, so the smallest instance at k=3 is K_4 with two edges gone at
one vertex: a unit triangle x, y, z and a pivot p at distance 1 from x and
1/sqrt(3) from y. The triangle takes all three colours, p differs from x, so p
shares with y or z. Forced, with neither leg forced alone, and y on the
classical circle.

That configuration is rigid, and its rigidity is the problem: |p - x| = 1 with
|p - y|^2 = 1/3 gives |p - z|^2 = (7 +- sqrt(33))/6, and the rotation closing
the second leg needs the square root of an element whose Galois conjugate is
negative -- outside every real multiquadratic field.

The fix is to stop hanging p off the triangle. In a connected patch of the
triangular lattice a 3-colouring is unique up to permutation, so every lattice
point of x's class carries x's colour in every colouring. Hanging p off any of
those gives the same forcing and frees its position: p sits on the unit circle
about w and the 1/sqrt(3) circle about y, and |p - z|^2 then varies with w.

So: sweep w over x's colour class, and keep the p whose second leg the field
can actually close.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.field import Field
from hn.generate import hex_ball
from hn.geometry import Point
from hn.graph import build_graph
from hn.mixed import (blocks_two_targets, circle_intersections,
                      joining_rotation_exists, two_orbit_block)
from hn.multispindle import cross_blocks
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "3"))
RADIUS = int(os.environ.get("HN_RADIUS", "5"))
# The field decides how many placements can be named at all. Over
# Q(sqrt3, sqrt11) only four of thirty survive; each generator added is
# another radical the intersections may use.
GENS = tuple(int(x) for x in os.environ.get("HN_GENS", "3,5,7,11").split(","))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/hn")
THIRD = Fraction(1, 3)


def main() -> None:
    pts = list(hex_ball(RADIUS, Field(GENS)))
    g = build_graph(pts)
    field = g.vertices[0].x.field
    print(f"hex ball r={RADIUS}: {g}  over {field}", flush=True)
    one = field.rational(1)
    third = field.rational(THIRD)

    # a unit triangle, and x's colour class: the lattice points whose squared
    # distance from x is divisible by 3, which is what shares its colour
    x = g.vertices[0]
    tri = None
    for i in range(g.n):
        for j in g.adj[i]:
            for k2 in g.adj[i]:
                if k2 > j and k2 in g.adj[j]:
                    tri = (g.vertices[i], g.vertices[j], g.vertices[k2])
                    break
            if tri:
                break
        if tri:
            break
    x, y, z = tri
    same = [q for q in g.vertices
            if q.dist2(x).is_rational() and q.dist2(x).c[0].denominator == 1
            and q.dist2(x).c[0] % 3 == 0]
    print(f"  triangle found; {len(same)} points in x's colour class",
          flush=True)

    t0, tried, ok = time.time(), 0, 0
    for w in same:
        for p in circle_intersections(w, one, y, third):
            tried += 1
            if p in set(g.vertices):
                continue
            d2z = p.dist2(z)
            if not joining_rotation_exists(field, d2z):
                continue
            ok += 1
            g2 = build_graph(list(g.vertices) + [p])
            if not is_k_colorable(g2, K)[0]:
                continue
            idx = {q: i for i, q in enumerate(g2.vertices)}
            bp, a, b = idx[p], idx[y], idx[z]
            if b in g2.adj[bp] or a in g2.adj[bp]:
                continue
            st = SeparationTest(g2, K, bp, [a, b])
            try:
                if not st.run(subset=[a])[0] or not st.run(subset=[b])[0]:
                    continue                      # a leg is forced alone
                if st.run(subset=[a, b])[0]:
                    continue                      # the pair is not forced
            finally:
                st.close()
            print(f"  GENUINE pair: p at d^2 1/3 from y, {d2z} from z"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
            copies = two_orbit_block(g2, bp, a, b)
            if copies is None:
                print("    but the block cannot be named here", flush=True)
                continue
            assert blocks_two_targets(g2, bp, [a, b], copies)
            assert cross_blocks(g2, bp, [a, b], copies)
            union = {q: None for q in g2.vertices}
            for r in copies:
                f = r.about(p)
                for q in g2.vertices:
                    union[f(q)] = None
            u = build_graph(list(union))
            good = is_k_colorable(u, K)[0]
            print(f"    union: {u}  {K}-colourable: {good}"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
            if not good:
                save_certificate(u, f"{OUT}/genuine_core2_no{K}col.json", k=K,
                                 claim=f"not {K}-colourable by the two-orbit "
                                 f"block on a genuine core of two")
                print(f"  *** chi >= {K + 1}, genuine core of two ***",
                      flush=True)
                return
    print(f"  {tried} placements named, {ok} with a closeable second leg"
          f"  [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
