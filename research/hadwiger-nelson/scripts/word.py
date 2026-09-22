"""Spindle words: build a graph as a chain of spindlings and measure it.

A spindle letter is a pair (centre, r^2).  Applying it to a point set H means
H u rho(H) where rho is the rotation about `centre` by 2*arcsin(1/(2r)), which
puts every point of H on the radius-r circle about `centre` at distance exactly
1 from its own image.

The cost of a letter is the radical sqrt((4r^2-1)) needs; the cost of a word is
the compositum.  de Grey's word is [(v_in_S, 3), (origin, 4), ((-2,0), 16)] and
costs Q(sqrt3,sqrt5,sqrt7,sqrt11).  This script runs single letters that cost
nothing at all on top of Q(sqrt3,sqrt11), which his construction never used.
"""
import sys, time
from fractions import Fraction as Fr

sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, required_radical
from hn.field import Field
from hn.graph import build_graph
from hn.coloring import is_k_colorable

F = Field((3, 5, 7, 11))          # ambient; the letters below stay in Q(3,11)
ORI = Point(F.zero(), F.zero())

def spindle(pts, centre, d2):
    rot = rotation_joining(d2, F).about(centre)
    seen, out = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen:
                seen.add(q); out.append(q)
    return out

def hinge(pts, centre, d2):
    return sum(1 for p in pts if (p - centre).norm2() == d2)

LETTERS = [
    ("origin", ORI, Fr(3)),      # the Moser rotation, free (sqrt11)
    ("origin", ORI, Fr(5, 9)),   # free (sqrt11), hinge 12 -- never used
    ("origin", ORI, Fr(7, 3)),   # free (sqrt3),  hinge  9
    ("origin", ORI, Fr(7)),      # free (sqrt3),  hinge  4
    ("origin", ORI, Fr(13, 3)),  # free (sqrt3),  hinge  4
    ("origin", ORI, Fr(4)),      # de Grey's own: costs sqrt15
]

Sa = build_Sa(F)
print(f"hub Sa: {len(Sa)} points", flush=True)
for name, c, d2 in LETTERS:
    t0 = time.time()
    h = hinge(Sa, c, d2)
    pts = spindle(Sa, c, d2)
    g = build_graph(pts)
    t1 = time.time()
    four, _ = is_k_colorable(g, 4, timeout=600)
    t2 = time.time()
    tag = "4-COLOURABLE" if four else ("NOT 4-colourable  <<<<<" if four is False else "timeout")
    print(f"  r^2={str(d2):<6} radical {required_radical(d2):<3} hinge {h:>3} "
          f"-> n={g.n:<5} m={sum(len(a) for a in g.adj)//2:<6} "
          f"[build {t1-t0:.0f}s, sat {t2-t1:.0f}s]  {tag}", flush=True)
