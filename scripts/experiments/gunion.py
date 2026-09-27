"""G unioned with rotated copies: the cheap half of the Sa experiment.

Unioning Sa with one rotated copy dropped rho from 7 to 5.  Measuring rho on
G costs a hard UNSAT, but PRESSURE does not -- it is a small optimisation over
the pivot's neighbourhood, and it is the quantity the whole core argument
needs to exceed 2.

So: G with copies rotated through the half-Moser angle and its powers, and the
pressure at the best pivot of each.  If unioning does for G at five colours
what it did for Sa at four, this is where it shows.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure

FLD = degrey.DEGREY_FIELD
HALF = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
                FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))
FULL = Rotation(FLD.rational(Fraction(5, 6)),
                FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))

pts = list(degrey.build_G(as_graph=False))
print(f"G: {len(pts)} points", flush=True)
cur = list(pts)
for label, rot in (("alone", None), ("+half", HALF), ("+full", FULL),
                   ("+half^3", lambda p: HALF(HALF(HALF(p))))):
    if rot is not None:
        seen = set(cur)
        for p in list(cur):
            q = rot(p)
            if q not in seen:
                seen.add(q)
                cur.append(q)
    g = build_graph(cur)
    deg = g.degrees()
    hub = max(range(g.n), key=lambda v: deg[v])
    t = time.time()
    rel = ColourRelations(g, 5)
    pr = pressure(rel, hub)
    # and the second-best pivot, in case the hub is not the rigid one
    others = sorted(range(g.n), key=lambda v: -deg[v])[1:4]
    prs = [pressure(rel, v) for v in others]
    print(f"  {label:9}: n={g.n:5} m={g.m:6}  hub degree {deg[hub]:3}, "
          f"pressure {pr}, next three {prs}"
          + ("   *** ABOVE 2 ***" if max([pr] + prs) > 2 else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
