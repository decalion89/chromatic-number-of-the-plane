"""The true minimum rho, by budget, not by greedy deletion.

Greedy deletion returns a MINIMAL forcing set -- one no single vertex can
leave -- which need not be the smallest.  The budget decision does return the
minimum: sweep the budget down, and the last value that still admits a forcing
set is rho.

This corrects two numbers recorded earlier from greedy deletion: rho(Sa,4) = 7
and rho(Sa u rot(Sa), 4) = 5.  Both are minimal-set sizes, and the budget run
already found a forcing SIX on Sa, so the first is wrong as a minimum.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
sys.path.insert(0, "/tmp/hn")
from is63 import decide

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))

pts = list(degrey.build_Sa())
for label in ("Sa", "Sa u rot(Sa)"):
    if label != "Sa":
        seen = set(pts)
        for p in list(pts):
            q = ROT(p)
            if q not in seen:
                seen.add(q)
                pts.append(q)
    g = build_graph(pts)
    print(f"\n{label}: {g.n} vertices (floor is rho >= 4)", flush=True)
    for budget in (7, 6, 5, 4, 3):
        ok, info = decide(g, 4, budget, report=10 ** 9)
        print(f"  budget {budget}: forcing set exists = {ok}", flush=True)
        if ok is False:
            print(f"  => rho({label}, 4) = {budget + 1}", flush=True)
            break
