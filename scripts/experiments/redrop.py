"""Re-measure the 7 -> 5 drop with minimums instead of minimal sets.

The drop was reported from greedy deletion, which returns a MINIMAL forcing
set -- one no single vertex can leave -- not a smallest one.  The budget
decision has since found a forcing six on Sa, then a forcing five.  If Sa
alone is already 5 and the union is also 5, there is no drop at all and the
headline goes with it.

So both graphs, both by decision, budget swept down to the floor rho >= k.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from rhotight import run

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
    print(f"\n{label}: {g.n} vertices, floor rho >= 4", flush=True)
    value = None
    for budget in (6, 5, 4):
        v = run(g, 4, budget, rounds=200000, report=10 ** 9)
        print(f"  budget {budget}: forcing set exists = {v}", flush=True)
        if v is True:
            value = budget
        elif v is False:
            print(f"  => rho({label}, 4) = {value}", flush=True)
            break
