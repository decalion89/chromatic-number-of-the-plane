"""The same push, at G's class 4/9 -- the instance CDCL has not answered.

Four solvers have been on "G at five colours with all 1558 pairs at squared
distance 4/9 forbidden" for over three hours without a verdict.  The
calibrated TabuCol plateaus at seven conflicts out of 9435 constraints, which
is closer to the line than any other instance here except Sa with
{4/9, 16/9} at five.

Seven is not a verdict.  If a run ever reaches zero the instance is
satisfiable and three hours of CDCL become moot; if a hundred runs do not,
that is evidence and stays evidence.  Either way it costs minutes rather than
hours, which is the point.

Balls too, smallest first: a ball that fails is a witness, and every ball that
colours is another instance the search has actually solved rather than
plateaued on.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point
from tabucol import instance, tabucol

t0 = time.time()
P = build_G(F, as_graph=False)
n, E = instance(P, {Fr(4, 9)})
print(f"G at five, class 4/9: {n} pts, {len(E)} constraints"
      f"  [{time.time()-t0:.0f}s]", flush=True)
got = Counter()
for sd in range(80):
    b, it = tabucol(n, E, 5, seed=2000 + sd, iters=800_000)
    got[b] += 1
    if b == 0:
        print(f"*** SATISFIABLE -- seed {2000+sd} reached zero after {it} "
              f"moves.  G does NOT carry the weak property on 4/9 at five, "
              f"and the three-hour CDCL run is moot  [{time.time()-t0:.0f}s]",
              flush=True)
        break
    if sd % 10 == 9:
        print(f"   {sd+1} runs, best {min(got)}, "
              f"distribution {dict(sorted(got.items()))}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
else:
    print(f"80 runs, none reached zero, best {min(got)}, "
          f"distribution {dict(sorted(got.items()))}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
