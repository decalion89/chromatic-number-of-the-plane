"""Push the calibrated search hard at the one instance sitting at five.

TabuCol, calibrated in both directions, plateaus at 21 conflicts on Sa with
all four carrying classes forbidden at five colours, and at FIVE on Sa with
just {4/9, 16/9}.  For reference it plateaus at 37 on an instance known to be
unsatisfiable, and reaches zero in two or three seconds on ones known to be
satisfiable.

Five conflicts is not a verdict.  A barely-unsatisfiable instance plateaus low
and so does a hard satisfiable one that the search nearly solves, and the two
look identical from outside.  What separates them is persistence: if any run
ever reaches zero the instance is satisfiable and the matter is closed, while
a hundred runs that do not are evidence and stay evidence.

So: many seeds, long runs, on that instance alone, and on the four-class one
behind it.  Report the best reached and how often, not an opinion.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from tabucol import instance, tabucol

t0 = time.time()
Sa = build_Sa(F)
JOBS = [("Sa, {4/9,16/9}, k=5", {Fr(4, 9), Fr(16, 9)}, 5),
        ("Sa, all four, k=5", {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}, 5)]
for label, cs, k in JOBS:
    n, E = instance(Sa, cs)
    got = Counter()
    hit = False
    for sd in range(60):
        b, it = tabucol(n, E, k, seed=1000 + sd, iters=600_000)
        got[b] += 1
        if b == 0:
            print(f"*** {label}: SATISFIABLE -- seed {1000+sd} reached zero "
                  f"after {it} moves  [{time.time()-t0:.0f}s]", flush=True)
            hit = True
            break
        if sd % 10 == 9:
            print(f"   {label}: {sd+1} runs, best {min(got)}, "
                  f"distribution {dict(sorted(got.items()))}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
    if not hit:
        print(f"{label}: 60 runs, none reached zero, best {min(got)}, "
              f"distribution {dict(sorted(got.items()))}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
