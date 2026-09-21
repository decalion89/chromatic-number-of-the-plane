"""More moves, not more seeds -- which is what settled G.

Sixty runs at 600k moves each on Sa at five with {4/9, 16/9} forbidden never
reached zero, and best-of-one was hit eighteen times.  That reads as evidence
for unsatisfiability until you look at how the G instance fell: the focused
run there found its colouring at move 375898 of an 800k budget, and the
earlier 400k runs had not.  The budget was the variable, not the seed count.

Sa's instance is smaller but denser in the relevant sense -- 2523 constraints
on 397 points against G's 9435 on 1581 -- and a tabu search on a small tight
instance needs longer, not more restarts, because each restart throws away the
structure it has found.

So: far longer runs, far fewer of them.  Three million moves apiece.  If one
lands, the instance is satisfiable and sixty short runs were the wrong
experiment; if none does, the evidence is worth more than it was.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from tabucol import instance, tabucol

t0 = time.time()
Sa = build_Sa(F)
n, E = instance(Sa, {Fr(4, 9), Fr(16, 9)})
print(f"Sa at five, {{4/9, 16/9}}: {n} pts, {len(E)} constraints"
      f"  [{time.time()-t0:.0f}s]", flush=True)
got = Counter()
for sd in range(12):
    b, it = tabucol(n, E, 5, seed=9000 + sd, iters=3_000_000)
    got[b] += 1
    if b == 0:
        print(f"*** SATISFIABLE -- seed {9000+sd} reached zero after {it} "
              f"moves  [{time.time()-t0:.0f}s]", flush=True)
        break
    print(f"   seed {9000+sd}: best {b} after {it} moves; so far "
          f"{dict(sorted(got.items()))}  [{time.time()-t0:.0f}s]", flush=True)
else:
    print(f"12 deep runs of 3M moves, none reached zero, best {min(got)}, "
          f"{dict(sorted(got.items()))}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
