"""The exact minimum, by MaxSAT rather than by greedy.

Greedy covering carries a ln(n) factor, so its 25 says little about the truth.
The covering condition is a plain CNF -- one clause per point of PG(5,5),
listing the directions whose hyperplane contains it -- and minimising the
number of chosen directions under it is exactly what an unweighted MaxSAT
solver does.  The answer is then the real size of the smallest step set over
Q(zeta_7) that no coset 5-colouring survives.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
from hn.homcol import has_homomorphism

S = [tuple(v) for v in json.load(open(
    "/tmp/hn/blocking.json"))]
D = 6
pts = [v for v in itertools.product(range(5), repeat=D)
       if any(v) and next(x for x in v if x) == 1]
print(f"{len(S)} directions, {len(pts)} projective points", flush=True)

hit = [[] for _ in pts]
for j, d in enumerate(S):
    r = tuple(x % 5 for x in d)
    for i, p in enumerate(pts):
        if sum(a * b for a, b in zip(r, p)) % 5 == 0:
            hit[i].append(j + 1)

w = WCNF()
for i, cl in enumerate(hit):
    if not cl:
        print(f"point {pts[i]} is covered by NOTHING -- these 87 cannot block")
        sys.exit()
    w.append(cl)
for j in range(len(S)):
    w.append([-(j + 1)], weight=1)

with RC2(w) as rc2:
    model = rc2.compute()
    chosen = [S[j] for j in range(len(S)) if model[j] > 0]

print(f"\nEXACT MINIMUM: {len(chosen)} directions "
      f"(greedy found 25; the counting floor is 5, the quadric floor 7)")
phi, _ = has_homomorphism(list(chosen), 5)
print("independent SAT screen: "
      + ("STILL 5-COLOURABLE (bug)" if phi else "NO homomorphism -- blocks"))
for d in chosen:
    print("   ", d)
json.dump([list(d) for d in chosen], open(
    "/tmp/hn/minblock.json", "w"))
