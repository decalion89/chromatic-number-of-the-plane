"""Why n = 8 lets the necklace through, and whether growth can close it.

The screen found a homomorphism to Z/8 whose Cayley graph needs only four
colours.  Growth adds directions, which both constrains phi further and
enlarges the connection set, so it ought to help -- UNLESS the failure is
structural, with every direction landing in a fixed pattern mod 8 no matter
how many are added.  Cay(Z/8, S) is bipartite as soon as S is all-odd, and no
amount of growth escapes that.

So the failing phi is printed with its connection set, which says at once
which case this is.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.homcol import periodic_screen, cayley_chromatic

with open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/neckdirs.pkl",
          "rb") as fh:
    dirs = [tuple(d) for d in pickle.load(fh)]
print(f"{len(dirs)} directions", flush=True)

r = periodic_screen(dirs, 8, cap=5, rounds=200)
print("verdict:", {k: v for k, v in r.items() if k != "phi"}, flush=True)
if r.get("colourable"):
    c = r["phi"]
    print("phi =", c, flush=True)
    S = sorted({sum(x * y for x, y in zip(c, d)) % 8 for d in dirs})
    print("connection set S =", S, flush=True)
    print("all odd?", all(s % 2 for s in S),
          "  |S| =", len(S), "of 7 possible", flush=True)
    print("chi(Cay(Z/8, S)) =", cayley_chromatic(8, S, cap=8), flush=True)
    # what would it take?  S must miss enough of Z/8 to stay 5-colourable
    for miss in range(1, 8):
        full = [s for s in range(1, 8) if s != miss]
        print(f"  if S were {full}: chi = {cayley_chromatic(8, full, cap=8)}",
              flush=True)
