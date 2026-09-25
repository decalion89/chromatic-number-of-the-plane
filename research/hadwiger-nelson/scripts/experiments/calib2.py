"""Does the counting threshold actually predict the cases we know?

The heuristic says a direction set of rank r over F_5 blocks once it carries
more than r.log5/log(5/4) = 7.21 r independent lines, because a random
functional survives L lines with probability (4/5)^L and there are 5^r of
them.  It is only a heuristic -- the directions are structured -- so it is
worth checking against every verdict already established.

If it calls all of them right, the prediction it makes for the sixteen-rotation
set (880 lines against a 693 threshold, so BLOCKS) is worth the wait.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.homcol import (denominator_29_directions, edge_vectors,
                       has_homomorphism, _rank_q, _rank_mod, blocks_at)
from hn.degrey import build_G
from hn.graph import build_graph

t0 = time.time()


def check(v, name, n=5):
    r = _rank_q(v)
    lines = len({tuple(x % n for x in w) for w in v}) // 2
    need = math.ceil(r * math.log(n) / math.log(n / (n - 1)))
    pred = "blocks" if lines > need else "escapes"
    real = "blocks" if blocks_at(v, n) else "escapes"
    mark = "AGREES" if pred == real else "*** WRONG ***"
    print(f"  {name}: rank {r}, {lines} lines, threshold {need} -> "
          f"predicts {pred}, truth {real}  {mark}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return pred == real


ok = []
ok.append(check(denominator_29_directions(3), "denominator-29 set"))
G = build_graph(build_G(as_graph=False))
D = edge_vectors(G)
ok.append(check(D, "de Grey's G"))
for n in (2, 3, 4):
    check(D, f"de Grey's G at n = {n}", n)
print(f"\n{sum(ok)} of {len(ok)} at the gate  [{time.time()-t0:.0f}s]",
      flush=True)
