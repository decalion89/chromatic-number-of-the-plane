"""Grow the tight set, because size is where the pigeonhole starts paying.

G carries a set of seven points whose palette is 4 out of 5: no proper
5-colouring of G shows all five colours on them.  Four hundred and fifty
random seven-sets, a hundred and fifty of them drawn from the very ball that
contains it, all show five.  So the set is real and it is rare.

It is not yet useful.  A set of m points held to t colours forces a
monochromatic class of size ceil(m/t), and the generic bound is ceil(m/k), so
the constraint only buys something when ceil(m/t) > ceil(m/k).  At m = 7,
t = 4, k = 5 both are 2, and a monochromatic pair among seven points is free.

The threshold is nearby.  Nine points at palette 4 give 3 against 2, and so
does ten; seven points at palette 3 give 3 against 2 as well.  Either forces a
monochromatic TRIPLE, which is a statement no colouring argument gets for
free.

So: keep the palette bound where it is and add points, one at a time, taking
whichever addition the bound survives.  Every step is verified exactly, and
the bound travels upward, so whatever survives here survives in every
supergraph of G.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
S = [7, 107, 109, 269, 406, 664, 668]
p0 = ring_palette_bound(cls, n * k, S, k)
print(f"starting set: {len(S)} points, palette {p0}, target <= {TARGET}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
assert p0 <= TARGET


def gain(m, t):
    return -(-m // t) > -(-m // k)


cur = list(S)
pool = set()
for v in cur:
    pool |= adj[v]
    for u in adj[v]:
        pool |= adj[u]
        for w in adj[u]:
            pool |= adj[w]
pool -= set(cur)
print(f"candidate pool: {len(pool)} points within three hops"
      f"  [{time.time()-t0:.0f}s]", flush=True)
stalled = False
while not stalled:
    stalled = True
    order = sorted(pool)
    for v in order:
        t = cur + [v]
        if ring_palette_bound(cls, n * k, t, k) <= TARGET:
            cur = t
            pool.discard(v)
            klass = -(-len(cur) // TARGET)
            note = (f"   <<< PIGEONHOLE GAIN: forces a monochromatic "
                    f"{klass}" if gain(len(cur), TARGET) else "")
            print(f"   grew to {len(cur)} points, palette still "
                  f"<= {TARGET}{note}  [{time.time()-t0:.0f}s]", flush=True)
            with open(SC + f"grown_{k}_{TARGET}.pkl", "wb") as f:
                pickle.dump((cur, [P[i] for i in cur]), f)
            stalled = False
            break
final = ring_palette_bound(cls, n * k, cur, k)
print(f"\nmaximal at {len(cur)} points, palette {final} of {k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"forces a monochromatic class of size "
      f"{-(-len(cur)//final)}; generic would be {-(-len(cur)//k)}", flush=True)
print(f"set: {cur}", flush=True)
print("DONE", flush=True)
