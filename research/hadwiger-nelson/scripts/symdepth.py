"""The palette of every ring of the symmetrised G, on the whole graph.

Ga -- G closed under the order-twelve group -- is the first carrier here that
is tight and symmetric at once: it contains G so it needs five colours, it
colours with five so it needs no more, and every point sits on a complete ring
about the origin.  That is de Grey's situation at four, one level up.

And it is small enough that no ball is needed.  A palette bound proved on the
whole graph is the strongest form of the statement: not "this ball cannot show
five colours", which leaves the rest of the graph unexamined, but "Ga cannot",
which then holds in every supergraph of Ga as well.

Measured for each ring about the origin and about the busiest other centres,
twice: the ring alone, and the ring together with its centre, which is the set
de Grey's lemma is actually about.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 6
PKL = sys.argv[3] if len(sys.argv) > 3 else "symG.pkl"
t0 = time.time()
with open(f"/tmp/claude-0/-home-user-darwin-50/"
          f"aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/{PKL}", "rb") as f:
    P = pickle.load(f)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
deg = {v: len(adj[v]) for v in range(n)}
print(f"{PKL}: {n} points, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
dm, D2 = b.dim, b.D * b.D
oi = next((i for i, p in enumerate(P)
           if not any(p.x.c) and not any(p.y.c)), 0)
centres = [oi] + [v for v in sorted(range(n), key=lambda v: -deg[v])
                  if v != oi][:NC - 1]
best = []
for ci in centres:
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        if off != ci:
            D = Fr(int(sq[off, 0]), D2)
            if 0 < D <= 4:
                grp[D].append(int(off))
    for D, mem in sorted(grp.items()):
        if len(mem) < k:
            continue
        pr = ring_palette_bound(cls, n * k, mem, k)
        pc = ring_palette_bound(cls, n * k, mem + [ci], k)
        flag = ""
        if pc < k:
            flag = "   <<< CAPPED"
            best.append((pc, ci, D, len(mem)))
        print(f"  centre {ci} (deg {deg[ci]}), D={D}, {len(mem)} pts: "
              f"ring palette {pr}, centre+ring {pc}{flag}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\ncapped centre+rings: {len(best)}", flush=True)
for e in sorted(best):
    print(f"   depth {e[0]} of {k}: centre {e[1]}, D={e[2]}, {e[3]} pts",
          flush=True)
print("DONE", flush=True)
