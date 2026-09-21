"""Palette depth of every ring of the pivot-symmetrised G, on the whole graph.

Gp contains G so it needs five colours; it colours with five so it needs no
more; and it is closed under the order-twelve group about the pivot, so every
point sits on a complete ring there.  That is the first object in this work
that is tight and symmetric at the same time, which is exactly the position
de Grey was in at four colours.

The question his lemma answers is not whether a ring shows all k colours.  It
is how few colours the ring and its centre are forced to share between them --
two out of four, in his case.  So both numbers are measured here, on the whole
graph rather than on a ball, which is the strongest form: a bound proved on Gp
holds in every supergraph of Gp.

A ring at distance one from its centre is capped for free, since the centre is
adjacent to all of it; those are reported and ignored.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = pickle.load(open(SC + "pivG.pkl", "rb"))
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
print(f"Gp: {n} points, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
PIV = Point(K.rational(-2), K.zero())
pi = next(i for i, p in enumerate(P) if p == PIV)
centres = [pi] + [v for v in sorted(range(n), key=lambda v: -deg[v])
                  if v != pi][:NC - 1]
dm, D2 = b.dim, b.D * b.D
found = []
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
    print(f"-- centre {ci} (deg {deg[ci]}): {len(grp)} rings"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for D, mem in sorted(grp.items(), key=lambda t: -len(t[1])):
        if len(mem) < k:
            continue
        pr = ring_palette_bound(cls, n * k, mem, k)
        pc = ring_palette_bound(cls, n * k, mem + [ci], k)
        note = ""
        if D == 1:
            note = "  (free: the centre touches the whole ring)"
        elif pc < k:
            note = "   <<< CAPPED, NOT FREE"
            found.append((pc, ci, D, len(mem)))
        elif pr < k:
            note = "   <<< ring capped, centre+ring not"
        print(f"   D={D}, {len(mem)} pts: ring {pr}, centre+ring {pc}{note}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nnon-trivial capped centre+rings: {len(found)}", flush=True)
for e in sorted(found):
    print(f"   depth {e[0]} of {k}: centre {e[1]}, D={e[2]}, {e[3]} pts",
          flush=True)
print("DONE", flush=True)
