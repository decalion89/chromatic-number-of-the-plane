"""Sa at four colours against G at five, exhaustively: every centre, every
radius, on the whole graph.

These are the two objects the comparison should always have been between.  Sa
is de Grey's symmetric seed closure, 397 points, 4-chromatic -- and his lemma
lives inside it.  G is 1581 points and 5-chromatic, the same kind of object
one level up.  Everything else measured here has been a universe or a ball,
where the answer depends on how much was included.

Both are small enough to do completely: every vertex as a centre, every
same-distance class about it with at least k points, no radius filter, and the
palette bound taken on the whole graph so that what comes out holds in every
supergraph.  Two numbers per ring -- the ring alone, and the ring with its
centre, which is the set de Grey's lemma is about -- plus whether the field
can bite that radius, since a cap on an unbiteable ring cannot be used.

If the count at four is large and the count at five is zero, that is the gap
stated on the two objects that matter, with nothing left to the choice of
universe.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

SEED = sys.argv[1] if len(sys.argv) > 1 else "Sa"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
builders = {"Sa": build_Sa, "Y": build_Y,
            "G": lambda f: build_G(f, as_graph=False)}
P = list(builders[SEED](K))
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{SEED}: {n} points, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])


def biteable(D):
    if D < Fr(1, 4):
        return None
    ct = 1 - Fr(1, 2) / D
    s2 = 1 - ct * ct
    for rr in CLASSES:
        q = s2 / rr
        a, c = q.numerator, q.denominator
        ra, rc = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc * rc == c:
            return f"cos={ct}" + (f", sin={Fr(ra,rc)}*sqrt{rr}" if rr != 1
                                  else "")
    return None


dm, D2 = b.dim, b.D * b.D
seen_ring = set()
caps, bitecaps, tested = [], [], 0
for ci in range(n):
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        if off != ci:
            grp[Fr(int(sq[off, 0]), D2)].append(int(off))
    for D, mem in grp.items():
        if len(mem) < k or D == 1:
            continue
        key = (ci, D)
        if key in seen_ring:
            continue
        seen_ring.add(key)
        tested += 1
        pr = ring_palette_bound(cls, n * k, mem, k)
        if pr >= k:
            continue
        pc = ring_palette_bound(cls, n * k, mem + [ci], k)
        bt = biteable(D)
        caps.append((pc, pr, ci, D, len(mem), bt))
        if bt:
            bitecaps.append(caps[-1])
        print(f"  capped: centre {ci}, D={D} (rho={float(D)**.5:.3f}), "
              f"{len(mem)} pts: ring {pr}, centre+ring {pc}, "
              f"{bt or 'not biteable in K'}  [{time.time()-t0:.0f}s]",
              flush=True)
    if ci % 200 == 0:
        print(f"   .. centre {ci}/{n}, {tested} rings tested, "
              f"{len(caps)} capped  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{SEED} at k={k}: {tested} rings over {n} centres, "
      f"{len(caps)} capped, {len(bitecaps)} of those biteable in K"
      f"  [{time.time()-t0:.0f}s]", flush=True)
deep = [c for c in caps if c[0] <= k - 2]
print(f"at de Grey's depth (centre+ring <= {k-2}): {len(deep)}", flush=True)
for c in sorted(deep)[:20]:
    print(f"   centre+ring {c[0]}, ring {c[1]}, centre {c[2]}, D={c[3]}, "
          f"{c[4]} pts, {c[5] or 'not biteable'}", flush=True)
print("DONE", flush=True)
