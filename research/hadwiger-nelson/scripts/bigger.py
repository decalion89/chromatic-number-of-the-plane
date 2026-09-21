"""Close under the largest dihedral group the field actually contains.

Sa is the twelve-element dihedral closure of a 39-point seed: six rotations by
sixty degrees and a reflection.  But a rotation matrix needs BOTH its cosine
and its sine in the field, and de Grey's field holds more than sixty degrees.
cos 30 = sqrt3/2 and sin 30 = 1/2 are both there, so the rotation of order
TWELVE is available and the dihedral group of order twenty-four with it --
twice what the construction used.  (Seventy-two degrees is not available, for
all that cos 72 is: its sine needs sqrt(10 + 2 sqrt5), which is quartic and
not multiquadratic, so no five-fold symmetry lives in this field.)

Twice the group means twice the ring.  Sa's D = 4 ring is a hexagon, six
points and three antipodal pairs; under the order-twelve rotation it becomes a
twelve-gon.  Since the cap -- the ring takes at most k - 2 colours -- is what
de Grey's construction actually rests on, and since a longer ring is a
stronger thing to cap, this is the first place to look that he did not.

Measured for both groups side by side, at four colours where the answer is
known and at five where it is not.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_S, build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, doubly_usable_ring, ring_palette_bound
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

t0 = time.time()
half = K.rational(Fr(1, 2))
rot60 = Rotation(half, K.sqrt(3) * half)
rot30 = Rotation(K.sqrt(3) * half, half)
assert (rot30.cos * rot30.cos + rot30.sin * rot30.sin) == K.rational(1)
ZERO = Point(K.zero(), K.zero())


def closure(seed, rot, order):
    seen, out = set(), []
    for p in seed:
        for base in (p, Point(p.x, -p.y)):
            q = base
            for _ in range(order):
                if q not in seen:
                    seen.add(q)
                    out.append(q)
                q = rot(q)
    return out


S = build_S(K)
print(f"seed S: {len(S)} points  [{time.time()-t0:.0f}s]", flush=True)


def report(P, label, ks=(4, 5)):
    b = IntBasis.covering(P)
    r = b.rows(P)
    head = b.overflow_headroom(r)
    assert head < 1.0, f"{label}: headroom {head}"
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    zi = P.index(ZERO)
    d = r - r[zi]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), d2)
        if v:
            grp[v].append(int(off))
    pops = {str(D): len(m) for D, m in
            sorted(grp.items(), key=lambda t: -len(t[1]))[:8]}
    print(f"\n{label}: {n} pts, {len(E)} edges ({len(E)/n:.2f}/v); "
          f"rings about the origin {pops}  [{time.time()-t0:.0f}s]",
          flush=True)
    for k in ks:
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, c in E:
            for col in range(k):
                cls.append([-(1 + a * k + col), -(1 + c * k + col)])
        s0 = Solver(name="cd15", bootstrap_with=cls)
        if not s0.solve():
            s0.delete()
            print(f"   k={k}: *** NOT {k}-COLOURABLE -- chi >= {k+1} ***",
                  flush=True)
            continue
        s0.delete()
        caps = []
        for D, mem in sorted(grp.items()):
            if D == 1 or len(mem) < 4 or not closable_distance(D):
                continue
            bnd = ring_palette_bound(cls, n * k, mem, k)
            caps.append((str(D), len(mem), bnd))
        best = [c for c in caps if c[2] < k]
        print(f"   k={k}: colours; ring palettes "
              f"{[(D, m, bd) for D, m, bd in caps]}", flush=True)
        if best:
            print(f"   k={k}: *** CAPPED RINGS {best} ***", flush=True)


report(closure(S, rot60, 6), "D6 closure (de Grey's Sa)")
report(closure(S, rot30, 12), "D12 closure (the field's largest)")
print("DONE", flush=True)
