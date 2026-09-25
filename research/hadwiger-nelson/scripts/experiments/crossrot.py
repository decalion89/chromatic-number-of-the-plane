"""Which rotation joins Sa to its own image most strongly?

De Grey chose the D = 4 bite for a structural reason -- it moves each point of
the hexagon by exactly one, threading the ring -- and it buys almost no
contact: Sa and rho(Sa) share ONE vertex and SIX edges.  That was enough at
four colours because Sa's statement about its ring is a palette cap, which is
strong enough to survive on six edges.  At five colours nothing is capped, so
the question is whether contact can be bought instead.

Every rotation about the origin commutes with the sixty-degree rotation, so
every union Sa u rho(Sa) stays dihedrally symmetric exactly as Sa is, and the
ring structure survives whichever rho is used.  So the choice is free, and it
can be made on cross-edge count rather than on the ring.

Enumerated over every rational squared distance the field can join, both
directions, with the exact edge finder and the overflow certificate checked.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, doubly_usable_ring

t0 = time.time()
Sa = build_Sa(K)
n0 = len(Sa)
b0 = IntBasis.covering(Sa)
r0 = b0.rows(Sa)
E0 = len(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b0, r0)))
print(f"Sa: {n0} pts, {E0} edges  [{time.time()-t0:.0f}s]", flush=True)

Ds = []
for den in range(1, 37):
    for num in range(1, 40 * den + 1):
        D = Fr(num, den)
        if D in Ds or D == 1:
            continue
        if closable_distance(D):
            Ds.append(D)
Ds = sorted(set(Ds))
print(f"{len(Ds)} rational rings the field can join  [{time.time()-t0:.0f}s]",
      flush=True)

rows = []
for D in Ds:
    rot = rotation_joining(D, K)
    for sgn, name in ((+1, "+"), (-1, "-")):
        rr = rot if sgn > 0 else Rotation(rot.cos, -rot.sin)
        img = [rr(p) for p in Sa]
        seen = set(Sa)
        U = list(Sa) + [q for q in img if q not in seen]
        b = IntBasis.covering(U)
        r = b.rows(U)
        if b.overflow_headroom(r) >= 1.0:
            continue
        E = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
        shared = len(Sa) + len(img) - len(U)
        cross = len(E) - 2 * E0 + 0   # edges inside each copy are E0 each,
        # minus any edge counted twice because both endpoints are shared
        rows.append((cross, len(U), shared, D, name))
    if len(rows) % 40 == 0 and rows:
        print(f"   {len(rows)} rotations tried, best so far "
              f"{max(rows)[0]} cross edges  [{time.time()-t0:.0f}s]",
              flush=True)
rows.sort(reverse=True)
print(f"\n{len(rows)} rotations measured  [{time.time()-t0:.0f}s]", flush=True)
print(f"{'cross':>7} {'points':>7} {'shared':>7}  ring      dir  usable",
      flush=True)
for cross, npts, shared, D, name in rows[:25]:
    print(f"{cross:7d} {npts:7d} {shared:7d}  D={str(D):8s} {name}   "
          f"{'doubly' if doubly_usable_ring(D) else 'spindle'}", flush=True)
print("\nde Grey's choice:", flush=True)
for cross, npts, shared, D, name in rows:
    if D == 4:
        print(f"{cross:7d} {npts:7d} {shared:7d}  D={str(D):8s} {name}",
              flush=True)
print("DONE", flush=True)
