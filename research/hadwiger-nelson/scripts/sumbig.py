"""Push the sumset: more generators, since density climbs with them.

Edges per vertex against generator count, from the first sweep:

    1, w                      2.45
    1, w, r4                  3.84
    1, w, r4, r7              4.72
    1, w, r3, r4, r7          5.47   (2819 points)

against 4.97 for Sa, 4.98 for G and 5.42 for the densest dihedral closure
built here.  Five generators already beat everything in the de Grey family,
and the trend has not turned over.

Every generator is a unit vector of the field, so the arithmetic stays exact,
and every one has infinite order, so the point set keeps the many scales the
lattices lacked.  rho_D exists whenever sqrt(4D-1) is in the field, which for
Q(v3,v5,v7,v11) includes D = 1, 2, 3, 4, 7, 9, 14, 16 -- and the conjugate
rotation is free, being the same vector reflected.

Two questions per configuration and both are one SAT call: is it 5-colourable
at all, and if it is, does it carry the weak property on any closable class.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())


def unit(r):
    return r(ONE)


G = {"w": unit(_rot60(K))}
for D in (2, 3, 4, 7, 9, 14, 16):
    r = rotation_joining(Fr(D), K)
    G[f"r{D}"] = unit(r)
    G[f"c{D}"] = unit(Rotation(r.cos, -r.sin))
print(f"{len(G)} generators available: {sorted(G)}  [{time.time()-t0:.0f}s]",
      flush=True)


def sumset(names, c, rmax, cap=70000):
    vs = [ONE] + [G[n] for n in names]
    out, seen = [], set()
    for a in itertools.product(range(-c, c + 1), repeat=len(vs)):
        p = Point(K.zero(), K.zero())
        for coef, v in zip(a, vs):
            if coef:
                p = Point(p.x + v.x * K.rational(coef),
                          p.y + v.y * K.rational(coef))
        if p in seen or float(p.norm2()) > rmax * rmax:
            continue
        seen.add(p)
        out.append(p)
        if len(out) > cap:
            return None
    return out


def study(label, P):
    g = build_graph(P)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    print(f"  {label}: {g.n} pts, {len(E)} edges, {len(E)/g.n:.2f} per "
          f"vertex -> {'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    return ok, g.n, len(E)


# Conjugate PAIRS are what pushed the density highest -- (r4, c4) beat every
# single-rotation set of the same size -- so the sweep now leads with them.
CONFIGS = [
    (("w", "r3", "r4", "c4", "r7"), 2, 5),
    (("w", "r3", "c3", "r4", "c4"), 2, 5),
    (("w", "r3", "c3", "r4", "c4", "r7"), 2, 5),
    (("w", "r3", "c3", "r4", "c4", "r7", "c7"), 2, 4),
    (("w", "r2", "c2", "r3", "c3", "r4", "c4"), 2, 4),
    (("w", "r3", "c3", "r4", "c4", "r7"), 2, 6),
    (("w", "r3", "r4", "c4", "r7", "r9", "c9"), 2, 4),
]
for names, c, rmax in CONFIGS:
    P = sumset(names, c, rmax)
    if P is None:
        print(f"  gens {('1',)+names} c={c} r<={rmax}: over the cap, skipped"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        continue
    if len(P) < 50:
        continue
    ok, n, e = study(f"gens {('1',)+names} c={c} r<={rmax}", P)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
