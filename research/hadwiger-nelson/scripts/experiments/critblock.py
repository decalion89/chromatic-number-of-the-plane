"""Can a critical graph block?

Every 6-chromatic unit-distance graph is blocked, and U = G u {w.G} shows the
gate is passable while keeping chi >= 5.  But in every construction so far the
blocking dies in the critical core: strip to a 4-critical subgraph and one
cheap spindle is left.

So the question is whether a critical graph can block at all.  Two facts pin
the small end down.

The Moser spindle is UNIQUE.  Its two arms are rhombi with tips at distance
sqrt3, joined by a unit edge, so the arm rotation r satisfies |1 - r|^2 = 1/3,
hence r + rbar = 5/3 and r = (5 +- sqrt-11)/6 -- no freedom at all.  And
blocking is invariant under a global rotation, since phi -> phi . u is a
bijection between the homomorphisms out of M and those out of uM.  So one test
settles every 7-vertex 4-critical unit-distance graph at once.

Past that, blocking needs directions: a set of hyperplanes covering PG(r-1,5)
needs at least six of them (a pencil through a fixed codimension-2 space is
the extremal cover), so a blocked graph has at least six edge directions.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from field24 import (K, D, e_zero, e_of, e_add, e_sub, e_mul, e_norm2,
                     ONE, Z6, RHO)
from hn.homcol import has_homomorphism
from pysat.solvers import Solver
import cmath


def ivecs(elts):
    den = 1
    for v in elts:
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
    out = []
    for v in elts:
        w = tuple(int(x * den) for x in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        w = tuple(t // g for t in w) if g > 1 else w
        if w not in out:
            out.append(w)
    return out


def zof(q):
    za = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.a))
    zb = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.b))
    return za + zb * cmath.sqrt(-11)


def udg(P):
    zs = [zof(q) for q in P]
    buck = {}
    for i, z in enumerate(zs):
        buck.setdefault((round(z.real * 1e6), round(z.imag * 1e6)), []).append(i)
    E = []
    R = 1e6
    for i, z in enumerate(zs):
        cx, cy = z.real * R, z.imag * R
        for dx in range(-int(R) - 1, int(R) + 2, 1) if False else ():
            pass
        for gx in range(int(cx - R) // 1, 0):
            pass
        for key, idxs in buck.items():
            pass
        break
    E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
         if abs(abs(zs[i] - zs[j]) - 1) < 1e-7
         and e_norm2(e_sub(P[j], P[i])) == 1]
    return E


def chrom(P, E, cap=5):
    for k in (3, 4, 5, 6):
        if k > cap:
            return None
        cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return None


def dirs_of(P, E, keep=None):
    out = []
    for a, b in E:
        if keep is not None and (a not in keep or b not in keep):
            continue
        d = e_sub(P[b], P[a])
        out.append(tuple(d.a) + tuple(d.b))
        d = e_sub(P[a], P[b])
        out.append(tuple(d.a) + tuple(d.b))
    return ivecs(out)


rh = [e_zero(), ONE, Z6, e_add(ONE, Z6)]
sp, seen = [], set()
for p in list(rh) + [e_mul(RHO, q) for q in rh[1:]]:
    if p not in seen:
        seen.add(p); sp.append(p)
Esp = udg(sp)
dv = dirs_of(sp, Esp)
phi, _ = has_homomorphism(dv, 5)
print(f"Moser spindle: {len(sp)} points, {len(Esp)} edges, "
      f"{len(dv)} projective directions", flush=True)
print("  " + ("BLOCKS" if phi is None else "admits a coset colouring")
      + "  -- and the spindle is unique, so this settles every 7-vertex "
        "4-critical unit-distance graph", flush=True)
