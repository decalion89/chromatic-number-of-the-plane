"""All three rotations of the class, at one pivot instead of six.

The symmetry needs the whole orbit of pivots, but the CONTRADICTION needs only
one: the class {a, b, c} shares a colour in every 4-colouring, and a rotation
about a that puts an image of one member at distance 1 from another closes it
by itself.  Using one pivot costs a quarter of the points and keeps all the
force, at the price of the group -- which is the right trade when the object
has to be tested at five colours, where every colouring of a 39313-vertex
graph costs many minutes.

The class supports three independent rotations about a:

    two free-angle ones, from |rho(b) - c| = 1, radical 23
    one classical spindle at |b - c|^2 = 64/9, radical 247

so the union is four copies of the carrier rather than nineteen.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 23, 247))
print(f"field {F.gens}, dimension {F.dim}", flush=True)
rot60 = _rot60(F); Sa = build_Sa(F)
def orb(p, refl=False):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if refl:
        for q in list(out):
            r = Point(q.x, -q.y)
            if r not in out: out.append(r)
    return out
t0 = time.time()
seen, U = set(Sa), list(Sa)
for w in orb(Sa[199], True):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
A_, B_, C_ = 726, 1526, 2730
print(f"carrier n={n}; class {A_},{B_},{C_}   [{time.time()-t0:.0f}s]", flush=True)
u = U[B_] - U[A_]; v = U[C_] - U[A_]
P = u.x * v.x + u.y * v.y
Q = u.x * v.y - u.y * v.x
R = F.rational((Fr(64, 3) + Fr(256, 9) - 1) / 2)
D = F.sqrt(23) * F.rational(Fr(13, 18))
inv = F.rational(Fr(27, 16384))
rots = []
for sgn in (1, -1):
    co = (P * R + Q * D * F.rational(sgn)) * inv
    si = (Q * R - P * D * F.rational(sgn)) * inv
    if co * co + si * si == F.rational(1):
        r = Rotation(co, si)
        if (r.about(U[A_])(U[B_]) - U[C_]).norm2() == F.rational(1):
            rots.append(("free-angle", r))
rots.append(("classical at 64/9", rotation_joining(Fr(64, 9), F)))
print(f"  {len(rots)} rotations: {[t for t, _ in rots]}   [{time.time()-t0:.0f}s]",
      flush=True)
seen2, V = set(U), list(U)
for tag, r in rots:
    rot = r.about(U[A_])
    for p in U:
        z = rot(p)
        if z not in seen2: seen2.add(z); V.append(z)
    print(f"    after {tag}: {len(V)} points   [{time.time()-t0:.0f}s]", flush=True)
gv = build_graph(V); nv = gv.n; mv = sum(len(x) for x in gv.adj) // 2
print(f"  union: n={nv} m={mv} deg={2.0*mv/nv:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
tri = gv.find_clique(3)
for KK in (4, 5, 6):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(nv)]
    for x, y in gv.edges():
        for c in range(KK):
            cc.append([-XX(x, c), -XX(y, c)])
    for i, x in enumerate(tri):
        cc.append([XX(x, i)])
        for c in range(KK):
            if c != i: cc.append([-XX(x, c)])
    t1 = time.time()
    sv = Solver(name="cd19", bootstrap_with=cc); ok = sv.solve(); sv.delete()
    print(f"  {KK}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if ok:
        if KK >= 6: print("  *** chi = 6 ***", flush=True)
        break
    print(f"  *** refuses {KK} ***", flush=True)
    if KK == 5:
        json.dump({"field_generators": list(F.gens), "n": nv, "m": mv,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gv.vertices]},
                  open(HN_DIR + "/data/six_candidate.json", "w"))
        print("  *** written data/six_candidate.json ***", flush=True)
