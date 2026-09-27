"""The free-angle spindle, which a forced class of three makes possible.

Forcing is transitive, so (a,b) and (a,c) forced-equal makes {a,b,c} a single
class.  An earlier measurement found every class in the SEQUENTIAL chain had
exactly two members and concluded the generalisation below was unavailable.
That was true of that chain and is false in general: on the symmetric carrier
from Sa[199] the two rare forced distances share their pivot, so

    (726, 1526) forced at d^2 = 64/3     and     (726, 2730) at d^2 = 256/9

put three vertices in one class.  And a class of three admits a spindle whose
angle is not constrained to make 4d^2 - 1 a square.  Any rotation about a
carrying b to distance exactly 1 from c will do, because

    c(rho(b)) = c(rho(a)) = c(a) = c(c)   and   rho(b) ~ c

so the only condition is that the circles meet:  | |a-b| - |a-c| | <= 1 <=
|a-b| + |a-c|.  Solving |rho(b) - c|^2 = r1 + r2 - 2 sqrt(r1 r2) cos t = 1 with
r1 = 64/3 and r2 = 256/9 gives

    cos t = 439 sqrt3 / 768      sin t = 13 sqrt69 / 768      69 = 3 * 23

and 578163 + 11661 = 589824 = 768^2 exactly, so this is a genuine rotation and
it needs no radical the 64/3 spindle did not already need.  The union is two
copies of a 3025-point carrier rather than the 39313 points that spending both
orbits the classical way costs.
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

F = Field((3, 11, 23))
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
print(f"carrier n={n} m={sum(len(a) for a in g.adj)//2}   [{time.time()-t0:.0f}s]",
      flush=True)
A_, B_, C_ = 726, 1526, 2730
d1 = (U[A_] - U[B_]).norm2(); d2 = (U[A_] - U[C_]).norm2()
print(f"  |a-b|^2 = {d1}   |a-c|^2 = {d2}", flush=True)
assert d1.is_rational() and Fr(d1.c[0]) == Fr(64, 3)
assert d2.is_rational() and Fr(d2.c[0]) == Fr(256, 9)
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
for (x, y, tag) in ((A_, B_, "a,b"), (A_, C_, "a,c"), (B_, C_, "b,c")):
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[X(x, 0), X(y, 1)]); s.delete()
    print(f"  {tag} can differ: {ok}  -> forced equal: {not ok}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    assert not ok, f"{tag} is not forced"
# The angle is NOT the target angle between the two vectors: rotating by that
# would only be right if they started aligned.  With u = b - a and v = c - a,
#
#     |rho(u) - v|^2 = r1 + r2 - 2 <rho(u), v>
#     <rho(u), v>    = cos t (u.v) + sin t (u x v)
#
# so cos t and sin t satisfy one linear equation, P cos t + Q sin t = R with
# P = u.v, Q = u x v, R = (r1 + r2 - 1)/2, together with cos^2 + sin^2 = 1.
# Since P^2 + Q^2 = r1 r2, the two solutions are
#
#     cos t = (P R +- Q D) / (r1 r2),   sin t = (Q R -+ P D) / (r1 r2)
#     D = sqrt(r1 r2 - R^2)
#
# and here r1 r2 - R^2 = 16384/27 - (439/18)^2 = 3887/324, whose square root is
# 13 sqrt23 / 18 -- rational times a radical the field already has.
u = U[B_] - U[A_]
v = U[C_] - U[A_]
P = u.x * v.x + u.y * v.y
Q = u.x * v.y - u.y * v.x
R = F.rational((Fr(64, 3) + Fr(256, 9) - 1) / 2)
D = F.sqrt(23) * F.rational(Fr(13, 18))
inv = F.rational(Fr(27, 16384))
built = None
for sgn in (1, -1):
    co = (P * R + Q * D * F.rational(sgn)) * inv
    si = (Q * R - P * D * F.rational(sgn)) * inv
    if co * co + si * si == F.rational(1):
        rot = Rotation(co, si)
        img = rot.about(U[A_])(U[B_])
        d = (img - U[C_]).norm2()
        print(f"  sign {sgn:+d}: cos^2+sin^2 = 1, |rho(b)-c|^2 = {d}", flush=True)
        if d == F.rational(1): built = rot
    else:
        print(f"  sign {sgn:+d}: cos^2+sin^2 = {co*co+si*si} -- not a rotation",
              flush=True)
assert built is not None, "neither sign put rho(b) one unit from c"
rot = built

seen2, V = set(U), list(U)
for w in orb(U[A_]):
    r = rot.about(w)
    for p in U:
        z = r(p)
        if z not in seen2: seen2.add(z); V.append(z)
gv = build_graph(V); nv = gv.n; mv = sum(len(x) for x in gv.adj) // 2
print(f"\n  free-angle spindle over 6 pivots: n={nv} m={mv} deg={2.0*mv/nv:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
tri = gv.find_clique(3)
for KK in (4, 5):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(nv)]
    for u, v in gv.edges():
        for c in range(KK):
            cc.append([-XX(u, c), -XX(v, c)])
    for i, v in enumerate(tri):
        cc.append([XX(v, i)])
        for c in range(KK):
            if c != i: cc.append([-XX(v, c)])
    t1 = time.time()
    sv = Solver(name="cd19", bootstrap_with=cc); ok = sv.solve(); sv.delete()
    print(f"  {KK}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if ok: break
    print(f"  *** refuses {KK} ***", flush=True)
    if KK == 4:
        json.dump({"field_generators": list(F.gens), "n": nv, "m": mv,
                   "cos": [[c.numerator, c.denominator] for c in rot.cos.c],
                   "sin": [[c.numerator, c.denominator] for c in rot.sin.c],
                   "class": [A_, B_, C_], "d2": ["64/3", "256/9"],
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gv.vertices]},
                  open(HN_DIR + "/data/free_angle.json", "w"))
        print("  written data/free_angle.json", flush=True)
