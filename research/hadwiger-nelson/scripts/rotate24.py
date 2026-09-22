"""Does rotating the tight 24-point graph buy circular chromatic number?

The tightness synthesis says caps live where a carrier is tight, and the
smallest tight carrier here is 24 points with chi_c = 4.  De Grey's step --
take a structure and rotate copies so their constraints conflict -- moved
chi_c from 4 (on Y) to 4.5 (on G) in this project's measurements.  Doing it
to 397 points and then again is out of reach; doing it to 24 is not.

The field contains cos 30 = sqrt3/2 and sin 30 = 1/2 exactly, so rotation by
multiples of 30 degrees about the origin is available with no adjunction and
no rounding.  Take the union of the graph with its rotations, measure chi_c
of each union, and see whether the step buys anything at all at this size.
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
t0 = time.time()
base = pickle.load(open(SC + "tightmin.pkl", "rb"))
half = K.rational(Fr(1, 2))
r3h = K.sqrt(3) * half
# rotation by 30 degrees, exactly
C, S = r3h, half


def rot(pt, c, s):
    return Point(pt.x * c - pt.y * s, pt.x * s + pt.y * c)


def compose(m):
    c, s = K.one(), K.zero()
    for _ in range(m):
        c, s = c * C - s * S, c * S + s * C
    return c, s


def dedupe(pts):
    seen, out = set(), []
    for p in pts:
        k = (tuple(p.x.c), tuple(p.y.c))
        if k not in seen:
            seen.add(k)
            out.append(p)
    return out


def chi_c(pts, hi=5):
    b = IntBasis.covering(pts)
    r = b.rows(pts)
    if b.overflow_headroom(r) >= 1.0:
        return None, None, None
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(pts)
    cands = sorted({Fr(a, q) for q in range(1, 10)
                    for a in range(2, hi * q + 1)
                    if Fr(a, q).denominator == q and 2 <= Fr(a, q) <= hi})
    for rr in cands:
        p, q = rr.numerator, rr.denominator
        if p * n > 400000:
            continue
        cl = [[1 + v * p + j for j in range(p)] for v in range(n)]
        for a, c in E:
            for j in range(p):
                for d in range(-(q - 1), q):
                    cl.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
        cl += [[-(1 + j)] for j in range(1, p)] + [[1]]
        s = Solver(name="cd15", bootstrap_with=cl)
        ok = s.solve()
        s.delete()
        if ok:
            return rr, n, len(E)
    return None, n, len(E)


print(f"base: {len(base)} points; chi_c = ?", flush=True)
r0, n0, m0 = chi_c(base)
print(f"   {n0} pts, {m0} edges, chi_c = {r0}  [{time.time()-t0:.0f}s]",
      flush=True)
acc = list(base)
for m in range(1, 12):
    c, s = compose(m)
    acc = dedupe(acc + [rot(p, c, s) for p in base])
    rr, n, mm = chi_c(acc)
    print(f"   + rotation {30*m:3d} deg: {n:5d} pts, {mm:6d} edges, "
          f"chi_c = {rr if rr else '> 5'}  [{time.time()-t0:.0f}s]",
          flush=True)
    if rr is None or float(rr) > 4.0:
        pickle.dump(acc, open(SC + f"rot24_{m}.pkl", "wb"))
        if rr is None or float(rr) > 4.5:
            print("   *** the step moved it past 4.5 ***", flush=True)
print("DONE", flush=True)
