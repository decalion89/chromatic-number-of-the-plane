"""The real step: spindle the tight graph, do not merely rotate it.

Rotating the 24-point tight core by every multiple of 30 degrees multiplied
its edges eightfold and left chi_c at 4.  That is because a rotation by a
pretty angle is not de Grey's step.  His step picks the angle so that a
CHOSEN pair lands exactly one unit from its own image, which creates a new
unit edge between the copy and the original and ties the two together.  The
angle is dictated by the geometry, not by symmetry.

For a pivot a and a point b at distance d, the image of b is a unit from b
when 2d sin(theta/2) = 1, that is

        cos theta = 1 - 1/(2 d^2).

d^2 lies in the field, so cos theta does too; the step exists exactly when
1 - cos^2(theta) is a SQUARE in the field, which is the condition this
project derived earlier for the bite's radius.  Every admissible (pivot,
point) pair is enumerated, the union taken, and chi_c measured.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver
from sqrtK import Roots

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
t0 = time.time()
base = pickle.load(open(SC + "tightmin.pkl", "rb"))
R = Roots(K)
one, two = K.one(), K.rational(2)


def rot(pt, piv, c, s):
    dx, dy = pt.x - piv.x, pt.y - piv.y
    return Point(piv.x + dx * c - dy * s, piv.y + dx * s + dy * c)


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
        return "overflow", len(pts), 0
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(pts)
    cands = sorted({Fr(a, q) for q in range(1, 10)
                    for a in range(2, hi * q + 1)
                    if Fr(a, q).denominator == q and 3 <= Fr(a, q) <= hi})
    for rr in cands:
        p, q = rr.numerator, rr.denominator
        if p * n > 500000:
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
    return "> 5", n, len(E)


print(f"base: {len(base)} points  [{time.time()-t0:.0f}s]", flush=True)
r0, n0, m0 = chi_c(base)
print(f"   chi_c = {r0} on {n0} pts, {m0} edges  [{time.time()-t0:.0f}s]",
      flush=True)
angles = []
for i, a in enumerate(base):
    for j, bpt in enumerate(base):
        if i == j:
            continue
        d2 = a.dist2(bpt)
        if float(d2) < 0.25:
            continue                      # needs 2d >= 1 for the chord
        cos = one - one / (two * d2)
        sin2 = one - cos * cos
        sn = R.sqrt(sin2)
        if sn is None:
            continue
        angles.append((i, cos, sn, float(d2)))
seen = set()
uniq = []
for i, c, s, d2 in angles:
    k = (i, tuple(c.c), tuple(s.c))
    if k not in seen:
        seen.add(k)
        uniq.append((i, c, s, d2))
print(f"{len(uniq)} admissible spindles (pivot, angle) out of "
      f"{len(base)*(len(base)-1)} pairs  [{time.time()-t0:.0f}s]", flush=True)
# One spindle never moves chi_c -- de Grey's construction is many copies,
# not one.  Accumulate instead: spindle the graph, spindle the result, and
# keep going, measuring chi_c at every step so the moment it moves (or
# stops moving) is visible rather than inferred.
acc = list(base)
cur = r0
for k, (i, c, s, d2) in enumerate(uniq):
    piv = base[i] if i < len(base) else acc[i % len(acc)]
    cop = [rot(p, piv, c, s) for p in acc]
    un = dedupe(acc + cop)
    if len(un) > 3000:
        print(f"   stopping at {len(un)} points, too large to measure",
              flush=True)
        break
    rr, n, m = chi_c(un)
    if rr == "overflow":
        print(f"   step {k}: overflow at {n} points, skipping", flush=True)
        continue
    moved = "   <<< MOVED" if rr != cur else ""
    print(f"   step {k} (pivot v{i}, d^2={d2:.3f}): {n:5d} pts, {m:6d} "
          f"edges, chi_c = {rr}{moved}  [{time.time()-t0:.0f}s]", flush=True)
    acc, cur = un, rr
    pickle.dump(acc, open(SC + "spun24.pkl", "wb"))
    if rr == "> 5":
        print("   *** chi_c > 5: that is chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
