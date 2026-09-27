"""The second floor: U2 = U1 u u.U1, and the pair test in each.

Sa has 81 rotations of K that bite it; U1 = Sa u u1.Sa has over a thousand.
That is the right direction -- de Grey's tower is S, Sa, Y, G, not one union
-- and the cheapest decisive probe is the one that settled the first floor:
his own pair (2,0),(-2,0), which sits in every one of these and which the
K-closable filter skips because D = 16 needs sqrt-15.

One query per union, ordered by how hard the rotation bites.  A forced pair
here is a 5-chromatic unit-distance graph over Q(m, sqrt-3, sqrt-11, sqrt-15),
whose residue degree at 5 is still at least 3 and which therefore can still
block, unlike anything de Grey's own field supports.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/hn/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.homcol import closable_over
from pysat.solvers import Solver

t0 = time.time()
KCL = (1, -3, -11, 33)
with open("/tmp/hn/ktower.pkl", "rb") as fh:
    rawU, rawR = pickle.load(fh)


def flat(x):
    out = []
    for part in x:
        out.extend(part)
    return tuple(out)


def unflat(t):
    return tuple(tuple(t[3 * i:3 * i + 3]) for i in range(4))


U1 = [unflat(v) for v in rawU]
ROT = [unflat(v) for v in rawR]
P2, M2 = krat(2), krat(-2)
u1set = set(U1)
assert P2 in u1set and M2 in u1set
print(f"U1: {len(U1)} points, {len(ROT)} rotations that bite it  "
      f"[{time.time()-t0:.0f}s]", flush=True)


def edges_of(pts, zs):
    cell = {}
    for i, z in enumerate(zs):
        cell.setdefault((int(z.real // 1), int(z.imag // 1)), []).append(i)
    out = []
    for i, z in enumerate(zs):
        cx, cy = int(z.real // 1), int(z.imag // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i or abs(abs(z - zs[j]) - 1) > 1e-6:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        out.append((i, j))
    return out


seen_rot, results = set(), []
for ri, u in enumerate(ROT):
    orb = []
    v = u
    for _ in range(6):
        orb.append(flat(v))
        v = kmul(v, Z6)
    key = min(orb)
    if key in seen_rot:
        continue
    seen_rot.add(key)
    img = [kmul(u, q) for q in U1]
    pts, seen2 = list(U1), set(U1)
    for q in img:
        if q not in seen2:
            seen2.add(q)
            pts.append(q)
    zs = [zof(q) for q in pts]
    E = edges_of(pts, zs)
    idx = {q: i for i, q in enumerate(pts)}
    a, b = idx[P2], idx[M2]
    n = len(pts)
    cls = [[1 + w * 4 + c for c in range(4)] for w in range(n)]
    for x, y in E:
        for c in range(4):
            cls.append([-(1 + x * 4 + c), -(1 + y * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            print(f"  *** rotation {ri}: {n} pts NOT 4-COLOURABLE ***",
                  flush=True)
            with open("/tmp/hn/kchi5.pkl", "wb") as fh:
                pickle.dump(([flat(q) for q in pts], E), fh)
            break
        sv.conf_budget(400000)
        r = sv.solve_limited(assumptions=[1 + a * 4, -(1 + b * 4)])
    tag = {False: "FORCED", True: "free", None: "over budget"}[r]
    results.append((n, len(E), tag))
    if r is not True:
        print(f"  *** rotation {ri}: {n} pts, {len(E)} edges -- "
              f"(2,0),(-2,0) {tag} ***  [{time.time()-t0:.0f}s]", flush=True)
        if r is False:
            with open("/tmp/hn/ktower2.pkl", "wb") as fh:
                pickle.dump((ri, [flat(q) for q in pts], E, a, b), fh)
            break
    if len(results) % 25 == 0:
        print(f"  ... {len(results)} unions tried, all free so far "
              f"(last {n} pts, {len(E)} edges)  [{time.time()-t0:.0f}s]",
              flush=True)
print(f"\n{len(results)} distinct unions, "
      f"{sum(1 for _, _, t in results if t == 'FORCED')} forced  "
      f"[{time.time()-t0:.0f}s]", flush=True)
