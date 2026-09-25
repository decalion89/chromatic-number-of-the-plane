"""Test de Grey's own pair in the K unions -- the one the filter skips.

The unions are built on HIS Sa, so (2,0) and (-2,0) are in every one of them.
But D = 16 is not closable over K -- it needs sqrt-15 -- so the pair scan
never asked about the very pair that Y forces.  It should have: if some
u.Sa forces it where rho_4.Sa does, the spindle is available one quadratic
step away, over

    K(sqrt-15) = Q(m, sqrt-3, sqrt-11, sqrt-15),

and adjoining a quadratic to Q(m) cannot drop the residue degree at 5 below 3
-- the prime with f = 3 either stays inert, splits into primes still of f = 3,
or ramifies with f = 3 -- so that field can still block.  The directions would
carry m through u, which is what blocking needs and what his own field can
never supply.

So: 81 unions, one query each.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.degrey import S_POINTS
from pysat.solvers import Solver

t0 = time.time()
SQ3 = ksub(kmul(krat(2), Z6), K1)


def to_k(xs, ys):
    p, q = Fr(xs.get(1, 0)), Fr(xs.get(33, 0))
    r, t = Fr(ys.get(3, 0)), Fr(ys.get(11, 0))
    return kadd(kadd(krat(p), kmul(krat(r), SQ3)),
                kmul(ksub(krat(t), kmul(krat(q), SQ3)), S))


Sa, seenp = [], set()
for z in [to_k(x, y) for x, y in S_POINTS]:
    for base in (z, kconj(z)):
        w = base
        for _ in range(6):
            if w not in seenp:
                seenp.add(w)
                Sa.append(w)
            w = kmul(Z6, w)
P2, M2 = krat(2), krat(-2)
assert P2 in seenp and M2 in seenp, "(2,0) and (-2,0) must be in Sa"
assert knorm2(ksub(P2, M2)) == krat(16)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/ksolve.pkl", "rb") as fh:
    raw = pickle.load(fh)
ROT = [tuple(tuple(v[3 * i:3 * i + 3]) for i in range(4)) for v in raw]
print(f"Sa {len(Sa)} points, {len(ROT)} rotations; testing the pair "
      f"(2,0),(-2,0) at D = 16 in each  [{time.time()-t0:.0f}s]", flush=True)


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


for ri, u in enumerate(ROT):
    pts, seen2 = list(Sa), set(Sa)
    for q in Sa:
        z = kmul(u, q)
        if z not in seen2:
            seen2.add(z)
            pts.append(z)
    zs = [zof(q) for q in pts]
    E = edges_of(pts, zs)
    idx = {q: i for i, q in enumerate(pts)}
    a, b = idx[P2], idx[M2]
    n = len(pts)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for x, y in E:
        for c in range(4):
            cls.append([-(1 + x * 4 + c), -(1 + y * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            print(f"  *** rotation {ri}: NOT 4-COLOURABLE ***", flush=True)
            break
        sv.conf_budget(200000)
        r = sv.solve_limited(assumptions=[1 + a * 4, -(1 + b * 4)])
    tag = {False: "FORCED", True: "free", None: "over budget"}[r]
    print(f"  rotation {ri}: {n} pts, {len(E)} edges -- (2,0),(-2,0) {tag}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if r is False:
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/kfour.pkl", "wb") as fh:
            pickle.dump((ri, [flat(q) for q in pts], E, a, b), fh)
        print("  *** spindle it over K(sqrt-15) ***", flush=True)
        break
