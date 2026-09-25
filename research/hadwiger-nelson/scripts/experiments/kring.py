"""Sa u rho_p(Sa) over K, scanned for a forced pair -- with a budget this time.

The earlier run of this used a list comprehension for the hits, so nothing
printed until every pair was done, and it sat for minutes with no output.  A
silent stall on this scan means one of two things and they look identical from
outside: a slow satisfiable query, or a real unsatisfiability proof -- which is
the prize.  So: a conflict budget on every query, progress as it goes, and the
ones that blow the budget revisited at the end without it.

K closes D = 1, 3, 7, 19, 25, ... so the unions available are rho_3 = the Moser
rotation and rho_7 = (13 + 3 sqrt-3)/14.  A forced pair at any K-closable
distance ends it: spindling it gives a 5-chromatic unit-distance graph over the
one field that also blocks at every modulus up to five.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.degrey import S_POINTS
from hn.homcol import closable_over
from pysat.solvers import Solver

t0 = time.time()
KCL = (1, -3, -11, 33)
SQ3 = ksub(kmul(krat(2), Z6), K1)
assert kmul(SQ3, SQ3) == krat(-3)


def to_k(xs, ys):
    p, q = Fr(xs.get(1, 0)), Fr(xs.get(33, 0))
    r, t = Fr(ys.get(3, 0)), Fr(ys.get(11, 0))
    return kadd(kadd(krat(p), kmul(krat(r), SQ3)),
                kmul(ksub(krat(t), kmul(krat(q), SQ3)), S))


Sa, seen = [], set()
for z in [to_k(x, y) for x, y in S_POINTS]:
    for base in (z, kconj(z)):
        q = base
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                Sa.append(q)
            q = kmul(Z6, q)
zsA = [zof(p) for p in Sa]
print(f"Sa over K: {len(Sa)} points  [{time.time()-t0:.0f}s]", flush=True)

ROT = {Fr(3): RHO,
       Fr(7): kmul(krat(Fr(1, 14)), kadd(krat(13), kmul(krat(3), SQ3)))}
for D, r in ROT.items():
    assert knorm2(r) == K1 and knorm2(ksub(K1, r)) == krat(Fr(1) / D)


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
                    if j <= i:
                        continue
                    if abs(abs(z - zs[j]) - 1) > 1e-7:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        out.append((i, j))
    return out


def about(rho, piv, q):
    return kadd(piv, kmul(rho, ksub(q, piv)))


def scan(pts, zs, name, budget=20000):
    E = edges_of(pts, zs)
    n = len(pts)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** {name}: {n} pts {len(E)} edges NOT 4-COLOURABLE ***",
              flush=True)
        sv.delete()
        return "chi5", (pts, E)
    cand = []
    for i in range(n):
        for j in range(i + 1, n):
            v = abs(zs[i] - zs[j]) ** 2
            if v > 40:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_over(D, KCL):
                continue
            cand.append((i, j, D))
    hits, hard = [], []
    for i, j, D in cand:
        sv.conf_budget(budget)
        r = sv.solve_limited(assumptions=[1 + i * 4, -(1 + j * 4)])
        if r is False:
            hits.append((i, j, D))
            print(f"  *** {name}: FORCED PAIR {i},{j} at D = {D} ***  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
        elif r is None:
            hard.append((i, j, D))
    for i, j, D in hard:
        if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
            hits.append((i, j, D))
            print(f"  *** {name}: FORCED PAIR (full run) {i},{j} at D = {D} "
                  f"***  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    return ("forced" if hits else None), (n, len(E), len(cand), len(hard))


tag, info = scan(Sa, zsA, "Sa")
print(f"  Sa alone: {info}  [{time.time()-t0:.0f}s]", flush=True)

for D in (Fr(3), Fr(7)):
    ring = Counter()
    fd = float(D)
    for i in range(len(Sa)):
        for j in range(len(Sa)):
            if i != j and abs(abs(zsA[i] - zsA[j]) ** 2 - fd) < 1e-7:
                ring[i] += 1
    # Sa is invariant under the 12-element dihedral group, so pivots in one
    # orbit give isomorphic unions -- a rotation g conjugates rho_p to
    # rho_{gp}, a reflection conjugates it to the inverse, and both rotations
    # are scanned anyway.  361 pivots become about thirty.
    pos = {p: i for i, p in enumerate(Sa)}
    reps, done = [], set()
    for i, _ in ring.most_common():
        if i in done:
            continue
        orb, q = set(), Sa[i]
        for base in (q, kconj(q)):
            r = base
            for _ in range(6):
                if r in pos:
                    orb.add(pos[r])
                r = kmul(Z6, r)
        done |= orb
        reps.append((i, ring[i], len(orb)))
    order = [(i, sz) for i, sz, _ in reps]
    print(f"\nD = {D}: {len(ring)} pivots in {len(order)} dihedral orbits, "
          f"largest ring {order[0][1] if order else 0}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for rho in (ROT[D], kconj(ROT[D])):
        for k, (pi_, sz) in enumerate(order):
            piv = Sa[pi_]
            pts, seen2 = list(Sa), set(Sa)
            for q in Sa:
                z = about(rho, piv, q)
                if z not in seen2:
                    seen2.add(z)
                    pts.append(z)
            zs = [zof(q) for q in pts]
            tag, info = scan(pts, zs, f"D={D} p={pi_} r={sz}")
            if tag:
                import pickle
                with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432"
                          "-5848-a506-39c59179b415/scratchpad/kforced.pkl",
                          "wb") as fh:
                    pickle.dump((tag, [flat(q) for q in pts]), fh)
                sys.exit(0)
            print(f"  p={pi_} r={sz}: {info}  [{time.time()-t0:.0f}s]",
                  flush=True)
