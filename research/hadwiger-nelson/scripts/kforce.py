"""The test that actually decides the last step.

Reading de Grey's final move backwards: G = Ya u Yb is two copies of Y turned
about p = (-2,0) by pi/2 +- arcsin(1/8), and the angle between them closes
distance 4 into 1.  For that to be a contradiction the two copies must carry
the SAME colour at the two images of one vertex, so what Y has to supply is
not a vague sphere condition but a pair:

    p, q at distance d, monochromatic in EVERY proper 4-colouring.

Given such a pair, and a rotation rho about p closing d (2d sin(theta/2) = 1),
the graph H u rho_p(H) is 5-chromatic outright: rho fixes p, so both copies
read c(q) = c(p), while q and rho(q) are one apart.  No sphere, no case split.

So the search is: over all pairs at a distance K can close, is any pair forced
equal?  Fixing c(a) = 0 is free (colours are interchangeable), so the query
"is there a 4-colouring with c(a) = 0 and c(b) != 0" is a single assumption
pair against one incremental solver, and UNSAT is the prize.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from pysat.solvers import Solver

t0 = time.time()


def squarefree(n):
    if n == 0:
        return 0
    s, d = (-1 if n < 0 else 1), abs(n)
    f = 2
    while f * f <= d:
        e = 0
        while d % f == 0:
            d //= f
            e += 1
        if e % 2:
            s *= f
        f += 1
    return s * d


def closable(D):
    r = Fr(1) - 4 * Fr(D)
    if r == 0:
        return False
    return squarefree(r.numerator * r.denominator) in (1, -3, -11, 33)


rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
X = close(seed, [r for _, r in rots])
print(f"X: {len(X)} points  [{time.time()-t0:.0f}s]", flush=True)


def edges_of(pts):
    zs = [zof(q) for q in pts]
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
    return out, zs


def scan(pts, name):
    E, zs = edges_of(pts)
    n = len(pts)
    print(f"{name}: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]",
          flush=True)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** {name} IS ALREADY 5-CHROMATIC ***", flush=True)
        return E, None
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            v = abs(zs[i] - zs[j]) ** 2
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable(D):
                continue
            pairs.append((i, j, D))
    print(f"  {len(pairs)} pairs at a closable distance  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    hits = []
    for k, (i, j, D) in enumerate(pairs):
        if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
            hits.append((i, j, D))
            print(f"  *** FORCED EQUAL: {i},{j} at D = {D}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
        if k and k % 500 == 0:
            print(f"  ... {k}/{len(pairs)}  [{time.time()-t0:.0f}s]",
                  flush=True)
    if not hits:
        print(f"  no forced pair in {name}  [{time.time()-t0:.0f}s]",
              flush=True)
    sv.delete()
    return E, hits


scan(X, "X")
