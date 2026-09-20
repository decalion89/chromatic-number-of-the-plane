"""Build the graph, do not just multiply its directions.

The claim to check: G u rho_p(G), rho a rotation of irrational chord living in
Q(m)(sqrt3, sqrt5, sqrt7, sqrt11), is 5-chromatic AND blocks at 2, 3, 4 and 5.
The first half is free -- it contains G -- but it is verified anyway, and the
second half is redone on the graph's own edge set, cross edges included, with
the directions put in a Z-basis of the lattice they generate rather than in the
ambient Z^d.

m is a root of x^3 - 10x^2 + 26x - 11, irreducible mod 5, so 5 has residue
degree 3 in Q(m) -- the condition the barrier theorem says a sixth colour
needs and no multiquadratic field can meet.
"""
import sys, time, pickle
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.homcol import has_homomorphism
from pysat.solvers import Solver
exec(open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/hnfcheck.py").read()
     .split("print(\"re-checking")[0].split("t0 = time.time()")[1])

t0 = time.time()
CUB = (Fr(11), Fr(-26), Fr(10))


def cz():
    return (F.zero(), F.zero(), F.zero())


def crat(q):
    return (F.rational(Fr(q)), F.zero(), F.zero())


def cemb(e):
    return (e, F.zero(), F.zero())


def cm():
    return (F.zero(), F.rational(1), F.zero())


def cadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def csub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cmul(a, b):
    r = [F.zero()] * 5
    for i in range(3):
        for j in range(3):
            if a[i] == F.zero() or b[j] == F.zero():
                continue
            r[i + j] = r[i + j] + a[i] * b[j]
    for d in (4, 3):
        c = r[d]
        if c == F.zero():
            continue
        r[d] = F.zero()
        for k in range(3):
            r[d - 3 + k] = r[d - 3 + k] + c * F.rational(CUB[k])
    return tuple(r[:3])


def cinv(x):
    """Solve for y with x.y = 1 by linear algebra over F."""
    cols = [cmul(x, e) for e in ((F.rational(1), F.zero(), F.zero()),
                                 (F.zero(), F.rational(1), F.zero()),
                                 (F.zero(), F.zero(), F.rational(1)))]
    A = [[cols[j][i] for j in range(3)] + [F.rational(1 if i == 0 else 0)]
         for i in range(3)]
    for c in range(3):
        p = next(t for t in range(c, 3) if A[t][c] != F.zero())
        A[c], A[p] = A[p], A[c]
        sc = A[c][c].inverse()
        A[c] = [v * sc for v in A[c]]
        for t in range(3):
            if t != c and A[t][c] != F.zero():
                f = A[t][c]
                A[t] = [u - f * v for u, v in zip(A[t], A[c])]
    return (A[0][3], A[1][3], A[2][3])


def cfloat(x):
    M = 8.3166247903554            # a real root of the cubic, to 13 places
    for _ in range(80):
        M = M - (M ** 3 - 10 * M ** 2 + 26 * M - 11) / (3 * M ** 2 - 20 * M + 26)
    return float(x[0]) + float(x[1]) * M + float(x[2]) * M * M


t = cm()
t2 = cmul(t, t)
den = cinv(cadd(crat(1), t2))
COS = cmul(csub(crat(1), t2), den)
SIN = cmul(cmul(crat(2), t), den)
assert cadd(cmul(COS, COS), cmul(SIN, SIN)) == crat(1), "not a rotation"
print(f"rotation verified: cos^2 + sin^2 = 1, cos ~ {cfloat(COS):.6f}, "
      f"sin ~ {cfloat(SIN):.6f}  [{time.time()-t0:.0f}s]", flush=True)
chord = csub(crat(2), cmul(crat(2), COS))
print(f"  chord = {cfloat(chord):.6f}, rational? "
      f"{chord[1] == F.zero() and chord[2] == F.zero()}", flush=True)

P = build_G(F, as_graph=False)
pts = [(cemb(p.x), cemb(p.y)) for p in P]
piv = pts[0]
rot = []
for x, y in pts:
    dx, dy = csub(x, piv[0]), csub(y, piv[1])
    rot.append((cadd(piv[0], csub(cmul(COS, dx), cmul(SIN, dy))),
                cadd(piv[1], cadd(cmul(SIN, dx), cmul(COS, dy)))))
allp, seen = [], set()
for q in pts + rot:
    if q not in seen:
        seen.add(q)
        allp.append(q)
print(f"G u rho_p(G): {len(allp)} points (G has {len(pts)})  "
      f"[{time.time()-t0:.0f}s]", flush=True)

zs = [(cfloat(x), cfloat(y)) for x, y in allp]
cell = {}
for i, (a, b) in enumerate(zs):
    cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
cands = []
for i, (a, b) in enumerate(zs):
    cx, cy = int(a // 1), int(b // 1)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            for j in cell.get((cx + da, cy + db), ()):
                if j > i and abs((a - zs[j][0]) ** 2
                                 + (b - zs[j][1]) ** 2 - 1) < 1e-7:
                    cands.append((i, j))
print(f"  {len(cands)} candidate edges by float  [{time.time()-t0:.0f}s]",
      flush=True)
E, vecs = [], set()
for i, j in cands:
    dx = csub(allp[j][0], allp[i][0])
    dy = csub(allp[j][1], allp[i][1])
    if cadd(cmul(dx, dx), cmul(dy, dy)) == crat(1):
        E.append((i, j))
        vecs.add(dx + dy)
        vecs.add(tuple(-c for c in dx) + tuple(-c for c in dy))
print(f"  {len(E)} edges confirmed exactly, {len(vecs)} signed directions  "
      f"[{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(allp))]
for a, b in E:
    for c in range(4):
        cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
with Solver(name="cd19", bootstrap_with=cls) as sv:
    four = sv.solve()
cls5 = [[1 + v * 5 + c for c in range(5)] for v in range(len(allp))]
for a, b in E:
    for c in range(5):
        cls5.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
with Solver(name="cd19", bootstrap_with=cls5) as sv:
    five = sv.solve()
print(f"  4-colourable {four}, 5-colourable {five}  -> chi = "
      f"{5 if (not four and five) else '?'}  [{time.time()-t0:.0f}s]",
      flush=True)

rows = []
for v in vecs:
    r = []
    for comp in v:
        r.extend(Fr(c) for c in comp.c)
    rows.append(r)
dn = 1
for r in rows:
    for q in r:
        dn = dn * q.denominator // gcd(dn, q.denominator)
ints = {tuple(int(q * dn) for q in r) for r in rows}
content = 0
for v in ints:
    for x in v:
        content = gcd(content, abs(x))
if content > 1:
    ints = {tuple(x // content for x in v) for v in ints}
ints = sorted(ints)
print(f"  {len(ints)} integer directions, dim {len(ints[0])}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
red, rank = in_basis(ints, len(ints[0]))
print(f"  lattice rank {rank}, {len(red)} distinct in the Z-basis  "
      f"[{time.time()-t0:.0f}s]", flush=True)
for name, VV in (("Z^d", ints), ("the lattice", red)):
    bl = [n for n in (2, 3, 4, 5) if has_homomorphism(VV, n)[0] is None]
    print(f"  over {name}: blocks at {bl}  [{time.time()-t0:.0f}s]",
          flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/blockG.pkl", "wb") as fh:
    pickle.dump({"points": len(allp), "edges": E, "dirs": red}, fh)
