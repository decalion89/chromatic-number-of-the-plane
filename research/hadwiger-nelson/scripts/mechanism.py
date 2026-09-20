"""What exactly does Y force?

de Grey's last step is G = Ya u Yb, two copies of Y turned about p = (-2, 0)
by pi/2 +- arcsin(1/8).  The angle between the copies is 2 arcsin(1/8), so a
point at distance 4 from p has its two images at distance 2*4*(1/8) = 1.
That is the whole mechanism, and it says what Y must be forcing:

    in every proper 4-colouring of Y, some vertex at distance 4 from p
    carries p's own colour.

If that is right the test is one SAT call: add, for every v with |v - p| = 4,
the constraint c(v) != c(p), and ask for a 4-colouring.  UNSAT means Y forces.

This matters because the forcing distance is what has to be ported to another
field, and 4 is not sqrt3.  The chord a rotation needs to close a pair at
distance d is c = 1/d^2, and c is available over F exactly when 3c(4 - c) is
a square there.  For d = 4 that is 189/256, needing sqrt21 -- de Grey has it
(sqrt3, sqrt7), K does not.  For d = sqrt3 it is 11/3, needing sqrt33 -- K
has it, de Grey does not.  So the two fields force at different distances,
and the port has to follow the arithmetic, not the picture.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_Y, build_Sa
from hn.geometry import DEGREY_FIELD, Point
from pysat.solvers import Solver

t0 = time.time()
F = DEGREY_FIELD
Y = build_Y(F)
idx = {p: i for i, p in enumerate(Y)}
print(f"Y: {len(Y)} vertices  [{time.time()-t0:.0f}s]", flush=True)

zs = [(float(p.x), float(p.y)) for p in Y]
cell = {}
for i, (a, b) in enumerate(zs):
    cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
E = []
for i, (a, b) in enumerate(zs):
    cx, cy = int(a // 1), int(b // 1)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            for j in cell.get((cx + da, cy + db), ()):
                if j <= i:
                    continue
                if abs((a - zs[j][0]) ** 2 + (b - zs[j][1]) ** 2 - 1) > 1e-7:
                    continue
                d = Y[i] - Y[j]
                if d.x * d.x + d.y * d.y == F.rational(1):
                    E.append((i, j))
print(f"  {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

P = Point(F.rational(-2), F.zero())
assert P in idx, "the pivot is not a vertex of Y"
pi_ = idx[P]
print(f"  pivot (-2,0) is vertex {pi_}", flush=True)

base = [[1 + v * 4 + c for c in range(4)] for v in range(len(Y))]
for a, b in E:
    for c in range(4):
        base.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
with Solver(name="cd19", bootstrap_with=base) as sv:
    print(f"  Y alone 4-colourable? {sv.solve()}  [{time.time()-t0:.0f}s]",
          flush=True)

# every squared distance from the pivot that Y actually realises, found by
# float and then confirmed exactly
from collections import Counter
cand = Counter()
for q in Y:
    d = q - P
    v = float(d.x) ** 2 + float(d.y) ** 2
    r = Fr(round(v * 144), 144)
    if abs(float(r) - v) < 1e-9:
        cand[r] += 1
rat = []
for r, n in sorted(cand.items()):
    ring = [i for i, q in enumerate(Y)
            if (q - P).x ** 2 + (q - P).y ** 2 == F.rational(r)]
    if ring:
        rat.append((r, len(ring)))
print(f"  {len(rat)} rational squared distances from the pivot: "
      f"{[(str(r), n) for r, n in rat]}", flush=True)

for target, _n in rat:
    ring = [i for i, q in enumerate(Y)
            if (q - P).x ** 2 + (q - P).y ** 2 == F.rational(target)]
    cls = list(base)
    for v in ring:
        for c in range(4):                      # c(v) != c(p)
            cls.append([-(1 + v * 4 + c), -(1 + pi_ * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        ok = sv.solve()
    verdict = "FORCES" if not ok else "does not force"
    print(f"  d^2 = {target}: {len(ring)} vertices on the sphere -- "
          f"{verdict}  [{time.time()-t0:.0f}s]", flush=True)
