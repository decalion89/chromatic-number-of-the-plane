# Referee check 7: (a) direct computation of S_N^theta for N = 25, 125 (all 625 / 15625 unit cells, all
# conditions |j| <= k applied to the square h + m + B_s) compared with the lifted computation;
# (b) theta just above 17/56: only main components at levels 1..10;
# (c) the main type-c polygon equals N h + P_k exactly (vertex sets), k <= 4.
import itertools, math
from fractions import Fraction as Fr
from gauss import G, RHO, H
from levels import level1, lift, cut, classify, HALF, run

def canon(P):
    return frozenset(P)

def direct(theta, k):
    s = HALF - theta; N = 5 ** k
    out = []
    for m1 in range(N):
        for m2 in range(N):
            cx, cy = HALF + m1, HALF + m2
            pieces = [[(cx - s, cy - s), (cx + s, cy - s), (cx + s, cy + s), (cx - s, cy + s)]]
            for j in range(-k, k + 1):
                pieces = [Q for P in pieces for Q in cut(P, j, s)]
            out.extend(pieces)
    return out

for theta in (Fr(3001611, 10000000), Fr(3036, 10001)):
    s = HALF - theta
    polys = level1(s)
    for k in (1, 2, 3):
        if k > 1:
            polys = lift(polys, k, s)
        d = direct(theta, k)
        print(f"theta={theta} k={k}: lifted {len(polys)} polygons, direct {len(d)} polygons, identical sets: {set(map(canon, polys)) == set(map(canon, d))}")

theta = Fr(17, 56) + Fr(1, 10 ** 9)
out, _ = run(theta, 10, verbose=False)
print("theta = 17/56 + 1e-9:", [(k, n, c['X']) for k, n, c, dg in out], "(level, polygons, extras)")

# (c) main polygon = N h + P_k
def Pk_vertices(k, s):
    hp = []
    for j in range(-k, k + 1):
        r = RHO ** j
        # conj(x) r = (x1 - i x2)(p + i q) = (x1 p + x2 q) + i (x1 q - x2 p)
        for (a, b) in ((r.re, r.im), (r.im, -r.re)):
            hp.append((a, b, s)); hp.append((-a, -b, s))
    pts = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(hp, 2):
        det = a1 * b2 - a2 * b1
        if det == 0: continue
        x = (c1 * b2 - c2 * b1) / det; y = (a1 * c2 - a2 * c1) / det
        if all(a * x + b * y <= c for a, b, c in hp):
            pts.add((x, y))
    return pts
for theta in (Fr(1, 3) - Fr(1, 100), Fr(3001611, 10000000)):
    s = HALF - theta
    polys = level1(s)
    for k in range(1, 5):
        if k > 1:
            polys = lift(polys, k, s)
        N = 5 ** k
        C = [P for P in polys if classify(P, k, s) == 'C']
        assert len(C) == 1
        shifted = {(v[0] - N * HALF - N * round((v[0] - N * HALF) / N), v[1] - N * HALF - N * round((v[1] - N * HALF) / N)) for v in C[0]}
        print(f"theta={theta} k={k}: type-c polygon - N h has {len(shifted)} vertices; equals P_k: {shifted == Pk_vertices(k, s)}")
