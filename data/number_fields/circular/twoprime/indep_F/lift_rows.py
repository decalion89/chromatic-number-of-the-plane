"""Referee program, second method (exact polygons, lifting by rows).

S^r(k,1) = S^r(k,0)  intersected with the conditions of the rows l = +1 and l = -1 (rotations rho^j sigma^{+-1}).
Step 1: S^r(k,0) modulo M = 5^k, cell by cell: the square [a+r, a+1-r] x [b+r, b+1-r] (a, b mod M) is split by the
        strips of the functionals of the rotations rho^j (|j| <= k).
Step 2: every piece is translated by M*(u + v i), 0 <= u, v < 13, and split by the strips of the functionals of
        the rotations rho^j sigma^{+-1} (|j| <= k).  The survivors are the components of S^r(k,1) mod 13M.
Exact arithmetic (Fractions); polygons in vertex form, clipped by closed half-planes (degenerate pieces kept).
Output: number of components and their index vectors (in the functional order of analyze_components.py), written
to a file for comparison with the vertex enumeration.
usage: python3 lift_rows.py k rn rd outfile"""
import sys, time
from fractions import Fraction as Fr

k, rn, rd, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
r = Fr(rn, rd)
M = 5 ** k
N = 13 * M
rots = [(A, B) for A in range(-N, N + 1) for B in range(-N, N + 1) if A * A + B * B == N * N]


def canon(a, b):
    return (-a, -b) if (a < 0 or (a == 0 and b < 0)) else (a, b)


row0, row1 = [], []
for (A, B) in rots:
    for (a, b) in ((A, B), (B, -A)):
        f = canon(a, b)
        tgt = row0 if (A % 13 == 0 and B % 13 == 0) else row1   # 13 | A, B  <=>  M*gamma in Z[i]  (row l = 0)
        if f not in row0 and f not in row1:
            tgt.append(f)
allf = sorted(row0 + row1)
print(f"k = {k}, N = {N}, r = {r}: {len(rots)} rotations; functionals: row 0: {len(row0)}, rows +-1: {len(row1)}")
g1 = [f for f in row0 if f in ((N, 0), (0, N))]
assert len(g1) == 2
rest0 = [f for f in row0 if f not in g1]


def fl(q):
    return q.numerator // q.denominator


def clip(poly, a, b, c):
    """{ a x + b y >= c } intersected with the convex polygon (list of vertices)."""
    s = [a * p[0] + b * p[1] - c for p in poly]
    n = len(poly)
    res = []
    for i in range(n):
        p, sp = poly[i], s[i]
        q, sq = poly[(i + 1) % n], s[(i + 1) % n]
        if sp >= 0:
            res.append(p)
        if (sp > 0 and sq < 0) or (sp < 0 and sq > 0):
            t = sp / (sp - sq)
            res.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    out = []
    for p in res:
        if not out or out[-1] != p:
            out.append(p)
    while len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def split(poly, f):
    """pieces of poly in the strips  n + r <= (A x + B y)/N <= n + 1 - r, with their n."""
    A, B = f
    vals = [Fr(A * p[0] + B * p[1], N) for p in poly]
    vmin, vmax = min(vals), max(vals)
    res = []
    for n in range(fl(vmin - (1 - r)) - 1, fl(vmax - r) + 2):
        P = clip(poly, Fr(A, N), Fr(B, N), n + r)
        if P:
            P = clip(P, Fr(-A, N), Fr(-B, N), -(n + 1 - r))
        if P:
            res.append((n, P))
    return res


t0 = time.time()
# step 1: S^r(k,0) modulo M
base = []
for a in range(M):
    for b in range(M):
        sq = [(a + r, b + r), (a + 1 - r, b + r), (a + 1 - r, b + 1 - r), (a + r, b + 1 - r)]
        pieces = [sq]
        for f in rest0:
            nxt = []
            for P in pieces:
                nxt.extend(P2 for (_, P2) in split(P, f))
            pieces = nxt
            if not pieces:
                break
        base.extend(pieces)
print(f"step 1: S^r({k},0) mod {M}: {len(base)} pieces ({time.time() - t0:.1f}s)")
# step 2
final = []
for P in base:
    for u in range(13):
        for v in range(13):
            Q = [(p[0] + M * u, p[1] + M * v) for p in P]
            pieces = [Q]
            for f in row1:
                nxt = []
                for P1 in pieces:
                    nxt.extend(P2 for (_, P2) in split(P1, f))
                pieces = nxt
                if not pieces:
                    break
            final.extend(pieces)
print(f"step 2: S^r({k},1) mod {N}: {len(final)} pieces ({time.time() - t0:.1f}s)")


def ivec(p):
    return tuple(fl(Fr(f[0] * p[0] + f[1] * p[1], N)) for f in allf)


vecs = set()
for P in final:
    cx = sum(p[0] for p in P) / len(P)
    cy = sum(p[1] for p in P) / len(P)
    vecs.add(ivec((cx, cy)))
print(f"distinct index vectors: {len(vecs)}; vertex counts: {sorted(len(P) for P in final)}")
with open(out, 'w') as fh:
    for v in sorted(vecs):
        fh.write(' '.join(map(str, v)) + '\n')
