"""Compare the index vectors of the components found by vertex_enum (file 1) and by lift_rows.py (file 2)."""
import sys
from fractions import Fraction as Fr
ve, lr = sys.argv[1], sys.argv[2]
N = None; pts = []
for line in open(ve):
    if line.startswith('# N'):
        N = int(line.split()[2])
    elif not line.startswith('#') and line.strip():
        X, Y, D = map(int, line.split()); pts.append((Fr(N * X, D), Fr(N * Y, D)))
rots = [(A, B) for A in range(-N, N + 1) for B in range(-N, N + 1) if A * A + B * B == N * N]
funs = set()
for (A, B) in rots:
    for (a, b) in ((A, B), (B, -A)):
        funs.add((-a, -b) if (a < 0 or (a == 0 and b < 0)) else (a, b))
funs = sorted(funs)
fl = lambda q: q.numerator // q.denominator
V1 = set(tuple(fl((f[0] * p[0] + f[1] * p[1]) / N) for f in funs) for p in pts)
V2 = set(tuple(map(int, l.split())) for l in open(lr) if l.strip())
print(f"vertex enumeration: {len(V1)} components; row lifting: {len(V2)} components; identical sets: {V1 == V2}")
