"""Is the 3-colouring saved by a stalled growth run a function of the point modulo mM (M = lattice of the
directions)? usage: colstruct.py STATE.json [mmax]"""
import sys, json
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
s = json.load(open(sys.argv[1])); mmax = int(sys.argv[2]) if len(sys.argv) > 2 else 40
U = [tuple(u) for u in s["U"]]; P = [tuple(p) for p in s["P"]]; col = s["col"]
H = hermite_normal_form(Matrix(U).T)
B = Matrix.hstack(*[H[:, j] for j in range(H.shape[1]) if any(H[:, j])])
Binv = B.inv()
p0 = Matrix(P[0])
C = []
for p in P:
    v = Binv * (Matrix(p) - p0)
    assert all(x.is_integer for x in v), "point not in p0 + M"
    C.append(tuple(int(x) for x in v))
print(f"{len(P)} points, {len(U)} directions, rank {B.shape[1]}; colours used {sorted(set(col))}")
for m in range(2, mmax + 1):
    f, ok = {}, True
    for c, k in zip(C, col):
        key = tuple(x % m for x in c)
        if f.setdefault(key, k) != k:
            ok = False; break
    if ok:
        print(f"  the colouring is a function of the point modulo {m}M ({len(f)} classes used)")
        break
else:
    print(f"  not a function of the point modulo mM for any m <= {mmax}")
