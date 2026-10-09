from fractions import Fraction as Fr
from nf import Field, unit_from, lconj
from limit_test import feasible

def quad(d):
    return Field([[-d, 0, 1]])

for d, ns in [(11, [1, 19]), (23, [1]), (35, [1, 55]), (7, [1, 2, 3, 5]), (3, [1, 2, 4]), (47, [1]), (59, [1, 19])]:
    Fd = quad(d)
    sq = Fd.gen(0)
    V = [(Fd.one(), {})]
    for n in ns:
        p = Fd.scal(Fr(n, d), sq)          # p = n / sqrt(d)
        u = unit_from(Fd, p)
        V.append(u)
        if n == 1:
            V.append(lconj(Fd, u))
    print(d, ns, feasible(Fd, V, want=True), flush=True)
