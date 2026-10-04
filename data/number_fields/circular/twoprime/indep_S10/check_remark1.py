"""Remark (1) of Section 10: the type-7 index vectors of the (2,1)-window, restricted to G(1,1) and translated by
65 Z[i], are among the 8 type-7 (EXTRA) vectors of the (1,1)-window, whose dual certificates give least margin <= 2/7."""
from fractions import Fraction as Fr


def cdiv(x, y):
    n = y[0] * y[0] + y[1] * y[1]
    return (Fr(x[0] * y[0] + x[1] * y[1], n), Fr(x[1] * y[0] - x[0] * y[1], n))


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z, e = cdiv((1, 0), z), -e
    for _ in range(e):
        w = cmul(w, z)
    return w


RHO, SIG = cdiv((2, 1), (2, -1)), cdiv((3, 2), (3, -2))


def fun(j, l, e):
    g = cmul(cpow(RHO, j), cpow(SIG, l))
    return (g[0], g[1]) if e == "Re" else (g[1], -g[0])


def load(fn):
    lines = open(fn).read().splitlines()
    fl_line = [ln for ln in lines if ln.startswith("# functionals")][0]
    keys = [tuple(x.strip("()").split(",")) for x in fl_line.split(":", 1)[1].split()]
    keys = [(int(j), int(l), e) for (j, l, e) in keys]
    out = []
    for ln in lines:
        if ln.startswith("component"):
            head, lab = ln.split("|", 1)
            v = [int(x) for x in head.split("indices")[1].split()]
            out.append((dict(zip(keys, v)), lab.strip().split()[0]))
    return out


w11 = load("run/prop7_certificates.txt")
w21 = load("run/cert_K2M1.txt")
K11 = [(j, l, e) for j in (-1, 0, 1) for l in (-1, 0, 1) for e in ("Re", "Im")]
base11 = {tuple(d[k] for k in K11): lab for d, lab in w11}
for d, lab in w21:
    a, b = d[(0, 0, "Re")], -d[(0, 0, "Im")] - 1
    m = (Fr(65 * (a // 65)), Fr(65 * (b // 65)))      # translate by -65 m
    v = []
    for k in K11:
        f = fun(*k)
        s = f[0] * m[0] + f[1] * m[1]
        assert s.denominator == 1
        v.append(d[k] - int(s))
    print(lab, "->", base11.get(tuple(v), "not among the 13 (1,1)-vectors"))
