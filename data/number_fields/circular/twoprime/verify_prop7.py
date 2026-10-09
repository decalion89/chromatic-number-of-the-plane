"""Independent exact verifier for prop7_certificates.txt (no code shared with probe2d.py / lpexact.py).

Rebuilds the nine rotations gamma_{j,l} = ((2+i)/(2-i))^j ((3+2i)/(3-2i))^l (j, l in {-1, 0, 1}) from Gaussian
integers, and checks:
  (1) every listed main component: the type point T + 65 m has exactly the listed strip indices at r = 7/25
      (T = 65(1+i)/2 or (65/3)(a + bi));
  (2) every listed extra component: the certificate lam >= 0 on three constraints
        '+' : t <= f(c) - n,     '-' : t <= n + 1 - f(c)
      has sum lam = 1, cancelling linear parts, and constant sum lam*const = 2/7, so that no point of the
      component has all margins > 2/7 (kappa <= 2/7); and the listed 7-point has all margins >= 2/7 and the
      listed indices (kappa = 2/7);
  (3) the 13 index vectors are distinct, 5 main and 8 extra.
Completeness (that these 13 are all the components) is the content of the two enumerations probe2d.py and
indep2d.py, which agree."""
from fractions import Fraction as Fr


def gdiv(a, b):
    # (a0 + a1 i) / (b0 + b1 i)
    n = b[0] * b[0] + b[1] * b[1]
    return (Fr(a[0] * b[0] + a[1] * b[1], n), Fr(a[1] * b[0] - a[0] * b[1], n))


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z = gdiv((1, 0), z)
        e = -e
    for _ in range(e):
        w = gmul(w, z)
    return w


RHO = gdiv((2, 1), (2, -1))
SIG = gdiv((3, 2), (3, -2))
assert RHO == (Fr(3, 5), Fr(4, 5)) and SIG == (Fr(5, 13), Fr(12, 13))


def functional(j, l, e):
    g = gmul(gpow(RHO, j), gpow(SIG, l))
    # conj(c) g with c = c1 + i c2: (c1 - i c2)(g0 + i g1) = (c1 g0 + c2 g1) + i (c1 g1 - c2 g0)
    return (g[0], g[1]) if e == 'Re' else (g[1], -g[0])


def fl(q):
    return q.numerator // q.denominator


r = Fr(7, 25)
lines = open("prop7_certificates.txt").read().splitlines()
funs = [tuple(x.strip('()').split(',')) for x in lines[1].split(':', 1)[1].split()]
funs = [(int(j), int(l), e) for (j, l, e) in funs]
F = {key: functional(*key) for key in funs}
vecs = set()
nmain = nextra = 0
for ln in lines[2:]:
    if not ln.startswith("component"):
        continue
    head, rest = ln.split("|", 1)
    iv = [int(x) for x in head.split("indices")[1].split()]
    assert tuple(iv) not in vecs
    vecs.add(tuple(iv))
    n_of = dict(zip(funs, iv))
    rest = rest.strip()
    if rest.startswith("MAIN"):
        parts = rest.split()
        name, mx, my = parts[1], int(parts[4]), int(parts[5])
        T = (Fr(65, 2), Fr(65, 2)) if name == 'C' else (Fr(65 * int(name[1]), 3), Fr(65 * int(name[2]), 3))
        T = (T[0] + 65 * mx, T[1] + 65 * my)
        for key in funs:
            v = F[key][0] * T[0] + F[key][1] * T[1]
            n = n_of[key]
            assert n + r <= v <= n + 1 - r, (name, key)
        nmain += 1
        continue
    assert rest.startswith("EXTRA")
    fields = [x.strip() for x in rest.split(";")]
    assert fields[0] == "EXTRA kappa = 2/7"
    sp = fields[2].split("(")[1].rstrip(")").split(",")
    seven = (Fr(sp[0].strip()), Fr(sp[1].strip()))
    certs = [fields[3][len("cert "):]] + fields[4:]
    lin = [Fr(0), Fr(0)]
    const = Fr(0)
    tot = Fr(0)
    for c in certs:
        j, l, e, sgn, lam = c.split()
        key = (int(j), int(l), e)
        lam = Fr(lam)
        assert lam >= 0
        f = F[key]
        n = n_of[key]
        if sgn == '+':
            lin[0] += lam * f[0]; lin[1] += lam * f[1]; const += lam * (-n)
        else:
            lin[0] -= lam * f[0]; lin[1] -= lam * f[1]; const += lam * (n + 1)
        tot += lam
    assert tot == 1 and lin == [0, 0] and const == Fr(2, 7), (tot, lin, const)
    # the 7-point: same indices, all margins >= 2/7, and it is a 7-torsion point N(a+bi)/7
    assert (7 * seven[0] / 65).denominator == 1 and (7 * seven[1] / 65).denominator == 1
    mn = Fr(1)
    for key in funs:
        v = F[key][0] * seven[0] + F[key][1] * seven[1]
        n = n_of[key]
        mn = min(mn, v - n, n + 1 - v)
    assert mn == Fr(2, 7), mn
    nextra += 1
assert (nmain, nextra) == (5, 8)
print(f"verified: {nmain} main components (type points with the listed indices), {nextra} extra components with "
      f"kappa = 2/7 exactly (dual certificates and 7-torsion points)")
