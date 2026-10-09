"""q3_11.py [CERT.json.gz]: 27 explicit unit vectors of Q(sqrt3, sqrt11)^2 with kappa = 1/4 (exact, integers and
fractions only).

U = { R60^j RA^k RG^l (1, 0) : j, k, l in {-1, 0, 1} }, where R60 is the rotation by 60 degrees, RA the rotation
with cosine 5/6 and sine sqrt11/6 (the angle of the Moser spindle) and RG the rotation with cosine 11/14 and sine
5 sqrt3/14.  These 27 vectors are pairwise distinct up to sign; an element of the field is written over the basis
(1, sqrt3, sqrt11, sqrt33), and a vector as 8 integers over the common denominator 84.

Prints: whether the vectors are those of the certificate (up to sign), and the least distance to Z of the character
of theta311.json, x = ((a0 + a1 s3 + a2 s11 + a3 s33)/84, (b0 + ...)/84) -> theta . (a0, ..., b3) mod 1, over the 27
vectors: exactly 1/4, so kappa(U) >= 1/4.  With the certificate (no character maps U into the open interval
(1/4, 3/4); check_open.py, check_open_indep.py), kappa(U) = 1/4."""
import gzip, json, os, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def fmul(a, b):
    """product in Q(sqrt3, sqrt11) over the basis (1, s3, s11, s33): s3 s11 = s33, s3 s33 = 3 s11, s11 s33 = 11 s3"""
    return (a[0] * b[0] + 3 * a[1] * b[1] + 11 * a[2] * b[2] + 33 * a[3] * b[3],
            a[0] * b[1] + a[1] * b[0] + 11 * (a[2] * b[3] + a[3] * b[2]),
            a[0] * b[2] + a[2] * b[0] + 3 * (a[1] * b[3] + a[3] * b[1]),
            a[0] * b[3] + a[3] * b[0] + a[1] * b[2] + a[2] * b[1])


def fadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def fneg(a):
    return tuple(-x for x in a)


def rotate(c, s, v, k):
    """apply the rotation (cos c, sin s) k times (k may be negative) to the vector v = (x, y)."""
    if k < 0:
        s, k = fneg(s), -k
    x, y = v
    for _ in range(k):
        x, y = fadd(fmul(c, x), fneg(fmul(s, y))), fadd(fmul(s, x), fmul(c, y))
    return x, y


ONE = (Fr(1), Fr(0), Fr(0), Fr(0))
ZERO = (Fr(0),) * 4
R60 = ((Fr(1, 2), Fr(0), Fr(0), Fr(0)), (Fr(0), Fr(1, 2), Fr(0), Fr(0)))          # cos 1/2, sin sqrt3/2
RA = ((Fr(5, 6), Fr(0), Fr(0), Fr(0)), (Fr(0), Fr(0), Fr(1, 6), Fr(0)))           # cos 5/6, sin sqrt11/6
RG = ((Fr(11, 14), Fr(0), Fr(0), Fr(0)), (Fr(0), Fr(5, 14), Fr(0), Fr(0)))        # cos 11/14, sin 5 sqrt3/14
D = 84


def units_311():
    out = []
    for j in (-1, 0, 1):
        for k in (-1, 0, 1):
            for l in (-1, 0, 1):
                v = rotate(*R60, (ONE, ZERO), j)
                v = rotate(*RA, v, k)
                v = rotate(*RG, v, l)
                n = fadd(fmul(v[0], v[0]), fmul(v[1], v[1]))
                assert n == ONE, "not a unit vector"
                w = [x * D for x in v[0] + v[1]]
                assert all(x.denominator == 1 for x in w)
                out.append([int(x) for x in w])
    pairs = {frozenset({tuple(u), tuple(-x for x in u)}) for u in out}
    assert len(pairs) == 27, "two vectors are equal or opposite"
    return out


def check(cert_path=None):
    U = units_311()
    print(f"27 unit vectors of Q(sqrt3, sqrt11)^2 (denominator {D}), pairwise distinct up to sign")
    if cert_path is None:
        return
    with (gzip.open(cert_path, "rt") if cert_path.endswith(".gz") else open(cert_path)) as fh:
        C = json.load(fh)
    same = ({frozenset({tuple(u), tuple(-x for x in u)}) for u in C["units"]} ==
            {frozenset({tuple(u), tuple(-x for x in u)}) for u in U}) and C["D"] == D and C["d"] == "3,11"
    print(f"the certificate's units are these vectors up to sign: {same}")
    T = json.load(open(os.path.join(HERE, "theta311.json")))
    th = [Fr(x) for x in T["theta"]]
    assert len(th) == 8 and T["D"] == D

    def dist(t):
        t = t - (t.numerator // t.denominator)
        return min(t, 1 - t)
    m = min(dist(sum(a * b for a, b in zip(th, u))) for u in U)
    print(f"theta311 = {[str(x) for x in th]}: least distance to Z over the 27 vectors: {m}")
    return same and m == Fr(1, 4)


if __name__ == "__main__":
    ok = check(sys.argv[1] if len(sys.argv) > 1 else None)
    sys.exit(0 if ok in (None, True) else 1)
