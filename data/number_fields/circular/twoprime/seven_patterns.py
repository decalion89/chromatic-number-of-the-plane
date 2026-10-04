"""Digit patterns of the 7-adic characters of A = Z[1/5][i] (exact).  For eta with 7*eta in E7 (8 classes), the
character a -> Re(conj(eta) a~) (a~ = a mod 7) has values V_j = conj(eta) rho^j mod Z[i] (rho^j taken mod 7, a
norm-one element); write V_j = h + zeta_j with zeta_j in [-3/14, 3/14]^2 and digits nu_j = zeta_j - rho zeta_{j-1}.
Asserts that the digits lie in Lambda_+ = ((2+i)/5) Z[i] and that the values have period 8, and prints:
  * whether no 7-pattern has a zero digit (False: zeros alternate with non-zero digits);
  * whether some 7-pattern contains a q-triple (nu, -i nu, -nu) (False);
  * whether three consecutive digits determine the pair (eta, position mod 8) (False, and not needed: the same
    triple occurs for (eta, j) and (eta rho^-1, j - 1));
  * that the 8 patterns are the 8 shifts of the pattern S of (2,3), that three consecutive digits of S determine the
    shift mod 8, and that every three consecutive digits of S contain a zero and a non-zero digit.
The last three checks were in the run that produced seven_patterns.txt but missing from the stored copy of this
program; they were written again when the folder was assembled, and the output is unchanged.  The referee's
indep_F/f49_and_digits.py checks the same facts and that two consecutive digits determine the pattern among all 13
(c, 4 q, 8 seven)."""
from fractions import Fraction as Fr
p = 7
E7 = [(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)]
inv5 = pow(5, -1, p)
rho7 = ((3 * inv5) % p, (4 * inv5) % p)
def mul7(x, y): return ((x[0]*y[0] - x[1]*y[1]) % p, (x[0]*y[1] + x[1]*y[0]) % p)
RHO = (Fr(3, 5), Fr(4, 5))
def cmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def sub(a, b): return (a[0]-b[0], a[1]-b[1])
def pattern(e, J=16):
    zs = []
    r = (1, 0)
    for j in range(J):
        # conj(eta) * rho^j with eta = e/7: Re = (e0 r0 + e1 r1)/7, Im = (e0 r1 - e1 r0)/7  (mod 1)
        re = Fr((e[0]*r[0] + e[1]*r[1]) % p, p); im = Fr((e[0]*r[1] - e[1]*r[0]) % p, p)
        zs.append((re - Fr(1, 2), im - Fr(1, 2)))
        r = mul7(r, rho7)
    digits = [sub(zs[j], cmul(RHO, zs[j-1])) for j in range(1, J)]
    return zs, digits
def in_lambda(nu):
    # nu in ((2+i)/5) Z[i]  <=>  nu * 5/(2+i) = nu (2-i) in Z[i]
    w = cmul(nu, (Fr(2), Fr(-1)))
    return w[0].denominator == 1 and w[1].denominator == 1
pats = {}
for e in E7:
    zs, ds = pattern(e)
    assert all(abs(z[0]) <= Fr(3, 14) and abs(z[1]) <= Fr(3, 14) for z in zs)
    assert all(in_lambda(d) for d in ds)
    assert zs[8] == zs[0]
    pats[e] = ds
    print(e, "digits (x5):", [(int(5*d[0]), int(5*d[1])) for d in ds[:8]])
zero = (Fr(0), Fr(0))
nonzero = all(d != zero for ds in pats.values() for d in ds)
print("no zero digit in any 7-pattern:", nonzero)
def is_qtriple(a, b, c):
    mi = lambda z: (z[1], -z[0])            # -i z
    return a != zero and b == mi(a) and c == mi(b)
qt = any(is_qtriple(ds[j], ds[j+1], ds[j+2]) for ds in pats.values() for j in range(len(ds) - 2))
print("some 7-pattern contains a q-triple:", qt)
triples = {}
ok = True
for e, ds in pats.items():
    for j in range(len(ds) - 2):
        key = (ds[j], ds[j+1], ds[j+2])
        # the character seen from position j: the class conj-shift; record (e, j mod 8)
        if key in triples and triples[key] != (e, j % 8):
            ok = False
        triples.setdefault(key, (e, j % 8))
print("three consecutive digits determine (eta, position mod 8):", ok, "; distinct triples:", len(triples))
S = pats[(2, 3)]
shifts = {e: [t for t in range(8) if all(pats[e][j] == S[(j + t) % 8] for j in range(8))] for e in E7}
print("each 7-pattern is a shift of S (pattern of (2,3)):", shifts,
      all(len(v) == 1 for v in shifts.values()) and sorted(v[0] for v in shifts.values()) == list(range(8)))
tri = [(S[t], S[(t + 1) % 8], S[(t + 2) % 8]) for t in range(8)]
print("three consecutive digits of S determine the shift mod 8:", len(set(tri)) == 8)
print("every 3 consecutive digits of S contain a zero and a non-zero digit:",
      all(any(d == zero for d in x) and any(d != zero for d in x) for x in tri))
