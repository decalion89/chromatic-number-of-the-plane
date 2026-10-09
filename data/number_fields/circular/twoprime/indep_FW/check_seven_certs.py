"""Referee check (own code) of the eight value-2/7 certificates of Proposition 9 (prop:window7) as used by the
draft's "Type 7 does not occur" step.

Conventions derived from scratch (not taken from the author's checker):
  * a complex number z = x + iy is a pair of Fractions;  <c, g> := Re(conj(c) g) = c_x g_x + c_y g_y.
  * functional (j, l, Re): f(c) = Re(conj(c) rho^j sigma^l)          -> rotation g = rho^j sigma^l
    functional (j, l, Im): f(c) = Im(conj(c) rho^j sigma^l)
                                 = Re(conj(c) * (-i) rho^j sigma^l)  -> rotation g = -i rho^j sigma^l
    (proved below numerically as an identity on random points as well).
  * strip index n = floor(f(c)); lower margin ("+") = f(c) - n, upper margin ("-") = n + 1 - f(c).
  * draft: eps = +1 for a lower margin, -1 for an upper margin; the linear part of the margin is eps * <c, g>.
Checks, for each of the 8 certificate lines:
  (1) all lambda > 0 and sum lambda = 1;
  (2) sum lambda_i eps_i g_i = 0 exactly in Q(i)   (cancellation of the linear parts);
  (3) the constant of sum lambda_i psi_i equals 2/7 (with the strip indices of the line);
  (4) at the listed 7-point c0: floors of all 30 functionals equal the index vector; all 60 margins >= 2/7;
      the three certificate margins equal exactly 2/7;
  (5) c0 = (325/7)(a+bi) with 3(a+bi) mod 7 in A'_7, A'_7 recomputed from its definition in the paper;
  (6) u_i := eps_i g_i (times rho^j v in the draft) satisfy <c0, u_i> == 2/7 mod 1 (value of the character);
      integer multiplicities n_i (proportional to lambda_i), their sum (walk length), divisibility by 7;
  (7) the index vectors coincide with the SEVEN lines of cert_K2M1.txt; also every EXTRA line there has value 1/4
      with positive multipliers summing to 1 and cancelling linear parts (so they are excluded at theta = 2/7 > 1/4).
"""
from fractions import Fraction as Fr
from math import lcm, floor
import os, random, itertools

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def cm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def conj(a):
    return (a[0], -a[1])

def cpow(z, e):
    r = (Fr(1), Fr(0))
    if e < 0:
        z = conj(z)          # |z| = 1, so z^{-1} = conj(z)
        e = -e
    for _ in range(e):
        r = cm(r, z)
    return r

RHO = (Fr(3, 5), Fr(4, 5))
SIG = (Fr(5, 13), Fr(12, 13))
assert RHO[0] ** 2 + RHO[1] ** 2 == 1 and SIG[0] ** 2 + SIG[1] ** 2 == 1
MI = (Fr(0), Fr(-1))  # -i

def rot(j, l, part):
    g = cm(cpow(RHO, j), cpow(SIG, l))
    if part == "Re":
        return g
    assert part == "Im"
    return cm(MI, g)

def ip(c, g):  # Re(conj(c) g)
    return c[0] * g[0] + c[1] * g[1]

def fval(c, j, l, part):
    """the functional computed directly from its definition (Re or Im of conj(c) rho^j sigma^l)"""
    g = cm(cpow(RHO, j), cpow(SIG, l))
    w = cm(conj(c), g)
    return w[0] if part == "Re" else w[1]

# sanity: Im(conj(c) g) == Re(conj(c) (-i g)) on random rational points
random.seed(1)
for _ in range(200):
    c = (Fr(random.randint(-10 ** 6, 10 ** 6), random.randint(1, 999)), Fr(random.randint(-10 ** 6, 10 ** 6), random.randint(1, 999)))
    j, l = random.randint(-2, 2), random.randint(-1, 1)
    for part in ("Re", "Im"):
        assert fval(c, j, l, part) == ip(c, rot(j, l, part))

def ffloor(q):
    return q.numerator // q.denominator

# A'_7 from the definition: e in F_49 with Re(conj(e) z) in {2,3,4,5} for all z in mu_8 (norm-one elements of F_49)
F49 = [(a, b) for a in range(7) for b in range(7)]
def m7(a, b):
    return ((a[0] * b[0] - a[1] * b[1]) % 7, (a[0] * b[1] + a[1] * b[0]) % 7)
mu8 = [z for z in F49 if (z[0] * z[0] + z[1] * z[1]) % 7 == 1]
assert len(mu8) == 8
A7 = sorted(e for e in F49 if all(m7(((e[0]) % 7, (-e[1]) % 7), z)[0] in (2, 3, 4, 5) for z in mu8))
assert A7 == [(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)], A7

def parse_keys(line):
    toks = line.split(":", 1)[1].split()
    out = []
    for t in toks:
        j, l, part = t.strip("()").split(",")
        out.append((int(j), int(l), part))
    return out

sev = open(os.path.join(D, "cert_K2M1_seven.txt")).read().splitlines()
main = open(os.path.join(D, "cert_K2M1.txt")).read().splitlines()
keys = parse_keys(sev[1])
assert keys == parse_keys(main[1]) and len(keys) == 30
assert sorted(set((j, l) for j, l, _ in keys)) == [(j, l) for j in range(-2, 3) for l in range(-1, 2)]

def parse_cert(s):
    terms = []
    for t in s.split(";"):
        j, l, part, sgn, lam = t.split()
        terms.append(((int(j), int(l), part), sgn, Fr(lam)))
    return terms

def check_cert(vec, terms, value):
    n = dict(zip(keys, vec))
    lams = [lam for _, _, lam in terms]
    assert all(lam > 0 for lam in lams), "nonpositive multiplier"
    assert sum(lams) == 1
    lin = (Fr(0), Fr(0))
    const = Fr(0)
    us = []
    for k, sgn, lam in terms:
        g = rot(*k)
        eps = 1 if sgn == "+" else -1
        lin = (lin[0] + lam * eps * g[0], lin[1] + lam * eps * g[1])
        const += lam * (-n[k] if sgn == "+" else n[k] + 1)
        us.append((eps, g))
    assert lin == (0, 0), lin
    assert const == value, const
    return us, lams

report = []
seven_vecs = {}
for ln in sev[2:]:
    head, label, cert = [s.strip() for s in ln.split("|")]
    vec = tuple(int(t) for t in head.split()[2:])
    assert len(vec) == 30
    x, y = (Fr(t) for t in label.split()[1:3])
    c0 = (x, y)
    # (5) the point is a type-7 point
    a, b = x * 7 / 325, y * 7 / 325
    assert a.denominator == 1 and b.denominator == 1
    a, b = int(a), int(b)
    assert (3 * a % 7, 3 * b % 7) in A7, (a, b)
    # (4) floors and margins at the 7-point, all 30 functionals and their negatives (60 rotations i^a g)
    n = dict(zip(keys, vec))
    for k in keys:
        v = fval(c0, *k)
        fl = ffloor(v)
        assert fl == n[k], (k, v, n[k])
        assert v - fl >= Fr(2, 7) and fl + 1 - v >= Fr(2, 7)
    # the other 30 rotations (-g) have values -v: margins are the same pair swapped
    assert cert.startswith("kappa = 2/7 ; cert ")
    terms = parse_cert(cert[len("kappa = 2/7 ; cert "):])
    us, lams = check_cert(vec, terms, Fr(2, 7))
    for (k, sgn, lam) in terms:
        v = fval(c0, *k)
        fl = ffloor(v)
        marg = v - fl if sgn == "+" else fl + 1 - v
        assert marg == Fr(2, 7), (k, sgn, marg)
    # (6) tight vectors u = eps g: value of the character Re(conj(c0) u) == 2/7 mod 1
    for eps, g in us:
        val = eps * ip(c0, g)
        assert (val - Fr(2, 7)).denominator == 1, val
    L = lcm(*(lam.denominator for lam in lams))
    mult = [int(lam * L) for lam in lams]
    tot = (sum(m * e * g[0] for m, (e, g) in zip(mult, us)), sum(m * e * g[1] for m, (e, g) in zip(mult, us)))
    assert tot == (0, 0)
    assert sum(mult) % 7 == 0
    seven_vecs[vec] = c0
    report.append((f"(325/7)({a}+{b}i)", [f"{j},{l},{p},{s},{lam}" for (j, l, p), s, lam in terms], mult, sum(mult)))

# (7) cross-check with cert_K2M1.txt
main_seven = {}
n_extra = 0
n_main = 0
for ln in main[2:]:
    if not ln.startswith("component"):
        continue
    head, rest = ln.split("|", 1)
    vec = tuple(int(t) for t in head.split()[2:])
    rest = rest.strip()
    if rest.startswith("SEVEN"):
        main_seven[vec] = tuple(Fr(t) for t in rest.split()[1:3])
    elif rest.startswith("EXTRA"):
        n_extra += 1
        assert rest.startswith("EXTRA kappa = 1/4 ; cert ")
        terms = parse_cert(rest[len("EXTRA kappa = 1/4 ; cert "):])
        check_cert(vec, terms, Fr(1, 4))
    elif rest.startswith("MAIN"):
        n_main += 1
    else:
        raise ValueError(rest)
assert main_seven == seven_vecs
print("A'_7 =", A7)
print(f"cert_K2M1.txt: {len(main_seven)} SEVEN, {n_extra} EXTRA (all value 1/4, lambda>0, sum 1, linear parts cancel),"
      f" {n_main} MAIN lines; total {len(main_seven) + n_extra + n_main}")
for r in report:
    print("7-point", r[0], "| cert", r[1], "| multiplicities", r[2], "| walk length", r[3])
print("ALL CHECKS PASSED: 8 certificates, lambda > 0, sum 1, sum lambda_i eps_i gamma_i = 0 in Q(i), constant 2/7,"
      " tight margins exactly 2/7 at the 7-points, all 60 margins >= 2/7 there, walk lengths divisible by 7")
