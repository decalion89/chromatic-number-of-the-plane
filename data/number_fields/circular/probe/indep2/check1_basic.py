# Referee check 1: congruences and identities of Lemma 8, Lemma 10/11 arithmetic,
# the 12-gon P and the hexagon Z((1+i)/3) of Lemma 9, the 5-adic points c_k of Remark (2).
from fractions import Fraction as Fr
import math, itertools
from gauss import G, I, RHO, H, modZ, mod1, dist_int, frac_center

ok_all = True
def check(cond, msg):
    global ok_all
    if not cond:
        ok_all = False
        print("FAIL:", msg)

# ---------- Lemma 8(i): G_N = { i^a rho^j : |j| <= k } ----------
for k in range(1, 4):
    N = 5 ** k
    sols = set()
    for a in range(-N, N + 1):
        b2 = N * N - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            for bb in {b, -b}:
                sols.add(G(Fr(a, N), Fr(bb, N)))
    claimed = set()
    for j in range(-k, k + 1):
        for a in range(4):
            claimed.add((I ** a) * (RHO ** j))
    check(sols == claimed, f"G_N k={k}")
    print(f"Lemma 8(i) k={k}: |G_N|={len(sols)} equals {{i^a rho^j}}: {sols == claimed}")

# ---------- Lemma 8(ii), (iii) ----------
for k in range(1, 13):
    N = 5 ** k
    for j in range(-k, k + 1):
        Nr = (RHO ** j) * N
        check(Nr.is_gint(), "N rho^j integral")
        # (ii)
        v = (H * N).conj() * (RHO ** j)
        check(modZ(v - H) == G(0), f"(ii) k={k} j={j}")
        # (iii) N rho^j = (-1)^k (-i)^j mod 3
        target = ((-1) ** k) * ((-I) ** j)
        d = Nr - target
        check(d.is_gint() and d.re % 3 == 0 and d.im % 3 == 0, f"(iii) cong k={k} j={j}")
print("Lemma 8(ii),(iii) congruences checked for k<=12, |j|<=k")

def eta(eps, k, j):
    """eta_j(eps) at level N=5^k: conj(N eps) rho^j in h + eta/6 + Z[i]."""
    N = 5 ** k
    v = (eps * N).conj() * (RHO ** j)
    w = modZ(v - H) * 6
    check(w.re in (1, -1) and w.im in (1, -1), f"eta not in {{+-1+-i}} k={k} j={j}")
    return w

EQ = [G(Fr(a, 3), Fr(b, 3)) for a in (1, 2) for b in (1, 2)]
for k in range(1, 9):
    for e in EQ:
        for j in range(-k + 2, k + 1):
            check(eta(e, k, j - 2) == -eta(e, k, j), "eta_{j-2} = -eta_j")
        # rotation: eta(i eps) = -i eta(eps)
        for j in range(-k, k + 1):
            check(eta(I * e, k, j) == (-I) * eta(e, k, j), "eta(i eps) = -i eta(eps)")
            # conjugation: eta_{-j}(conj eps) = conj eta_j(eps)
            check(eta(e.conj(), k, -j) == eta(e, k, j).conj(), "eta conj")
        # level-N etas for |j|<=1 are those of +-eps at level 5
        sgn = (-1) ** (k - 1)
        for j in (-1, 0, 1):
            check(eta(e, k, j) == eta(e * sgn, 1, j), "eta level N vs level 5")
e0 = G(Fr(1, 3), Fr(1, 3))
print("eta at level 5 for (1+i)/3: j=-1,0,1 ->", [eta(e0, 1, j) for j in (-1, 0, 1)])
check([eta(e0, 1, j) for j in (-1, 0, 1)] == [G(1, 1), G(1, -1), G(-1, -1)], "eta table (1+i)/3")
# multiplication by i permutes the four classes of E_q cyclically
cls = lambda z: (mod1(z.re), mod1(z.im))
orb = [cls(e0)]
z = e0
for _ in range(4):
    z = I * z
    orb.append(cls(z))
print("orbit of (1+i)/3 under i (classes mod Z[i]):", orb)
check(len(set(orb[:4])) == 4 and orb[4] == orb[0], "E_q cyclic")

# ---------- Lambda_+ / Z[i] classes and nu_pm ----------
Lp = G(2, 1) / 5
reps = set()
for t in range(5):
    reps.add(modZ(Lp * t))
print("Lambda_+/Z[i] classes (centred reps):", sorted(reps, key=lambda g: (g.re, g.im)))
for k in range(1, 7):
    N = 5 ** k
    for a in range(5):
        for b in range(5):
            lam = G(a, b)
            nup = lam.conj() * (RHO ** k) * Fr(N, 5)
            num = lam.conj() * (RHO ** (-k)) * Fr(N, 5)
            check((nup * G(2, -1)).is_gint(), "nu_+ in Lambda_+")   # Lambda_+ = Z[i]/(2-i)
            check((num * G(2, 1)).is_gint(), "nu_- in Lambda_-")
            d_p = (lam / G(2, 1)).is_gint()
            d_m = (lam / G(2, -1)).is_gint()
            check(nup.is_gint() == d_p, "nu_+ integral iff (2+i)|lam")
            check(num.is_gint() == d_m, "nu_- integral iff (2-i)|lam")
print("nu_pm(lambda) facts checked for k<=6, all lambda mod 5")

# ---------- identities of the base case and of Lemma 10/11 ----------
def zr(z): return z * RHO
def zri(z): return z * (RHO ** -1)
for z in (G(1), G(0, 1)):
    v1 = 7 * zr(z).re + 25 * zri(z).re - 24 * zr(z).im
    v2 = 25 * zr(z).im - 24 * zri(z).re + 7 * zri(z).im
    check(v1 == 0 and v2 == 0, "base-case identities")
print("7(zr)_1+25(zr^-1)_1-24(zr)_2 = 0 and 25(zr)_2-24(zr^-1)_1+7(zr^-1)_2 = 0: verified on R-basis")
check(RHO ** -2 == G(Fr(-7, 25), Fr(-24, 25)), "rho^-2")
# 7*(7u1-24u2)+24*(24u1+7u2) = 625 u1, 7*25+24*25 = 775
check(7 * 7 + 24 * 24 == 625 and -7 * 24 + 24 * 7 == 0 and 7 * 25 + 24 * 25 == 775, "625/775")
# lambda orbits, nu reps
for lam, want_p, want_m in [(G(1), G(Fr(-2, 5), Fr(-1, 5)), G(Fr(-2, 5), Fr(1, 5))),
                           (G(2), G(Fr(1, 5), Fr(-2, 5)), G(Fr(1, 5), Fr(2, 5))),
                           (G(2, 2), G(Fr(-1, 5), Fr(2, 5)), G(Fr(-2, 5), Fr(1, 5)))]:
    nup = modZ(lam.conj() * RHO); num = modZ(lam.conj() * RHO ** -1)
    check(nup == want_p and num == want_m, f"nu reps lam={lam}")
    print(f"lambda={lam}: nu_+ = {nup}, nu_- = {num}")
# orbits of the 16 classes under multiplication by i
classes16 = []
for a in range(5):
    for b in range(5):
        lam = G(a, b)
        if not (lam / G(2, 1)).is_gint() and not (lam / G(2, -1)).is_gint():
            classes16.append((a, b))
def red5(z): return (int(z.re) % 5, int(z.im) % 5)
orbits = []
seen = set()
for c in classes16:
    if c in seen: continue
    o = []
    z = G(*c)
    for _ in range(4):
        o.append(red5(z)); seen.add(red5(z)); z = I * z
    orbits.append(o)
print("16 classes, orbits under i:", orbits)
reps_given = [(1, 0), (2, 0), (2, 2), (1, 1)]
check(len(classes16) == 16 and len(orbits) == 4, "16 classes / 4 orbits")
check(all(any(r in o for o in orbits) for r in reps_given) and
      len({next(idx for idx, o in enumerate(orbits) if r in o) for r in reps_given}) == 4, "reps hit 4 orbits")
# lambda = 1+i: 5h + lambda - 5*(2+2i)/3 = (1+i)/6
check(H * 5 + G(1, 1) - G(Fr(10, 3), Fr(10, 3)) == G(Fr(1, 6), Fr(1, 6)), "lambda=1+i shift")
# type-q lift: alpha' - 5 alpha = 0 mod 3 with alpha = 2 alpha' mod 3
for ap in (1, 2):
    al = (2 * ap) % 3
    check((ap - 5 * al) % 3 == 0 and al in (1, 2), "alpha lift")
# (2+i)^(2k+1) = 2+i mod 5
for k in range(0, 60):
    d = G(2, 1) ** (2 * k + 1) - G(2, 1)
    check(d.re % 5 == 0 and d.im % 5 == 0, f"(2+i)^(2k+1) k={k}")
print("(2+i)^(2k+1) = 2+i mod 5 for k<60")
# numerical constants
s0 = Fr(11, 56)
check(s0 * (1 + math.sqrt(2)) < 0.6, "s(1+sqrt2)<3/5")
check(Fr(1, 2) - Fr(17, 56) == Fr(11, 56) and Fr(1, 3) - Fr(17, 56) == Fr(5, 168), "17/56 <-> 11/56, 5/168")
check(Fr(1, 6) + Fr(5, 168) == Fr(11, 56), "s = 1/6 + delta at threshold")
print("sqrt2*(11/56+1/6) =", math.sqrt(2) * (11 / 56 + 1 / 6), "< 0.52;  1/6+0.52 =", 1 / 6 + 0.52, "< 1-11/56 =", 1 - 11 / 56)

# ---------- polygons P and Z ----------
def vertices(halfplanes):
    """halfplanes: list of (a,b,c) meaning a*x+b*y <= c (Fractions). Exact vertex enumeration."""
    pts = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(halfplanes, 2):
        det = a1 * b2 - a2 * b1
        if det == 0: continue
        x = (c1 * b2 - c2 * b1) / det
        y = (a1 * c2 - a2 * c1) / det
        if all(a * x + b * y <= c for a, b, c in halfplanes):
            pts.add((x, y))
    cx = sum(p[0] for p in pts) / len(pts); cy = sum(p[1] for p in pts) / len(pts)
    return sorted(pts, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))

def form(w_of_x):
    """w_of_x: function x(G)->G (R-linear). Return (re-coeffs, im-coeffs) as ((a,b),(a,b)) in x=(x1,x2)."""
    w1 = w_of_x(G(1)); w2 = w_of_x(G(0, 1))
    return (w1.re, w2.re), (w1.im, w2.im)

# P_1 with s=1: |(conj(x) rho^j)_e| <= 1, j in {-1,0,1}
hp = []
for j in (-1, 0, 1):
    fr, fi = form(lambda x, j=j: x.conj() * (RHO ** j))
    for (a, b) in (fr, fi):
        hp.append((a, b, Fr(1))); hp.append((-a, -b, Fr(1)))
VP = vertices(hp)
claimed = set()
for s1 in (1, -1):
    for s2 in (1, -1):
        claimed |= {(Fr(s1), Fr(s2, 3)), (Fr(s1, 3), Fr(s2)), (Fr(5 * s1, 7), Fr(5 * s2, 7))}
print("P: %d vertices" % len(VP), [(str(x), str(y)) for x, y in VP])
check(set(VP) == claimed, "12-gon P vertices")
print("P vertices equal the claimed 12:", set(VP) == claimed, "; max |x|^2 =", max(x * x + y * y for x, y in VP))
# check the example: x = 1 + i/3
x = G(1, Fr(1, 3))
print("x=1+i/3: conj(x)rho =", x.conj() * RHO, " conj(x)rho^-1 =", x.conj() * RHO ** -1)
check(x.conj() * RHO == G(Fr(13, 15), Fr(3, 5)) and x.conj() * RHO ** -1 == G(Fr(1, 3), -1), "example vertex")

# Z((1+i)/3) at level 5: eta_{j,e} (conj(y) rho^j)_e <= 1
hz = []
for j in (-1, 0, 1):
    et = eta(e0, 1, j)
    fr, fi = form(lambda y, j=j: y.conj() * (RHO ** j))
    hz.append((et.re * fr[0], et.re * fr[1], Fr(1)))
    hz.append((et.im * fi[0], et.im * fi[1], Fr(1)))
VZ = vertices(hz)
claimedZ = {(Fr(-5, 4), Fr(0)), (Fr(0), Fr(-5, 4)), (Fr(-5, 7), Fr(-5, 7)), (Fr(1), Fr(-1, 2)), (Fr(1), Fr(1)), (Fr(-1, 2), Fr(1))}
print("Z: %d vertices" % len(VZ), [(str(x), str(y)) for x, y in VZ])
check(set(VZ) == claimedZ, "hexagon Z vertices")
print("Z vertices equal the claimed 6:", set(VZ) == claimedZ, "; max |y|^2 =", max(x * x + y * y for x, y in VZ))
# Z(i eps) = i Z(eps)
for e in EQ:
    hze = []
    for j in (-1, 0, 1):
        et = eta(e, 1, j)
        fr, fi = form(lambda y, j=j: y.conj() * (RHO ** j))
        hze.append((et.re * fr[0], et.re * fr[1], Fr(1)))
        hze.append((et.im * fi[0], et.im * fi[1], Fr(1)))
    V = vertices(hze)
    # find rotation a with V = i^a Z(e0)
    found = None
    for a in range(4):
        rot = {((I ** a) * G(x, y)).re for x, y in VZ}
        rotset = {(((I ** a) * G(x, y)).re, ((I ** a) * G(x, y)).im) for x, y in VZ}
        if rotset == set(V):
            found = a
    print(f"Z({e}) = i^{found} Z((1+i)/3)")
    check(found is not None, "Z rotation")

# ---------- Remark (2): the 5-adic points c_k ----------
def in_S(c, k, theta):
    """exact: min over gamma in G_N of ||Re(conj(c) gamma)||, and membership in S_N^theta."""
    m = None
    for j in range(-k, k + 1):
        w = c.conj() * (RHO ** j)
        for t in (w.re, w.im):
            d = dist_int(t)
            m = d if m is None or d < m else m
    return m >= theta, m
for k in range(1, 16):
    N = 5 ** k
    ck = H * N + (RHO ** k) / 5 - G(2, -1) * (5 ** (k - 1))
    inside, m = in_S(ck, k, Fr(3, 10))
    check(inside and m == Fr(3, 10), f"c_k k={k}")
    # j<k: conj(c_k) rho^j = h + rho^{j-k}/5 mod Z[i]; j=k: h-(1+i)/5
    for j in range(-k, k):
        check(modZ(ck.conj() * RHO ** j - H - (RHO ** (j - k)) / 5) == G(0), "c_k pattern j<k")
    check(modZ(ck.conj() * RHO ** k - H + G(1, 1) / 5) == G(0), "c_k pattern j=k")
    # distances: to N E_c and N E_q (exact min over nearby lattice points)
    dc = None
    base = ck - H * N
    for mu_r in range(-3, 4):
        for mu_i in range(-3, 4):
            # nearest multiples: shift by rounding first
            r0 = round(base.re / N); i0 = round(base.im / N)
            d = base - G(N * (r0 + mu_r), N * (i0 + mu_i))
            dc = d.norm() if dc is None or d.norm() < dc else dc
    dq = None
    for e in EQ:
        b2 = ck - e * N
        r0 = round(b2.re / N); i0 = round(b2.im / N)
        for mu_r in range(-2, 3):
            for mu_i in range(-2, 3):
                d = b2 - G(N * (r0 + mu_r), N * (i0 + mu_i))
                dq = d.norm() if dq is None or d.norm() < dq else dq
    lbc = math.sqrt(5) * 5 ** (k - 1) - 0.2
    lbq = 7 / 6 * 5 ** (k - 1) - 0.2
    check(math.sqrt(dc) >= lbc - 1e-9 and math.sqrt(dq) >= lbq - 1e-9, f"c_k distances k={k}")
    if k <= 4 or k == 15:
        print(f"c_{k}: in S^(3/10), min ||.|| = {m}; dist to N E_c = {math.sqrt(dc):.4f} >= {lbc:.4f}; to N E_q = {math.sqrt(dq):.4f} >= {lbq:.4f}")
print("c_k checked exactly for k <= 15")
print("ALL OK" if ok_all else "SOME CHECK FAILED")
