"""
Helper module of check_ie_indep.py: exact arithmetic, the embedded 23-point set, rigorous enclosures
of J0, and the individual checks.  See check_ie_indep.py for what is checked and which classical
facts are assumed.  Standard library only, except numpy for the exact check over all 2^23 subsets.
"""
import sys, os, json, time, math, random
from fractions import Fraction as F
from math import isqrt, gcd

def log(*a):
    print(*a, flush=True)

# =====================================================================================
# 1. Exact arithmetic in K = Q(sqrt3, sqrt11), basis (1, sqrt3, sqrt11, sqrt33)
#    (1, sqrt3, sqrt11, sqrt33 are linearly independent over Q, so equality is coefficientwise)
# =====================================================================================
KZERO = (F(0), F(0), F(0), F(0))
KONE = (F(1), F(0), F(0), F(0))

def kadd(x, y): return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])
def ksub(x, y): return (x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3])
def kmul(x, y):
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return (a1 * a2 + 3 * b1 * b2 + 11 * c1 * c2 + 33 * d1 * d2,      # 1
            a1 * b2 + b1 * a2 + 11 * (c1 * d2 + d1 * c2),              # sqrt3   (sqrt11 sqrt33 = 11 sqrt3)
            a1 * c2 + c1 * a2 + 3 * (b1 * d2 + d1 * b2),               # sqrt11  (sqrt3 sqrt33 = 3 sqrt11)
            a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2)                     # sqrt33  (sqrt3 sqrt11 = sqrt33)
def kstr(x):
    names = ['', '*sqrt3', '*sqrt11', '*sqrt33']
    parts = [str(c) + nm for c, nm in zip(x, names) if c != 0]
    return ' + '.join(parts) if parts else '0'
def csub(z, w): return (ksub(z[0], w[0]), ksub(z[1], w[1]))
def cabs2(z): return kadd(kmul(z[0], z[0]), kmul(z[1], z[1]))

# =====================================================================================
# 2. The 23 points of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki, in the order of their list X.
#    Row: real part (c_1, c_sqrt3, c_sqrt11, c_sqrt33), imaginary part (c_1, c_sqrt3, c_sqrt11, c_sqrt33).
# =====================================================================================
_POINTS = [
    ('0', '0', '0', '0',            '0', '0', '0', '0'),       # 1: 0
    ('1', '0', '0', '0',            '0', '0', '0', '0'),       # 2: 1
    ('1/2', '0', '0', '0',          '0', '1/2', '0', '0'),     # 3: 1/2 + (sqrt3/2) i
    ('3/2', '0', '0', '0',          '0', '1/2', '0', '0'),     # 4: 3/2 + (sqrt3/2) i
    ('5/6', '0', '0', '0',          '0', '0', '1/6', '0'),     # 5: 5/6 + (sqrt11/6) i
    ('5/12', '0', '0', '-1/12',     '0', '5/12', '1/12', '0'),
    ('5/4', '0', '0', '-1/12',      '0', '5/12', '1/4', '0'),
    ('13/12', '0', '0', '1/12',     '0', '1/12', '-1/12', '0'),
    ('1/4', '0', '0', '1/12',       '0', '1/12', '-1/4', '0'),
    ('13/12', '0', '0', '-1/12',    '0', '-1/12', '-1/12', '0'),
    ('2/3', '0', '0', '0',          '0', '1/2', '-1/6', '0'),
    ('7/12', '0', '0', '-1/12',     '0', '5/12', '-1/12', '0'),
    ('5/6', '0', '0', '-1/6',       '0', '-1/6', '1/6', '0'),
    ('7/12', '0', '0', '1/12',      '0', '7/12', '-1/12', '0'),
    ('1/12', '0', '0', '-1/12',     '0', '-1/12', '-1/12', '0'),
    ('2/3', '0', '0', '-1/6',       '0', '1/3', '-1/6', '0'),
    ('19/12', '0', '0', '-1/12',    '0', '5/12', '-1/12', '0'),
    ('1/3', '0', '0', '0',          '0', '1/2', '1/6', '0'),
    ('3/4', '0', '0', '-1/12',      '0', '-1/12', '1/4', '0'),
    ('1', '0', '0', '-1/6',         '0', '-1/6', '0', '0'),
    ('7/6', '0', '0', '0',          '0', '0', '-1/6', '0'),
    ('4/3', '0', '0', '-1/6',       '0', '1/3', '1/6', '0'),
    ('19/12', '0', '0', '-1/4',     '0', '1/4', '-1/12', '0'),  # 23: -sqrt33/4 + 19/12 + (sqrt3/4 - sqrt11/12) i
]
POINTS = [((F(r[0]), F(r[1]), F(r[2]), F(r[3])), (F(r[4]), F(r[5]), F(r[6]), F(r[7]))) for r in _POINTS]
NN = len(POINTS)

# The published multiset of the 253 squared distances: (a, d, multiplicity) for a + d*sqrt33.
_TALLY = [
    ('1', '0', 47), ('3', '0', 9), ('3/2', '1/6', 15), ('3/2', '-1/6', 18), ('7/6', '-1/6', 16),
    ('2', '-1/3', 16), ('5/3', '0', 15), ('1/3', '0', 23), ('10/3', '-1/3', 3), ('7/6', '1/6', 11),
    ('29/6', '-5/6', 1), ('4/3', '0', 11), ('8/3', '-1/3', 8), ('5/2', '-1/6', 11), ('7/2', '-1/2', 7),
    ('10/3', '1/3', 2), ('2', '1/3', 8), ('5/2', '1/6', 8), ('7/3', '0', 8), ('5', '0', 1),
    ('8/3', '1/3', 3), ('4', '-1/3', 3), ('9/2', '-1/2', 3), ('16/3', '-1/3', 1), ('35/6', '-5/6', 1),
    ('19/6', '-1/2', 3), ('5', '-2/3', 1),
]

class Geometry:
    """Exact squared distances, the unit-distance graph and the distinct distances (with assertions)."""
    def __init__(self):
        self.sq = {}
        for i in range(NN):
            for j in range(i + 1, NN):
                self.sq[(i + 1, j + 1)] = cabs2(csub(POINTS[i], POINTS[j]))
        self.edges = sorted(p for p, s in self.sq.items() if s == KONE)
        self.distinct = {}
        for s in self.sq.values():
            self.distinct[s] = self.distinct.get(s, 0) + 1
        tally = {(F(a), F(0), F(0), F(d)): c for a, d, c in _TALLY}
        assert all(s != KZERO for s in self.sq.values()), 'two points coincide'
        assert len(self.edges) == 47, 'transcription check: expected 47 unit pairs'
        assert len(self.distinct) == 27, 'transcription check: expected 27 distinct distances'
        assert self.distinct == tally, 'transcription check: distance multiset differs from the published one'
    def d2(self, i, j):
        return self.sq[(i, j)] if i < j else self.sq[(j, i)]

# =====================================================================================
# 3. Fixed-point interval arithmetic on integers: a number N stands for N / 2^PB; every operation
#    rounds outward, so [lo, hi] always encloses the exact value.
# =====================================================================================
PB = 200
ONE = 1 << PB
ONE2 = 1 << (2 * PB)
TINY = 1 << (PB - 160)

def cdiv(a, b):                 # ceil(a/b) for b > 0
    return -((-a) // b)
def fx_floor(q): return (q.numerator << PB) // q.denominator
def fx_ceil(q): return cdiv(q.numerator << PB, q.denominator)
def fx_to_F(N): return F(N, ONE)
def imul(alo, ahi, blo, bhi):
    p1, p2, p3, p4 = alo * blo, alo * bhi, ahi * blo, ahi * bhi
    return min(p1, p2, p3, p4) >> PB, -((-max(p1, p2, p3, p4)) >> PB)
def iqmul_lo(q, lo, hi):        # lower bound of q*[lo,hi] for a rational q
    p, d = q.numerator, q.denominator
    return (p * lo) // d if p >= 0 else (p * hi) // d
def iqmul_hi(q, lo, hi):
    p, d = q.numerator, q.denominator
    return cdiv(p * hi, d) if p >= 0 else cdiv(p * lo, d)
def isqrt_lo(N): return isqrt(N << PB)
def isqrt_hi(N): return isqrt(N << PB) + 1
def inv_sqrt_up(q):             # rational upper bound of 1/sqrt(q), q > 0
    return F(ONE, isqrt((q.numerator << (2 * PB)) // q.denominator))
def sqrt_up(q):                 # rational upper bound of sqrt(q), q >= 0
    return F(isqrt(cdiv(q.numerator << (2 * PB), q.denominator)) + 1, ONE)

# pi from Machin's formula pi = 16 atan(1/5) - 4 atan(1/239), alternating-series bounds
def _atan_inv(m, n):
    s = F(0)
    for k in range(n):
        s += F((-1) ** k, (2 * k + 1) * m ** (2 * k + 1))
    nxt = F(1, (2 * n + 1) * m ** (2 * n + 1))
    return (s, s + nxt) if n % 2 == 0 else (s - nxt, s)
_a5 = _atan_inv(5, 120); _a239 = _atan_inv(239, 40)
PI_LO_Q = 16 * _a5[0] - 4 * _a239[1]
PI_HI_Q = 16 * _a5[1] - 4 * _a239[0]
assert F(314159265358979323846264338327950288419716939937510, 10 ** 50) < PI_LO_Q < PI_HI_Q \
    < F(314159265358979323846264338327950288419716939937511, 10 ** 50)
PI_LO = fx_floor(PI_LO_Q); PI_HI = fx_ceil(PI_HI_Q)

SQ_ENC = {}
for _n in (3, 11, 33):
    _lo = isqrt(_n << (2 * PB)); SQ_ENC[_n] = (_lo, _lo + 1)

def k_enclose(x):
    """fixed-point enclosure of a + b sqrt3 + c sqrt11 + d sqrt33"""
    lo, hi = fx_floor(x[0]), fx_ceil(x[0])
    for c, n in ((x[1], 3), (x[2], 11), (x[3], 33)):
        if c == 0: continue
        Lo, Hi = SQ_ENC[n]
        lo += iqmul_lo(c, Lo, Hi); hi += iqmul_hi(c, Lo, Hi)
    return lo, hi

# =====================================================================================
# 4. Rigorous enclosures of J0 on [Xlo, Xhi] (fixed point, 0 <= Xlo <= Xhi)
# =====================================================================================
NEVAL = {'series': 0, 'hankel': 0}
OPTS = {'series_only': False}

def j0_series(Xlo, Xhi, extra=0):
    """J0(x) = sum_k (-1)^k t_k, t_k = y^k/(k!)^2, y = x^2/4, evaluated at the exact point Xlo with
    PB+extra bits, every t_k bracketed by directed rounding; once t_k <= 2^-160 and q = y/(k+2)^2 <= 1/2
    the tail is at most t_{k+1}/(1-q) <= 2 t_{k+1}.  Widened by (Xhi - Xlo) since |J0'| <= 1."""
    W = PB + extra
    WONE = 1 << W
    tiny = 1 << (W - 160)
    Xl = Xlo << extra
    Ylo = (Xl * Xl) >> (W + 2)
    Yhi = -((-(Xl * Xl)) >> (W + 2))
    tlo = thi = WONE
    slo = shi = WONE
    k = 0
    while True:
        k += 1
        kk = k * k
        tlo = ((tlo * Ylo) >> W) // kk
        thi = cdiv(-((-(thi * Yhi)) >> W), kk)
        if k & 1:
            slo -= thi; shi -= tlo
        else:
            slo += tlo; shi += thi
        if thi <= tiny and 2 * Yhi <= ((k + 2) * (k + 2)) << W:
            break
    tn = cdiv(-((-(thi * Yhi)) >> W), (k + 1) * (k + 1))
    tail = 2 * tn + 1
    slo -= tail; shi += tail
    NEVAL['series'] += 1
    wid = Xhi - Xlo + 1
    return (slo >> extra) - wid, -((-shi) >> extra) + wid

def _sincos_point(C):
    """sin, cos of the exact fixed-point number C (|C| < 1): Taylor polynomials with Lagrange remainder."""
    a = C if C >= 0 else -C
    neg = C < 0
    mlo = mhi = ONE
    clo = chi = ONE
    slo = shi = 0
    n = 0
    while True:
        n += 1
        mlo = ((mlo * a) >> PB) // n
        mhi = cdiv(-((-(mhi * a)) >> PB), n)
        if n & 1:
            sgn = -1 if ((n - 1) // 2) % 2 else 1
            if neg: sgn = -sgn
            if sgn > 0: slo += mlo; shi += mhi
            else: slo -= mhi; shi -= mlo
        else:
            sgn = -1 if (n // 2) % 2 else 1
            if sgn > 0: clo += mlo; chi += mhi
            else: clo -= mhi; chi -= mlo
        if mhi <= TINY:
            break
    rem = cdiv(-((-(mhi * a)) >> PB), n + 1) + 1
    return (slo - rem, shi + rem), (clo - rem, chi + rem)

def j0_hankel(Xlo, Xhi):
    """x > 30.  From Hankel's integral
         H0^(1)(x) = sqrt(2/(pi x)) e^{iw} (1/Gamma(1/2)) int_0^oo e^{-u} u^{-1/2} (1 + iu/(2x))^{-1/2} du,
    w = x - pi/4, Taylor's formula with integral remainder for (1+z)^{-1/2} (|1 + s z| >= 1 for z on
    the imaginary axis) gives  J0(x) = sqrt(2/(pi x)) (P cos w - Q sin w) + err  with
    P = g0 - g2 + g4 - ..., -Q = g1 - g3 + ..., g_k = |a_k| x^-k, a_k = (-1)^k (1*3*...*(2k-1))^2/(k! 8^k),
    terms k < n used and |err| <= sqrt(2/(pi x)) g_n."""
    Zlo = ONE2 // Xhi
    Zhi = cdiv(ONE2, Xlo)
    xf = Xlo / ONE
    glo = [ONE]; ghi = [ONE]
    k = 0
    while True:
        k += 1
        m = (2 * k - 1) ** 2
        glo.append(((glo[-1] * m * Zlo) >> PB) // (8 * k))
        ghi.append(cdiv(-((-(ghi[-1] * m * Zhi)) >> PB), 8 * k))
        if k >= 3 and (ghi[-1] <= TINY or (2 * k + 1) ** 2 > 0.9 * 8 * (k + 1) * xf):
            break
    K = k - 1
    Plo = Phi = Qlo = Qhi = 0
    for i in range(K + 1):
        if i % 2 == 0:
            if (i // 2) % 2 == 0: Plo += glo[i]; Phi += ghi[i]
            else: Plo -= ghi[i]; Phi -= glo[i]
        else:
            if ((i - 1) // 2) % 2 == 0: Qlo += glo[i]; Qhi += ghi[i]
            else: Qlo -= ghi[i]; Qhi -= glo[i]
    gnext = ghi[K + 1]
    s2lo = (1 << (3 * PB + 1)) // (PI_HI * Xhi)
    s2hi = cdiv(1 << (3 * PB + 1), PI_LO * Xlo)
    Slo, Shi = isqrt_lo(s2lo), isqrt_hi(s2hi)
    q = int(round((xf - math.pi / 4) / (math.pi / 2)))
    k4 = 2 * q + 1
    thlo = Xlo - cdiv(k4 * PI_HI, 4)
    thhi = Xhi - (k4 * PI_LO) // 4
    (sl, sh), (cl, ch) = _sincos_point(thlo)
    wid = thhi - thlo
    sl -= wid; sh += wid; cl -= wid; ch += wid
    r = q % 4
    if r == 0: cw, sw = (cl, ch), (sl, sh)
    elif r == 1: cw, sw = (-sh, -sl), (cl, ch)
    elif r == 2: cw, sw = (-ch, -cl), (-sh, -sl)
    else: cw, sw = (sl, sh), (-ch, -cl)
    a1 = imul(Plo, Phi, cw[0], cw[1])
    a2 = imul(Qlo, Qhi, sw[0], sw[1])
    jlo, jhi = imul(Slo, Shi, a1[0] + a2[0], a1[1] + a2[1])
    err = -((-(Shi * gnext)) >> PB) + 1
    NEVAL['hankel'] += 1
    return jlo - err, jhi + err

X_SWITCH = 30 << PB
def j0_enc(Xlo, Xhi):
    if Xlo <= X_SWITCH:
        return j0_series(Xlo, Xhi)
    if OPTS['series_only']:
        return j0_series(Xlo, Xhi, extra=int(1.4427 * (Xhi >> PB)) + 80)
    return j0_hankel(Xlo, Xhi)

# reference values J0(x) (50 significant digits) used only by the self-test
_J0_REF = [
    ('0.5', '0.93846980724081290422840467359971262556892679709682'),
    ('1', '0.76519768655796655144971752610266322090927428975533'),
    ('2.404825557695773', '-1.2011950073676861231254998810877307352097857488789e-16'),
    ('3.8317059702075123', '-0.40275939570255297209600218642716515386976423230379'),
    ('10', '-0.24593576445134833519776086248532875382960007282657'),
    ('29.99', '-0.087551353531462972776129474905939240650240442590076'),
    ('30.01', '-0.085176372734292647939675373469499682986945933631956'),
    ('31.7', '0.12399787757698112276972378067904618440464338251814'),
    ('45.5', '0.088176093155092095249092814321580254772781253198244'),
    ('100', '0.019985850304223122424228390950848990680633578859028'),
    ('117.06', '-0.073701057764123245202255719983504186531646259288251'),
    ('500.25', '-0.03563094809285174482593410255118398620662619468178'),
    ('999.81', '0.025233671145099068513931404562416274086772769053751'),
    ('1700.3', '-0.019276488474784176748668181677852432874982777544108'),
]

def selftest_j0():
    """Consistency tests of the J0 code (not part of the proof): tabulated values must lie in the
    enclosures, and the large-argument enclosures must overlap high-precision power-series enclosures."""
    bad = 0
    tol = F(1, 10 ** 45)
    for xs, vs in _J0_REF:
        x, v = F(xs), F(vs)
        lo, hi = j0_enc(fx_floor(x), fx_ceil(x))
        if not (F(lo, ONE) - tol <= v <= F(hi, ONE) + tol):
            bad += 1
            log('   J0 self-test: value at x = %s outside the enclosure' % xs)
    rnd = random.Random(12345)
    xs = [30 + 15 * k / 40 + 0.0123 for k in range(40)] + [rnd.uniform(30, 1800) for _ in range(20)]
    for x in xs:
        X = int(x * 2 ** 40) << (PB - 40)
        a = j0_series(X, X, extra=int(1.4427 * x) + 80)
        b = j0_hankel(X, X)
        if max(a[0], b[0]) > min(a[1], b[1]):
            bad += 1
            log('   J0 self-test: series and Hankel enclosures disagree at x = %r' % x)
    log('J0 self-test: %d tabulated values, %d series/Hankel comparisons, %d failures'
        % (len(_J0_REF), len(xs), bad))
    if bad:
        raise RuntimeError('J0 self-test failed')

def prove_j0_lower():
    """J0(x) >= -0.4028 for all x >= 0: on [0,4] a grid of step 2^-8 with min(J0(a), J0(b)) - h^2/8
    (|J0''| <= 1); for x >= 4, |J0(x)| <= sqrt(2/(pi x)) <= sqrt(1/(2 pi)) < 0.4028."""
    bound = F(-1007, 2500)
    GB = 8
    vals = [F(j0_enc(k << (PB - GB), k << (PB - GB))[0], ONE) for k in range(4 * 2 ** GB + 1)]
    worst = min(min(vals[k], vals[k + 1]) - F(1, 8 * 4 ** GB) for k in range(len(vals) - 1))
    ok = worst >= bound and F(1) / (2 * PI_LO_Q) < bound * bound
    log('J0 >= -0.4028: certified minimum on [0,4] >= %.9f; for x >= 4 |J0| <= %.6f' %
        (float(worst), math.sqrt(float(F(1) / (2 * PI_LO_Q)))))
    return ok

# =====================================================================================
# 5. Certificates
# =====================================================================================
def num(v, what):
    if isinstance(v, bool) or isinstance(v, float):
        raise ValueError('%s: non-exact number %r' % (what, v))
    if isinstance(v, int):
        return F(v)
    if isinstance(v, str):
        return F(v.strip())
    raise ValueError('%s: bad number %r' % (what, v))

def load_cert(path):
    """Returns a dict with 'kind' in {'B', 'A'}; A certificates get zeta (0 if absent) and Bfile (or None)."""
    c = json.load(open(path))
    if not isinstance(c, dict):
        raise ValueError('not a JSON object')
    side = c.get('side', 'A')
    cert = {'name': os.path.basename(path)}
    if side == 'B':
        cert['kind'] = 'B'
        for k in ('psi0', 'psiT'):
            cert[k] = num(c[k], k)
        raw = c['pi']
        if isinstance(raw, list):
            cert['pi'] = None
            cert['pi_list'] = []
            for ent in raw:
                if not (isinstance(ent, list) and len(ent) == 3):
                    raise ValueError('pi list entry %r is not [i, j, value]' % (ent,))
                cert['pi_list'].append((int(ent[0]), int(ent[1]), num(ent[2], 'pi')))
        else:
            cert['pi'] = num(raw, 'pi')
            cert['pi_list'] = None
    elif side == 'A':
        cert['kind'] = 'A'
        for k in ('D', 'd0', 'mu', 'psi0', 'psiT', 'eta', 'tau'):
            cert[k] = num(c[k], k)
        cert['zeta'] = num(c['zeta'], 'zeta') if 'zeta' in c else F(0)
        cert['Bfile'] = c.get('B')
        if cert['Bfile'] is not None and not isinstance(cert['Bfile'], str):
            raise ValueError('key B must be a file name')
    else:
        raise ValueError('unknown side %r' % (side,))
    cert['lam'] = [num(x, 'lam') for x in c['lam']]
    cert['u'] = [(int(i), int(j), num(v, 'u')) for i, j, v in c['u']]
    cert['wc'] = [([int(a) for a in I], [int(b) for b in J], num(v, 'wc')) for I, J, v in c['wc']]
    return cert

def structure(cert, geo, verdict, tag):
    """lam length, indices, u not on unit pairs, exact congruence of every wc entry; returns 'malformed'."""
    ok, malformed, msgs = True, False, []
    if len(cert['lam']) != NN:
        ok = False; malformed = True; msgs.append('lam has %d entries' % len(cert['lam']))
    unit = set(geo.edges)
    for (i, j, v) in cert['u']:
        if not (1 <= i < j <= NN):
            ok = False; malformed = True; msgs.append('bad u pair %r' % ((i, j),)); continue
        if (i, j) in unit:
            ok = False; msgs.append('u pair %r is a unit pair' % ((i, j),))
    npairs = 0
    for idx, (I, J, v) in enumerate(cert['wc']):
        if len(I) != len(J) or len(I) == 0:
            ok = False; malformed = True; msgs.append('wc[%d]: lengths %d, %d' % (idx, len(I), len(J))); continue
        if any(not (1 <= a <= NN) for a in I + J):
            ok = False; malformed = True; msgs.append('wc[%d]: index out of range' % idx); continue
        if len(set(I)) != len(I) or len(set(J)) != len(J):
            ok = False; malformed = True; msgs.append('wc[%d]: repeated index' % idx); continue
        for a in range(len(I)):
            for b in range(a + 1, len(I)):
                npairs += 1
                if geo.d2(I[a], I[b]) != geo.d2(J[a], J[b]):
                    ok = False
                    msgs.append('wc[%d] not congruent: |x%d-x%d|^2 = %s but |x%d-x%d|^2 = %s' %
                                (idx, I[a], I[b], kstr(geo.d2(I[a], I[b])), J[a], J[b], kstr(geo.d2(J[a], J[b]))))
    for m in msgs[:10]:
        log('   ' + m)
    if len(msgs) > 10:
        log('    ... (%d messages)' % len(msgs))
    log('%s: %d lam, %d u, %d wc entries; %d pairwise distance equalities checked'
        % (cert['name'], len(cert['lam']), len(cert['u']), len(cert['wc']), npairs))
    verdict('(%s) indices, u not on unit pairs, wc congruences (exact)' % tag, ok)
    return malformed

def common_denominator(vals):
    L = 1
    for q in vals:
        L = L * q.denominator // gcd(L, q.denominator)
    return L

def resolve_pi(cert, geo):
    """{edge: pi_e} for the 47 unit pairs (scalar pi: same value on every edge); None if malformed."""
    if cert['pi_list'] is None:
        pi = cert['pi']
        return {e: pi for e in geo.edges}, pi >= 0, ['pi is a single rational %s' % pi]
    E = set(geo.edges)
    pe, msgs, malformed = {}, [], False
    for (i, j, v) in cert['pi_list']:
        if not (1 <= i < j <= NN):
            malformed = True; msgs.append('pi entry %r: indices must satisfy 1 <= i < j <= 23' % ((i, j),)); continue
        if (i, j) not in E:
            malformed = True; msgs.append('pi entry %r is not a unit-distance pair' % ((i, j),)); continue
        if (i, j) in pe:
            malformed = True; msgs.append('pi entry %r appears twice' % ((i, j),)); continue
        pe[(i, j)] = v
    missing = sorted(E - set(pe))
    if missing:
        malformed = True; msgs.append('unit pairs without a pi entry: %s' % missing)
    nonneg = all(v >= 0 for v in pe.values())
    if not malformed:
        msgs.insert(0, 'pi is a list of 47 per-edge values, min %.6g, max %.6g, sum %.10g'
                    % (float(min(pe.values())), float(max(pe.values())), float(sum(pe.values(), F(0)))))
    return (None if malformed else pe), nonneg, msgs

# =====================================================================================
# 6. (A): independent sets of the unit-distance graph, exact
# =====================================================================================
def check_atoms(cert, geo):
    n = NN
    adj = [0] * n
    for (i, j) in geo.edges:
        adj[i - 1] |= 1 << (j - 1); adj[j - 1] |= 1 << (i - 1)
    indep = []
    def rec(k, S):
        if k == n:
            indep.append(S); return
        rec(k + 1, S)
        if not (adj[k] & S):
            rec(k + 1, S | (1 << k))
    rec(0, 0)
    # cross-check of the enumeration: all sets distinct and independent, and their number equals the
    # count from the recursion I(G) = I(G - v) + I(G - N[v]) on vertex masks
    memo = {}
    def count(mask):
        if mask == 0:
            return 1
        if mask in memo:
            return memo[mask]
        v = (mask & -mask).bit_length() - 1
        r = count(mask & ~(1 << v)) + count(mask & ~(1 << v) & ~adj[v])
        memo[mask] = r
        return r
    enum_ok = (len(set(indep)) == len(indep) == count((1 << n) - 1)
               and all(not (adj[k] & S) for S in indep for k in range(n) if S >> k & 1))
    lam, U, WC, psiT = cert['lam'], cert['u'], cert['wc'], cert['psiT']
    L = common_denominator([psiT] + lam + [v for _, _, v in U] + [v for _, _, v in WC])
    def toint(q):
        r = q * L; assert r.denominator == 1; return r.numerator
    lam_i = [toint(x) for x in lam]
    u_m = [((1 << (i - 1)) | (1 << (j - 1)), toint(v)) for i, j, v in U]
    wc_m = []
    for I, J, v in WC:
        mI = 0
        for a in I: mI |= 1 << (a - 1)
        mJ = 0
        for b in J: mJ |= 1 << (b - 1)
        wc_m.append((mI, mJ, toint(v)))
    Qempty = None; minQ = None; argmin = None; nneg = 0
    for S in indep:
        q = toint(psiT)
        for k in range(n):
            if S >> k & 1: q += lam_i[k]
        for m, v in u_m:
            if S & m == m: q += v
        for mI, mJ, v in wc_m:
            if S & mI == mI: q += v
            if S & mJ == mJ: q -= v
        if S == 0:
            Qempty = F(q, L); continue
        if q < 0: nneg += 1
        if minQ is None or q < minQ:
            minQ, argmin = q, S
    return dict(nsets=len(indep), enum_ok=enum_ok, Qempty=Qempty, minQ=F(minQ, L), nneg=nneg,
                argmin=[k + 1 for k in range(n) if argmin >> k & 1])

# =====================================================================================
# 7. (AB): all 2^23 subsets, exact (int64 over the common denominator, with an overflow bound)
# =====================================================================================
def check_all_subsets(cert, geo, pi_e, crosscheck=False):
    import numpy as np
    n = NN
    lam, U, WC, psiT = cert['lam'], cert['u'], cert['wc'], cert['psiT']
    L = common_denominator([psiT] + list(pi_e.values()) + lam + [v for _, _, v in U] + [v for _, _, v in WC])
    def toint(q):
        r = q * L; assert r.denominator == 1; return r.numerator
    terms = [(0, toint(psiT))]
    terms += [(1 << k, toint(lam[k])) for k in range(n)]
    terms += [((1 << (i - 1)) | (1 << (j - 1)), toint(v)) for i, j, v in U]
    terms += [((1 << (i - 1)) | (1 << (j - 1)), toint(pi_e[(i, j)])) for (i, j) in geo.edges]
    for I, J, v in WC:
        mI = sum(1 << (a - 1) for a in I); mJ = sum(1 << (b - 1) for b in J)
        vi = toint(v)
        terms.append((mI, vi)); terms.append((mJ, -vi))
    f = {}
    for m, v in terms:
        f[m] = f.get(m, 0) + v
    S_abs = sum(abs(v) for m, v in terms)
    log('(AB) common denominator %d%s; sum of |terms| x denominator < 2^%d (int64 safe below 2^62)'
        % (L, ' = 2^%d' % (L.bit_length() - 1) if L & (L - 1) == 0 else '', S_abs.bit_length()))
    if S_abs >= 1 << 62:
        return False, 'overflow risk: data too large for exact int64 arithmetic', None
    # subset-sum transform g(T) = sum_{M subset T} f(M); every partial sum is bounded by S_abs
    g = np.zeros(1 << n, dtype=np.int64)
    for m, v in f.items():
        g[m] = v
    for i in range(n):
        g2 = g.reshape(-1, 2, 1 << i)
        g2[:, 1, :] += g2[:, 0, :]
    gmin = int(g.min()); garg = int(g.argmin()); nneg = int(np.count_nonzero(g < 0))
    agree = None
    if crosscheck:
        idx = np.arange(1 << n, dtype=np.int64)
        bits = [((idx >> i) & 1).astype(bool) for i in range(n)]
        del idx
        Q = np.zeros(1 << n, dtype=np.int64)
        ind = np.empty(1 << n, dtype=bool)
        for m, v in f.items():
            if v == 0: continue
            if m == 0:
                Q += np.int64(v); continue
            ks = [i for i in range(n) if m >> i & 1]
            np.copyto(ind, bits[ks[0]])
            for i in ks[1:]:
                np.logical_and(ind, bits[i], out=ind)
            np.add(Q, np.int64(v), out=Q, where=ind)
        agree = bool(np.array_equal(Q, g))
        del Q, bits, ind
        log('(AB) direct recomputation of all values: identical = %s' % agree)
    # exact recomputation from the definition at sample subsets (independent of the table f)
    rnd = random.Random(2024)
    sample = [0, (1 << n) - 1, garg] + [1 << k for k in range(n)] + [rnd.getrandbits(n) for _ in range(1000)]
    edge_masks = [((1 << (i - 1)) | (1 << (j - 1)), (i, j)) for (i, j) in geo.edges]
    bad = 0
    for T in sample:
        q = psiT + sum((lam[k] for k in range(n) if T >> k & 1), F(0))
        q += sum((v for i, j, v in U if (T >> (i - 1) & 1) and (T >> (j - 1) & 1)), F(0))
        q += sum((pi_e[e] for m, e in edge_masks if T & m == m), F(0))
        for I, J, v in WC:
            inI = all(T >> (a - 1) & 1 for a in I)
            inJ = all(T >> (b - 1) & 1 for b in J)
            if inI != inJ:
                q += v if inI else -v
        if q != F(int(g[T]), L):
            bad += 1
    del g
    minQ = F(gmin, L)
    log('(AB) %d subsets; exact recomputation at %d of them: %d mismatches; min Q_B = %s (~%.3e) at T = %s; '
        '#negative = %d' % (1 << n, len(sample), bad, minQ, float(minQ),
                            [k + 1 for k in range(n) if garg >> k & 1], nneg))
    return (nneg == 0 and bad == 0 and agree is not False), '(#T with Q_B < 0: %d)' % nneg, minQ

# =====================================================================================
# 8. Sweeps on [0, T]
# =====================================================================================
GRID = 40               # grid points are multiples of 2^-40
INV_SQRT2_UP = F(70711, 100000)
CTX = {}

def _point(A):
    """Enclosures at t = A 2^-GRID of Phi_s(t) = psi0 J0(t) + C0 + sum_g w_g J0(r_g t) (w_g = -V_g) and of
    W(t) = 529 - H(t) = sum_d 2 m_d (1 - J0(r_d t)) (only when needed)."""
    X = A << (PB - GRID)
    J = {}
    for key, Rlo, Rhi in CTX['dists']:
        J[key] = j0_enc(X, X) if key == KONE else j0_enc((Rlo * A) >> GRID, cdiv(Rhi * A, 1 << GRID))
    jlo, jhi = J[KONE]
    lo = CTX['C0lo'] + iqmul_lo(CTX['psi0'], jlo, jhi)
    hi = CTX['C0hi'] + iqmul_hi(CTX['psi0'], jlo, jhi)
    for key, w in CTX['ucoef']:
        lo += iqmul_lo(w, *J[key]); hi += iqmul_hi(w, *J[key])
    Wlo = Whi = 0
    for key, m2 in CTX['wcoef']:
        jlo, jhi = J[key]
        Wlo += m2 * (ONE - jhi); Whi += m2 * (ONE - jlo)
    return lo, hi, Wlo, Whi

def _e_up(A):
    if A == 0:
        return CTX['mu_up']
    D = CTX['D']
    v = cdiv((1 << (2 * PB + 1 + GRID)) * D.denominator, PI_LO * D.numerator * A)
    return min(CTX['mu_up'], isqrt_hi(v))

def _e_lo(A):
    D = CTX['D']
    v = ((1 << (2 * PB + 1 + GRID)) * D.denominator) // (PI_HI * D.numerator * A)
    return min(CTX['mu_lo'], isqrt_lo(v))

def _F_upper(A):
    lo, hi, Wlo, Whi = _point(A)
    if not CTX['hasE']:
        return hi
    el = _e_lo(A)
    return hi - iqmul_lo(CTX['eta'], el, el) - ((iqmul_lo(CTX['zetaL'], el, el) * max(Wlo, 0)) >> PB)

def _fail_status(A):
    cands = ([A] if (A > 0 or not CTX['hasE']) else []) + [A + (1 << k) for k in range(51)]
    up, Av = min((_F_upper(Av), Av) for Av in cands)
    if up < ONE:
        return ('sweep stopped at t = %.12g; VIOLATION PROVEN: the function is <= 1 - %.6e at t = %.12g'
                % (A / 2 ** GRID, (ONE - up) / ONE, Av / 2 ** GRID))
    return 'certified lower bound at t = %.12g does not reach 1 (cannot certify)' % (A / 2 ** GRID)

def sweep_chunk(job):
    """Adaptive sweep of [A0, A1].  On [a, b], E >= e(D t) and K = zeta L_B E give
       F(t) >= G(t) - eta E,  G = Phi_s - K W   (W >= 0, eta, zeta L_B >= 0),
       chord:     G(t) >= min(G(a), G(b)) - (M2 + K M2W) h^2/8,
       lipschitz: G(t) >= (G(a) + G(b))/2 - (M1 + K M1W) h/2."""
    A0, A1, mode = job
    for k in NEVAL: NEVAL[k] = 0
    M1, M2, M1W, M2W, hasE = CTX['M1'], CTX['M2'], CTX['M1W'], CTX['M2W'], CTX['hasE']
    A, pa = A0, _point(A0)
    steps = nevals = 0
    min_margin = None; min_at = None
    while A < A1:
        if hasE:
            E = _e_up(A)
            etaE = iqmul_hi(CTX['eta'], 0, E)
            K = iqmul_hi(CTX['zetaL'], 0, E)
        else:
            E = etaE = K = 0
        Ga = pa[0] - cdiv(K * max(pa[3], 0), ONE)
        slack = Ga - etaE - ONE
        if slack <= 0:
            return dict(ok=False, status=_fail_status(A))
        Kf = K / ONE
        if mode == 'chord':
            h = math.sqrt(2.0 * (slack / ONE) / (float(M2) + Kf * float(M2W)))
        else:
            h = (slack / ONE) / (float(M1) + Kf * float(M1W))
        H = max(1, min(int(h * (1 << GRID)), 1 << GRID, A1 - A))
        while True:
            B = A + H
            pb = _point(B)
            nevals += 1
            Gb = pb[0] - cdiv(K * max(pb[3], 0), ONE)
            if mode == 'chord':
                hh = cdiv(H * H << PB, 8 << (2 * GRID))
                loss = iqmul_hi(M2, 0, hh) + iqmul_hi(M2W, 0, cdiv(K * hh, ONE))
                lb = min(Ga, Gb) - loss - etaE
            else:
                hh = cdiv(H << PB, 2 << GRID)
                loss = iqmul_hi(M1, 0, hh) + iqmul_hi(M1W, 0, cdiv(K * hh, ONE))
                lb = (Ga + Gb) // 2 - loss - etaE
            if lb >= ONE:
                break
            H //= 2
            if H == 0:
                return dict(ok=False, status='step size underflow; ' + _fail_status(A))
        if min_margin is None or lb - ONE < min_margin:
            min_margin, min_at = lb - ONE, (A, B)
        A, pa = B, pb
        steps += 1
    return dict(ok=True, steps=steps, nevals=nevals, min_margin=min_margin, min_at=min_at, neval=dict(NEVAL))

def sweep(T, mode, workers=1):
    t0 = time.time()
    TU = int(T * (1 << GRID))
    assert F(TU, 1 << GRID) == T
    if workers <= 1:
        res = [sweep_chunk((0, TU, mode))]
    else:
        nch = 8 * workers
        bounds = sorted(set([0] + [int(TU * (k / nch) ** 0.8) for k in range(1, nch)] + [TU]))
        import multiprocessing as mp
        with mp.get_context('fork').Pool(workers) as pool:
            res = pool.map(sweep_chunk, [(bounds[k], bounds[k + 1], mode) for k in range(len(bounds) - 1)],
                           chunksize=1)
    if not all(r['ok'] for r in res):
        return False, '; '.join(r['status'] for r in res if not r['ok'])
    mm = min(res, key=lambda r: r['min_margin'])
    return True, ('T = %s, %d steps, %d J0 enclosures per grid point (%d by power series, %d by Hankel '
                  'expansion); smallest certified lower bound of (function - 1) on an interval: %.4e; %.1f s'
                  % (T, sum(r['steps'] for r in res), len(CTX['dists']),
                     sum(r['neval']['series'] for r in res), sum(r['neval']['hankel'] for r in res),
                     mm['min_margin'] / ONE, time.time() - t0))

def setup_ctx(cert, geo, withW, hasE, zetaL=F(0)):
    groups = {}
    for i, j, v in cert['u']:
        s = geo.d2(i, j)
        groups[s] = groups.get(s, F(0)) + v          # u terms with exactly equal distance are merged
    need = set(groups) | {KONE}
    if withW:
        need |= set(geo.distinct)
    enc = {}
    for s in need:
        slo, shi = k_enclose(s)
        assert slo > 0
        enc[s] = (slo, shi, isqrt_lo(slo), isqrt_hi(shi))
    C0 = 2 * cert.get('d0', F(0)) * cert.get('tau', F(0)) - sum(cert['lam'], F(0))
    M2 = abs(cert['psi0']) + sum((abs(V) * fx_to_F(enc[s][1]) for s, V in groups.items()), F(0))
    M1 = INV_SQRT2_UP * (abs(cert['psi0']) + sum((abs(V) * fx_to_F(enc[s][3]) for s, V in groups.items()), F(0)))
    M2W = M1W = F(0)
    if withW:
        M2W = sum((2 * m * fx_to_F(enc[s][1]) for s, m in geo.distinct.items()), F(0))
        M1W = INV_SQRT2_UP * sum((2 * m * fx_to_F(enc[s][3]) for s, m in geo.distinct.items()), F(0))
    CTX.clear()
    CTX.update(psi0=cert['psi0'], C0lo=fx_floor(C0), C0hi=fx_ceil(C0),
               dists=[(s, enc[s][2], enc[s][3]) for s in need], ucoef=[(s, -V) for s, V in groups.items()],
               wcoef=[(s, 2 * m) for s, m in geo.distinct.items()] if withW else [],
               hasE=hasE, M1=M1, M2=M2, M1W=M1W, M2W=M2W)
    if hasE:
        CTX.update(D=cert['D'], mu_up=fx_ceil(cert['mu']), mu_lo=fx_floor(cert['mu']), eta=cert['eta'], zetaL=zetaL)
    sqrt2pi = sqrt_up(F(2) / PI_LO_Q)
    Csum = abs(cert['psi0']) + sum((abs(V) * inv_sqrt_up(fx_to_F(enc[s][2])) for s, V in groups.items()), F(0))
    return groups, C0, sqrt2pi, Csum

def choose_T(Phiinf, Ct, verdict, label):
    """tail: function >= Phiinf - Ct/sqrt(t); returns T with a certified bound >= 1 beyond it."""
    if Phiinf <= 1:
        verdict('%s tail: limit %.6g > 1' % (label, float(Phiinf)), False)
        return None
    Tneed = (Ct / (Phiinf - 1)) ** 2
    T = F(math.ceil(float(Tneed) * 1.02 / 10) * 10)
    tail = Phiinf - Ct * inv_sqrt_up(T)
    verdict('%s tail t >= T = %s: function >= %.8g (T >= %.2f needed)' % (label, T, float(tail), float(Tneed)),
            tail >= 1)
    return T if tail >= 1 else None

# =====================================================================================
# 9. Side B and side A
# =====================================================================================
def run_B(cert, geo, opts, verdict):
    log('\n--- side B: %s ---' % cert['name'])
    malformed = structure(cert, geo, verdict, 'GB')
    pi_e, nonneg, msgs = resolve_pi(cert, geo)
    for m in msgs:
        log('   ' + m)
    if cert['pi_list'] is not None:
        verdict('(GB) pi list covers exactly the 47 unit pairs', pi_e is not None)
    verdict('(GB) pi_e >= 0', nonneg)
    if malformed or pi_e is None:
        verdict('(GB) certificate well formed', False)
        return None
    okAB, info, minQB = check_all_subsets(cert, geo, pi_e, crosscheck=opts['ab_crosscheck'])
    verdict('(AB) Q_B(T) >= 0 for all 2^23 subsets T (exact)', okAB, info)
    sumlam = sum(cert['lam'], F(0))
    PhiB0 = cert['psi0'] - sumlam - sum((v for _, _, v in cert['u']), F(0))
    log('Phi_B(0) = %s (~%.10g); Phi_B(oo) = %.10g' % (PhiB0, float(PhiB0), float(-sumlam)))
    groups, C0, sqrt2pi, Csum = setup_ctx(cert, geo, withW=False, hasE=False)
    T = choose_T(-sumlam, sqrt2pi * Csum, verdict, '(PB)')
    if T is not None:
        ok, info = sweep(T, 'chord', opts['workers'])
        verdict('(PB) Phi_B(t) >= 1 on [0, T] (t = 0 included)', ok, info)
        if opts['lipschitz']:
            ok, info = sweep(T, 'lipschitz', opts['workers'])
            verdict('(PB) redundant sweep with |J1| <= 1/sqrt2', ok, info)
    LB = max(F(0), cert['psi0'] + sum(pi_e.values(), F(0)))
    log('psiT_B = %s (~%.12g); L_B = max(0, psi0 + sum pi_e) = %s (~%.12g)'
        % (cert['psiT'], float(cert['psiT']), LB, float(LB)))
    return dict(psiT=cert['psiT'], L=LB)

def run_A(cert, Bres, geo, opts, verdict):
    hit = cert['Bfile'] is not None
    log('\n--- side A: %s (%s) ---' % (cert['name'], 'with B certificate' if hit else 'without B certificate'))
    for k in ('D', 'd0', 'mu', 'eta', 'tau', 'zeta'):
        log('   %-4s = %s (~%.10g)' % (k, cert[k], float(cert[k])))
    malformed = structure(cert, geo, verdict, 'G')
    D, d0, mu, eta, tau, zeta = (cert[k] for k in ('D', 'd0', 'mu', 'eta', 'tau', 'zeta'))
    verdict('(G) eta >= 0, tau >= 0, zeta >= 0, 0 < d0 < 1, D > 0, mu >= 0.4028',
            eta >= 0 and tau >= 0 and zeta >= 0 and 0 < d0 < 1 and D > 0 and mu >= F('0.4028'))
    if not hit and zeta != 0:
        verdict('(G) zeta = 0 when no B certificate is given', False)
    if malformed:
        verdict('(G) certificate well formed', False)
        return None
    sumlam = sum(cert['lam'], F(0)); sumu = sum((v for _, _, v in cert['u']), F(0))
    Phi0 = cert['psi0'] - sumlam - sumu + eta + tau * (2 * d0 - 1)
    verdict('(Z) Phi0 = %s >= 1 (exact; Phi0 - 1 ~ %.4e)' % (Phi0, float(Phi0 - 1)), Phi0 >= 1)
    r = check_atoms(cert, geo)
    verdict('(A) %d independent sets enumerated (cross-checked by a second count)' % r['nsets'], r['enum_ok'])
    verdict('(A) Q(S) >= 0 for every nonempty independent set (exact; min %s ~ %.4e at S = %s)'
            % (r['minQ'], float(r['minQ']), r['argmin']), r['nneg'] == 0, '(#violations %d)' % r['nneg'])
    verdict('(A) Q(empty) = psiT >= zeta (psiT - zeta = %s ~ %.4e)' % (cert['psiT'] - zeta, float(cert['psiT'] - zeta)),
            r['Qempty'] == cert['psiT'] and cert['psiT'] >= zeta)
    verdict('(P) J0(x) >= -0.4028 for all x >= 0', prove_j0_lower())
    LB, psiTB = (Bres['L'], Bres['psiT']) if hit and Bres else (F(0), F(0))
    if hit and Bres is None:
        verdict('(P) a valid B certificate is required', False)
        return None
    if not (D > 0 and eta >= 0 and mu >= 0 and zeta >= 0):
        verdict('(P) not attempted (needs D > 0, eta >= 0, mu >= 0, zeta >= 0)', False)
        return dict(bound=cert['psiT'] + tau * d0 * d0 - zeta * (1 - psiTB), D=D)
    zetaL = zeta * LB
    withW = zetaL > 0
    groups, C0, sqrt2pi, Csum = setup_ctx(cert, geo, withW=withW, hasE=True, zetaL=zetaL)
    if hit:
        log('zeta * L_B = %.10g' % float(zetaL))
    Ct = sqrt2pi * (Csum + (eta + 529 * zetaL) * inv_sqrt_up(D))
    T = choose_T(C0, Ct, verdict, '(P)')
    if T is not None:
        ok, info = sweep(T, 'chord', opts['workers'])
        verdict('(P) the function is >= 1 on (0, T]', ok, info)
        if opts['lipschitz']:
            ok, info = sweep(T, 'lipschitz', opts['workers'])
            verdict('(P) redundant sweep with |J1| <= 1/sqrt2', ok, info)
    bound = cert['psiT'] + tau * d0 * d0 - zeta * (1 - psiTB)
    return dict(bound=bound, D=D)
