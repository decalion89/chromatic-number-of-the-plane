"""
Helper module of check_m1_indep.py: exact arithmetic in Q(sqrt3, sqrt11), certificate validation,
the exact check over all independent sets, rigorous enclosures of J0 and the sweep for Phi.
See check_m1_indep.py for what is checked and which classical facts are assumed.
"""
import sys, os, json, time, math, random
from fractions import Fraction as F
from math import isqrt, gcd

def log(*a):
    print(*a, flush=True)

# =====================================================================================
# 1. Exact arithmetic in K = Q(sqrt3, sqrt11), basis (1, sqrt3, sqrt11, sqrt33)
#    (a Q-basis of K, so two elements are equal iff their four coefficients are equal)
# =====================================================================================
KZERO = (F(0), F(0), F(0), F(0))
KONE = (F(1), F(0), F(0), F(0))

def kadd(x, y): return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])
def ksub(x, y): return (x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3])
def kmul(x, y):
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return (a1 * a2 + 3 * b1 * b2 + 11 * c1 * c2 + 33 * d1 * d2,
            a1 * b2 + b1 * a2 + 11 * (c1 * d2 + d1 * c2),
            a1 * c2 + c1 * a2 + 3 * (b1 * d2 + d1 * b2),
            a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2)
def kstr(x):
    names = ['', '*sqrt3', '*sqrt11', '*sqrt33']
    parts = [str(c) + nm for c, nm in zip(x, names) if c != 0]
    return ' + '.join(parts) if parts else '0'
def csub(z, w): return (ksub(z[0], w[0]), ksub(z[1], w[1]))
def cabs2(z): return kadd(kmul(z[0], z[0]), kmul(z[1], z[1]))

# The published 23-point set (only used to report whether a certificate uses exactly this set).
_X23 = [
    ('0', '0', '0', '0', '0', '0', '0', '0'), ('1', '0', '0', '0', '0', '0', '0', '0'),
    ('1/2', '0', '0', '0', '0', '1/2', '0', '0'), ('3/2', '0', '0', '0', '0', '1/2', '0', '0'),
    ('5/6', '0', '0', '0', '0', '0', '1/6', '0'), ('5/12', '0', '0', '-1/12', '0', '5/12', '1/12', '0'),
    ('5/4', '0', '0', '-1/12', '0', '5/12', '1/4', '0'), ('13/12', '0', '0', '1/12', '0', '1/12', '-1/12', '0'),
    ('1/4', '0', '0', '1/12', '0', '1/12', '-1/4', '0'), ('13/12', '0', '0', '-1/12', '0', '-1/12', '-1/12', '0'),
    ('2/3', '0', '0', '0', '0', '1/2', '-1/6', '0'), ('7/12', '0', '0', '-1/12', '0', '5/12', '-1/12', '0'),
    ('5/6', '0', '0', '-1/6', '0', '-1/6', '1/6', '0'), ('7/12', '0', '0', '1/12', '0', '7/12', '-1/12', '0'),
    ('1/12', '0', '0', '-1/12', '0', '-1/12', '-1/12', '0'), ('2/3', '0', '0', '-1/6', '0', '1/3', '-1/6', '0'),
    ('19/12', '0', '0', '-1/12', '0', '5/12', '-1/12', '0'), ('1/3', '0', '0', '0', '0', '1/2', '1/6', '0'),
    ('3/4', '0', '0', '-1/12', '0', '-1/12', '1/4', '0'), ('1', '0', '0', '-1/6', '0', '-1/6', '0', '0'),
    ('7/6', '0', '0', '0', '0', '0', '-1/6', '0'), ('4/3', '0', '0', '-1/6', '0', '1/3', '1/6', '0'),
    ('19/12', '0', '0', '-1/4', '0', '1/4', '-1/12', '0'),
]
X23 = [((F(r[0]), F(r[1]), F(r[2]), F(r[3])), (F(r[4]), F(r[5]), F(r[6]), F(r[7]))) for r in _X23]

# =====================================================================================
# 2. Fixed-point interval arithmetic on integers (N stands for N / 2^PB; outward rounding)
# =====================================================================================
PB = 200
ONE = 1 << PB
ONE2 = 1 << (2 * PB)
TINY = 1 << (PB - 160)

def cdiv(a, b):
    return -((-a) // b)
def fx_floor(q): return (q.numerator << PB) // q.denominator
def fx_ceil(q): return cdiv(q.numerator << PB, q.denominator)
def fx_to_F(N): return F(N, ONE)
def imul(alo, ahi, blo, bhi):
    p1, p2, p3, p4 = alo * blo, alo * bhi, ahi * blo, ahi * bhi
    return min(p1, p2, p3, p4) >> PB, -((-max(p1, p2, p3, p4)) >> PB)
def iqmul_lo(q, lo, hi):
    p, d = q.numerator, q.denominator
    return (p * lo) // d if p >= 0 else (p * hi) // d
def iqmul_hi(q, lo, hi):
    p, d = q.numerator, q.denominator
    return cdiv(p * hi, d) if p >= 0 else cdiv(p * lo, d)
def isqrt_lo(N): return isqrt(N << PB)
def isqrt_hi(N): return isqrt(N << PB) + 1
def inv_sqrt_up(q): return F(ONE, isqrt((q.numerator << (2 * PB)) // q.denominator))
def sqrt_up(q): return F(isqrt(cdiv(q.numerator << (2 * PB), q.denominator)) + 1, ONE)

def _atan_inv(m, n):
    s = F(0)
    for k in range(n):
        s += F((-1) ** k, (2 * k + 1) * m ** (2 * k + 1))
    nxt = F(1, (2 * n + 1) * m ** (2 * n + 1))
    return (s, s + nxt) if n % 2 == 0 else (s - nxt, s)
_a5 = _atan_inv(5, 120); _a239 = _atan_inv(239, 40)
PI_LO_Q = 16 * _a5[0] - 4 * _a239[1]          # Machin: pi = 16 atan(1/5) - 4 atan(1/239)
PI_HI_Q = 16 * _a5[1] - 4 * _a239[0]
assert F(314159265358979323846264338327950288419716939937510, 10 ** 50) < PI_LO_Q < PI_HI_Q \
    < F(314159265358979323846264338327950288419716939937511, 10 ** 50)
PI_LO = fx_floor(PI_LO_Q); PI_HI = fx_ceil(PI_HI_Q)

SQ_ENC = {}
for _n in (3, 11, 33):
    _lo = isqrt(_n << (2 * PB)); SQ_ENC[_n] = (_lo, _lo + 1)

def k_enclose(x):
    lo, hi = fx_floor(x[0]), fx_ceil(x[0])
    for c, n in ((x[1], 3), (x[2], 11), (x[3], 33)):
        if c == 0: continue
        Lo, Hi = SQ_ENC[n]
        lo += iqmul_lo(c, Lo, Hi); hi += iqmul_hi(c, Lo, Hi)
    return lo, hi

# =====================================================================================
# 3. Rigorous enclosures of J0 on [Xlo, Xhi] (fixed point, 0 <= Xlo <= Xhi)
# =====================================================================================
NEVAL = {'series': 0, 'hankel': 0}
OPTS = {'series_only': False}

def j0_series(Xlo, Xhi, extra=0):
    """Power series at the exact point Xlo with PB+extra bits and directed rounding; tail <= 2 t_{k+1}
    once t_k <= 2^-160 and y/(k+2)^2 <= 1/2; widened by Xhi - Xlo (|J0'| <= 1)."""
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
    """x > 30: J0(x) = sqrt(2/(pi x)) (P cos w - Q sin w) + err, w = x - pi/4, P = g0 - g2 + ...,
    -Q = g1 - g3 + ..., g_k = |a_k| x^-k, a_k = (-1)^k (1*3*...*(2k-1))^2/(k! 8^k), |err| <=
    sqrt(2/(pi x)) g_n for the first omitted index n (from Hankel's integral for H0^(1))."""
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

_J0_REF = [   # J0(x) to 50 significant digits (self-test only)
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
    bad = 0
    tol = F(1, 10 ** 45)
    for xs, vs in _J0_REF:
        x, v = F(xs), F(vs)
        lo, hi = j0_enc(fx_floor(x), fx_ceil(x))
        if not (F(lo, ONE) - tol <= v <= F(hi, ONE) + tol):
            bad += 1
            log('   J0 self-test: tabulated value at x = %s outside the enclosure' % xs)
    rnd = random.Random(12345)
    xs = [30 + 15 * k / 40 + 0.0123 for k in range(40)] + [rnd.uniform(30, 1800) for _ in range(20)]
    for x in xs:
        X = int(x * 2 ** 40) << (PB - 40)
        a = j0_series(X, X, extra=int(1.4427 * x) + 80)
        b = j0_hankel(X, X)
        if max(a[0], b[0]) > min(a[1], b[1]):
            bad += 1
            log('   J0 self-test: series and large-argument enclosures disagree at x = %r' % x)
    log('J0 self-test: %d tabulated values, %d series/large-argument comparisons, %d failures'
        % (len(_J0_REF), len(xs), bad))
    if bad:
        raise RuntimeError('J0 self-test failed')

# =====================================================================================
# 4. Certificate parsing and the exact structural checks (G)
# =====================================================================================
def num(v, what):
    if isinstance(v, bool) or isinstance(v, float):
        raise ValueError('%s: non-exact number %r' % (what, v))
    if isinstance(v, int):
        return F(v)
    if isinstance(v, str):
        return F(v.strip())
    raise ValueError('%s: bad number %r' % (what, v))

def index(v, what):
    if isinstance(v, bool) or not isinstance(v, int):
        raise ValueError('%s: index %r is not an integer' % (what, v))
    return v

def load_cert(path):
    c = json.load(open(path))
    if not isinstance(c, dict):
        raise ValueError('certificate is not a JSON object')
    if c.get('side') != 'M1':
        raise ValueError('"side" must be "M1"')
    cert = {'name': os.path.basename(path)}
    pts = []
    for k, p in enumerate(c['points']):
        if not (isinstance(p, list) and len(p) == 2 and all(isinstance(q, list) and len(q) == 4 for q in p)):
            raise ValueError('point %d is not [[a,b,c,d],[a,b,c,d]]' % (k + 1))
        pts.append((tuple(num(x, 'point %d' % (k + 1)) for x in p[0]), tuple(num(x, 'point %d' % (k + 1)) for x in p[1])))
    cert['points'] = pts
    cert['psi0'] = num(c['psi0'], 'psi0')
    cert['psiT'] = num(c['psiT'], 'psiT')
    cert['lam'] = [num(x, 'lam') for x in c['lam']]
    cert['u'] = []
    for e in c['u']:
        if not (isinstance(e, list) and len(e) == 3):
            raise ValueError('u entry %r is not [i, j, value]' % (e,))
        cert['u'].append((index(e[0], 'u'), index(e[1], 'u'), num(e[2], 'u')))
    cert['wc'] = []
    for e in c['wc']:
        if not (isinstance(e, list) and len(e) == 3 and isinstance(e[0], list) and isinstance(e[1], list)):
            raise ValueError('wc entry %r is not [I, J, value]' % (e,))
        cert['wc'].append(([index(a, 'wc') for a in e[0]], [index(b, 'wc') for b in e[1]], num(e[2], 'wc')))
    return cert

class Geometry:
    def __init__(self, pts):
        self.n = len(pts)
        self.sq = {}
        for i in range(self.n):
            for j in range(i + 1, self.n):
                self.sq[(i + 1, j + 1)] = cabs2(csub(pts[i], pts[j]))
        self.coincide = sorted(p for p, s in self.sq.items() if s == KZERO)
        self.edges = sorted(p for p, s in self.sq.items() if s == KONE)
        self.distinct = {}
        for s in self.sq.values():
            self.distinct[s] = self.distinct.get(s, 0) + 1
    def d2(self, i, j):
        return self.sq[(i, j)] if i < j else self.sq[(j, i)]

def check_structure(cert, geo, verdict):
    """(G): returns True if the data are well formed enough for the remaining checks."""
    n = geo.n
    msgs = []
    ok_points = n >= 1 and not geo.coincide
    if geo.coincide:
        msgs.append('coinciding points: %s' % geo.coincide[:5])
    ok_lam = len(cert['lam']) == n
    if not ok_lam:
        msgs.append('lam has %d entries for %d points' % (len(cert['lam']), n))
    unit = set(geo.edges)
    ok_u, seen = True, set()
    for (i, j, v) in cert['u']:
        if not (1 <= i <= n and 1 <= j <= n):
            ok_u = False; msgs.append('u pair %r: index out of range 1..%d' % ((i, j), n)); continue
        if not i < j:
            ok_u = False; msgs.append('u pair %r: needs i < j' % ((i, j),)); continue
        if (i, j) in seen:
            ok_u = False; msgs.append('u pair %r appears twice' % ((i, j),)); continue
        seen.add((i, j))
        if (i, j) in unit:
            ok_u = False; msgs.append('u pair %r is at distance 1' % ((i, j),))
    ok_wc, npairs, sizes = True, 0, {}
    for idx, (I, J, v) in enumerate(cert['wc']):
        if len(I) != len(J) or len(I) == 0:
            ok_wc = False; msgs.append('wc[%d]: |I| = %d, |J| = %d' % (idx, len(I), len(J))); continue
        if any(not (1 <= a <= n) for a in I + J):
            ok_wc = False; msgs.append('wc[%d]: index out of range 1..%d in %s / %s' % (idx, n, I, J)); continue
        if len(set(I)) != len(I) or len(set(J)) != len(J):
            ok_wc = False; msgs.append('wc[%d]: repeated index' % idx); continue
        sizes[len(I)] = sizes.get(len(I), 0) + 1
        for a in range(len(I)):
            for b in range(a + 1, len(I)):
                npairs += 1
                if geo.d2(I[a], I[b]) != geo.d2(J[a], J[b]):
                    ok_wc = False
                    msgs.append('wc[%d] not congruent: |x%d-x%d|^2 = %s but |x%d-x%d|^2 = %s'
                                % (idx, I[a], I[b], kstr(geo.d2(I[a], I[b])), J[a], J[b], kstr(geo.d2(J[a], J[b]))))
    for m in msgs[:12]:
        log('   ' + m)
    if len(msgs) > 12:
        log('   ... (%d messages)' % len(msgs))
    log('n = %d points, %d unit-distance pairs, %d distinct distances; %d u entries; %d wc entries (sizes %s), '
        '%d pairwise distance equalities checked' % (n, len(geo.edges), len(geo.distinct), len(cert['u']),
                                                     len(cert['wc']), dict(sorted(sizes.items())), npairs))
    verdict('(G) at least one point, no two points equal', ok_points)
    verdict('(G) one lam per point', ok_lam)
    verdict('(G) u: indices in 1..n, i < j, no repeated pair, no pair at distance 1', ok_u)
    verdict('(G) wc: |I| = |J| >= 1, indices in 1..n and distinct, congruent in order (exact)', ok_wc)
    return ok_points and ok_lam and ok_u and ok_wc

# =====================================================================================
# 5. (A): exact check over all independent sets (numpy for the set enumeration and the transform)
# =====================================================================================
def common_denominator(vals):
    L = 1
    for q in vals:
        L = L * q.denominator // gcd(L, q.denominator)
    return L

def check_independent_sets(cert, geo, max_sets=30_000_000):
    import numpy as np
    n = geo.n
    if n > 62:
        return dict(ok=False, msg='more than 62 points are not supported (vertex masks are int64)')
    adj = [0] * n
    for (i, j) in geo.edges:
        adj[i - 1] |= 1 << (j - 1); adj[j - 1] |= 1 << (i - 1)
    # number of independent sets by the recursion I(G) = I(G - v) + I(G - N[v]) (memoised on vertex masks)
    memo = {}
    sys.setrecursionlimit(max(1000, 10 * n + 100))
    def count(mask):
        if mask == 0:
            return 1
        r = memo.get(mask)
        if r is None:
            v = (mask & -mask).bit_length() - 1
            rest = mask & ~(1 << v)
            r = count(rest) + count(rest & ~adj[v])
            memo[mask] = r
        return r
    total = count((1 << n) - 1)
    log('(A) number of independent sets (including the empty set): %d' % total)
    if total > max_sets:
        return dict(ok=False, msg='too many independent sets (%d > %d) for this checker' % (total, max_sets))
    # enumeration, vertex by vertex: every independent set arises exactly once
    sets = np.zeros(1, dtype=np.int64)
    adj_np = [np.int64(a) for a in adj]
    for v in range(n):
        add = sets[(sets & adj_np[v]) == 0] | np.int64(1 << v)
        sets = np.concatenate([sets, add])
    sets.sort()
    enum_ok = (len(sets) == total and bool(np.all(np.diff(sets) > 0)))
    # coefficient f(M) of [M subset S] for every vertex mask M
    lam, U, WC, psiT = cert['lam'], cert['u'], cert['wc'], cert['psiT']
    L = common_denominator([psiT] + lam + [v for _, _, v in U] + [v for _, _, v in WC])
    def toint(q):
        r = q * L; assert r.denominator == 1; return r.numerator
    terms = [(0, toint(psiT))] + [(1 << k, toint(lam[k])) for k in range(n)]
    terms += [((1 << (i - 1)) | (1 << (j - 1)), toint(v)) for i, j, v in U]
    for I, J, v in WC:
        vi = toint(v)
        terms.append((sum(1 << (a - 1) for a in I), vi)); terms.append((sum(1 << (b - 1) for b in J), -vi))
    f = {}
    for m, v in terms:
        f[m] = f.get(m, 0) + v
    def independent(m):
        return all(not (adj[k] & m) for k in range(n) if m >> k & 1)
    dead = sum(1 for m in f if not independent(m))       # masks never contained in an independent set
    S_abs = sum(abs(v) for m, v in terms)
    use_int64 = S_abs < (1 << 62)
    log('(A) common denominator %d%s; sum of |terms| x denominator < 2^%d -> %s arithmetic; '
        '%d coefficient masks are not independent (they never contribute)'
        % (L, ' = 2^%d' % (L.bit_length() - 1) if L & (L - 1) == 0 else '', S_abs.bit_length(),
           'int64' if use_int64 else 'Python-integer', dead))
    g = np.zeros(len(sets), dtype=np.int64 if use_int64 else object)
    if not use_int64:
        g[:] = 0
    for m, v in f.items():
        if v != 0 and independent(m):
            pos = int(np.searchsorted(sets, m))
            assert sets[pos] == m
            g[pos] += v
    # subset-sum transform on the family of independent sets (closed under taking subsets):
    # after the pass for bit i, g[S] = sum of f(M) over M subset S agreeing with S above bit i
    for i in range(n):
        bit = np.int64(1 << i)
        sel = np.nonzero(sets & bit)[0]
        if len(sel) == 0:
            continue
        par = np.searchsorted(sets, sets[sel] ^ bit)
        g[sel] += g[par]
    gmin = min(g) if not use_int64 else int(g.min())
    garg = int(np.argmin(g)) if use_int64 else min(range(len(g)), key=lambda k: g[k])
    nneg = int(np.count_nonzero(g < 0)) if use_int64 else sum(1 for x in g if x < 0)
    # exact recomputation from the definition at sample sets (independent of the table f)
    rnd = random.Random(7)
    sample = sorted(set([0, garg] + [rnd.randrange(len(sets)) for _ in range(min(1500, len(sets)))]))
    bad = 0
    for k in sample:
        S = int(sets[k])
        q = psiT + sum((lam[i] for i in range(n) if S >> i & 1), F(0))
        q += sum((v for i, j, v in U if (S >> (i - 1) & 1) and (S >> (j - 1) & 1)), F(0))
        for I, J, v in WC:
            inI = all(S >> (a - 1) & 1 for a in I)
            inJ = all(S >> (b - 1) & 1 for b in J)
            if inI != inJ:
                q += v if inI else -v
        if q != F(int(g[k]), L):
            bad += 1
    minQ = F(int(gmin), L)
    Sstar = [i + 1 for i in range(n) if int(sets[garg]) >> i & 1]
    log('(A) exact recomputation from the definition at %d sets: %d mismatches' % (len(sample), bad))
    return dict(ok=(nneg == 0 and bad == 0 and enum_ok), nsets=len(sets), enum_ok=enum_ok, minQ=minQ, argmin=Sstar,
                nneg=nneg, bad=bad, Qempty=F(int(g[0]), L))

# =====================================================================================
# 6. (P): Phi(t) >= 1 for all t >= 0
# =====================================================================================
GRID = 40
INV_SQRT2_UP = F(70711, 100000)
CTX = {}

def _phis(A):
    X = A << (PB - GRID)
    J1 = j0_enc(X, X)
    lo = CTX['C0lo'] + iqmul_lo(CTX['psi0'], *J1)
    hi = CTX['C0hi'] + iqmul_hi(CTX['psi0'], *J1)
    for w, Rlo, Rhi in CTX['groups']:
        J = j0_enc((Rlo * A) >> GRID, cdiv(Rhi * A, 1 << GRID))
        lo += iqmul_lo(w, *J); hi += iqmul_hi(w, *J)
    return lo, hi

def _fail_status(A):
    up, Av = min((_phis(Av)[1], Av) for Av in [A] + [A + (1 << k) for k in range(51)])
    if up < ONE:
        return 'sweep stopped at t = %.12g; VIOLATION PROVEN: Phi <= 1 - %.6e at t = %.12g' % (
            A / 2 ** GRID, (ONE - up) / ONE, Av / 2 ** GRID)
    return 'certified lower bound at t = %.12g does not reach 1 (cannot certify)' % (A / 2 ** GRID)

def sweep_chunk(job):
    """chord: Phi(t) >= min(Phi(a), Phi(b)) - M2 h^2/8;  lipschitz: Phi(t) >= (Phi(a)+Phi(b))/2 - M1 h/2."""
    A0, A1, mode = job
    for k in NEVAL: NEVAL[k] = 0
    M1, M2 = CTX['M1'], CTX['M2']
    A, La = A0, _phis(A0)[0]
    steps = 0
    min_margin = None; min_at = None
    while A < A1:
        slack = La - ONE
        if slack <= 0:
            return dict(ok=False, status=_fail_status(A))
        h = math.sqrt(2.0 * (slack / ONE) / float(M2)) if mode == 'chord' else (slack / ONE) / float(M1)
        H = max(1, min(int(h * (1 << GRID)), 1 << GRID, A1 - A))
        while True:
            B = A + H
            Lb = _phis(B)[0]
            if mode == 'chord':
                lb = min(La, Lb) - iqmul_hi(M2, 0, cdiv(H * H << PB, 8 << (2 * GRID)))
            else:
                lb = (La + Lb) // 2 - iqmul_hi(M1, 0, cdiv(H << PB, 2 << GRID))
            if lb >= ONE:
                break
            H //= 2
            if H == 0:
                return dict(ok=False, status='step size underflow; ' + _fail_status(A))
        if min_margin is None or lb - ONE < min_margin:
            min_margin, min_at = lb - ONE, (A, B)
        A, La = B, Lb
        steps += 1
    return dict(ok=True, steps=steps, min_margin=min_margin, min_at=min_at, neval=dict(NEVAL))

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
        return False, '; '.join(r['status'] for r in res if not r['ok']), None
    mm = min(res, key=lambda r: r['min_margin'])
    margin = F(mm['min_margin'], ONE)
    return True, ('T = %s, %d steps, %d J0 enclosures per grid point (%d by power series, %d by large-argument '
                  'expansion); smallest certified lower bound of Phi - 1 on a step: %.4e on [%.6g, %.6g]; %.1f s'
                  % (T, sum(r['steps'] for r in res), 1 + len(CTX['groups']),
                     sum(r['neval']['series'] for r in res), sum(r['neval']['hankel'] for r in res),
                     float(margin), mm['min_at'][0] / 2 ** GRID, mm['min_at'][1] / 2 ** GRID,
                     time.time() - t0)), margin

def check_phi(cert, geo, opts, verdict):
    groups = {}
    for i, j, v in cert['u']:
        s = geo.d2(i, j)
        groups[s] = groups.get(s, F(0)) + v          # u terms with exactly equal distance are merged
    enc = {}
    for s in groups:
        slo, shi = k_enclose(s)
        assert slo > 0
        enc[s] = (slo, shi, isqrt_lo(slo), isqrt_hi(shi))
    sumlam = sum(cert['lam'], F(0))
    C0 = -sumlam
    Phi0 = cert['psi0'] - sumlam - sum((v for _, _, v in cert['u']), F(0))
    M2 = abs(cert['psi0']) + sum((abs(V) * fx_to_F(enc[s][1]) for s, V in groups.items()), F(0))
    M1 = INV_SQRT2_UP * (abs(cert['psi0']) + sum((abs(V) * fx_to_F(enc[s][3]) for s, V in groups.items()), F(0)))
    CTX.clear()
    CTX.update(psi0=cert['psi0'], C0lo=fx_floor(C0), C0hi=fx_ceil(C0), M1=M1, M2=M2,
               groups=[(-V, enc[s][2], enc[s][3]) for s, V in groups.items()])
    log('(P) Phi(0) = %s (~%.10g); Phi(oo) = -sum(lam) = %.10g; %d u entries in %d distance classes; '
        'M2 = %.6g' % (Phi0, float(Phi0), float(-sumlam), len(cert['u']), len(groups), float(M2)))
    if -sumlam <= 1:
        verdict('(P) Phi(oo) = -sum(lam) > 1', False)
        return None
    Ct = sqrt_up(F(2) / PI_LO_Q) * (abs(cert['psi0'])
                                     + sum((abs(V) * inv_sqrt_up(fx_to_F(enc[s][2])) for s, V in groups.items()), F(0)))
    Tneed = (Ct / (-sumlam - 1)) ** 2
    T = F(math.ceil(float(Tneed) * 1.02 / 10) * 10)
    tail = -sumlam - Ct * inv_sqrt_up(T)
    verdict('(P) tail: Phi(t) >= %.8g for t >= T = %s (T >= %.2f needed; |J0(x)| <= sqrt(2/(pi x)))'
            % (float(tail), T, float(Tneed)), tail >= 1)
    if tail < 1:
        return None
    ok, info, margin = sweep(T, 'chord', opts['workers'])
    verdict('(P) Phi(t) >= 1 on [0, T] (t = 0 included)', ok, info)
    if ok and opts['lipschitz']:
        ok2, info2, _ = sweep(T, 'lipschitz', opts['workers'])
        verdict('(P) redundant sweep with |J1| <= 1/sqrt2', ok2, info2)
    return dict(margin=margin, tail=tail, T=T, Phi0=Phi0) if ok else None
