"""Rigorous check of a certificate for m_1(R^3) <= B: the inclusion-exclusion bound of Ambrus, Csiszarik, Matolcsi,
Varga and Zsamboki, carried to R^3 (kernel Om(x) = sin(x)/x) and combined with the atom kappa(0) = delta^2 of the
spectral measure, for a finite configuration X whose points are q/sqrt2 with q in Q^3 (exact rational arithmetic, Arb).
Certificate (side 'M1R3', scale '1/sqrt2'): points ([x, y, z], rationals as strings), psi0, psiT, m, lam (one per
point), u ([i, j, value], 1-based, pairs not at distance 1), wc ([I, J, value]: congruent in this order). Checks
  (G) the exact unit-distance graph (squared distances |q - q'|^2/2 are rational); indices in range; no u on a unit
      pair; every wc entry congruent (all squared distances equal; a distance-preserving bijection between finite subsets
      of R^3 extends to an isometry of R^3);
  (A) Q(S) = psiT + sum_{i in S} lam_i + sum_{i<j in S} u_ij + sum wc([I in S] - [J in S]) >= 0 for every independent
      set S, including the empty set (exact integers; count cross-checked by I(G) = I(G - v) + I(G - N[v]));
  (P) Phi(t) = psi0 Om(t) - sum lam - sum u_ij Om(r_ij t) >= m for all t > 0 (Arb boxes on (0, T]; |Om(x)| <= 1/x
      beyond T);
  (B) c = Phi(0) - m = psi0 - sum lam - sum u - m > 0.
Then for every periodic measurable 1-avoiding A in R^3 of density delta: m delta + c delta^2 <= psiT (the average of
the patterns of A on congruent copies of X over translations and Haar measure on O(3), reflections included; the atom
kappa(0) = delta^2 of the spectral measure), so delta is at most the positive root B of c x^2 + m x - psiT; by
periodic approximation (Oliveira Filho and Vallentin 2010, Section 2) the same bound holds for the upper density of
every measurable 1-avoiding set. Facts about Om used: Om(x) = int_0^1 cos(xs) ds, so |Om| <= 1, |Om''| <= 1/3;
Om >= -0.2173 (proved here); |Om(x)| <= 1/x. No assert statements (runs under python3 -O). Exit code 0 and a last line
CERTIFIED: ... when every check passes; exit code 1 and REJECTED: ... otherwise. Needs numpy and python-flint.
usage: python3 check_m1r3.py CERT.json"""
import sys, json, time
from fractions import Fraction as F
from math import lcm, isqrt
import numpy as np
from flint import arb, ctx
ctx.prec = 96

def fail(msg):
    print(f"REJECTED: {msg}", flush=True)
    sys.exit(1)
def A(x):
    x = F(x); return arb(x.numerator) / arb(x.denominator)
def lo(x): return x.lower()
def hi(x): return x.upper()
def exact(x):
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return F(m * 2 ** e) if e >= 0 else F(m, 2 ** -e)
MIN = F(-2173, 10000)
def om_ball(x):
    """Arb ball for Om(x) = sin(x)/x, x an arb (may contain 0 only as the point 0)"""
    return x.sinc()
def om_box(a, b, use_global=True):
    """Arb lower and upper bounds for Om on [a, b] (Fractions, 0 <= a <= b): Taylor at the midpoint (|Om''| <= 1/3),
    the envelope 1/a, the global range [-0.2173, 1]"""
    mid = (a + b) / 2; r = (b - a) / 2
    am, ar = A(mid), A(r)
    v = am.sinc()
    if mid > 0:
        s, c = am.sin(), am.cos(); d1 = (am * c - s) / (am * am)
    else:
        d1 = arb(0)
    w = abs(d1) * ar + ar * ar / 6
    L = lo(v - w); U = hi(v + w)
    if a > 0:
        env = hi(arb(1) / A(a))
        if env < -L: L = -env
        if env < U: U = env
    if use_global:
        if lo(A(MIN)) > L: L = lo(A(MIN))
        if U > 1: U = arb(1)
    return L, U
def check_ommin():
    """Om(x) >= -0.2173 for x >= 0: boxes on [0, 5], and 1/x < 0.2173 for x >= 5"""
    if not hi(arb(1) / arb(5)) < -A(MIN): return False
    for i in range(5000):
        L, _ = om_box(F(i, 1000), F(i + 1, 1000), use_global=False)
        if not lo(L) >= A(MIN): return False
    return True
def rational(x, what):
    if isinstance(x, bool) or not isinstance(x, (int, str)): fail(f'{what}: rationals must be integers or strings')
    try: return F(x)
    except (ValueError, ZeroDivisionError): fail(f'{what}: not a rational: {x!r}')
try:
    C = json.load(open(sys.argv[1]))
except (OSError, ValueError, IndexError) as e:
    fail(f'cannot read the certificate ({e})')
if not isinstance(C, dict) or C.get('side') != 'M1R3' or C.get('scale') != '1/sqrt2':
    fail("the certificate must have side 'M1R3' and scale '1/sqrt2'")
if not isinstance(C.get('points'), list): fail('no points')
pts = []
for p in C['points']:
    if not (isinstance(p, list) and len(p) == 3): fail('a point is not [x, y, z]')
    pts.append(tuple(rational(c, 'point') for c in p))
n = len(pts)
if not 1 <= n <= 60: fail('between 1 and 60 points')
def index(a):
    if not (isinstance(a, int) and not isinstance(a, bool) and 1 <= a <= n): fail(f'index out of range: {a!r}')
    return a - 1
d2 = {}
for i in range(n):
    for j in range(i + 1, n):
        d2[(i, j)] = d2[(j, i)] = sum((a - b) ** 2 for a, b in zip(pts[i], pts[j])) / 2
        if d2[(i, j)] == 0: fail(f'points {i + 1} and {j + 1} coincide')
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if d2[(i, j)] == 1]
psi0, psiT, m = rational(C.get('psi0'), 'psi0'), rational(C.get('psiT'), 'psiT'), rational(C.get('m'), 'm')
if not isinstance(C.get('lam'), list) or len(C['lam']) != n: fail('one lam per point')
lam = [rational(x, 'lam') for x in C['lam']]
U = {}
for e in C.get('u', []):
    if not (isinstance(e, list) and len(e) == 3): fail('a u entry is not [i, j, value]')
    a, b = index(e[0]), index(e[1])
    if not a < b: fail('u pairs must have i < j')
    if (a, b) in U: fail('a u pair is listed twice')
    if d2[(a, b)] == 1: fail(f'u on the unit pair {a + 1}, {b + 1}')
    U[(a, b)] = rational(e[2], 'u')
WC = []
for e in C.get('wc', []):
    if not (isinstance(e, list) and len(e) == 3 and isinstance(e[0], list) and isinstance(e[1], list)):
        fail('a wc entry is not [I, J, value]')
    I = [index(a) for a in e[0]]; J = [index(b) for b in e[1]]
    if not (1 <= len(I) == len(J) and len(set(I)) == len(I) and len(set(J)) == len(J)): fail(f'bad congruence entry {e[:2]}')
    for a in range(len(I)):
        for b in range(a + 1, len(I)):
            if d2[(I[a], I[b])] != d2[(J[a], J[b])]: fail(f'not congruent: {e[0]}, {e[1]}')
    WC.append((I, J, rational(e[2], 'wc')))
print(f"(G) {n} points, {len(edges)} unit pairs, {len(WC)} congruences checked exactly", flush=True)

# ---- (A)
nbm = [0] * n
for i, j in edges: nbm[i] |= 1 << j; nbm[j] |= 1 << i
atoms = []
def rec(start, cur, mask):
    atoms.append(cur)
    for v in range(start, n):
        if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nbm[v])
sys.setrecursionlimit(10000)
rec(0, 0, 0)
memo = {}
def count(alive):
    if alive == 0: return 1
    if alive in memo: return memo[alive]
    v = alive.bit_length() - 1
    r = count(alive & ~(1 << v)) + count(alive & ~(1 << v) & ~nbm[v])
    memo[alive] = r; return r
if count((1 << n) - 1) != len(atoms) or len(set(atoms)) != len(atoms) or atoms[0] != 0:
    fail('the enumeration of the independent sets is inconsistent')
if any(a & nbm[i] for a in atoms for i in range(n) if a >> i & 1): fail('an enumerated set is not independent')
den = 1
for x in [psiT] + lam + list(U.values()) + [v for _, _, v in WC]: den = lcm(den, x.denominator)
num = lambda x: int(x * den)
total = abs(num(psiT)) + sum(abs(num(x)) for x in lam) + sum(abs(num(x)) for x in U.values()) \
    + 2 * sum(abs(num(v)) for _, _, v in WC)
dt = np.int64 if total < 2 ** 62 else object
M = np.array(atoms, dtype=np.int64); Q = np.full(len(M), num(psiT), dtype=dt)
for i in range(n):
    if lam[i]: Q += num(lam[i]) * ((M >> i) & 1).astype(dt)
for (i, j), v in U.items():
    if v: Q += num(v) * (((M >> i) & 1) & ((M >> j) & 1)).astype(dt)
for I, J, v in WC:
    mi = sum(1 << a for a in I); mj = sum(1 << b for b in J)
    Q += num(v) * (((M & mi) == mi).astype(dt) - ((M & mj) == mj).astype(dt))
k = int(np.argmin(Q)); qmin = int(Q[k])
print(f"(A) {len(atoms)} independent sets (count cross-checked): min Q = {float(F(qmin, den)):.3e} "
      f"(Q(empty) = psiT = {float(psiT):.3e})", flush=True)
if qmin < 0:
    fail(f'Q(S) = {F(qmin, den)} < 0 for the independent set S = {[i + 1 for i in range(n) if atoms[k] >> i & 1]}')
del Q, M

# ---- (B)
c = psi0 - sum(lam) - sum(U.values()) - m
if not c > 0: fail(f'c = Phi(0) - m = {c} is not positive')
print(f"(B) c = Phi(0) - m = {float(c):.9f} > 0, m = {float(m):.9f}", flush=True)

# ---- (P)
if not check_ommin(): fail('Om >= -0.2173 could not be certified')
G = {}
for (i, j), v in U.items():
    if v: G[d2[(i, j)]] = G.get(d2[(i, j)], F(0)) + v
def enclose(q):
    r = A(q).sqrt(); l, u = exact(r.lower()), exact(r.upper())
    if not 0 < l <= u: fail('a distance could not be enclosed')
    return l, u
terms = [(A(-v), enclose(q), q, -v >= 0) for q, v in G.items() if v]
c0 = A(-sum(lam))
def lower(t0, t1):
    Lo, Uo = om_box(t0, t1)
    v = c0 + (A(psi0) * Lo if psi0 >= 0 else A(psi0) * Uo)
    for cf, (rl, rh), _, pos in terms:
        Lr, Ur = om_box(rl * t0, rh * t1); v = v + (cf * Lr if pos else cf * Ur)     # sign of the exact rational
    return lo(v)
def phi_ball(t):
    at = A(t); v = c0 + A(psi0) * at.sinc()
    for cf, _, q, _ in terms: v = v + cf * (A(q).sqrt() * at).sinc()
    return v
T = F(100)
while True:
    tl = c0 - abs(A(psi0)) / A(T)
    for cf, (rl, rh), _, _ in terms: tl = tl - abs(cf) / (A(rl) * A(T))
    if lo(tl) >= A(m): break
    T *= 2
    if T > 10 ** 9: fail('the tail bound fails: -sum lam is not above m')
stack = [(F(0), T)]; nb = 0; t_start = time.time(); Am = A(m)
while stack:
    t0, t1 = stack.pop(); nb += 1
    if lower(t0, t1) >= Am: continue
    if t1 - t0 <= F(1, 10 ** 15):
        b = phi_ball(t0)
        if hi(b) < Am: fail(f'Phi({float(t0)}) < m (proved)')
        fail(f'Phi >= m could not be certified near t = {float(t0)}')
    mid = (t0 + t1) / 2; stack.append((mid, t1)); stack.append((t0, mid))
print(f"(P) Phi >= m on [0, {T}] ({nb} boxes, {time.time() - t_start:.0f}s) and on [{T}, oo) by |Om(x)| <= 1/x", flush=True)
if psiT < 0: fail('psiT < 0')
# a rational Bq with Bq > 0, Bq >= -m/(2c) and c Bq^2 + m Bq - psiT >= 0 is at least the positive root
Bf = (-float(m) + (float(m) ** 2 + 4 * float(c) * float(psiT)) ** 0.5) / (2 * float(c))
k = 10 ** 8; Bq = F(int(Bf * k) + 1, k)
while not (c * Bq * Bq + m * Bq - psiT >= 0 and Bq >= -m / (2 * c)): Bq += F(1, k)
bi = int(Bq * k)
print(f"CERTIFIED: every measurable set in R^3 avoiding the distance 1 has upper density at most {bi // k}.{bi % k:08d}"
      f" = 1/{1 / float(Bq):.5f} "
      f"(m delta + c delta^2 <= psiT with m = {float(m):.9f}, c = {float(c):.9f}, psiT = {float(psiT):.3e})")
