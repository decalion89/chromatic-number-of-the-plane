"""Rigorous check of a certificate for m_1(R^2) <= psiT: the inclusion-exclusion bound of Ambrus, Csiszarik, Matolcsi,
Varga and Zsamboki for a finite configuration X with coordinates in Q(sqrt3, sqrt11), in exact arithmetic and Arb.
Certificate (side 'M1'): points (real and imaginary parts as [a, b, c, d] for a + b sqrt3 + c sqrt11 + d sqrt33,
rationals as strings), psi0, psiT, lam (one value per point), u ([i, j, value], 1-based, pairs not at distance 1),
wc ([I, J, value]: lists of distinct point indices, congruent in this order). Checks
  (G) the exact unit-distance graph of X; every index in 1..n; no u on a unit pair; every wc entry congruent
      (all squared distances equal, exactly);
  (A) Q(S) = psiT + sum_{i in S} lam_i + sum_{i<j in S} u_ij + sum wc([I in S] - [J in S]) >= 0 for every independent
      set S of the unit-distance graph, including the empty set (exact integers; the number of independent sets is
      cross-checked by the recursion I(G) = I(G - v) + I(G - N[v]));
  (P) Phi(t) = psi0 J0(t) - sum lam - sum u_ij J0(r_ij t) >= 1 for all t >= 0 (Arb boxes on [0, T], and
      |J0(x)| <= sqrt(2/(pi x)) beyond T).
Then every measurable planar set with no two points at distance 1 has upper density at most psiT.
Facts about J0 used: J0(x) = (1/pi) int_0^pi cos(x sin s) ds, so |J0| <= 1 and |J0''| <= 1; J0 >= -0.4028 (proved
here); |J0(x)| <= sqrt(2/(pi x)) for x > 0 (from Nicholson's formula, DLMF 10.9.30; Watson, Bessel Functions, 13.74).
The checks do not use assert, so they also run under python3 -O. Exit code 0 and a last line CERTIFIED: ... when every
check passes, exit code 1 and REJECTED: ... otherwise. Self-contained: needs numpy and python-flint (Arb).
usage: python3 check_m1.py CERT.json"""
import sys, json, time
from fractions import Fraction as F
from math import lcm
import numpy as np
from flint import arb, ctx
ctx.prec = 96
PI = arb.pi()

def fail(msg):
    print(f"REJECTED: {msg}", flush=True)
    sys.exit(1)

# ---- exact arithmetic in K = Q(sqrt3, sqrt11), basis 1, sqrt3, sqrt11, sqrt33
class K:
    __slots__ = ('v',)
    def __init__(self, v): self.v = tuple(F(x) for x in v)
    def __add__(s, o): return K(tuple(a + b for a, b in zip(s.v, o.v)))
    def __sub__(s, o): return K(tuple(a - b for a, b in zip(s.v, o.v)))
    def __mul__(s, o):
        a, b, c, d = s.v; e, f, g, h = o.v
        return K((a*e + 3*b*f + 11*c*g + 33*d*h, a*f + b*e + 11*c*h + 11*d*g,
                  a*g + c*e + 3*b*h + 3*d*f, a*h + d*e + b*g + c*f))
    def __eq__(s, o): return s.v == o.v
    def __hash__(s): return hash(s.v)
ZERO = K((0, 0, 0, 0)); ONE = K((1, 0, 0, 0))
def dist2(p, q):
    dx = p[0] - q[0]; dy = p[1] - q[1]
    return dx * dx + dy * dy

# ---- Arb helpers
def A(x):
    x = F(x); return arb(x.numerator) / arb(x.denominator)
def lo(x): return x.lower()
def hi(x): return x.upper()
M0 = F(-4028, 10000)
def j0_box(a, b, use_global=True):
    """Arb lower and upper bounds for J0 on [a, b] (Fractions, 0 <= a <= b): Taylor at the midpoint with |J0''| <= 1,
    the envelope sqrt(2/(pi a)) and the global range [-0.4028, 1]"""
    m = (a + b) / 2; r = (b - a) / 2
    am, ar = A(m), A(r)
    j0m = am.bessel_j(0); j1m = am.bessel_j(1)
    w = abs(j1m) * ar + ar * ar / 2
    L = lo(j0m - w); U = hi(j0m + w)
    if a > 0:
        env = hi((A(2) / (PI * A(a))).sqrt())
        if env < -L: L = -env
        if env < U: U = env
    if use_global:
        if lo(A(M0)) > L: L = lo(A(M0))
        if U > 1: U = arb(1)
    return L, U
def check_j0min():
    """J0(x) >= -0.4028 for x >= 0: interval enclosures on [0, 4], envelope sqrt(2/(pi x)) < 0.4028 on [4, oo)"""
    if not hi((A(2) / (PI * A(4))).sqrt()) < -A(M0): return False
    for i in range(4000):
        L, _ = j0_box(F(i, 1000), F(i + 1, 1000), use_global=False)
        if not lo(L) >= A(M0): return False
    return True
def K_arb(x):
    a, b, c, d = x.v
    return A(a) + A(b) * arb(3).sqrt() + A(c) * arb(11).sqrt() + A(d) * arb(33).sqrt()
def exact(x):
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return F(m * 2 ** e) if e >= 0 else F(m, 2 ** -e)
def enclose(q):
    """rational lower and upper bounds for sqrt(q), q in K positive"""
    r = K_arb(q).sqrt(); l, u = exact(r.lower()), exact(r.upper())
    if not 0 < l <= u: fail('a distance could not be enclosed')
    return l, u
env = lambda x: (A(2) / (PI * x)).sqrt()

# ---- reading the certificate
def rational(x, what):
    if isinstance(x, bool) or not isinstance(x, (int, str)): fail(f'{what}: rationals must be integers or strings')
    try: return F(x)
    except (ValueError, ZeroDivisionError): fail(f'{what}: not a rational: {x!r}')
try:
    C = json.load(open(sys.argv[1]))
except (OSError, ValueError, IndexError) as e:
    fail(f'cannot read the certificate ({e})')
if not isinstance(C, dict) or C.get('side') != 'M1': fail("the certificate must have side 'M1'")
def kpt(p):
    if not (isinstance(p, list) and len(p) == 2 and all(isinstance(c, list) and len(c) == 4 for c in p)): fail('bad point')
    return (K(tuple(rational(x, 'point') for x in p[0])), K(tuple(rational(x, 'point') for x in p[1])))
if not isinstance(C.get('points'), list): fail('no points')
pts = [kpt(p) for p in C['points']]; n = len(pts)
if not 1 <= n <= 40: fail('between 1 and 40 points')
def index(a):
    """a 1-based point index, checked to lie in 1..n"""
    if not (isinstance(a, int) and not isinstance(a, bool) and 1 <= a <= n): fail(f'index out of range: {a!r}')
    return a - 1
d2 = {}
for i in range(n):
    for j in range(i + 1, n):
        d2[(i, j)] = d2[(j, i)] = dist2(pts[i], pts[j])
        if d2[(i, j)] == ZERO: fail(f'points {i + 1} and {j + 1} coincide')
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if d2[(i, j)] == ONE]
psi0, psiT = rational(C.get('psi0'), 'psi0'), rational(C.get('psiT'), 'psiT')
if not isinstance(C.get('lam'), list) or len(C['lam']) != n: fail('one lam per point')
lam = [rational(x, 'lam') for x in C['lam']]
U = {}
for e in C.get('u', []):
    if not (isinstance(e, list) and len(e) == 3): fail('a u entry is not [i, j, value]')
    a, b = index(e[0]), index(e[1])
    if not a < b: fail('u pairs must have i < j')
    if (a, b) in U: fail('a u pair is listed twice')
    if d2[(a, b)] == ONE: fail(f'u on the unit pair {a + 1}, {b + 1}')
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
rec(0, 0, 0)
memo = {}
def count(alive):
    """number of independent sets of the graph induced on the vertex set alive (a bitmask)"""
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
dt = np.int64 if total < 2 ** 62 else object                   # Python integers when int64 could overflow
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
      f"(Q(empty) = psiT = {float(psiT):.10f}){'' if dt is np.int64 else ' [Python integers]'}", flush=True)
if qmin < 0:
    S = [i + 1 for i in range(n) if atoms[k] >> i & 1]
    fail(f'Q(S) = {F(qmin, den)} < 0 for the independent set S = {S}')
del Q, M

# ---- (P)
if not check_j0min(): fail('J0 >= -0.4028 could not be certified')
G = {}
for (i, j), v in U.items():
    if v: G[d2[(i, j)]] = G.get(d2[(i, j)], F(0)) + v
terms = [(A(-v), enclose(q), q) for q, v in G.items() if v]
c0 = A(-sum(lam))
def lower(t0, t1):
    Lj, Uj = j0_box(t0, t1)
    v = c0 + (A(psi0) * Lj if psi0 >= 0 else A(psi0) * Uj)
    for cf, (rl, rh), _ in terms:
        Lr, Ur = j0_box(rl * t0, rh * t1); v = v + (cf * Lr if lo(cf) >= 0 else cf * Ur)
    return lo(v)
def phi_ball(t):
    """an Arb ball containing Phi(t), for a rational t"""
    at = A(t); v = c0 + A(psi0) * at.bessel_j(0)
    for cf, _, q in terms:
        v = v + cf * (K_arb(q).sqrt() * at).bessel_j(0)
    return v
T = F(100)
while True:
    tl = c0 - abs(A(psi0)) * env(A(T))
    for cf, (rl, rh), _ in terms: tl = tl - abs(cf) * env(A(rl) * A(T))
    if lo(tl) >= 1: break
    T *= 2
    if T > 10 ** 9: fail('the tail bound fails: -sum lam is not above 1')
stack = [(F(0), T)]; nb = 0; t_start = time.time()
while stack:
    t0, t1 = stack.pop(); nb += 1
    if lower(t0, t1) >= 1: continue
    if t1 - t0 <= F(1, 10 ** 15):
        b = phi_ball(t0)
        if hi(b) < 1: fail(f'Phi({float(t0)}) < 1 (proved: Phi is in [{float(lo(b)):.12f}, {float(hi(b)):.12f}])')
        fail(f'Phi >= 1 could not be certified near t = {float(t0)} (Phi there in [{float(lo(b)):.12f}, {float(hi(b)):.12f}])')
    m = (t0 + t1) / 2; stack.append((m, t1)); stack.append((t0, m))
print(f"(P) Phi >= 1 on [0, {T}] ({nb} boxes, {time.time() - t_start:.0f}s) and on [{T}, oo) by the envelope", flush=True)
c = -((-psiT * 10 ** 8) // 1)
print(f"CERTIFIED: every measurable planar set avoiding the distance 1 has upper density at most psiT = "
      f"{c // 10 ** 8}.{c % 10 ** 8:08d} (rounded up) = 1/{1 / float(psiT):.5f}")
