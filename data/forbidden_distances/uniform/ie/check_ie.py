"""Rigorous check of the inclusion-exclusion certificates (Propositions 19 and 21 and Lemma 20 of the paper), built on the 23 points
X23 of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki (x23.py), in exact arithmetic of Q(sqrt3, sqrt11) and Arb.

A-side certificate (side 'A'): D, d0, mu, psi0, psiT, eta, tau, zeta, lam (23 values), u ([i, j, value], 1-based,
non-unit pairs), wc ([I, J, value]: index lists congruent in this order), and, when zeta > 0, 'B': the file of a
B-side certificate. Checks
  (G) the unit pairs of X23 (47, exact), every congruence (all squared distances equal, exactly), no u on a unit pair,
      eta, tau, zeta >= 0, mu >= 0.4028 >= -min J0 (proved by Arb on [0, 4] and the envelope beyond);
  (A) Q(S) = psiT + sum_S lam + sum_{pairs in S} u + sum wc([I in S] - [J in S]) >= 0 for every nonempty independent
      set S of the unit-distance graph of X23, and Q(empty) = psiT >= zeta (exact);
  (Z) Phi0 = psi0 - sum lam - sum u + eta + tau (2 d0 - 1) >= 1 (exact);
  (P) Phi(t) = psi0 J0(t) - sum lam - sum u J0(r t) - eta e(Dt) + 2 d0 tau - zeta L_B e(Dt) (529 - H(t)) >= 1 for all
      t > 0, with e(y) = min(mu, sqrt(2/(pi y))) and H(t) = sum_{i,j} J0(r_ij t) (Arb boxes on (0, T], envelope beyond).
  Then m_{1,d}(R^2) <= psiT + tau d0^2 - zeta (1 - psiT_B) for every d >= D.
B-side certificate (side 'B'): psi0, psiT, pi (one number, or [i, j, value] for each of the 47 unit pairs), lam, u, wc.
  Checks the congruences, Q_B(T) = psiT + sum_T lam + sum_{non-unit pairs in T} u + sum_{unit pairs e in T} pi_e
  + sum wc([I in T] - [J in T]) >= 0 for ALL 2^23 subsets T (exact integers),
  and Phi_B(t) = psi0 J0(t) - sum lam - sum u J0(r t) >= 1 for all t >= 0. Then a set B (or a family of sets) whose
  averaged autocorrelation at distance 1 is eps has density at most psiT_B + L_B eps, L_B = max(0, psi0 + sum pi_e).
Every point index must lie in 1..23. The checker also asks 0 < d0 < 1 and pi_e >= 0, which the paper does not need.
usage: python3 check_ie.py CERT.json [k]        (k: the number of colours to certify, default from the bound)"""
import sys, os, json, time
from fractions import Fraction as F
from math import lcm
import numpy as np
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.dirname(HERE))
from kfield import K, ONE, dist2
from x23 import P
from uniform import A, lo, hi, j0_box, check_j0min
ctx.prec = 96
PI = arb.pi()
n = len(P)
pts = [(K(re), K(im)) for re, im in P]
d2 = {(i, j): dist2(pts[i], pts[j]) for i in range(n) for j in range(i + 1, n)}
edges = [(i, j) for (i, j), v in d2.items() if v == ONE]
assert len(edges) == 47, len(edges)
def K_arb(x):
    a, b, c, d = x.v
    return A(a) + A(b) * arb(3).sqrt() + A(c) * arb(11).sqrt() + A(d) * arb(33).sqrt()
def exact(x):
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return F(m * 2 ** e) if e >= 0 else F(m, 2 ** -e)
def enclose(q):
    """rational lower and upper bounds for sqrt(q), q in K positive"""
    r = K_arb(q).sqrt(); l, u = exact(r.lower()), exact(r.upper()); assert 0 < l <= u; return l, u
def index(a):
    """a 1-based point index, checked to lie in 1..23 (an index outside would make [I in S] vacuous)"""
    assert isinstance(a, int) and not isinstance(a, bool) and 1 <= a <= n, ('index out of range', a)
    return a - 1
def parse(C):
    lam = [F(x) for x in C['lam']]; assert len(lam) == n
    U = {(index(i), index(j)): F(v) for i, j, v in C['u']}
    assert len(U) == len(C['u']), 'a pair is listed twice'
    WC = [([index(a) for a in I], [index(b) for b in J], F(v)) for I, J, v in C['wc']]
    assert all(0 <= i < j < n for i, j in U) and not any(d2[e] == ONE for e in U), 'u on a unit pair'
    for I, J, v in WC:
        assert 1 <= len(I) == len(J) and len(set(I)) == len(I) and len(set(J)) == len(J)
        for a in range(len(I)):
            for b in range(a + 1, len(I)):
                p, q = sorted((I[a], I[b])); r, s = sorted((J[a], J[b]))
                assert d2[(p, q)] == d2[(r, s)], ('not congruent', I, J)
    return lam, U, WC
def grouped(U):
    G = {}
    for (i, j), v in U.items():
        if v: G[d2[(i, j)]] = G.get(d2[(i, j)], F(0)) + v
    return [(A(-v), enclose(q)) for q, v in G.items() if v]
env = lambda x: (A(2) / (PI * x)).sqrt()
def sweep(lower, T, label):
    stack = [(F(0), T)]; nbox = 0; t_start = time.time()
    while stack:
        t0, t1 = stack.pop(); nbox += 1
        lb = lower(t0, t1)
        if lb >= 1: continue
        assert t1 - t0 > F(1, 10 ** 15), f'{label} < 1 near t = {float(t0)} (lower bound {float(lb)})'
        m = (t0 + t1) / 2; stack.append((m, t1)); stack.append((t0, m))
    return nbox, time.time() - t_start
def check_B(path):
    C = json.load(open(path)); assert C['side'] == 'B'
    psi0, psiT = F(C['psi0']), F(C['psiT'])
    if isinstance(C['pi'], list):                          # one penalty per unit pair
        PIe = {(index(i), index(j)): F(v) for i, j, v in C['pi']}
        assert len(PIe) == len(C['pi']) == len(edges) and set(PIe) == set(edges), \
            'the penalties must be given on exactly the 47 unit pairs'
    else:
        PIe = {e: F(C['pi']) for e in edges}
    assert all(v >= 0 for v in PIe.values())
    lam, U, WC = parse(C)
    L = 1
    for x in [psiT] + list(PIe.values()) + lam + list(U.values()) + [v for _, _, v in WC]: L = lcm(L, x.denominator)
    num = lambda x: int(x * L)
    total = abs(num(psiT)) + sum(abs(num(x)) for x in PIe.values()) + sum(abs(num(x)) for x in lam) \
        + sum(abs(num(x)) for x in U.values()) + 2 * sum(abs(num(v)) for _, _, v in WC)
    assert total < 2 ** 62, 'the integers would overflow'
    T = np.arange(1 << n, dtype=np.int64); Q = np.full(1 << n, num(psiT), dtype=np.int64)
    for i in range(n): Q += num(lam[i]) * ((T >> i) & 1)
    for (i, j), v in U.items(): Q += num(v) * (((T >> i) & 1) & ((T >> j) & 1))
    for (i, j), v in PIe.items(): Q += num(v) * (((T >> i) & 1) & ((T >> j) & 1))
    for I, J, v in WC:
        mi = sum(1 << a for a in I); mj = sum(1 << b for b in J)
        Q += num(v) * (((T & mi) == mi).astype(np.int64) - ((T & mj) == mj).astype(np.int64))
    qmin = int(Q.min()); del Q, T
    print(f"  (B) all {1 << n} subsets of X23: min Q_B = {F(qmin, L)}", flush=True); assert qmin >= 0
    terms = grouped(U); c0 = A(-sum(lam))
    def lower(t0, t1):
        Lj, Uj = j0_box(t0, t1)
        v = c0 + (A(psi0) * Lj if psi0 >= 0 else A(psi0) * Uj)
        for cf, (rl, rh) in terms:
            Lr, Ur = j0_box(rl * t0, rh * t1); v = v + (cf * Lr if lo(cf) >= 0 else cf * Ur)
        return lo(v)
    TB = F(100)
    while True:
        tl = c0 - abs(A(psi0)) * env(A(TB))
        for cf, (rl, rh) in terms: tl = tl - abs(cf) * env(A(rl) * A(TB))
        if lo(tl) >= 1: break
        TB *= 2; assert TB < 10 ** 9
    nb, sec = sweep(lower, TB, 'Phi_B')
    LB = max(F(0), psi0 + sum(PIe.values()))
    print(f"  (B) Phi_B >= 1 on [0, {TB}] ({nb} boxes, {sec:.0f}s) and beyond by the envelope; "
          f"psiT_B = {float(psiT):.8f}, L_B = {float(LB):.4f}", flush=True)
    return psiT, LB
assert check_j0min()
C = json.load(open(sys.argv[1]))
if C['side'] == 'B':
    psiTB, LB = check_B(sys.argv[1]); print(f"B-side certificate {sys.argv[1]}: VALID"); sys.exit(0)
D, d0, mu = F(C['D']), F(C['d0']), F(C['mu'])
psi0, psiT, eta, tau, zeta = F(C['psi0']), F(C['psiT']), F(C['eta']), F(C['tau']), F(C.get('zeta', '0'))
assert eta >= 0 and tau >= 0 and zeta >= 0 and D > 0 and 0 < d0 < 1 and mu >= F(4028, 10000)
lam, U, WC = parse(C)
print(f"(G) {len(edges)} exact unit edges, {len(WC)} congruences checked exactly, J0 >= -0.4028 proved", flush=True)
psiTB = LB = F(0)
if zeta:
    psiTB, LB = check_B(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), C['B']))
nbm = [0] * n
for i, j in edges: nbm[i] |= 1 << j; nbm[j] |= 1 << i
atoms = []
def rec(start, cur, mask):
    atoms.append(cur)
    for v in range(start, n):
        if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nbm[v])
rec(0, 0, 0); assert atoms[0] == 0 and len(atoms) == 13552
La = 1
for x in [psiT] + lam + list(U.values()) + [v for _, _, v in WC]: La = lcm(La, x.denominator)
M = np.array(atoms, dtype=np.int64); Qa = [int(psiT * La)] * len(atoms)
def add(mask, val):
    iv = int(val * La)
    for a in np.flatnonzero((M & mask) == mask): Qa[a] += iv
for i in range(n):
    if lam[i]: add(1 << i, lam[i])
for (i, j), v in U.items():
    if v: add((1 << i) | (1 << j), v)
for I, J, v in WC:
    if v: add(sum(1 << i for i in I), v); add(sum(1 << j for j in J), -v)
qmin = min(Qa[1:])
print(f"(A) {len(atoms)} independent atoms: min Q over nonempty atoms = {float(F(qmin, La)):.3e}, "
      f"Q(empty) - zeta = {float(psiT - zeta):.3e}", flush=True)
assert qmin >= 0 and psiT >= zeta
Phi0 = psi0 - sum(lam) - sum(U.values()) + eta + tau * (2 * d0 - 1)
print(f"(Z) Phi0 = {float(Phi0):.9f}", flush=True); assert Phi0 >= 1
terms = grouped(U)
mult = {}
for q in d2.values(): mult[q] = mult.get(q, 0) + 1
Hterms = [(m, (F(1), F(1)) if q == ONE else enclose(q)) for q, m in mult.items()]
c0 = A(-sum(lam) + 2 * d0 * tau); zL = A(zeta * LB)
def e_up(t0):
    if t0 <= 0: return A(mu)
    s = (A(2) / (PI * A(D) * A(t0))).sqrt()
    return s if hi(s) < A(mu) else A(mu)
def lower(t0, t1):
    Lj, Uj = j0_box(t0, t1); e0 = e_up(t0)
    v = c0 + (A(psi0) * Lj if psi0 >= 0 else A(psi0) * Uj) - A(eta) * e0
    for cf, (rl, rh) in terms:
        Lr, Ur = j0_box(rl * t0, rh * t1); v = v + (cf * Lr if lo(cf) >= 0 else cf * Ur)
    if zeta:
        Hlo = A(n)
        for m, (rl, rh) in Hterms:
            Lr, _ = j0_box(rl * t0, rh * t1); Hlo = Hlo + 2 * m * Lr
        gap = A(n * n) - Hlo
        if hi(gap) > n * n: gap = A(n * n)                   # 529 - H <= 529 because H >= 0
        v = v - zL * e0 * gap
    return lo(v)
T = F(100)
while True:
    tl = c0 - abs(A(psi0)) * env(A(T)) - A(eta) * e_up(T) - zL * e_up(T) * (n * n)
    for cf, (rl, rh) in terms: tl = tl - abs(cf) * env(A(rl) * A(T))
    if lo(tl) >= 1: break
    T *= 2; assert T < 10 ** 9, 'the tail bound fails'
nb, sec = sweep(lower, T, 'Phi')
print(f"(P) Phi >= 1 on (0, {T}] ({nb} boxes, {sec:.0f}s) and on [{T}, oo) by the envelope", flush=True)
B = psiT + tau * d0 * d0 - zeta * (1 - psiTB)
k = int(sys.argv[2]) if len(sys.argv) > 2 else int(1 / B) + 1
ok = B < F(1, k - 1)
c = -((-B * 10 ** 8) // 1)                                 # the bound rounded up at the 8th decimal
print(f"bound for every d >= {D}: {c // 10 ** 8}.{c % 10 ** 8:08d} (rounded up) = 1/{1 / float(B):.4f} -> "
      f"{'CERTIFIED: chi(R^2, {1, d}) >= %d for every real algebraic d >= %s' % (k, D) if ok else 'NOT below 1/%d' % (k - 1)}")
sys.exit(0 if ok else 1)
