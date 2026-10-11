#!/usr/bin/env python3
"""
Independent rigorous checker for "M1R3" certificates (unit-distance avoiding
sets in R^3, density bound m*d + c*d^2 <= psiT).

Usage:  python3 -I check_m1r3_indep.py CERT.json

Exit code 0 and last line "ACCEPTED: ..." if every check passes,
exit code 1 and last line "REJECTED: <reason>" otherwise.
No assert statements are used (safe under python3 -O).

Checks
  (G) exact geometry: squared distances |q-q'|^2/2 (rational), unit graph,
      distinct points, u-pairs not unit, wc tuples congruent.
  (A) Q(S) >= 0 for every independent set S of the unit graph (exact ints).
      Independent sets are listed by a DFS; their number is cross-checked by
      a memoised deletion recursion I(G)=I(G-v)+I(G-N[v]) and (n<=24) by a
      brute-force sweep over all 2^n subsets.
  (B) c = psi0 - sum lam - sum u - m > 0 and psiT >= 0 (exact).
  (P) Phi(t) = psi0*Om(t) - sum lam - sum u_ij Om(r_ij t) >= m on (0,inf),
      Om(x) = sin(x)/x, rigorous (see the section "Bounds for Omega").

-------------------------------------------------------------------------
Bounds for Omega (proofs)
-------------------------------------------------------------------------
For all real x:   Om(x) = sin(x)/x = int_0^1 cos(s x) ds        (Om(0)=1).
Differentiating under the integral (smooth integrand, compact range):
    Om^(k)(x) = int_0^1 s^k cos^(k)(s x) ds,
and |cos^(k)| <= 1, hence for all x and k >= 0
    (D_k)   |Om^(k)(x)| <= int_0^1 s^k ds = 1/(k+1).
In particular |Om'''| <= 1/4.
Also  Om(x) - 1 = int_0^1 (cos(sx) - 1) ds  and  |cos y - 1| <= y^2/2, so
    (Z)     |Om(x) - 1| <= x^2 int_0^1 s^2/2 ds = x^2/6.
For x > 0, |sin x| <= 1 gives
    (T)     |Om(x)| <= 1/x.

Write Phi(t) = psi0*Om(t) - L - sum_d U_d Om(r_d t), where the u-weights have
been summed exactly per distinct distance r_d (U_d = sum of u_ij with
r_ij = r_d).  Then Phi^(k)(t) = psi0 Om^(k)(t) - sum_d U_d r_d^k Om^(k)(r_d t)
and by (D_3)
    |Phi'''(t)| <= M3 := ( |psi0| + sum_d |U_d| r_d^3 ) / 4   for all t.
Taylor with Lagrange remainder around t0, |s| <= h:
    Phi(t0+s) >= Phi(t0) + Phi'(t0) s + Phi''(t0) s^2/2 - M3 h^3/6.
The quadratic part is bounded below rigorously over s in [-h,h]
(function quad_lower).  Phi(t0), Phi'(t0), Phi''(t0) are evaluated in
ball arithmetic (python-flint arb, rigorous enclosures of sin/cos/sqrt) via
    Om   = sin x / x
    Om'  = cos x / x - sin x / x^2
    Om'' = -sin x / x - 2 cos x / x^2 + 2 sin x / x^3 .
Intervals that fail are bisected (adaptive).

Near t = 0 (interval [0, t1]): by (Z),
    Phi(t) - Phi(0) = psi0 (Om(t)-1) - sum_d U_d (Om(r_d t)-1)
                    >= -( |psi0| + sum_d |U_d| r_d^2 ) t^2 / 6,
and Phi(0) = psi0 - L - sum u = c + m.  So Phi >= m on [0,t1] as soon as
c - K2 t1^2/6 > 0 with K2 = |psi0| + sum |U_d| r_d^2 (exact rational: r_d^2
is rational).

Tail t >= T: by (T),  Phi(t) >= -L - ( |psi0| + sum_d |U_d|/r_d ) / t,
which is >= m whenever -L-m > 0 and t >= T := (|psi0| + sum|U_d|/r_d)/(-L-m).
We use a rational upper bound for 1/r_d = sqrt(1/d2) (integer sqrt), so T
is a rational upper bound for the threshold.
"""

import sys
import json
import time
from fractions import Fraction
from math import gcd, isqrt


class Reject(Exception):
    pass


def fail(msg):
    raise Reject(msg)


def parse_q(s, what):
    if isinstance(s, bool) or not isinstance(s, (str, int)):
        fail("%s: expected rational string, got %r" % (what, s))
    try:
        return Fraction(s)
    except (ValueError, ZeroDivisionError):
        fail("%s: bad rational %r" % (what, s))


def lcm(a, b):
    return a // gcd(a, b) * b


def sqrt_upper(q, bits=200):
    """Rational upper bound for sqrt(q), q >= 0 rational."""
    if q < 0:
        fail("internal: sqrt of negative")
    a, b = q.numerator, q.denominator
    # sqrt(a/b) = sqrt(a*b)/b ; sqrt(a*b*4^k) <= isqrt(a*b*4^k)+1
    k = bits
    s = isqrt(a * b * (1 << (2 * k))) + 1
    r = Fraction(s, b << k)
    if r * r < q:  # paranoia, cannot happen
        fail("internal: sqrt_upper wrong")
    return r


# ---------------------------------------------------------------------------
# Loading and format checks
# ---------------------------------------------------------------------------
def load(path):
    try:
        with open(path) as f:
            d = json.load(f)
    except Exception as e:  # noqa
        fail("cannot read JSON: %s" % e)
    if not isinstance(d, dict):
        fail("top level is not an object")
    for k in ("side", "scale", "points", "psi0", "psiT", "m", "lam", "u", "wc"):
        if k not in d:
            fail("missing key %r" % k)
    if d["side"] != "M1R3":
        fail("side is %r, expected 'M1R3'" % d["side"])
    if d["scale"] != "1/sqrt2":
        fail("scale is %r, expected '1/sqrt2'" % d["scale"])
    pts = d["points"]
    if not isinstance(pts, list) or len(pts) == 0:
        fail("points must be a nonempty list")
    P = []
    for a, p in enumerate(pts):
        if not isinstance(p, list) or len(p) != 3:
            fail("point %d is not a 3-vector" % (a + 1))
        P.append(tuple(parse_q(x, "point %d" % (a + 1)) for x in p))
    n = len(P)
    if n > 62:
        fail("too many points for this checker (%d)" % n)
    psi0 = parse_q(d["psi0"], "psi0")
    psiT = parse_q(d["psiT"], "psiT")
    m = parse_q(d["m"], "m")
    if not isinstance(d["lam"], list) or len(d["lam"]) != n:
        fail("lam must have exactly n=%d entries" % n)
    lam = [parse_q(x, "lam[%d]" % (a + 1)) for a, x in enumerate(d["lam"])]
    U = {}
    if not isinstance(d["u"], list):
        fail("u must be a list")
    for e in d["u"]:
        if not isinstance(e, list) or len(e) != 3:
            fail("bad u entry %r" % (e,))
        i, j, v = e
        if not (isinstance(i, int) and isinstance(j, int)) or isinstance(i, bool) or isinstance(j, bool):
            fail("u indices must be integers: %r" % (e,))
        if not (1 <= i < j <= n):
            fail("u entry needs 1 <= i < j <= n: %r" % (e,))
        if (i - 1, j - 1) in U:
            fail("duplicate u pair (%d,%d)" % (i, j))
        U[(i - 1, j - 1)] = parse_q(v, "u(%d,%d)" % (i, j))
    WC = []
    if not isinstance(d["wc"], list):
        fail("wc must be a list")
    for e in d["wc"]:
        if not isinstance(e, list) or len(e) != 3:
            fail("bad wc entry %r" % (e,))
        I, J, v = e
        if not isinstance(I, list) or not isinstance(J, list):
            fail("wc index lists malformed: %r" % (e,))
        if len(I) != len(J) or len(I) == 0:
            fail("wc lists must be nonempty and of equal length: %r" % (e,))
        for x in I + J:
            if not isinstance(x, int) or isinstance(x, bool) or not (1 <= x <= n):
                fail("wc index out of range: %r" % (e,))
        if len(set(I)) != len(I) or len(set(J)) != len(J):
            fail("wc lists must have distinct indices: %r" % (e,))
        WC.append(([x - 1 for x in I], [x - 1 for x in J], parse_q(v, "wc weight")))
    return P, psi0, psiT, m, lam, U, WC


# ---------------------------------------------------------------------------
# (G) geometry
# ---------------------------------------------------------------------------
def geometry(P, U, WC):
    n = len(P)
    D2 = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            s = sum((P[i][k] - P[j][k]) ** 2 for k in range(3)) / 2
            if s == 0:
                fail("points %d and %d coincide" % (i + 1, j + 1))
            D2[i][j] = D2[j][i] = s
    adj = [0] * n
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if D2[i][j] == 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
                edges.append((i, j))
    for (i, j) in U:
        if D2[i][j] == 1:
            fail("u placed on unit-distance pair (%d,%d)" % (i + 1, j + 1))
    for idx, (I, J, w) in enumerate(WC):
        k = len(I)
        for a in range(k):
            for b in range(a + 1, k):
                if D2[I[a]][I[b]] != D2[J[a]][J[b]]:
                    fail("wc entry %d not congruent: |p%d-p%d|^2=%s vs |p%d-p%d|^2=%s"
                         % (idx + 1, I[a] + 1, I[b] + 1, D2[I[a]][I[b]],
                            J[a] + 1, J[b] + 1, D2[J[a]][J[b]]))
    return D2, adj, edges


# ---------------------------------------------------------------------------
# (A) independent sets
# ---------------------------------------------------------------------------
def indep_dfs(n, adj):
    """List all independent sets as bitmasks (vertices added in increasing
    order; each set produced exactly once)."""
    out = []
    stack = [(0, 0, 0)]  # (set, forbidden, next vertex)
    while stack:
        S, forb, v = stack.pop()
        out.append(S)
        for w in range(v, n):
            if not (forb >> w) & 1:
                stack.append((S | (1 << w), forb | adj[w] | ((1 << (w + 1)) - 1), w + 1))
    return out


def indep_count_rec(n, adj):
    memo = {}
    sys.setrecursionlimit(10000)

    def f(mask):
        if mask == 0:
            return 1
        r = memo.get(mask)
        if r is not None:
            return r
        v = mask.bit_length() - 1  # highest vertex
        rest = mask & ~(1 << v)
        r = f(rest) + f(rest & ~adj[v])
        memo[mask] = r
        return r

    return f((1 << n) - 1)


def indep_count_brute(n, edges):
    import numpy as np
    N = 1 << n
    ok = np.ones(N, dtype=bool)
    ar = np.arange(N, dtype=np.int64)
    for (i, j) in edges:
        pm = (1 << i) | (1 << j)
        ok &= (ar & pm) != pm
    return int(ok.sum())


def check_Q(n, adj, edges, psiT, lam, U, WC):
    import numpy as np
    sets = indep_dfs(n, adj)
    K = len(sets)
    if len(set(sets)) != K:
        fail("internal: DFS produced duplicate sets")
    for S in sets:
        for i in range(n):
            if (S >> i) & 1 and (adj[i] & S):
                fail("internal: DFS produced a non-independent set")
    K2 = indep_count_rec(n, adj)
    if K2 != K:
        fail("independent-set count mismatch: DFS %d vs recursion %d" % (K, K2))
    K3 = None
    if n <= 24:
        K3 = indep_count_brute(n, edges)
        if K3 != K:
            fail("independent-set count mismatch: DFS %d vs brute force %d" % (K, K3))

    # clear denominators
    allq = [psiT] + lam + list(U.values()) + [w for _, _, w in WC]
    D = 1
    for q in allq:
        D = lcm(D, q.denominator)
    iz = lambda q: (q * D).numerator  # exact since D*q is an integer
    for q in allq:
        if (q * D).denominator != 1:
            fail("internal: denominator clearing failed")
    total = sum(abs(iz(q)) for q in allq) + abs(iz(psiT)) * 2 + sum(abs(iz(w)) for _, _, w in WC)
    dt = np.int64 if total < (1 << 62) else object

    M = np.array(sets, dtype=np.int64)
    Q = np.full(K, iz(psiT), dtype=dt)
    for i in range(n):
        if lam[i] != 0:
            Q += np.where((M >> i) & 1 == 1, iz(lam[i]), 0).astype(dt)
    for (i, j), v in U.items():
        if v != 0:
            pm = (1 << i) | (1 << j)
            Q += np.where((M & pm) == pm, iz(v), 0).astype(dt)
    for I, J, w in WC:
        if w == 0:
            continue
        mi = sum(1 << x for x in I)
        mj = sum(1 << x for x in J)
        Q += (np.where((M & mi) == mi, iz(w), 0) - np.where((M & mj) == mj, iz(w), 0)).astype(dt)
    k = int(np.argmin(Q)) if dt is np.int64 else min(range(K), key=lambda a: Q[a])
    qmin = Fraction(int(Q[k]), D)
    if qmin < 0:
        Sl = [i + 1 for i in range(n) if (sets[k] >> i) & 1]
        fail("Q(S) < 0 for independent set S=%s : Q=%s" % (Sl, qmin))
    # independent spot check of a few sets with pure Fraction arithmetic
    for S in [sets[0], sets[k], sets[-1], sets[K // 2]]:
        q = psiT + sum(lam[i] for i in range(n) if (S >> i) & 1)
        q += sum(v for (i, j), v in U.items() if (S >> i) & 1 and (S >> j) & 1)
        for I, J, w in WC:
            q += w * ((1 if all((S >> x) & 1 for x in I) else 0) - (1 if all((S >> x) & 1 for x in J) else 0))
        idx = sets.index(S)
        if Fraction(int(Q[idx]), D) != q:
            fail("internal: vectorised Q disagrees with exact Fraction evaluation")
    return K, K2, K3, qmin, len(edges)


# ---------------------------------------------------------------------------
# (P) Phi(t) >= m on (0, inf)
# ---------------------------------------------------------------------------
def check_Phi(psi0, L, m, c, D2, U):
    import flint
    from flint import arb, fmpq
    flint.ctx.prec = 160

    # aggregate u by squared distance (exact)
    agg = {}
    for (i, j), v in U.items():
        agg[D2[i][j]] = agg.get(D2[i][j], Fraction(0)) + v
    agg = {d: v for d, v in agg.items() if v != 0}
    dist = sorted(agg)
    if any(d <= 0 for d in dist):
        fail("internal: nonpositive distance")

    # exact rational constants
    r_up = {d: sqrt_upper(d) for d in dist}            # >= r_d
    rinv_up = {d: sqrt_upper(1 / d) for d in dist}     # >= 1/r_d
    K2 = abs(psi0) + sum(abs(agg[d]) * d for d in dist)            # exact
    M3 = (abs(psi0) + sum(abs(agg[d]) * d * r_up[d] for d in dist)) / 4
    gap = -L - m
    if gap <= 0:
        fail("tail argument impossible: -sum(lam) - m = %s <= 0 (Phi(t)->-sum lam as t->inf)" % gap)
    T = (abs(psi0) + sum(abs(agg[d]) * rinv_up[d] for d in dist)) / gap
    # [0, t1]: need c - K2 t1^2/6 > 0
    t1 = Fraction(1, 1)
    while c - K2 * t1 * t1 / 6 <= 0:
        t1 /= 2
        if t1 < Fraction(1, 1 << 60):
            fail("cannot certify Phi >= m near t = 0")
    # round T up to a dyadic for nicer subdivision
    Tq = Fraction(-((-T.numerator * 1024) // T.denominator), 1024)
    if Tq < T or Tq <= t1:
        Tq = max(Tq, t1 + 1)

    def A(q):
        return arb(fmpq(q.numerator, q.denominator))

    psi0a = A(psi0)
    La = A(L)
    ma = A(m)
    M3a = A(M3)
    terms = [(A(agg[d]), A(d).sqrt()) for d in dist]

    def phi012(t0):
        t = A(t0)
        s, co = t.sin(), t.cos()
        it = 1 / t
        o0 = s * it
        o1 = (co - o0) * it
        o2 = -o0 - 2 * o1 * it
        f0 = psi0a * o0 - La
        f1 = psi0a * o1
        f2 = psi0a * o2
        for w, r in terms:
            x = r * t
            s, co = x.sin(), x.cos()
            ix = 1 / x
            o0 = s * ix
            o1 = (co - o0) * ix          # = cos/x - sin/x^2
            o2 = -o0 - 2 * o1 * ix       # = -sin/x - 2cos/x^2 + 2 sin/x^3
            f0 -= w * o0
            f1 -= w * r * o1
            f2 -= w * r * r * o2
        return f0, f1, f2

    def quad_lower(f0, f1, f2, h):
        """Rigorous lower bound (arb, use .lower()) for
        min_{|s|<=h} f0 + f1 s + f2 s^2/2, for all values in the balls."""
        ha = A(h)
        e1 = f0 + f1 * ha + f2 * ha * ha / 2
        e2 = f0 - f1 * ha + f2 * ha * ha / 2
        lo_end = e1.lower() if e1.lower() < e2.lower() else e2.lower()
        if f2 > 0:
            # convex: global minimum f0 - f1^2/(2 f2) is a valid lower bound;
            # if the vertex is certainly outside [-h,h] the min is at an end.
            vert = -f1 / f2
            if (vert > ha) or (vert < -ha):
                return lo_end
            v = f0 - f1 * f1 / (2 * f2)
            return v.lower()
        if f2 <= 0:
            # concave (or linear): min over the interval at an endpoint
            return lo_end
        # f2 ball straddles 0: f0 + f1 s + f2 s^2/2 >= f0 - |f1| h + min(0,f2) h^2/2
        b = f0 - abs(f1) * ha + (f2.lower() if f2.lower() < 0 else arb(0)) * ha * ha / 2
        return b.lower()

    stack = [(Tq - Fraction(k, 1) * (Tq - t1) / 64 - (Tq - t1) / 64,
              Tq - Fraction(k, 1) * (Tq - t1) / 64) for k in range(64)]
    nint = 0
    worst = None
    min_h = Fraction(1, 1 << 40)
    while stack:
        a, b = stack.pop()
        t0 = (a + b) / 2
        h = (b - a) / 2
        f0, f1, f2 = phi012(t0)
        lb = quad_lower(f0, f1, f2, h) - M3a * A(h) ** 3 / 6
        if (f0 - ma) < 0:
            fail("Phi(t0) < m at t0 = %s ~ %.12f (Phi(t0)-m in %s)"
                 % (t0, float(t0), (f0 - ma).str(10)))
        marg = lb - ma
        if marg > 0:
            nint += 1
            if worst is None or marg.lower() < worst[0]:
                worst = (marg.lower(), t0)
            continue
        if h < min_h:
            fail("Phi(t) >= m fails (or cannot be certified) near t = %.12f "
                 "(Phi(t)-m in %s)" % (float(t0), (f0 - ma).str(10)))
        stack.append((a, t0))
        stack.append((t0, b))
    return {
        "t1": t1, "T": T, "Tq": Tq, "intervals": nint,
        "min_margin": worst[0], "argmin": worst[1],
        "ndist": len(dist), "M3": M3,
    }


# ---------------------------------------------------------------------------
def upper_root(c, m, psiT, digits=8):
    """Smallest k/10^digits >= B, where B = positive root of c x^2 + m x - psiT
    (c > 0, psiT >= 0).  For x >= 0, f(x) >= 0  <=>  x >= B."""
    f = lambda x: c * x * x + m * x - psiT
    S = 10 ** digits
    # float guess, then exact correction
    disc = float(m) ** 2 + 4 * float(c) * float(psiT)
    g = (-float(m) + disc ** 0.5) / (2 * float(c))
    k = max(0, int(g * S) - 5)
    while f(Fraction(k, S)) >= 0 and k > 0:
        k -= 1
    while f(Fraction(k, S)) < 0:
        k += 1
    if f(Fraction(k, S)) < 0 or (k > 0 and f(Fraction(k - 1, S)) >= 0 and psiT > 0):
        fail("internal: root bracketing failed")
    return Fraction(k, S), k


def main(argv):
    if len(argv) != 2:
        print("usage: check_m1r3_indep.py CERT.json")
        return 1
    t_start = time.time()
    try:
        P, psi0, psiT, m, lam, U, WC = load(argv[1])
        n = len(P)
        print("loaded: n=%d points, %d u-pairs, %d wc entries" % (n, len(U), len(WC)))

        D2, adj, edges = geometry(P, U, WC)
        print("(G) ok: %d unit edges, points distinct, u off unit pairs, all %d wc congruent"
              % (len(edges), len(WC)))

        L = sum(lam, Fraction(0))
        su = sum(U.values(), Fraction(0))
        c = psi0 - L - su - m
        if psiT < 0:
            fail("psiT < 0")
        if c <= 0:
            fail("c = psi0 - sum lam - sum u - m = %s <= 0" % c)
        print("(B) ok: c = %s ~ %.12g > 0, psiT = %s ~ %.6g >= 0" % (c, float(c), psiT, float(psiT)))

        t = time.time()
        K, K2, K3, qmin, ne = check_Q(n, adj, edges, psiT, lam, U, WC)
        print("(A) ok: %d independent sets (recursion %d, brute force %s); min Q(S) = %s ~ %.6g [%.1fs]"
              % (K, K2, K3, qmin, float(qmin), time.time() - t))

        t = time.time()
        info = check_Phi(psi0, L, m, c, D2, U)
        print("(P) ok: Phi >= m on [0,%s] by t^2 bound, on [t1,%.4f] by %d Taylor intervals "
              "(min certified margin %.3e near t=%.4f), tail t >= %.4f by |Om(x)|<=1/x; "
              "%d distinct distances [%.1fs]"
              % (info["t1"], float(info["Tq"]), info["intervals"],
                 float(info["min_margin"].mid()), float(info["argmin"]),
                 float(info["T"]), info["ndist"], time.time() - t))

        Bq, k = upper_root(c, m, psiT)
        print("B <= %d/10^8 = %s (exact: c*x^2+m*x-psiT >= 0 at x=%s)"
              % (k, "%.8f" % (k / 1e8), "%.8f" % (k / 1e8)))
        print("ACCEPTED: upper density <= B <= %.8f  (total %.1fs)" % (k / 1e8, time.time() - t_start))
        return 0
    except Reject as e:
        print("REJECTED: %s" % e)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
