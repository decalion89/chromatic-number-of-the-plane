# Referee's own exact computation of S_N^theta, N = 5^k, level by level, as convex polygons mod N.
# c in S_N^theta  <=>  for |j| <= k, Re and Im of conj(c) rho^j lie in 1/2 + [-s, s] + Z  (s = 1/2 - theta).
# Level 1 starts from the 25 squares h + m + B_s (condition at rho^0), cut at rho^{+-1};
# level k: the 25 translates by (N/5){0..4}^2 of each level-(k-1) polygon, cut at rho^{+-k}.
# Every polygon is then classified: inside N h + P_k (type c), inside N eps + Y_k(eps) (type q), or extra.
import sys, time, math
from fractions import Fraction as Fr
from gauss import G, I, RHO, H, modZ

HALF = Fr(1, 2)

def forms(j):
    """Linear forms (Re, Im) of conj(c) rho^j for c=(a,b): returns ((p,q),(q,-p))."""
    r = RHO ** j
    return ((r.re, r.im), (r.im, -r.re))

def clip(poly, a, b, c):
    """Keep a*x + b*y <= c on a convex polygon (list of vertices, possibly degenerate)."""
    n = len(poly)
    if n == 0:
        return []
    vals = [a * p[0] + b * p[1] - c for p in poly]
    if all(v <= 0 for v in vals):
        return poly
    if all(v > 0 for v in vals):
        return []
    if n == 1:
        return [] if vals[0] > 0 else poly
    out = []
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        fp, fq = vals[i], vals[(i + 1) % n]
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    res = []
    for p in out:
        if p not in res:
            res.append(p)
    return res

def cut(poly, j, s):
    """Intersect with {Re, Im of conj(c) rho^j in 1/2 + [-s,s] + Z}. Returns list of polygons."""
    pieces = [poly]
    for (a, b) in forms(j):
        new = []
        for P in pieces:
            vals = [a * p[0] + b * p[1] for p in P]
            lo, hi = min(vals), max(vals)
            nlo = math.ceil(lo - HALF - s); nhi = math.floor(hi - HALF + s)
            for n in range(nlo, nhi + 1):
                Q = clip(P, a, b, n + HALF + s)
                Q = clip(Q, -a, -b, -(n + HALF - s))
                if Q:
                    new.append(Q)
        pieces = new
    return pieces

def level1(s):
    polys = []
    for m1 in range(5):
        for m2 in range(5):
            cx, cy = HALF + m1, HALF + m2
            sq = [(cx - s, cy - s), (cx + s, cy - s), (cx + s, cy + s), (cx - s, cy + s)]
            pieces = [sq]
            for j in (0, 1, -1):
                pieces = [Q for P in pieces for Q in cut(P, j, s)]
            polys.extend(pieces)
    return polys

def lift(polys, k, s):
    N = 5 ** k; step = N // 5
    out = []
    for P in polys:
        for a in range(5):
            for b in range(5):
                T = [(p[0] + a * step, p[1] + b * step) for p in P]
                pieces = [T]
                for j in (k, -k):
                    pieces = [Q for R in pieces for Q in cut(R, j, s)]
                out.extend(pieces)
    return out

EQ = [G(Fr(a, 3), Fr(b, 3)) for a in (1, 2) for b in (1, 2)]

def classify(P, k, s):
    """Return 'C', 'Q' or 'X' (extra)."""
    N = 5 ** k
    rj = {j: RHO ** j for j in range(-k, k + 1)}
    def reduce(v, base):
        x = G(v[0], v[1]) - base
        return G(x.re - N * round(x.re / N), x.im - N * round(x.im / N))
    def inPk(x):
        xc = x.conj()
        for j in range(-k, k + 1):
            w = xc * rj[j]
            if abs(w.re) > s or abs(w.im) > s:
                return False
        return True
    if all(inPk(reduce(v, H * N)) for v in P):
        return 'C'
    for e in EQ:
        base = e * N
        et = {j: modZ(base.conj() * rj[j] - H) * 6 for j in range(-k, k + 1)}
        def inY(y):
            yc = y.conj()
            for j in range(-k, k + 1):
                w = yc * rj[j]
                if abs(et[j].re / 6 + w.re) > s or abs(et[j].im / 6 + w.im) > s:
                    return False
            return True
        if all(inY(reduce(v, base)) for v in P):
            return 'Q'
    return 'X'

def run(theta, kmax, verbose=True, stop_when_clean=False):
    s = HALF - theta
    out = []
    polys = level1(s)
    for k in range(1, kmax + 1):
        t0 = time.time()
        if k > 1:
            polys = lift(polys, k, s)
        cls = [classify(P, k, s) for P in polys]
        cnt = {t: cls.count(t) for t in 'CQX'}
        degen = sum(1 for P in polys if len(P) <= 2)
        out.append((k, len(polys), cnt, degen))
        if verbose:
            print(f"theta={theta} level {k}: {len(polys)} polygons  C={cnt['C']} Q={cnt['Q']} extra={cnt['X']}  (degenerate: {degen})  [{time.time()-t0:.1f}s]", flush=True)
    return out, polys

if __name__ == "__main__":
    theta = Fr(sys.argv[1])
    kmax = int(sys.argv[2])
    run(theta, kmax)
