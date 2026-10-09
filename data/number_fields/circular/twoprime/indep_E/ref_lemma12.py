"""Referee E: exact test of the conclusion of Lemma 12 for small k.

S(k) = S^theta(k,1) = {C : Re(conj(C) gamma) in [theta,1-theta]+Z, gamma in G(k,1)}, period D_k = 13*5^k.
Level 1: direct clipping over all 65^2 cells (as ref_clip.py, re-done here with the generator-built rotations).
Level k >= 2: exact lifting. S(k) is contained in S(k-1) (G(k-1,1) in G(k,1)), whose period is D_k/5, so S(k)
modulo D_k is contained in the 25 translates by (D_k/5){0..4}^2 of the components of S(k-1); each translate is
clipped by the new functionals (gamma = rho^{+-k} sigma^l, |l| <= 1, Re and Im), branching over strip indices.
Optionally (mode 'direct2') level 2 is also computed by direct clipping over all 325^2 cells, as a cross-check.

For every component, every vertex C is tested exactly against Lemma 12's conclusion:
  C in N h + P_k + N Z[i],  or  C in N eps + Y_k(eps) + N Z[i] for some eps in E_q,   N = 5^k,
with P_k = {x : conj(x) rho^j in B_s, |j|<=k},  Y_k(eps) = {y : eta_j(eps)/6 + conj(y) rho^j in B_s, |j|<=k}.

Usage: python3 ref_lemma12.py NUM DEN KMAX [direct2]
"""
import sys
import time
from fractions import Fraction as Fr
from math import floor, ceil

from ref_rot import rotations_from_generators, gmul


def funcs_level(k):
    """functionals (p,q) over D = 13*5^k for G(k,1): one per +-pair, as dict key (j,l,re/im)."""
    D, rots = rotations_from_generators(k, 1)
    out = {}
    for (a, j, l), z in rots.items():
        if a == 0:
            out[(j, l, 0)] = z                      # Re(conj(c) gamma)
            out[(j, l, 1)] = gmul((0, -1), z)       # Re(conj(c) (-i gamma)) = Im(conj(c) gamma)
    return D, out


def fval(f, P, D):
    return (f[0] * P[0] + f[1] * P[1]) / D


def clip(poly, f, C, keep_le, D):
    if not poly:
        return []
    out = []
    m = len(poly)
    vals = []
    for P in poly:
        v = fval(f, P, D) - C
        vals.append(v if keep_le else -v)
    for k in range(m):
        P, Q = poly[k], poly[(k + 1) % m]
        vP, vQ = vals[k], vals[(k + 1) % m]
        if vP <= 0:
            out.append(P)
        if (vP < 0 < vQ) or (vQ < 0 < vP):
            t = vP / (vP - vQ)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    res = []
    for P in out:
        if not res or res[-1] != P:
            res.append(P)
    while len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


def refine(polys, flist, theta, D):
    """Clip each polygon by all functionals in flist (branching over strip indices); return leaves."""
    leaves = []
    stack = [(p, 0) for p in polys]
    while stack:
        poly, depth = stack.pop()
        if depth == len(flist):
            leaves.append(poly)
            continue
        f = flist[depth]
        vs = [fval(f, P, D) for P in poly]
        lo, hi = min(vs), max(vs)
        for n in range(ceil(lo - 1 + theta), floor(hi - theta) + 1):
            q2 = clip(clip(poly, f, n + theta, False, D), f, n + 1 - theta, True, D)
            if q2:
                stack.append((q2, depth + 1))
    return leaves


def level1_direct(theta):
    D, F = funcs_level(1)
    fx, fy = F[(0, 0, 0)], F[(0, 0, 1)]
    assert fx == (65, 0) and fy == (0, -65)  # -y in strips iff y in strips
    others = [F[key] for key in sorted(F) if key[:2] != (0, 0)]
    polys = []
    for a in range(D):
        for b in range(D):
            x0, x1, y0, y1 = a + theta, a + 1 - theta, b + theta, b + 1 - theta
            polys.append([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    return refine(polys, others, theta, D)


def level2_direct(theta):
    D, F = funcs_level(2)
    fx, fy = F[(0, 0, 0)], F[(0, 0, 1)]
    assert fx == (325, 0) and fy == (0, -325)
    # put the rho^{+-1} and sigma functionals first (cheap pruning), all others after
    others = [F[key] for key in sorted(F, key=lambda t: (abs(t[0]) + abs(t[1]), t)) if key[:2] != (0, 0)]
    leaves = []
    for a in range(D):
        polys = []
        for b in range(D):
            x0, x1, y0, y1 = a + theta, a + 1 - theta, b + theta, b + 1 - theta
            polys.append([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
        leaves.extend(refine(polys, others, theta, D))
    return leaves


def lift(comps, k, theta):
    """components of S(k-1) (mod D_{k-1}) -> components of S(k) (mod D_k)."""
    D, F = funcs_level(k)
    step = D // 5
    newf = [F[key] for key in sorted(F) if abs(key[0]) == k]
    assert len(newf) == 12
    polys = []
    for poly in comps:
        for a in range(5):
            for b in range(5):
                polys.append([(P[0] + a * step, P[1] + b * step) for P in poly])
    return refine(polys, newf, theta, D)


def extreme(poly):
    """remove collinear (non-extreme) vertices of a convex polygon given in cyclic order."""
    pts = list(poly)
    changed = True
    while changed and len(pts) > 2:
        changed = False
        for t in range(len(pts)):
            A, B, C = pts[t - 1], pts[t], pts[(t + 1) % len(pts)]
            cr = (B[0] - A[0]) * (C[1] - B[1]) - (B[1] - A[1]) * (C[0] - B[0])
            if cr == 0:
                pts.pop(t)
                changed = True
                break
    return pts


def in_box(z, s):
    return -s <= z[0] <= s and -s <= z[1] <= s


def cmul_conj_rho(x, j, k):
    """conj(x) * rho^j for complex x as (re, im) Fractions, rho = (3+4i)/5."""
    re, im = x[0], -x[1]
    r = (Fr(3, 5), Fr(4, 5)) if j >= 0 else (Fr(3, 5), Fr(-4, 5))
    for _ in range(abs(j)):
        re, im = re * r[0] - im * r[1], re * r[1] + im * r[0]
    return (re, im)


def frac_class(z):
    """z mod Z[i] as (re - floor, im - floor)."""
    return (z[0] - floor(z[0]), z[1] - floor(z[1]))


def reduce_mod(z, N):
    """representative of z modulo N Z[i] with coordinates in [-N/2, N/2)."""
    return tuple(t - N * floor((t + Fr(N, 2)) / N) for t in z)


def classify(poly, k, theta):
    """Return ('c',) or ('q', alpha, beta) if every vertex satisfies Lemma 12's conclusion with that type
    (and the same class mod N), else None."""
    N = 5 ** k
    s = Fr(1, 2) - theta
    h = (Fr(1, 2), Fr(1, 2))
    cands = [("c", (N * h[0], N * h[1]), None)]
    for al in (1, 2):
        for be in (1, 2):
            eps = (Fr(al, 3), Fr(be, 3))
            # eta_j(eps): conj(N eps) rho^j = h + eta_j/6 mod Z[i]
            etas = {}
            for j in range(-k, k + 1):
                w = frac_class(cmul_conj_rho((N * eps[0], N * eps[1]), j, k))
                et = (6 * (w[0] - h[0]), 6 * (w[1] - h[1]))
                assert et[0] in (1, -1) and et[1] in (1, -1), et
                etas[j] = et
            cands.append((("q", al, be), (N * eps[0], N * eps[1]), etas))
    for name, center, etas in cands:
        ok = True
        shifts = set()
        for P in poly:
            x = reduce_mod((P[0] - center[0], P[1] - center[1]), N)
            shifts.add(((P[0] - center[0] - x[0]) / N, (P[1] - center[1] - x[1]) / N))
            for j in range(-k, k + 1):
                z = cmul_conj_rho(x, j, k)
                if etas is not None:
                    z = (z[0] + Fr(etas[j][0], 6), z[1] + Fr(etas[j][1], 6))
                if not in_box(z, s):
                    ok = False
                    break
            if not ok:
                break
        if ok and len(shifts) == 1:
            return name, shifts.pop()
    return None


def report(level, comps, theta, t0):
    D = 13 * 5 ** level
    N = 5 ** level
    lines = []
    bad = 0
    kinds = {}
    pos13 = {}
    for poly in comps:
        cl = classify(poly, level, theta)
        if cl is None:
            bad += 1
            lines.append("   NOT of the form of Lemma 12: vertices %s" % [(float(P[0]), float(P[1])) for P in poly][:4])
            continue
        name, m = cl
        kinds[str(name)] = kinds.get(str(name), 0) + 1
        # position of the component's centre modulo 13N relative to 13N*E (is it at 13N*eps?)
        pos13.setdefault(str(name), set()).add((m[0] % 13, m[1] % 13))
    lines.insert(0, "level k=%d (N=%d, period %d): %d components mod %dZ[i]; failing Lemma 12's conclusion: %d ; types %s ; "
                 "N-multiples m mod 13 by type: %s ; %.1fs"
                 % (level, N, D, len(comps), D, bad, kinds, {k2: sorted(v) for k2, v in pos13.items()}, time.time() - t0))
    return "\n".join(lines)


def main():
    num, den, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    theta = Fr(num, den)
    mode = sys.argv[4] if len(sys.argv) > 4 else ""
    t0 = time.time()
    print("theta = %s = %.9f" % (theta, float(theta)))
    comps = level1_direct(theta)
    print(report(1, comps, theta, t0), flush=True)
    if mode == "direct2":
        t1 = time.time()
        d2 = level2_direct(theta)
        print("direct level 2:", report(2, d2, theta, t1), flush=True)
    for k in range(2, kmax + 1):
        t1 = time.time()
        comps = lift(comps, k, theta)
        print(report(k, comps, theta, t1), flush=True)
        if mode == "direct2" and k == 2:
            a = sorted(sorted(extreme(p)) for p in comps)
            b = sorted(sorted(extreme(p)) for p in d2)
            print("   level 2 lifting == direct clipping (same vertex sets):", a == b, flush=True)


if __name__ == "__main__":
    main()
