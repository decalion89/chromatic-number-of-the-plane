"""Certificates for the two-prime window W(K, M) = {rho^j sigma^l : |j| <= K, |l| <= M} (N = 5^K 13^M) at r0.

Writes  cert_K{K}M{M}.txt :  the functionals; every index vector of S^{r0}(K, M) (from probe2d.py, lifts 5 ... 5
then 13 ... 13) with its label:
    MAIN C / Q{a}{b} m1 m2    : the type point (N h or (N/3)(a+bi)) translated by N m has these indices;
    SEVEN a b                 : the 7-torsion point N(a+bi)/7 has these indices and all margins >= 2/7;
    EXTRA kappa ; cert ...    : exact LP dual certificate: lam >= 0, sum lam = 1, linear parts cancel,
                                sum lam*const = kappa (so no point of the component has all margins > kappa);
and  tree_K{K}M{M}.txt : a Farkas-certified branch tree over all N^2 cells (format of gen_tree.py) whose leaves are
exactly these index vectors.  Checked by verify_window.py (independent code).
usage: python3 certify_window.py K M r0"""
from fractions import Fraction as Fr
import sys, time
from probe2d import Probe, centroid, gamma, floor
from lpexact import lp_max

K, M, R0 = int(sys.argv[1]), int(sys.argv[2]), Fr(sys.argv[3])
N = 5 ** K * 13 ** M
LO, HI = R0, 1 - R0
PAIRS = [(j, l) for j in range(-K, K + 1) for l in range(-M, M + 1)]
FUN = []
for (j, l) in PAIRS:
    A, B = gamma(j, l)
    FUN.append(((j, l, 'Re'), (A, B)))
    FUN.append(((j, l, 'Im'), (B, -A)))
T_RE = [i for i, (key, f) in enumerate(FUN) if key == (0, 0, 'Re')][0]
T_IM = [i for i, (key, f) in enumerate(FUN) if key == (0, 0, 'Im')][0]
ORDER = [t for t in range(len(FUN)) if t not in (T_RE, T_IM)]


def ev(f, p):
    return f[0] * p[0] + f[1] * p[1]


def ivec(p):
    return tuple(floor(ev(f, p) - LO) for (_, f) in FUN)


def fmt(q):
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def components():
    P = Probe(R0)
    for _ in range(K):
        P.lift5()
    for _ in range(M):
        P.lift13()
    assert P.N == N
    return P.comps


def label(pt, iv):
    types = [('C', (Fr(N, 2), Fr(N, 2)))] + [(f'Q{a}{b}', (Fr(N * a, 3), Fr(N * b, 3))) for a in (1, 2) for b in (1, 2)]
    for name, T in types:
        m0 = ((floor(pt[0]) - floor(T[0])) // N, (floor(pt[1]) - floor(T[1])) // N)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                mm = (m0[0] + dx, m0[1] + dy)
                if ivec((T[0] + N * mm[0], T[1] + N * mm[1])) == iv:
                    return f"MAIN {name} {mm[0]} {mm[1]}"
    for a in range(7):
        for b in range(7):
            if (a, b) == (0, 0):
                continue
            base = (Fr(N * a, 7), Fr(N * b, 7))
            m0 = ((floor(pt[0]) - floor(base[0])) // N, (floor(pt[1]) - floor(base[1])) // N)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    q = (base[0] + N * (m0[0] + dx), base[1] + N * (m0[1] + dy))
                    if ivec(q) == iv:
                        return f"SEVEN {fmt(q[0])} {fmt(q[1])}"
    return None


def main():
    t0 = time.time()
    comps = components()
    lines = [f"# window K={K} M={M}, N={N}, r0={R0}; {len(comps)} components",
             "# functionals: " + " ".join(f"({j},{l},{e})" for ((j, l, e), _) in FUN)]
    targets = []
    count = {}
    for C in comps:
        pt = centroid(C)
        iv = ivec(pt)
        targets.append(iv)
        lab = label(pt, iv)
        head = "component indices " + " ".join(str(n) for n in iv) + " | "
        if lab is None:
            cons, meta = [], []
            for (key, f), n in zip(FUN, iv):
                cons.append((f[0], f[1], Fr(-n))); meta.append((key, '+'))
                cons.append((-f[0], -f[1], Fr(n + 1))); meta.append((key, '-'))
            t, x, y, lam = lp_max(cons)
            cert = " ; ".join(f"{meta[i][0][0]} {meta[i][0][1]} {meta[i][0][2]} {meta[i][1]} {fmt(lam[i])}" for i in sorted(lam))
            lab = f"EXTRA kappa = {fmt(t)} ; cert {cert}"
        lines.append(head + lab)
        kind = lab.split()[0] + ("" if not lab.startswith("EXTRA") else " " + lab.split()[3])
        count[kind] = count.get(kind, 0) + 1
    assert len(set(targets)) == len(targets)
    with open(f"cert_K{K}M{M}.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"components: {len(comps)}; labels {count} ({time.time()-t0:.0f}s)", flush=True)
    write_tree(set(targets))


def clip(poly, a, b, c):
    out = []
    n = len(poly)
    vals = [c - (a * p[0] + b * p[1]) for p in poly]
    if n == 1:
        return list(poly) if vals[0] >= 0 else []
    for i in range(n):
        P, Q = poly[i], poly[(i + 1) % n]
        fp, fq = vals[i], vals[(i + 1) % n]
        if fp >= 0:
            out.append(P)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    while len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


def farkas(poly, cons, new):
    (an, bn, cn), refn = new
    v = min(poly, key=lambda p: an * p[0] + bn * p[1])
    act = [(h, r) for (h, r) in cons if h[0] * v[0] + h[1] * v[1] == h[2]]
    for i in range(len(act)):
        (a1, b1, c1), r1 = act[i]
        for j in range(i, len(act)):
            (a2, b2, c2), r2 = act[j]
            det = a1 * b2 - a2 * b1
            if det == 0:
                if i == j and b1 * an - a1 * bn == 0:
                    lam = (-an / a1) if a1 != 0 else (-bn / b1)
                    if lam >= 0 and lam * c1 + cn < 0:
                        return [(r1, lam), (refn, Fr(1))]
                continue
            l1 = (-an * b2 + bn * a2) / det
            l2 = (-a1 * bn + b1 * an) / det
            if l1 >= 0 and l2 >= 0 and l1 * c1 + l2 * c2 + cn < 0:
                return [(r1, l1), (r2, l2), (refn, Fr(1))]
    raise RuntimeError("no Farkas certificate")


def write_tree(targets):
    t0 = time.time()
    stats = {'S': 0, 'D': 0, 'L': 0}
    leaves = set()
    out = open(f"tree_K{K}M{M}.txt", "w")
    out.write(f"# window K={K} M={M} N={N} r0={R0}; functional order: {' '.join(str(t) for t in ORDER)}\n")

    def rec(poly, cons, depth, idx, corners, lines):
        if depth == len(ORDER):
            iv = tuple(idx[t] for t in range(len(FUN)))
            assert iv in targets, iv
            lines.append("L")
            stats['L'] += 1
            leaves.add(iv)
            return
        t = ORDER[depth]
        a, b = FUN[t][1]
        vals = [a * p[0] + b * p[1] for p in corners]
        mn, mx = min(vals), max(vals)
        for n in [n for n in range(floor(mn - HI) - 1, floor(mx - LO) + 2) if n + HI >= mn and n + LO <= mx]:
            lo_c = ((-a, -b, -(n + LO)), f"{t} lo")
            hi_c = ((a, b, n + HI), f"{t} hi")
            p1 = clip(poly, *lo_c[0])
            if not p1:
                cert = farkas(poly, cons, lo_c)
                lines.append(f"D {n} " + " ; ".join(f"{ref} {fmt(l)}" for (ref, l) in cert)); stats['D'] += 1
                continue
            p2 = clip(p1, *hi_c[0])
            if not p2:
                cert = farkas(p1, cons + [lo_c], hi_c)
                lines.append(f"D {n} " + " ; ".join(f"{ref} {fmt(l)}" for (ref, l) in cert)); stats['D'] += 1
                continue
            lines.append(f"S {n}"); stats['S'] += 1
            idx2 = dict(idx); idx2[t] = n
            rec(p2, cons + [lo_c, hi_c], depth + 1, idx2, corners, lines)

    for a in range(N):
        for b in range(N):
            cons = [((Fr(1), Fr(0), a + HI), "x hi"), ((Fr(-1), Fr(0), -(a + LO)), "x lo"),
                    ((Fr(0), Fr(1), b + HI), "y hi"), ((Fr(0), Fr(-1), -(b + LO)), "y lo")]
            corners = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
            lines = [f"C {a} {b}"]
            rec(list(corners), cons, 0, {T_RE: a, T_IM: -b - 1}, corners, lines)
            out.write("\n".join(lines) + "\n")
    out.close()
    assert leaves == targets, (len(leaves), len(targets))
    print(f"tree: {stats}; leaves = the {len(targets)} index vectors ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
