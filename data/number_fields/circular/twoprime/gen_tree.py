"""Farkas-certified enumeration of S^r(1,1) (primes 5, 13; N = 65) at r = 7/25.

For every cell (a, b) in {0..64}^2 (the square a + r <= x <= a + 1 - r, b + r <= y <= b + 1 - r; these are the strips
of the two functionals x = Re(conj c) and -y = Im(conj c)), branch on the other 16 functionals f_t (fixed order) of
the nine rotations rho^j sigma^l (j, l in {-1,0,1}).  For functional t the candidate strip indices are all integers n
with [n + r, n + 1 - r] meeting the range of f_t on the CELL (the verifier recomputes it from the 4 corners).  Each
candidate is either
   D n ; ref lam ; ref lam ; ...   dead: exact Farkas certificate, constraints a x + b y <= c referenced as
        'x hi', 'x lo', 'y hi', 'y lo' (cell sides) or 't hi', 't lo' (strip sides of a fixed functional t, including
        the new one), with lam >= 0, sum lam (a, b) = (0, 0) and sum lam c < 0;
   S n   followed by the subtree for the next functional.
After the last functional a leaf 'L' is written; its index vector must be one of the 13 index vectors listed in
prop7_certificates.txt.  Output: tree_r7_25.txt; checked by verify_tree.py (independent code)."""
from fractions import Fraction as Fr
import time
from probe2d import gamma

R = Fr(7, 25)
LO, HI = R, 1 - R
FUN = []
for j in (-1, 0, 1):
    for l in (-1, 0, 1):
        A, B = gamma(j, l)
        FUN.append((A, B))       # Re(conj(c) g) = A x + B y
        FUN.append((B, -A))      # Im(conj(c) g) = B x - A y
T_RE, T_IM = 8, 9                # (j, l) = (0, 0)
assert FUN[T_RE] == (1, 0) and FUN[T_IM] == (0, -1)
ORDER = [t for t in range(18) if t not in (T_RE, T_IM)]


def floor(q):
    return q.numerator // q.denominator


def clip(poly, a, b, c):
    """part of the convex polygon with a x + b y <= c"""
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
    """poly: non-empty polygon = {c : all cons}; new = ((a, b, c), ref) with poly ∩ {a x + b y <= c} empty.
    Returns [(ref, lam), ...] with sum lam*(a,b) = 0 and sum lam*c < 0 (two old constraints and the new one)."""
    (an, bn, cn), refn = new
    # vertex of poly minimising an x + bn y
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
    raise RuntimeError("no Farkas certificate found")


def fmt(q):
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def main():
    targets = set()
    for ln in open("prop7_certificates.txt"):
        if ln.startswith("component"):
            targets.add(tuple(int(x) for x in ln.split("|")[0].split("indices")[1].split()))
    stats = {'S': 0, 'D': 0, 'L': 0}
    leaves = []
    t0 = time.time()
    out = open("tree_r7_25.txt", "w")
    out.write(f"# r = {R}; cells a, b in 0..64; functional order: {' '.join(str(t) for t in ORDER)}\n")

    def rec(poly, cons, depth, idx, corners, lines):
        if depth == len(ORDER):
            iv = tuple(idx[t] for t in range(18))
            assert iv in targets, iv
            lines.append("L")
            stats['L'] += 1
            leaves.append(iv)
            return
        t = ORDER[depth]
        a, b = FUN[t]
        vals = [a * p[0] + b * p[1] for p in corners]
        mn, mx = min(vals), max(vals)
        cands = [n for n in range(floor(mn - HI) - 1, floor(mx - LO) + 2) if n + HI >= mn and n + LO <= mx]
        for n in cands:
            lo_c = ((-a, -b, -(n + LO)), f"{t} lo")      # f >= n + r
            hi_c = ((a, b, n + HI), f"{t} hi")           # f <= n + 1 - r
            p1 = clip(poly, *lo_c[0])
            if not p1:
                cert = farkas(poly, cons, lo_c)
                lines.append(f"D {n} " + " ; ".join(f"{ref} {fmt(l)}" for (ref, l) in cert))
                stats['D'] += 1
                continue
            p2 = clip(p1, *hi_c[0])
            if not p2:
                cert = farkas(p1, cons + [lo_c], hi_c)
                lines.append(f"D {n} " + " ; ".join(f"{ref} {fmt(l)}" for (ref, l) in cert))
                stats['D'] += 1
                continue
            lines.append(f"S {n}")
            stats['S'] += 1
            idx2 = dict(idx)
            idx2[t] = n
            rec(p2, cons + [lo_c, hi_c], depth + 1, idx2, corners, lines)

    for a in range(65):
        for b in range(65):
            cons = [((Fr(1), Fr(0), a + HI), "x hi"), ((Fr(-1), Fr(0), -(a + LO)), "x lo"),
                    ((Fr(0), Fr(1), b + HI), "y hi"), ((Fr(0), Fr(-1), -(b + LO)), "y lo")]
            corners = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
            idx = {T_RE: a, T_IM: -b - 1}
            lines = [f"C {a} {b}"]
            rec(list(corners), cons, 0, idx, corners, lines)
            out.write("\n".join(lines) + "\n")
    out.close()
    print(f"tree written: {stats}; leaves {len(leaves)} ({len(set(leaves))} distinct index vectors); {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
