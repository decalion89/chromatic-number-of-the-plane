"""Referee E, method A: direct exact enumeration of S^theta(1,1) modulo 65 Z[i].

For each of the 65^2 unit cells (fixed by the functionals x and y, i.e. gamma = 1 and gamma = i),
the convex polygon {c : all 18 functionals in their strips} is computed by exact clipping of the
cell square (vertex representation, Fractions), branching over the strip index of each functional.
Degenerate pieces (segments, points) are kept: the strips are closed.

Usage: python3 ref_clip.py NUM DEN [outfile]     (theta = NUM/DEN)
"""
import sys
import time
from fractions import Fraction as Fr
from math import floor, ceil

from ref_rot import eighteen_functionals

D = 65


def fval(f, P):
    return (f[0] * P[0] + f[1] * P[1]) / D


def clip(poly, f, C, keep_le):
    """Keep the part of convex polygon poly (list of vertices, possibly degenerate) where
    f(P) <= C (keep_le True) or f(P) >= C (keep_le False).  Exact."""
    if not poly:
        return []
    out = []
    m = len(poly)
    vals = []
    for P in poly:
        v = fval(f, P) - C
        vals.append(v if keep_le else -v)
    for k in range(m):
        P, Q = poly[k], poly[(k + 1) % m]
        vP, vQ = vals[k], vals[(k + 1) % m]
        if vP <= 0:
            out.append(P)
        if (vP < 0 < vQ) or (vQ < 0 < vP):
            t = vP / (vP - vQ)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    # remove consecutive duplicates (and wrap-around duplicate)
    res = []
    for P in out:
        if not res or res[-1] != P:
            res.append(P)
    while len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


def run(theta, order=None, verbose=True):
    funcs = eighteen_functionals()
    fx, fy = (65, 0), (0, 65)
    others = [f for f in funcs if f not in (fx, fy)]
    assert len(others) == 16
    if order is not None:
        others = [others[k] for k in order]
    leaves = []
    nodes = 0
    for a in range(D):
        for b in range(D):
            x0, x1 = a + theta, a + 1 - theta
            y0, y1 = b + theta, b + 1 - theta
            sq = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            stack = [(sq, 0, (a, b))]
            while stack:
                poly, depth, idx = stack.pop()
                nodes += 1
                if depth == 16:
                    leaves.append((poly, idx))
                    continue
                f = others[depth]
                vs = [fval(f, P) for P in poly]
                lo, hi = min(vs), max(vs)
                # strips [n+theta, n+1-theta] meeting [lo, hi]
                nmin = ceil(lo - 1 + theta)
                nmax = floor(hi - theta)
                for n in range(nmin, nmax + 1):
                    q1 = clip(poly, f, n + theta, False)
                    q2 = clip(q1, f, n + 1 - theta, True)
                    if q2:
                        stack.append((q2, depth + 1, idx + (n,)))
    return funcs, [fx, fy] + others, leaves, nodes


def index_vector(funcs_order, P, theta):
    """strip indices of the point P (Fractions) for the functionals in funcs_order; None if P not in S^theta."""
    out = []
    for f in funcs_order:
        v = fval(f, P)
        n = floor(v)
        fr = v - n
        if not (theta <= fr <= 1 - theta):
            return None
        out.append(n)
    return tuple(out)


def main():
    num, den = int(sys.argv[1]), int(sys.argv[2])
    theta = Fr(num, den)
    outfile = sys.argv[3] if len(sys.argv) > 3 else None
    t0 = time.time()
    funcs, order, leaves, nodes = run(theta)
    t1 = time.time()
    lines = []
    lines.append("theta = %s = %.10f ; nodes visited %d ; leaves (index vectors mod 65Z[i]) %d ; time %.1fs"
                 % (theta, float(theta), nodes, len(leaves), t1 - t0))
    # type points and 7-torsion points, reduced into [0,65)^2
    h = (Fr(65, 2), Fr(65, 2))
    types = {"Tc": h}
    for al in (1, 2):
        for be in (1, 2):
            types["T%d%d" % (al, be)] = (Fr(65 * al, 3), Fr(65 * be, 3))
    tors = {}
    for aa in range(7):
        for bb in range(7):
            if (aa, bb) != (0, 0):
                tors["7t(%d,%d)" % (aa, bb)] = (Fr(65 * aa, 7), Fr(65 * bb, 7))
    special = dict(types)
    special.update(tors)
    special_iv = {}
    for name, P in special.items():
        iv = index_vector(order, P, theta)
        if iv is not None:
            special_iv.setdefault(iv, []).append(name)
    ntype = 0
    for k, (poly, idx) in enumerate(leaves):
        names = special_iv.get(idx, [])
        if any(nm.startswith("T") for nm in names):
            ntype += 1
        xs = [float(P[0]) for P in poly]
        ys = [float(P[1]) for P in poly]
        lines.append("leaf %2d: cell (%d,%d) nverts %d  x in [%.4f,%.4f] y in [%.4f,%.4f]  contains: %s"
                     % (k, idx[0], idx[1], len(poly), min(xs), max(xs), min(ys), max(ys), ",".join(names) or "-"))
    lines.append("leaves containing a type point: %d ; others: %d" % (ntype, len(leaves) - ntype))
    lines.append("special points in S^theta (name -> found as a leaf?):")
    leafset = set(idx for _, idx in leaves)
    for iv, names in special_iv.items():
        lines.append("  %s : %s" % (",".join(names), "leaf" if iv in leafset else "NOT A LEAF (error)"))
    txt = "\n".join(lines)
    print(txt)
    if outfile:
        with open(outfile, "w") as fh:
            fh.write(txt + "\n")
            # also store the index vectors (for the LP step), one per line
            fh.write("ORDER " + " ".join("%d,%d" % f for f in order) + "\n")
            for poly, idx in leaves:
                fh.write("IV " + " ".join(str(n) for n in idx) + "\n")


if __name__ == "__main__":
    main()
