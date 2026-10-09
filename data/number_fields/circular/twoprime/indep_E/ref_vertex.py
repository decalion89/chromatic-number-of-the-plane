"""Referee E, method B: vertex enumeration with exact integer arithmetic (numpy int64, overflow-checked).

S = {c : (p x + q y)/D in [theta, 1-theta] + Z for all functionals (p,q)}, theta = a/b, is invariant under D Z[i].
Every nonempty connected component (a convex polygon with fixed strip indices, inside one unit cell) has a vertex,
which is the intersection of two boundary lines  b(p x + q y) = K, K in {D(b n + a), D(b n + b - a)}.
We enumerate all such intersection points in the closed square [0, D]^2, keep those in S, and group them by
index vector (n_f = floor of each functional).  Components = distinct index vectors of points with x, y in [0, D).

Usage: python3 ref_vertex.py a b [window]   window in {65, cross65, 85, 221, 65minus:<k>} (default 65)
"""
import sys
import itertools
import numpy as np

from ref_rot import gmul, eighteen_functionals


def functionals_for(window):
    if window == "65":
        return 65, eighteen_functionals()
    if window == "cross65":
        # {1, rho^{+-1}, sigma^{+-1}} times units: 5 rotation classes -> 10 functionals, period 65
        rots = [(65, 0), (39, 52), (39, -52), (25, 60), (25, -60)]
        fs = set()
        for (p, q) in rots:
            for u in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
                z = gmul(u, (p, q))
                if z[0] > 0 or (z[0] == 0 and z[1] > 0):
                    fs.add(z)
        return 65, sorted(fs)
    if window in ("85", "221"):
        D = int(window)
        sols = [(p, q) for p in range(-D, D + 1) for q in range(-D, D + 1) if p * p + q * q == D * D]
        reps = sorted(set(v for v in sols if v[0] > 0 or (v[0] == 0 and v[1] > 0)))
        return D, reps
    raise ValueError(window)


def run(a, b, D, funcs, chunk=4_000_000):
    """Return dict index_vector -> list of points (as (X, Y, Den) integer triples, point = (X/Den, Y/Den))."""
    F = np.array(funcs, dtype=np.int64)  # shape (m, 2)
    m = len(funcs)
    # boundary constants per functional: lines meeting [0,D]^2
    Ks = []
    for (p, q) in funcs:
        vals = [0, p * D, q * D, (p + q) * D]  # D * value range of p x + q y on the square's corners, times 1
        lo, hi = min(vals) // D - 2, max(vals) // D + 2  # range of (p x + q y)/D, i.e. n range (generous)
        ks = []
        for n in range(lo, hi + 1):
            ks.append(D * (b * n + a))
            ks.append(D * (b * n + b - a))
        Ks.append(np.array(ks, dtype=np.int64))
    found = {}
    npts = 0
    nin = 0
    for i1, i2 in itertools.combinations(range(m), 2):
        p1, q1 = funcs[i1]
        p2, q2 = funcs[i2]
        det = p1 * q2 - p2 * q1
        if det == 0:
            continue
        K1 = Ks[i1][:, None]
        K2 = Ks[i2][None, :]
        # b(p1 x + q1 y) = K1, b(p2 x + q2 y) = K2  ->  x = (K1 q2 - K2 q1)/(b det), y = (p1 K2 - p2 K1)/(b det)
        X = (K1 * q2 - K2 * q1).ravel()
        Y = (p1 * K2 - p2 * K1).ravel()
        Den = b * det
        if Den < 0:
            X, Y, Den = -X, -Y, -Den
        # keep points in [0, D)^2 :  0 <= X < D*Den
        keep = (X >= 0) & (X < D * Den) & (Y >= 0) & (Y < D * Den)
        X, Y = X[keep], Y[keep]
        npts += X.size
        if X.size == 0:
            continue
        M = D * Den  # value of functional = (p X + q Y)/(D Den)
        # overflow guards (int64): |X|,|Y| < D*Den after filtering, |p|,|q| <= D
        assert int(Ks[i1].__abs__().max()) * D * 2 < 2 ** 62 and int(Ks[i2].__abs__().max()) * D * 2 < 2 ** 62
        assert 2 * D * (D * Den) < 2 ** 62
        assert b * M < 2 ** 62
        V = F[:, 0:1] * X[None, :] + F[:, 1:2] * Y[None, :]  # (m, npts)
        R = np.mod(V, M)
        ok = np.all((a * M <= b * R) & (b * R <= (b - a) * M), axis=0)
        if not ok.any():
            continue
        X, Y, V = X[ok], Y[ok], V[:, ok]
        nin += X.size
        Nidx = np.floor_divide(V, M)  # strip indices (floor), shape (m, k)
        for col in range(X.size):
            iv = tuple(int(t) for t in Nidx[:, col])
            found.setdefault(iv, set()).add((int(X[col]), int(Y[col]), int(Den)))
    return found, npts, nin


def main():
    a, b = int(sys.argv[1]), int(sys.argv[2])
    window = sys.argv[3] if len(sys.argv) > 3 else "65"
    D, funcs = functionals_for(window)
    found, npts, nin = run(a, b, D, funcs)
    print("window %s: D=%d, %d functionals, theta=%d/%d ; candidate points in [0,D)^2: %d ; in S: %d ; components: %d"
          % (window, D, len(funcs), a, b, npts, nin, len(found)))
    from fractions import Fraction as Fr
    for iv, pts in sorted(found.items(), key=lambda kv: min((Fr(x, d), Fr(y, d)) for x, y, d in kv[1])):
        pl = sorted(set((Fr(x, d), Fr(y, d)) for x, y, d in pts))
        cx = sum(P[0] for P in pl) / len(pl)
        cy = sum(P[1] for P in pl) / len(pl)
        print("  component: %3d vertices, centroid (%.4f, %.4f)" % (len(pl), float(cx), float(cy)))


if __name__ == "__main__":
    main()
