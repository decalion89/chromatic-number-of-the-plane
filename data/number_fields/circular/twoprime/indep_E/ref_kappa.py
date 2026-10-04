"""Referee E: exact kappa (max over the component of the minimal margin) for each index vector.

Index vectors are read from the output of ref_clip.py (lines 'ORDER ...' and 'IV ...').
Margins: mu+_f(c) = f(c) - n_f, mu-_f(c) = n_f + 1 - f(c), f(c) = (p x + q y)/65, 36 of them.
kappa = max_{x,y} min_f mu(c)  (an LP in x, y, t).
Primal: enumerate all triples of constraints mu_i(c) >= t, solve the 3x3 system exactly, keep feasible
points, take the max t (the feasible region is a pointed polyhedron and t <= 1/2, so the optimum is a vertex).
Dual: for all triples of margins, solve lambda >= 0, sum lambda = 1, sum lambda_i * grad(mu_i) = 0; then
min mu <= sum lambda_i mu_i = constant for every c.  The minimum of these constants is the dual optimum.
Both are reported and must be equal.
"""
import sys
import itertools
from fractions import Fraction as Fr


def solve3(A, rhs):
    """Exact Cramer solve of 3x3; returns None if singular."""
    def det3(M):
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    d = det3(A)
    if d == 0:
        return None
    out = []
    for col in range(3):
        M = [row[:] for row in A]
        for r in range(3):
            M[r][col] = rhs[r]
        out.append(Fr(det3(M), 1) / d if isinstance(d, int) else det3(M) / d)
    return out


def margins(order, iv, D=65):
    """list of (gx, gy, const): mu(c) = gx*x + gy*y + const."""
    ms = []
    for (p, q), n in zip(order, iv):
        ms.append((Fr(p, D), Fr(q, D), Fr(-n)))          # f - n
        ms.append((Fr(-p, D), Fr(-q, D), Fr(n + 1)))     # n + 1 - f
    return ms


def kappa_primal(ms):
    best = None
    arg = None
    for i, j, k in itertools.combinations(range(len(ms)), 3):
        # mu_r(c) = t  for r in (i,j,k):  gx x + gy y - t = -const
        A = [[ms[r][0], ms[r][1], Fr(-1)] for r in (i, j, k)]
        rhs = [-ms[r][2] for r in (i, j, k)]
        sol = solve3(A, rhs)
        if sol is None:
            continue
        x, y, t = sol
        if all(g0 * x + g1 * y + c0 >= t for (g0, g1, c0) in ms):
            if best is None or t > best:
                best, arg = t, (x, y)
    return best, arg


def kappa_dual(ms):
    best = None
    cert = None
    for i, j, k in itertools.combinations(range(len(ms)), 3):
        A = [[ms[r][0] for r in (i, j, k)], [ms[r][1] for r in (i, j, k)], [Fr(1)] * 3]
        sol = solve3(A, [Fr(0), Fr(0), Fr(1)])
        if sol is None or min(sol) < 0:
            continue
        # check linear parts cancel exactly
        assert sum(l * ms[r][0] for l, r in zip(sol, (i, j, k))) == 0
        assert sum(l * ms[r][1] for l, r in zip(sol, (i, j, k))) == 0
        val = sum(l * ms[r][2] for l, r in zip(sol, (i, j, k)))
        if best is None or val < best:
            best, cert = val, ((i, j, k), sol)
    # pairs: only opposite margins of one functional, value 1/2
    return best, cert


def main():
    fn = sys.argv[1]
    order = None
    ivs = []
    for line in open(fn):
        if line.startswith("ORDER "):
            order = [tuple(int(t) for t in w.split(",")) for w in line.split()[1:]]
        elif line.startswith("IV "):
            ivs.append(tuple(int(t) for t in line.split()[1:]))
    print("file %s: %d index vectors" % (fn, len(ivs)))
    for iv in ivs:
        ms = margins(order, iv)
        kp, arg = kappa_primal(ms)
        kd, cert = kappa_dual(ms)
        names = []
        x, y = arg
        # identify the maximiser
        for aa in range(7):
            for bb in range(7):
                if (x, y) == (Fr(65 * aa, 7), Fr(65 * bb, 7)):
                    names.append("65(%d+%di)/7" % (aa, bb))
        print("cell (%2d,%2d): kappa primal = %s, dual = %s, equal: %s ; maximiser (%s, %s) %s ; dual cert on margins %s lambda %s"
              % (iv[0], iv[1], kp, kd, kp == kd, x, y, " ".join(names), cert[0], [str(l) for l in cert[1]]))


if __name__ == "__main__":
    main()
