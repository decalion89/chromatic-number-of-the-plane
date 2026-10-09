"""Certificates for the two-prime window proposition (window |j|, |l| <= 1, primes 5 and 13, N = 65).

Enumerates S^r(1,1) at r = 7/25 (< 2/7) with probe2d.py (lift 5 then 13), and writes for every component:
  * its strip-index vector: for the 18 functionals f_{j,l,e}(c) = Re / Im (conj(c) rho^j sigma^l), j, l in {-1,0,1},
    the integer n with f(c) in [n + r, n + 1 - r] on the component;
  * main components: the type point T (65h or (65/3)(a + bi)) and m in Z[i] with T + 65 m having the same indices;
  * other components: an exact LP certificate that kappa <= 2/7 (three constraints t <= +-(f - n) + const with
    multipliers lam >= 0, sum lam = 1, the linear parts cancelling, and sum lam * const = 2/7), and the 7-torsion
    point of the component with all margins >= 2/7 (so kappa = 2/7 exactly).
Output: prop7_certificates.txt (plain text, exact rationals).  Checked by verify_prop7.py (independent code)."""
from fractions import Fraction as Fr
from probe2d import Probe, centroid, gamma, floor
from lpexact import lp_max

R0 = Fr(7, 25)
FUN = []                                   # (j, l, e, (a, b)) with f(c) = a c1 + b c2
for j in (-1, 0, 1):
    for l in (-1, 0, 1):
        A, B = gamma(j, l)
        FUN.append((j, l, 'Re', (A, B)))     # Re(conj(c) g) = A c1 + B c2
        FUN.append((j, l, 'Im', (B, -A)))    # Im(conj(c) g) = B c1 - A c2


def ev(f, p):
    return f[0] * p[0] + f[1] * p[1]


def index_vector(p, r=R0):
    return tuple(floor(ev(f, p) - r) for (_, _, _, f) in FUN)


def main():
    P = Probe(R0)
    P.lift5()
    P.lift13()
    N = P.N
    types = [('C', (Fr(N, 2), Fr(N, 2)))] + [(f'Q{a}{b}', (Fr(N * a, 3), Fr(N * b, 3))) for a in (1, 2) for b in (1, 2)]
    out = []
    out.append(f"# two-prime window (1,1), N = {N}, r = {R0}; {len(P.comps)} components")
    out.append("# functionals (j, l, Re/Im): " + " ".join(f"({j},{l},{e})" for (j, l, e, _) in FUN))
    nX = 0
    for idx, C in enumerate(P.comps):
        pt = centroid(C)
        iv = index_vector(pt)
        lab = None
        for name, T in types:
            # translate T by 65 m (m in Z[i]) towards the cell of pt, then test nearby translates
            m = ((floor(pt[0]) - floor(T[0])) // N, (floor(pt[1]) - floor(T[1])) // N)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    mm = (m[0] + dx, m[1] + dy)
                    TT = (T[0] + N * mm[0], T[1] + N * mm[1])
                    if index_vector(TT) == iv:
                        lab = (name, mm)
        line = f"component {idx}: indices {' '.join(str(n) for n in iv)}"
        if lab is not None:
            out.append(line + f" | MAIN {lab[0]} m = {lab[1][0]} {lab[1][1]}")
            continue
        nX += 1
        cons = []
        meta = []
        for (j, l, e, f), n in zip(FUN, iv):
            # t <= f(c) - n   and   t <= n + 1 - f(c)
            cons.append((f[0], f[1], Fr(-n)))
            meta.append((j, l, e, +1, n))
            cons.append((-f[0], -f[1], Fr(n + 1)))
            meta.append((j, l, e, -1, n))
        t, x, y, lam = lp_max(cons)
        assert t == Fr(2, 7), t
        cert = " ; ".join(f"{meta[i][0]} {meta[i][1]} {meta[i][2]} {'+' if meta[i][3] > 0 else '-'} {lam[i]}" for i in sorted(lam))
        # the 7-torsion point: N (a + bi)/7 nearest to the optimum
        sev = (Fr(round(7 * x / N)) * N / 7, Fr(round(7 * y / N)) * N / 7)
        out.append(line + f" | EXTRA kappa = {t} ; optimum ({x}, {y}) ; seven-point ({sev[0]}, {sev[1]}) ; cert {cert}")
    out.append(f"# {len(P.comps) - nX} main, {nX} extra")
    with open("prop7_certificates.txt", "w") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n".join(out[-3:]))


if __name__ == "__main__":
    main()
