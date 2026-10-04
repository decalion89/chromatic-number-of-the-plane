# Referee check 3: thresholds of the induction step of Proposition 6 with the full hypotheses
# x in P_{k-1} (type c), resp. y in Y_{k-1}(eps) (type q), in the coordinates u = conj(x) rho^k.
# theta_C(k) = 1/2 - min s such that some u has |(u rho^{-m})_e| <= s (1<=m<=2k-1) and
#              u + nu in B_s + n  for some nu in Lambda_+ \ Z[i] (centred rep) and integer shift n;
# theta_Q(k) = same with eta's: |eta_{j,e}/6 + (u rho^{j-k})_e| <= s (|j|<=k-1) and
#              eta_k/6 + u + nu in B_s + n.  All four eps in E_q are treated, no symmetry used.
# The rho^{-k} / Lambda_- side is checked as well (directly, not by conjugation).
import sys, time
from fractions import Fraction as Fr
from gauss import G, I, RHO, H, modZ
from exactlp import certified_min_s, brute_min_s

def rows(fn):
    w1 = fn(G(1)); w2 = fn(G(0, 1))
    return [(w1.re, w2.re), (w1.im, w2.im)]

def add_box(A, b, rws, off):
    for (a1, a2), o in zip(rws, off):
        A.append([a1, a2, Fr(-1)]); b.append(-o)
        A.append([-a1, -a2, Fr(-1)]); b.append(o)

def eta(eps, k, j):
    N = 5 ** k
    v = (eps * N).conj() * (RHO ** j)
    return modZ(v - H) * 6

NU = {+1: [modZ(G(2, 1) / 5 * t) for t in range(1, 5)], -1: [modZ(G(2, -1) / 5 * t) for t in range(1, 5)]}
SH = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
EQ = [G(Fr(a, 3), Fr(b, 3)) for a in (1, 2) for b in (1, 2)]

def theta_C(k, sign, solver):
    best = None
    rp = {m: rows(lambda u, m=m: u * RHO ** (-sign * m)) for m in range(1, 2 * k)}
    for nu in NU[sign]:
        for sh in SH:
            A = []; b = []
            for m in range(1, 2 * k):
                add_box(A, b, rp[m], (Fr(0), Fr(0)))
            add_box(A, b, rows(lambda u: u), (nu.re - sh[0], nu.im - sh[1]))
            v = solver(A, b)[0]
            if v is not None and (best is None or v < best):
                best = v
    return Fr(1, 2) - best

def theta_Q(k, sign, solver):
    best = None
    rp = {}
    for j in range(-(k - 1), k):
        # u = conj(y) rho^{sign k};  conj(y) rho^j = u rho^{j - sign k}
        rp[j] = rows(lambda u, j=j: u * RHO ** (j - sign * k))
    for eps in EQ:
        et = {j: eta(eps, k, j) for j in range(-k, k + 1)}
        ek = et[sign * k]
        for nu in NU[sign]:
            for sh in SH:
                A = []; b = []
                for j in range(-(k - 1), k):
                    add_box(A, b, rp[j], (et[j].re / 6, et[j].im / 6))
                add_box(A, b, rows(lambda u: u), (ek.re / 6 + nu.re - sh[0], ek.im / 6 + nu.im - sh[1]))
                v = solver(A, b)[0]
                if v is not None and (best is None or v < best):
                    best = v
    return Fr(1, 2) - best

if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    # validate the certified solver against brute force on k = 2, 3
    for k in (2, 3):
        a = theta_C(k, 1, certified_min_s); b_ = theta_C(k, 1, brute_min_s)
        c = theta_Q(k, 1, certified_min_s); d = theta_Q(k, 1, brute_min_s)
        print(f"validation k={k}: C certified {a} brute {b_}; Q certified {c} brute {d}; agree: {a == b_ and c == d}")
    prev = None
    for k in range(2, kmax + 1):
        t0 = time.time()
        tc = theta_C(k, 1, certified_min_s); tcm = theta_C(k, -1, certified_min_s)
        tq = theta_Q(k, 1, certified_min_s); tqm = theta_Q(k, -1, certified_min_s)
        mx = max(tc, tq, tcm, tqm)
        mono = "" if prev is None else ("  non-increasing" if mx <= prev else "  INCREASED!")
        prev = mx
        fmt = lambda f: f"{float(f):.9f}" + (f" ({f})" if f.denominator < 10 ** 8 else "")
        print(f"k={k:2d}: theta_C = {fmt(tc)}  theta_Q = {fmt(tq)}  (minus side equal: {tc == tcm and tq == tqm}){mono}  [{time.time() - t0:.1f}s]", flush=True)
