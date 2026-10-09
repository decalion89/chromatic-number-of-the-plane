"""Print the exact LP certificates behind the lemmas (C_2) and (Q_2), in readable form, for every nu in
Lambda_+ \\ Z[i] with |nu| <= 1.3 (C) and every q in Q_25 (Q), both signs; and check them (lp3.check_certificate)."""
from fractions import Fraction as F
from lp3 import lp3_max, check_certificate
from lemma_lp import lam, lin, pattern
from snr import rho_pow

def show_row(r, names):
    a1, a2, a3, b = r
    terms = []
    for coef, nm in zip((a1, a2, a3), names):
        if coef != 0:
            terms.append(f"({coef})*{nm}")
    return " + ".join(terms) + f" <= {b}"

k = 2
print("== (C_2): x in P_1^(s), conj(x) rho^{2 sign} + nu in Sq_s, nu in Lambda_sign \\ Z[i]; variables x1, x2, s")
for sign in (1, -1):
    res = []
    for nu in lam(sign):
        rows = []
        labels = []
        for j in range(-(k - 1), k):
            for e, (a, b) in enumerate(lin(j)):
                nm = ("Re" if e == 0 else "Im") + f"(conj(x) rho^{j})"
                rows.append((a, b, F(-1), F(0))); labels.append(f"{nm} <= s")
                rows.append((-a, -b, F(-1), F(0))); labels.append(f"-{nm} <= s")
        for e, (a, b) in enumerate(lin(sign * k)):
            nm = ("Re" if e == 0 else "Im") + f"(conj(x) rho^{sign*k}) + ({nu[e]})"
            rows.append((a, b, F(-1), -nu[e])); labels.append(f"{nm} <= s")
            rows.append((-a, -b, F(-1), nu[e])); labels.append(f"-[{nm}] <= s")
        cert = lp3_max((0, 0, -1), rows)
        check_certificate((0, 0, -1), rows, cert)
        res.append((-cert[0], nu, cert, rows, labels))
    res.sort(key=lambda t: t[0])
    print(f" sign {sign}: thresholds over nu: " + ", ".join(f"nu=({nu[0]},{nu[1]}): s>={t}" for t, nu, *_ in res[:6]) + " ...")
    t, nu, cert, rows, labels = res[0]
    val, z, tri, m = cert
    print(f"   minimum s = {t} at x = ({z[0]}, {z[1]}), nu = ({nu[0]}, {nu[1]}); certificate (multiplier * inequality):")
    for mi, i in zip(m, tri):
        print(f"      {mi} * [ {labels[i]} ]   i.e. {show_row(rows[i], ('x1','x2','s'))}")

print("\n== (Q_2): q in Q_25, y in Y_1(q), conj(q+y) rho^{2 sign} + nu in h + Sq_s + Z[i]; variables y1, y2, eps (r = 1/3 - eps)")
N = 25
for a0 in (1, 2):
    for b0 in (1, 2):
        q = (F(N * a0, 3), F(N * b0, 3))
        for sign in (1, -1):
            best = None
            for nu in lam(sign):
                rows = []; labels = []
                for j in range(-(k - 1), k):
                    p = pattern(q, j)
                    for e, (a, b) in enumerate(lin(j)):
                        nm = ("Re" if e == 0 else "Im") + f"(conj(y) rho^{j})"
                        rows.append((-a, -b, F(-1), p[e] - F(1, 3))); labels.append(f"{p[e]} + {nm} >= 1/3 - eps")
                        rows.append((a, b, F(-1), F(2, 3) - p[e])); labels.append(f"{p[e]} + {nm} <= 2/3 + eps")
                rows.append((F(0), F(0), F(-1), F(0))); labels.append("eps >= 0")
                pk = pattern(q, sign * k)
                for e, (a, b) in enumerate(lin(sign * k)):
                    c0 = pk[e] + nu[e] - F(1, 2)
                    nm = ("Re" if e == 0 else "Im") + f"(conj(y) rho^{sign*k})"
                    rows.append((a, b, F(-1), F(1, 6) - c0)); labels.append(f"{nm} + ({c0}) <= 1/6 + eps")
                    rows.append((-a, -b, F(-1), F(1, 6) + c0)); labels.append(f"-[{nm} + ({c0})] <= 1/6 + eps")
                cert = lp3_max((0, 0, -1), rows)
                assert cert is not None
                check_certificate((0, 0, -1), rows, cert)
                eps = -cert[0]
                if best is None or eps < best[0]:
                    best = (eps, nu, cert, rows, labels, pk)
            eps, nu, cert, rows, labels, pk = best
            val, z, tri, m = cert
            print(f" q = 25/3*({a0}+{b0}i), sign {sign}: values at q: j=-1: {pattern(q,-1)}, j=0: {pattern(q,0)}, j=1: {pattern(q,1)}, j={2*sign}: {pk}")
            print(f"   min eps = {eps} (r = {F(1,3)-eps}) at y = ({z[0]}, {z[1]}), nu = ({nu[0]}, {nu[1]}); certificate:")
            for mi, i in zip(m, tri):
                print(f"      {mi} * [ {labels[i]} ]")
