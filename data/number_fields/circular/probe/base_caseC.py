"""Base case k = 1, case C (both mu_+ and mu_- non-integral).  Variables (x1, x2, s); u = conj(x) rho,
u' = conj(x) conj(rho).  Conditions: conj(x) in Sq_s, u - nu in Sq_s, u' - nu' in Sq_s, for nu in the 4 nearest nonzero
points of Lambda_+ and nu' in the 4 nearest nonzero points of Lambda_-.  For each of the 16 pairs: the least s for
which the system is feasible (exact LP with certificate), and whether the corner x0 in {(+-1 +- i)/6} solves it at
s = 1/6."""
from fractions import Fraction as F
from lp3 import lp3_max, check_certificate
from snr import rho_pow, cmul, conj
NUP = [(F(2, 5), F(1, 5)), (F(-1, 5), F(2, 5)), (F(-2, 5), F(-1, 5)), (F(1, 5), F(-2, 5))]
NUM = [(a, -b) for (a, b) in NUP]
def lin(j):
    A, B = rho_pow(j)
    return [(A, B), (B, -A)]       # (Re, Im)(conj(x) rho^j)
def rows_for(nu, nup):
    rows = []; labels = []
    for e, (a, b) in enumerate(lin(0)):
        rows.append((a, b, F(-1), F(0))); labels.append(f"{'Re' if e==0 else 'Im'}(conj x) <= s")
        rows.append((-a, -b, F(-1), F(0))); labels.append(f"-{'Re' if e==0 else 'Im'}(conj x) <= s")
    for j, v in ((1, nu), (-1, nup)):
        for e, (a, b) in enumerate(lin(j)):
            nm = f"{'Re' if e==0 else 'Im'}(conj x rho^{j})"
            rows.append((a, b, F(-1), v[e])); labels.append(f"{nm} - ({v[e]}) <= s")
            rows.append((-a, -b, F(-1), -v[e])); labels.append(f"-[{nm} - ({v[e]})] <= s")
    return rows, labels
corners = [(F(a, 6), F(b, 6)) for a in (1, -1) for b in (1, -1)]
for nu in NUP:
    for nup in NUM:
        rows, labels = rows_for(nu, nup)
        cert = lp3_max((0, 0, -1), rows)
        check_certificate((0, 0, -1), rows, cert)
        smin = -cert[0]
        cs = [x0 for x0 in corners if all(r[0] * x0[0] + r[1] * x0[1] + r[2] * F(1, 6) <= r[3] for r in rows)]
        line = f"nu=({nu[0]},{nu[1]}) nu'=({nup[0]},{nup[1]}): min s = {smin} = {float(smin):.5f}; corners solving at s=1/6: {[(str(a), str(b)) for a, b in cs]}"
        if smin <= F(11, 56):
            print(line + "  <-- feasible below 11/56")
        else:
            val, z, tri, m = cert
            print(line + "; certificate: " + " + ".join(f"{mi}*[{labels[i]}]" for mi, i in zip(m, tri)))
