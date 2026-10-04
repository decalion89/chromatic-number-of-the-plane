"""Base case k = 1 at r = 17/56 (s = 11/56): S_5^(r) inside S_1^(r) = h + Sq_s + Z[i]; write c = c*_5 + x + lam,
x in Sq_s, lam in Z[i] mod 5.  For each lam: mu_+ = conj(lam) rho (in Lambda_+), mu_- = conj(lam) conj(rho)
(in Lambda_-); list the pieces {x : conj(x) rho^{+-1} + mu_{+-} in Sq_s + Z[i]}."""
from snr import *
import sys
r = F(sys.argv[1]) if len(sys.argv) > 1 else F(17, 56)
s = F(1, 2) - r
k = 1; N = 5
cst = (F(5, 2), F(5, 2))
fams = []
for j in (1, -1):
    A, B = rho_pow(j); fams += [(A, B), (B, -A)]
rows = []
for a in range(5):
    for b in range(5):
        lam = (F(a), F(b))
        # mu_+ = conj(lam) rho, mu_- = conj(lam) conj(rho)
        mup = cmul(conj(lam), RHO); mum = cmul(conj(lam), conj(RHO))
        intp = mup[0].denominator == 1 and mup[1].denominator == 1
        intm = mum[0].denominator == 1 and mum[1].denominator == 1
        sq = [(cst[0] + lam[0] + dx, cst[1] + lam[1] + dy) for dx, dy in ((-s, -s), (s, -s), (s, s), (-s, s))]
        pieces = [sq]
        for (fa, fb) in fams:
            nxt = []
            for P in pieces:
                nxt.extend(q for _, q in strip_split(P, fa, fb, r, 1 - r))
            pieces = nxt
        desc = []
        for P in pieces:
            P = clean(P)
            x = [(v[0] - cst[0] - lam[0], v[1] - cst[1] - lam[1]) for v in P]
            kap = component_kappa_fast(P, 1)[0]
            desc.append(f"{len(P)}-gon, kappa {kap}, x-vertices {[(str(v[0]), str(v[1])) for v in x][:4]}{'...' if len(x) > 4 else ''}")
        rows.append(((a, b), intp, intm, desc))
for (a, b), intp, intm, desc in rows:
    print(f"lam = {a}+{b}i: mu_+ in Z[i]: {intp}, mu_- in Z[i]: {intm}: " + ("; ".join(desc) if desc else "empty"))
