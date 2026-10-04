"""Base case k = 1, referee's version.  c = c*_5 + x + lam, x in Sq_s, lam = a + bi (a, b mod 5).
For each class: integrality of mu_+- = conj(lam) rho^(+-1), and the exact least s for which
  |x_e| <= s,  |Re/Im(conj(x) rho^(+-1)) + nu_+-| <= s   (nu_+- over ALL representatives of mu_+- mod Z[i] in a box)
is feasible (exact 3-variable LP by vertex enumeration).  Then the probe's printed certificates are re-evaluated."""
from fractions import Fraction as Fr
import re
from hgeom import funcs, gmul, gconj, RHO, fl
from lp3x import lp_min

def reps(mu, box=Fr(3, 2)):
    out = []
    base = (mu[0] - fl(mu[0]), mu[1] - fl(mu[1]))
    for m in range(-3, 4):
        for n in range(-3, 4):
            z = (base[0] + m, base[1] + n)
            if abs(z[0]) <= box and abs(z[1]) <= box:
                out.append(z)
    return out

def least_s(nup, num):
    cons = []
    for f in funcs(0):
        cons += [(f[0], f[1], Fr(-1), Fr(0)), (-f[0], -f[1], Fr(-1), Fr(0))]
    for j, nu in ((1, nup), (-1, num)):
        for e, f in enumerate(funcs(j)):
            cons += [(f[0], f[1], Fr(-1), -nu[e]), (-f[0], -f[1], Fr(-1), nu[e])]
    return lp_min((0, 0, 1), cons)

summary = {}
for a in range(5):
    for b in range(5):
        lam = (Fr(a), Fr(b))
        mup = gmul(gconj(lam), RHO)
        mum = gmul(gconj(lam), gconj(RHO))
        ip = mup[0].denominator == 1 and mup[1].denominator == 1
        im = mum[0].denominator == 1 and mum[1].denominator == 1
        best = None
        for nup in reps(mup):
            for num in reps(mum):
                res = least_s(nup, num)
                if res is not None and (best is None or res[0] < best[0]):
                    best = (res[0], res[1], nup, num)
        kind = "both integral" if ip and im else ("exactly one integral" if ip != im else "both non-integral")
        summary.setdefault(kind, []).append(best[0])
        print(f"lam = {a}+{b}i: {kind:22s} least s = {best[0]}  at x = ({best[1][0]}, {best[1][1]})")
for kind, v in summary.items():
    print(kind, len(v), "classes; least s values:", sorted(set(v)))

# re-evaluate the probe's printed certificates (base_caseC.txt) with the referee's own linear forms
import os
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "quadwt", "data", "number_fields", "circular", "probe", "base_caseC.txt")
def form(expr):
    # returns (c1, c2, cs, const) of the inequality  c1 x1 + c2 x2 + const <= cs * s  written as LHS <= s
    expr = expr.strip()
    m = re.fullmatch(r"(-?)(Re|Im)\(conj x\) <= s", expr)
    if m:
        sg = -1 if m.group(1) else 1
        f = funcs(0)[0 if m.group(2) == "Re" else 1]
        return (sg * f[0], sg * f[1], Fr(0))
    m = re.fullmatch(r"(-\[)?(Re|Im)\(conj x rho\^(-?1)\) - \((-?\d+/\d+)\)\]? <= s", expr)
    if m:
        sg = -1 if m.group(1) else 1
        f = funcs(int(m.group(3)))[0 if m.group(2) == "Re" else 1]
        v = Fr(m.group(4))
        return (sg * f[0], sg * f[1], -sg * v)
    raise ValueError(expr)

for line in open(P):
    if "certificate:" not in line:
        continue
    head, cert = line.split("certificate:")
    thr = Fr(re.search(r"min s = (\S+) =", head).group(1))
    terms = re.findall(r"(\d+/\d+|\d+)\*\[(.*?<= s)\]", cert)
    tot = [Fr(0), Fr(0), Fr(0)]; ms = Fr(0)
    for mult, ex in terms:
        mult = Fr(mult); c = form(ex)
        for t in range(3):
            tot[t] += mult * c[t]
        ms += mult
    ok = tot[0] == 0 and tot[1] == 0 and ms == 1 and tot[2] == thr
    print(f"certificate for {head.split(':')[0]}: {len(terms)} terms, sum gives {tot[2]} <= s  (threshold {thr}): {'OK' if ok else 'FAIL'}")
