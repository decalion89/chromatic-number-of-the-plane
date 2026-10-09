"""Checker for cert_K2M1_seven.txt (shares no code with seven_certificates_K2M1.py or certify_window.py).

The window G(2,1) has modulus 325 = (2+i)^2 (2-i)^2 (3+2i) (3-2i), and 325 rho^j sigma^l is the Gaussian integer
(2+i)^(2+j) (2-i)^(2-j) (3+2i)^(1+l) (3-2i)^(1-l), so every functional is (integer linear form)/325.  Checks:
  (1) the eight index vectors are exactly those labelled SEVEN in cert_K2M1.txt, with the same 7-points;
  (2) at its 7-point each vector is the vector of floors, and all 60 margins are at least 2/7;
  (3) each certificate: lam >= 0, sum lam = 1, the linear parts of the three margins cancel, and the constant is 2/7.
So the largest least margin over each of these eight components is exactly 2/7."""
from fractions import Fraction as Q

def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def power(z, e):
    r = (1, 0)
    for _ in range(e):
        r = mul(r, z)
    return r
def form(j, l, part):
    """integers (p, q) with Re/Im(conj(c) rho^j sigma^l) = (p x + q y)/325 for c = x + iy"""
    g = (1, 0)
    for z, e in (((2, 1), 2 + j), ((2, -1), 2 - j), ((3, 2), 1 + l), ((3, -2), 1 - l)):
        g = mul(g, power(z, e))
    return (g[0], g[1]) if part == "Re" else (g[1], -g[0])

def parse_keys(line):
    out = []
    for tok in line.split(":", 1)[1].split():
        j, l, part = tok.strip("()").split(",")
        out.append((int(j), int(l), part))
    return out

main = open("cert_K2M1.txt").read().splitlines()
sev = open("cert_K2M1_seven.txt").read().splitlines()
keys = parse_keys(main[1])
assert parse_keys(sev[1]) == keys and len(keys) == 30
assert sorted({(j, l) for j, l, _ in keys}) == [(j, l) for j in range(-2, 3) for l in range(-1, 2)]
FORM = {k: form(*k) for k in keys}
seven_main = {}
for ln in main[2:]:
    if "| SEVEN" in ln:
        vec = tuple(int(t) for t in ln.split("|")[0].split()[2:])
        pt = tuple(Q(t) for t in ln.split("SEVEN")[1].split())
        seven_main[vec] = pt
assert len(seven_main) == 8
seen = {}
for ln in sev[2:]:
    head, label, cert = [s.strip() for s in ln.split("|")]
    vec = tuple(int(t) for t in head.split()[2:])
    pt = tuple(Q(t) for t in label.split()[1:3])
    assert label.startswith("SEVEN") and len(vec) == 30
    seen[vec] = pt
    x, y = pt
    n = dict(zip(keys, vec))
    for k in keys:                                   # (2): floors and margins at the 7-point
        p, q = FORM[k]
        v = Q(p * x + q * y, 1) / 325
        fl = v.numerator // v.denominator
        assert fl == n[k] and v - fl >= Q(2, 7) and fl + 1 - v >= Q(2, 7), (k, v)
    assert cert.startswith("kappa = 2/7 ; cert ")
    terms = [t.split() for t in cert[len("kappa = 2/7 ; cert "):].split(";")]
    assert len(terms) == 3
    lx = ly = c0 = Q(0); s = Q(0)
    for j, l, part, sign, lam in terms:
        k = (int(j), int(l), part); lam = Q(lam)
        assert k in n and sign in "+-" and lam >= 0
        p, q = FORM[k]
        if sign == "+":                              # margin f(c) - n
            lx += lam * Q(p, 325); ly += lam * Q(q, 325); c0 -= lam * n[k]
        else:                                        # margin n + 1 - f(c)
            lx -= lam * Q(p, 325); ly -= lam * Q(q, 325); c0 += lam * (n[k] + 1)
        s += lam
    assert s == 1 and lx == 0 and ly == 0 and c0 == Q(2, 7), (s, lx, ly, c0)
assert seen == seven_main
print("verified: the 8 index vectors of family 7 in the window (2,1) have largest least margin exactly 2/7 "
      "(7-points with all margins >= 2/7; dual certificates with value 2/7)")
