"""Referee: exact post-processing of the output of vertex_enum.c.

* rebuilds the rotations (x^2 + y^2 = N^2) and the functionals itself and checks them against the C header;
* re-checks every printed point exactly (Fractions): it lies in [0,N)^2 and in S^r(N);
* groups the points by their vector of strip indices (one group = one component of S^r(N) modulo N Z[i]);
* labels each component: 'c' (contains N h), 'q' (contains (N/3)(a+bi), a,b in {1,2}), '7' (contains a point
  N(a+bi)/7), or 'other';
* computes kappa = max over the component of the smallest of its 2*nf margins EXACTLY: by enumerating the
  vertices of the 3-dimensional LP  max t  s.t.  margin_i(c) >= t  (all triples of constraints that can be tight,
  checked against all constraints), and gives a dual certificate: <= 3 tight margins with weights lam >= 0,
  sum lam = 1, whose combination has zero linear part and constant value kappa (so no point of the component has
  all margins > kappa).
usage: python3 analyze_components.py file.txt [labelfile]"""
import sys
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict

fn = sys.argv[1]
N = rn = rd = None
cfun = []
raw = []
for line in open(fn):
    if line.startswith('# N'):
        p = line.split()
        N, rn, rd = int(p[2]), int(p[4]), int(p[5])
    elif line.startswith('# F'):
        p = line.split()
        cfun.append((int(p[2]), int(p[3])))
    elif line.startswith('#') or not line.strip():
        continue
    else:
        raw.append(tuple(map(int, line.split())))
r = Fr(rn, rd)
lo, hi = r, 1 - r

# rotations and functionals, rebuilt here
rots = [(A, B) for A in range(-N, N + 1) for B in range(-N, N + 1) if A * A + B * B == N * N]
funs = []
for (A, B) in rots:
    for (a, b) in ((A, B), (B, -A)):
        if a < 0 or (a == 0 and b < 0):
            a, b = -a, -b
        if (a, b) not in funs:
            funs.append((a, b))
assert sorted(funs) == sorted(cfun), "functional list differs from the C program"
funs = sorted(funs)
nf = len(funs)
print(f"N = {N}, r = {r} = {float(r):.8f}: {len(rots)} rotations, {nf} functionals up to sign; "
      f"{len(raw)} points printed by vertex_enum")


def fl(q):
    return q.numerator // q.denominator


def val(f, p):
    return (f[0] * p[0] + f[1] * p[1]) / N


def in_S(p):
    for f in funs:
        v = val(f, p)
        fr = v - fl(v)
        if fr < lo or fr > hi:
            return False
    return True


def ivec(p):
    return tuple(fl(val(f, p)) for f in funs)


pts = set()
for (X, Y, D) in raw:
    p = (Fr(N * X, D), Fr(N * Y, D))
    assert 0 <= p[0] < N and 0 <= p[1] < N
    assert in_S(p), p
    pts.add(p)
comps = defaultdict(set)
for p in pts:
    comps[ivec(p)].add(p)
print(f"{len(pts)} distinct points, {len(comps)} components (distinct index vectors)")

# type points reduced to [0,N)^2
def red(q):
    return q - N * fl(q / N)


types = [('c', (red(Fr(N, 2)), red(Fr(N, 2))))]
for a in (1, 2):
    for b in (1, 2):
        types.append((f'q{a}{b}', (red(Fr(N * a, 3)), red(Fr(N * b, 3)))))
for a in range(7):
    for b in range(7):
        if (a, b) != (0, 0):
            types.append((f'7:{a},{b}', (red(Fr(N * a, 7)), red(Fr(N * b, 7)))))
tidx = {}
for name, T in types:
    if in_S(T):
        tidx.setdefault(ivec(T), []).append(name)
print("type points lying in S^r:", sum(len(v) for v in tidx.values()), sorted(n for v in tidx.values() for n in v))


def constraints(iv):
    """margins scaled by N:  a x + b y + d >= N t."""
    cs = []
    for f, n in zip(funs, iv):
        cs.append((f[0], f[1], -N * n, ('+', f)))          # N*(f(c) - n)
        cs.append((-f[0], -f[1], N * (n + 1), ('-', f)))   # N*(n + 1 - f(c))
    return cs


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def kappa(iv, verts):
    cs = constraints(iv)
    # U = min_i max_{polygon} margin_i  (an upper bound for kappa); a constraint whose minimum over the polygon
    # exceeds U cannot be tight at the optimum.
    def mval(c, p):
        return c[0] * p[0] + c[1] * p[1] + c[2]
    U = min(max(mval(c, p) for p in verts) for c in cs)
    cand = [c for c in cs if min(mval(c, p) for p in verts) <= U]
    best = None
    for (c1, c2, c3) in combinations(cand, 3):
        M = [[c1[0], c1[1], -1], [c2[0], c2[1], -1], [c3[0], c3[1], -1]]
        dd = det3(M)
        if dd == 0:
            continue
        rhs = [-c1[2], -c2[2], -c3[2]]
        Mx = [[rhs[i], M[i][1], M[i][2]] for i in range(3)]
        My = [[M[i][0], rhs[i], M[i][2]] for i in range(3)]
        Mt = [[M[i][0], M[i][1], rhs[i]] for i in range(3)]
        x, y, T = Fr(det3(Mx), dd), Fr(det3(My), dd), Fr(det3(Mt), dd)
        if best is not None and T <= best[0]:
            continue
        if all(c[0] * x + c[1] * y + c[2] >= T for c in cs):
            best = (T, (x, y))
    T, (x, y) = best
    kap = T / N
    # dual certificate among the tight constraints
    tight = [c for c in cs if c[0] * x + c[1] * y + c[2] == T]
    cert = None
    for k in (2, 3):
        for sub in combinations(tight, k):
            # find lam >= 0, sum lam = 1, sum lam (a,b) = 0
            if k == 2:
                (a1, b1), (a2, b2) = (sub[0][0], sub[0][1]), (sub[1][0], sub[1][1])
                if a1 * b2 - a2 * b1 == 0 and (a1 * a2 < 0 or b1 * b2 < 0) and (a1, b1) != (0, 0):
                    # parallel opposite: lam1 * |v1| = lam2 * |v2|
                    n1 = abs(a1) + abs(b1)
                    n2 = abs(a2) + abs(b2)
                    lam = [Fr(n2, n1 + n2), Fr(n1, n1 + n2)]
                    if lam[0] * a1 + lam[1] * a2 == 0 and lam[0] * b1 + lam[1] * b2 == 0:
                        cert = list(zip(lam, sub))
                        break
            else:
                A = [[sub[0][0], sub[1][0], sub[2][0]], [sub[0][1], sub[1][1], sub[2][1]], [1, 1, 1]]
                dd = det3(A)
                if dd == 0:
                    continue
                lam = []
                for col in range(3):
                    Ac = [row[:] for row in A]
                    for rr in range(3):
                        Ac[rr][col] = (0, 0, 1)[rr]
                    lam.append(Fr(det3(Ac), dd))
                if all(l >= 0 for l in lam):
                    cert = list(zip(lam, sub))
                    break
        if cert:
            break
    # verify the certificate: linear parts cancel, constant = kappa * N, lam >= 0, sum = 1
    ok = False
    if cert:
        sa = sum(l * c[0] for l, c in cert)
        sb = sum(l * c[1] for l, c in cert)
        sc = sum(l * c[2] for l, c in cert)
        ok = sa == 0 and sb == 0 and sc == T and sum(l for l, _ in cert) == 1 and all(l >= 0 for l, _ in cert)
    return kap, (x, y), cert, ok


rows = []
for iv, verts in comps.items():
    lab = tidx.get(iv, ['other'])
    kap, opt, cert, ok = kappa(iv, list(verts))
    rows.append((lab, kap, opt, cert, ok, len(verts), iv))
rows.sort(key=lambda t: (t[0][0][0], -t[1]))
from collections import Counter
cnt = Counter()
for lab, kap, opt, cert, ok, nv, iv in rows:
    fam = lab[0][0] if lab[0] != 'other' else 'other'
    cnt[fam] += 1
print("components by family:", dict(cnt))
for lab, kap, opt, cert, ok, nv, iv in rows:
    cs = "; ".join(f"{l}*mu{c[3][0]}{c[3][1]}" for l, c in cert) if cert else "NONE"
    print(f"  {','.join(lab):10s} vertices {nv:2d}  kappa = {kap} ({float(kap):.6f})  optimum ({opt[0]}, {opt[1]})"
          f"  dual cert ok: {ok}  [{cs}]")
kaps = Counter((lab[0][0] if lab[0] != 'other' else 'other', kap) for lab, kap, *_ in rows)
print("kappa summary (family, kappa): count ->", {f"{k[0]} {k[1]}": v for k, v in sorted(kaps.items())})
print("all dual certificates verified:", all(t[4] for t in rows))
