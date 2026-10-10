"""Rigorous check of a clique-strengthened LP bound (Oliveira Filho-Vallentin, JEMS 2010, Thm 1.1 and Sec. 3.1).

For a measurable planar set A avoiding the distances D, with f the radial average of the autocorrelation normalised by
f(0) = 1: f(d) = 0 for d in D, and f(|w_1 - p|) + ... + f(|w_m - p|) <= 1 for every point set {w_1..w_m} whose pairwise
distances all lie in D and every p (the sets {x in A : x + w_i - p in A} are disjoint subsets of A). Hence, if
   g(t) = z0 + sum_k z_k J0(d_k t) + sum_c y_c sum_i J0(r_ci t) >= 0 for all t >= 0, with y_c >= 0 and
   s = z0 + sum_k z_k + sum_c y_c m_c > 0,
then the upper density of A is at most (z0 + sum_c y_c) / s.
usage: python3 certify_lp2.py lp/NAME.json k [certificates/NAME.cert.json]
Exactness: pairwise clique distances are checked exactly with sympy (the minimal polynomial of d^2 - d_k^2 is X);
distances and radii are enclosed in Arb balls from their exact expressions; the coefficients are exact rationals.
On [0, T] the sweep evaluates g at points t in Arb and uses g >= g(t) - L (s - t) on [t, s], where
L = sum |coef| * radius / sqrt(2) bounds |g'| because |J_1(x)| <= 1/sqrt(2) for real x (DLMF 10.14.1, nu >= 1);
on [T, oo) it uses |J0(x)| <= sqrt(2/(pi x)) for x > 0 (Sonin-Polya).
"""
import sys, json, itertools
from fractions import Fraction as Fr
import sympy as sp
import flint
from flint import arb
from shapes import SHAPES

flint.ctx.prec = 96

def sym(e):
    return sp.sympify(e, locals={'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'pi': sp.pi})

def arb_of(e):
    """Arb enclosure of an exact sympy expression (rationals, pi, +, *, integer and half-integer powers, sin, cos)."""
    e = sp.sympify(e)
    if e.is_Rational:
        return arb(int(e.p)) / arb(int(e.q))
    if e is sp.pi:
        return arb.pi()
    if e.is_Add:
        return sum((arb_of(a) for a in e.args), arb(0))
    if e.is_Mul:
        out = arb(1)
        for a in e.args:
            out *= arb_of(a)
        return out
    if e.is_Pow and e.exp.is_Integer:
        return arb_of(e.base) ** int(e.exp)
    if e.is_Pow and e.exp.is_Rational and e.exp.q == 2:
        return arb_of(e.base).sqrt() ** int(e.exp.p)
    if isinstance(e, sp.sin):
        return arb_of(e.args[0]).sin()
    if isinstance(e, sp.cos):
        return arb_of(e.args[0]).cos()
    raise ValueError(f'no enclosure for {e}')

def main():
    J = json.load(open(sys.argv[1])); k = int(sys.argv[2])
    Dsym = [sym(e) for e in J['D'].split(';')]
    D = [arb_of(d) for d in Dsym]
    grid = Fr(J.get('grid', '1'))
    digits = 10**9
    z = [Fr(round(v * digits), digits) for v in J['z'][1:]]
    # distance terms
    coef_r = [(z[i], D[i]) for i in range(len(D))]
    ysum = Fr(0); mterm = Fr(0)
    for c in J['cliques']:
        y = Fr(round(c['y'] * digits), digits)
        if y <= 0:
            continue
        W = [(sym(x), sym(yv)) for x, yv in SHAPES[c['shape']]]
        for a, b in itertools.combinations(range(len(W)), 2):
            d2 = (W[a][0] - W[b][0])**2 + (W[a][1] - W[b][1])**2
            if not any(sp.minimal_polynomial(d2 - dd**2, sp.Symbol('X')) == sp.Symbol('X') for dd in Dsym):
                raise SystemExit(f"not a clique: {c['shape']}, squared distance {d2}")
        cx = sum(w[0] for w in W) / len(W); cy = sum(w[1] for w in W) / len(W)
        if 'ox' in c:
            px = cx + sp.Rational(c['ox']); py = cy + sp.Rational(c['oy'])
        else:
            px = cx + sp.Rational(c['i']) * sp.Rational(grid.numerator, grid.denominator)
            py = cy + sp.Rational(c['j']) * sp.Rational(grid.numerator, grid.denominator)
        for (wx, wy) in W:
            r2 = sp.expand((wx - px)**2 + (wy - py)**2)
            coef_r.append((y, arb_of(r2).sqrt()))
        ysum += y; mterm += y * len(W)
    # rigorous minimum of h(t) = sum coef J0(r t) on [0, T]
    sq2 = arb(2).sqrt()
    L = sum(abs(arb(cf.numerator) / arb(cf.denominator)) * r.upper() for cf, r in coef_r) / sq2
    zest = J['z'][0]
    # tail threshold: sum |coef| sqrt(2/(pi r T)) <= zest/3
    S = sum(abs(float(cf)) / float(r.lower().mid()) ** 0.5 for cf, r in coef_r)
    T = int((3 * S / zest) ** 2 * 2 / 3.14159) + 2
    # adaptive rigorous sweep of [0, T]: with slack eps added to z0, if h(t) + z0 >= L*step at a point t then
    # h + z0 >= 0 on [t, t + step] (Lipschitz). We take z0 = zest + eps and verify h(t) + z0 - L*step >= 0 at every step.
    eps = Fr(1, 10**4)
    z0 = Fr(zest).limit_denominator(10**9) + eps
    Z0 = arb(z0.numerator) / arb(z0.denominator)
    t = Fr(0); nsteps = 0; lo = None
    hmin = Fr(1, 10**15)
    while t < T:
        tm = arb(t.numerator) / arb(t.denominator)
        h = sum(arb(cf.numerator) / arb(cf.denominator) * (r * tm).bessel_j(0) for cf, r in coef_r)
        val = (h + Z0).lower()          # rigorous lower bound of g(t)
        gl = float(val.mid()) - float(val.rad())
        if gl <= 0:
            raise SystemExit(f'g(t) <= 0 at t = {float(t)} (lower bound {gl}); increase eps')
        step = Fr(gl / float(L.upper().mid()) * 0.999)   # the exact (dyadic) value of the float
        if step < hmin:
            raise SystemExit(f'step too small at t = {float(t)}')
        # rigorous check of the Lipschitz condition with this step
        if not (h + Z0 - L * (arb(step.numerator) / arb(step.denominator))).lower() >= 0:
            raise SystemExit(f'Lipschitz step not certified at t = {float(t)}')
        t += step; nsteps += 1
        if lo is None or val < lo:
            lo = val
    print(f'  sweep: {nsteps} steps on [0, {T}], min lower bound of g = {float(lo.mid()):.6f}')
    tail = sum(abs(arb(cf.numerator) / arb(cf.denominator)) * (arb(2) / (arb.pi() * r * T)).sqrt() for cf, r in coef_r)
    if not arb(z0.numerator) / arb(z0.denominator) - tail > 0:
        raise SystemExit('tail bound fails')
    cert = {'D': J['D'], 'T': T, 'z0': str(z0), 'z': [str(v) for v in z],
            'cliques': [{'shape': c['shape'], 'ox': c.get('ox', str(Fr(c.get('i', 0)) * grid)),
                         'oy': c.get('oy', str(Fr(c.get('j', 0)) * grid)), 'y': str(Fr(round(c['y'] * digits), digits))}
                        for c in J['cliques'] if round(c['y'] * digits) > 0]}
    out = sys.argv[3] if len(sys.argv) > 3 else sys.argv[1].replace('.json', '') + '.cert.json'
    json.dump(cert, open(out, 'w'), indent=1)
    s = z0 + sum(z) + mterm
    if not s > 0:
        raise SystemExit('s <= 0')
    B = (z0 + ysum) / s
    print(f"D = {J['D']}: float LP {J['lp']:.6f}; {len(J['cliques'])} clique placements; T = {T}, adaptive steps")
    print(f"  certified z0 = {float(z0):.8f}, sum y = {float(ysum):.8f}, s = {float(s):.8f}")
    print(f"  bound B = {float(B):.8f} -> 1/B = {float(1/B):.5f};  B < 1/{k}: {B < Fr(1, k)}  (exact B = {B})")
    if not B < Fr(1, k):
        raise SystemExit(1)

if __name__ == '__main__':
    main()
