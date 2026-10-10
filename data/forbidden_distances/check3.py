"""Third check of a certificate NAME.cert.json written by certify_lp2.py: rigorous, and independent of both other programs.

Method (after the independent checker of the second referee):
  * the clique vertices are built by hand in exact arithmetic of Q(sqrt 3), the pairwise distances are checked
    exactly there, and each radius r with r^2 = a + b sqrt 3 (a, b rational) becomes an Arb ball without any
    decimal evaluation; the distances in D are parsed by a small arithmetic parser (integers, + - * /, sqrt);
  * g(t) = z0 + sum_k z_k J0(d_k t) + sum_c y_c sum_i J0(r_ci t) is enclosed by Arb over whole t-intervals (ball
    input), bisecting until each enclosure is certainly >= 0, on [0, T] with the T of the certificate;
  * on [T, oo), |J0(x)| <= sqrt(2/(pi x)) for x > 0.
usage: python3 check3.py NAME.cert.json k
"""
import ast
import itertools
import json
import sys
from fractions import Fraction as Fr

import flint
from flint import arb

flint.ctx.prec = 120
S3 = arb(3).sqrt()


class Q3:
    """a + b sqrt(3) with a, b rational."""

    def __init__(self, a, b=0):
        self.a, self.b = Fr(a), Fr(b)

    def _c(self, o):
        return o if isinstance(o, Q3) else Q3(o)

    def __add__(self, o):
        o = self._c(o)
        return Q3(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q3(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-self._c(o))

    def __mul__(self, o):
        o = self._c(o)
        return Q3(self.a * o.a + 3 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        if o.b != 0:
            raise ValueError('division by an irrational number of Q(sqrt 3) is not needed here')
        return Q3(self.a / o.a, self.b / o.a)

    def __eq__(self, o):
        o = self._c(o)
        return self.a == o.a and self.b == o.b

    def arb(self):
        return A(self.a) + A(self.b) * S3


def A(q):
    q = Fr(q)
    return arb(q.numerator) / arb(q.denominator)


def rational_sqrt(q):
    """Square root of a nonnegative rational square, or None."""
    q = Fr(q)
    for v in (q.numerator, q.denominator):
        if v < 0 or round(v ** 0.5) ** 2 != v:
            return None
    return Fr(round(q.numerator ** 0.5), round(q.denominator ** 0.5))


def parse(expr, kind):
    """Value of an expression in integers, + - * /, unary minus and sqrt(.): an Arb ball (kind 'arb') or an exact
    element of Q(sqrt 3) (kind 'q3', where only sqrt(3) and square roots of rational squares are allowed)."""
    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return arb(n.value) if kind == 'arb' else Q3(n.value)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return -ev(n.operand)
        if isinstance(n, ast.BinOp) and type(n.op) in (ast.Add, ast.Sub, ast.Mult, ast.Div):
            a, b = ev(n.left), ev(n.right)
            return {ast.Add: lambda: a + b, ast.Sub: lambda: a - b, ast.Mult: lambda: a * b,
                    ast.Div: lambda: a / b}[type(n.op)]()
        if (isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == 'sqrt' and len(n.args) == 1
                and not n.keywords):
            x = ev(n.args[0])
            if kind == 'arb':
                return x.sqrt()
            if x == Q3(3):
                return Q3(0, 1)
            if x.b == 0 and rational_sqrt(x.a) is not None:
                return Q3(rational_sqrt(x.a))
            raise ValueError('square root outside Q(sqrt 3)')
        raise ValueError('unsupported expression: ' + ast.dump(n))
    return ev(ast.parse(expr, mode='eval'))


h = Fr(1, 2)
SHAPES_Q3 = {  # exact vertex coordinates in Q(sqrt 3)^2 (the same shapes as in shapes.py)
    'tri111': [(Q3(0), Q3(0)), (Q3(1), Q3(0)), (Q3(h), Q3(0, h))],
    'tri222': [(Q3(0), Q3(0)), (Q3(2), Q3(0)), (Q3(1), Q3(0, 1))],
    'rhomb1': [(Q3(0), Q3(0)), (Q3(1), Q3(0)), (Q3(h), Q3(0, h)), (Q3(h), Q3(0, -h))],
}


def main():
    C = json.load(open(sys.argv[1]))
    k = int(sys.argv[2])
    dexprs = C['D'].split(';')
    z0 = Fr(C['z0'])
    z = [Fr(v) for v in C['z']]
    if len(z) != len(dexprs):
        raise SystemExit('check3: one coefficient per distance expected')
    terms = [(z[i], parse(dexprs[i], 'arb')) for i in range(len(dexprs))]
    ysum = Fr(0)
    m = Fr(0)
    if C['cliques']:
        Dsq = [d * d for d in (parse(e, 'q3') for e in dexprs)]
    for c in C['cliques']:
        y = Fr(c['y'])
        if y < 0:
            raise SystemExit('check3: negative clique weight')
        W = SHAPES_Q3[c['shape']]
        for P, R in itertools.combinations(W, 2):
            d2 = (P[0] - R[0]) * (P[0] - R[0]) + (P[1] - R[1]) * (P[1] - R[1])
            if not any(d2 == dd for dd in Dsq):
                raise SystemExit(f"check3: not a clique: {c['shape']}")
        n = len(W)
        cx = sum((w[0] for w in W), Q3(0)) / n
        cy = sum((w[1] for w in W), Q3(0)) / n
        px, py = cx + Fr(c['ox']), cy + Fr(c['oy'])
        for wx, wy in W:
            r2 = (wx - px) * (wx - px) + (wy - py) * (wy - py)
            terms.append((y, r2.arb().sqrt()))
        ysum += y
        m += y * n
    Z0 = A(z0)
    coefs = [(A(cf), r) for cf, r in terms]
    T = Fr(C['T'])
    # [T, oo): |sum_i c_i J0(r_i t)| <= sum_i |c_i| sqrt(2/(pi r_i T)) < z0
    tail = sum((abs(cf) * (arb(2) / (arb.pi() * r * A(T))).sqrt() for cf, r in coefs), arb(0))
    if not Z0 - tail > 0:
        raise SystemExit('check3: tail bound fails')
    # [0, T]: enclose g over whole intervals, bisecting until the enclosure is certainly >= 0
    stack = [(Fr(i) * T / 64, Fr(i + 1) * T / 64) for i in range(64)]
    pieces = 0
    worst = None
    while stack:
        a, b = stack.pop()
        mid, rad = (a + b) / 2, (b - a) / 2
        t = A(mid) + arb(0, rad.numerator / rad.denominator * (1 + 1e-12) + 1e-40)
        g = Z0 + sum((cf * (r * t).bessel_j(0) for cf, r in coefs), arb(0))
        if g >= 0:
            pieces += 1
            lo = g.lower()
            worst = lo if worst is None or lo < worst else worst
            continue
        if b - a < Fr(1, 10**15):
            raise SystemExit(f'check3: cannot certify g >= 0 near t = {float(a)}')
        stack += [(a, mid), (mid, b)]
    s = z0 + sum(z) + m
    if not s > 0:
        raise SystemExit('check3: s <= 0')
    B = (z0 + ysum) / s
    print(f"check3: {pieces} intervals on [0, {C['T']}], min lower bound of g {float(worst.mid()):.3e}, "
          f"tail {float(tail.mid()):.3e} < z0 {float(z0):.6f}")
    print(f"check3: B = {float(B):.8f}, 1/B = {float(1 / B):.5f}, B < 1/{k}: {B < Fr(1, k)}")
    if not B < Fr(1, k):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
