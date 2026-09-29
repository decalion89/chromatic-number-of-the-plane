"""CNF encodings of 'G_13 has a proper 5-colouring', with optional symmetry breaking and redundant constraints.

Variables: x(v, c) = 1 + 5 v + c  (vertex v = 13 x + y of F_13^2, colour c = 0..4).  Auxiliary variables follow.
Clauses of the plain formula: at least one colour per vertex, at most one, and -x(a,c) | -x(b,c) per edge and colour.

Options used by the plan (soundness in plan_C.py's docstring; each option only restricts the colourings to a
normal form that every colouring can be brought to):
  --big0 s --cap s     class 0 has exactly s points (with --cap: every class has at most s)
  --dom0               class 0 is dominating (maximal independent)
  --domcap s           every class 1..4 with s points is dominating
  --lex0 L             class 0 is lex-leader under all 4 732 automorphisms, on the first L positions of lex_order()
  --lexcircles 9,11,.. circles first in lex_order() (default 0, N = 2..12, unit circle last)
  --vp 1234            value precedence of colours 1 < 2 < 3 < 4 along lex_order()
  --rigid34            case M* = 34: classes 0..3 have 34 points and are dominating, class 4 has 33;
                       value precedence 1 < 2 < 3
  --pin-edge           the baseline: x(0,0), x(e,1), value precedence
Exploratory only: --nf (lex-max normal form: sorted class sizes + Grundy), --maxmult, --units.
Plan D's case formulas:  F36 = --big0 36 --cap 36 --dom0 --lex0 25 --vp 1234 --domcap 36 (F35 likewise),
                         F34 = --rigid34 --lex0 25.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g13 import (Q, K, NV, POINTS, INDEX, UNITS, EDGES, ADJ, POW, GEN, add, sub, mul, conj, norm, circle,
                 stabiliser, automorphisms)


def X(v, c):
    return 1 + v * K + c


class CNF:
    def __init__(self):
        self.top = NV * K
        self.clauses = []

    def new(self, k=1):
        out = list(range(self.top + 1, self.top + k + 1))
        self.top += k
        return out

    def add(self, cl):
        self.clauses.append(list(cl))

    def text(self):
        return f"p cnf {self.top} {len(self.clauses)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in self.clauses)


def base(cnf, amo=True):
    for v in range(NV):
        cnf.add([X(v, c) for c in range(K)])
        if amo:
            for c in range(K):
                for d in range(c + 1, K):
                    cnf.add([-X(v, c), -X(v, d)])
    for a, b in EDGES:
        for c in range(K):
            cnf.add([-X(a, c), -X(b, c)])


def totalizer_atleast(cnf, lits, t):
    """at least t of lits (unary totalizer, outputs capped at t)"""
    if t <= 0:
        return
    if t > len(lits):
        cnf.add([])
        return
    def join(left, right):
        m = min(t, len(left) + len(right))
        out = cnf.new(m)
        for a in range(len(left) + 1):
            for b in range(len(right) + 1):
                j = a + b + 1
                if j > m:
                    continue
                cl = ([left[a]] if a < len(left) else []) + ([right[b]] if b < len(right) else [])
                cnf.add(cl + [-out[j - 1]])
        return out

    def build(items):
        if len(items) == 1:
            return items[0]
        mid = len(items) // 2
        return join(build(items[:mid]), build(items[mid:]))

    root = build([[l] for l in lits])
    cnf.add([root[t - 1]])


def unary_count(cnf, lits, t):
    """totalizer outputs o[0..t-1], o[j-1] <-> at least j of lits (both directions), capped at t"""
    def join(left, right):
        m = min(t, len(left) + len(right))
        out = cnf.new(m)
        for a in range(len(left) + 1):
            for b in range(len(right) + 1):
                j = a + b
                # upward: left >= a and right >= b  ->  out >= a + b
                if 1 <= j <= m:
                    cl = ([-left[a - 1]] if a > 0 else []) + ([-right[b - 1]] if b > 0 else [])
                    cnf.add(cl + [out[j - 1]])
                # downward: left < a+1 and right < b+1 -> out < a + b + 1
                j = a + b + 1
                if j <= m:
                    cl = ([left[a]] if a < len(left) else []) + ([right[b]] if b < len(right) else [])
                    cnf.add(cl + [-out[j - 1]])
        return out

    def build(items):
        if len(items) == 1:
            return items[0]
        mid = len(items) // 2
        return join(build(items[:mid]), build(items[mid:]))

    return build([[l] for l in lits])


def totalizer_atmost(cnf, lits, t):
    """at most t of lits: 'at least n - t of the negations'"""
    totalizer_atleast(cnf, [-l for l in lits], len(lits) - t)


def seq_atmost(cnf, lits, t):
    """at most t of lits, sequential counter (Sinz); s[i][j] = at least j+1 of lits[0..i]"""
    n = len(lits)
    if t >= n:
        return
    if t == 0:
        for l in lits:
            cnf.add([-l])
        return
    s = [cnf.new(t) for _ in range(n - 1)]
    cnf.add([-lits[0], s[0][0]])
    for j in range(1, t):
        cnf.add([-s[0][j]])
    for i in range(1, n - 1):
        cnf.add([-lits[i], s[i][0]])
        cnf.add([-s[i - 1][0], s[i][0]])
        for j in range(1, t):
            cnf.add([-lits[i], -s[i - 1][j - 1], s[i][j]])
            cnf.add([-s[i - 1][j], s[i][j]])
        cnf.add([-lits[i], -s[i - 1][t - 1]])
    cnf.add([-lits[n - 1], -s[n - 2][t - 1]])


def lex_order(unit_last=True, circles=None):
    """0, then the circles N = 2..12 (unit circle last), each as r, g r, ..., g^13 r; with circles, those first
    (in the given order), then the others ascending, the unit circle last"""
    out = [INDEX[(0, 0)]]
    cs = list(range(2, Q)) + [1] if unit_last else list(range(1, Q))
    if circles:
        cs = list(circles) + [c for c in cs if c not in circles]
    for c in cs:
        r = circle(c)[0]
        out += [INDEX[mul(p, r)] for p in POW]
    return out


def lex_leader_sets(cnf, lits, perms, order, L):
    """for each p: (lits[v] for v in order[:L]) >=_lex (lits[p[v]] for v in order[:L]); lits indexed by vertex.
    i.e. the set S = {v : lits[v]} is lex-greater or equal to p^-1(S)."""
    ordL = order[:L]
    for p in perms:
        prev = None
        pos = [v for v in ordL if p[v] != v]
        for k, v in enumerate(pos):
            x, y = lits[v], lits[p[v]]
            cnf.add(([-prev] if prev else []) + [x, -y])
            if k == len(pos) - 1:
                break
            e = cnf.new()[0]
            cnf.add(([-prev] if prev else []) + [x, e])
            cnf.add(([-prev] if prev else []) + [-y, e])
            prev = e


def value_precedence(cnf, colours, order):
    """colours = [c_1 < c_2 < ...]: colour c_{i+1} is used at a vertex only if c_i is used earlier in order.
    Encoded with y(i, c) = 'colour c used among order[:i+1]'."""
    for a, b in zip(colours, colours[1:]):
        # b at position i requires a at some position < i
        used = None
        for i, v in enumerate(order):
            if used is None:
                cnf.add([-X(v, b)])
            else:
                cnf.add([-X(v, b), used])
            nu = cnf.new()[0]
            # nu <-> used | x(v, a)
            if used is None:
                cnf.add([-nu, X(v, a)])
                cnf.add([nu, -X(v, a)])
            else:
                cnf.add([-nu, used, X(v, a)])
                cnf.add([nu, -used])
                cnf.add([nu, -X(v, a)])
            used = nu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--no-amo", action="store_true")
    ap.add_argument("--pin-edge", action="store_true", help="x(0,0), x(e,1), value precedence 2<3<4 (baseline)")
    ap.add_argument("--big0", type=int, default=0, help="class 0 has at least this many vertices")
    ap.add_argument("--cap", type=int, default=0, help="every class has at most this many vertices")
    ap.add_argument("--cap-others", type=int, default=0, help="classes 1..4 have at most this many vertices")
    ap.add_argument("--lex0", type=int, default=0, help="lex-leader of class 0 under Aut, prefix length")
    ap.add_argument("--lexcircles", default="", help="circles first in the lex order, e.g. 9,11,7,6")
    ap.add_argument("--lexstab", action="store_true", help="only the stabiliser of 0 (and pin 0 in class 0)")
    ap.add_argument("--dom0", action="store_true", help="class 0 dominating (maximal independent)")
    ap.add_argument("--vp", default="", help="value precedence on these colours, e.g. 1234")
    ap.add_argument("--nf", action="store_true",
                    help="lex-max normal form: |C0| >= ... >= |C4| and Grundy (v in C_j sees every C_i, i < j)")
    ap.add_argument("--domcap", type=int, default=0,
                    help="case M* = s: every class of colours 1..4 with s points is dominating")
    ap.add_argument("--rigid34", action="store_true",
                    help="case M* = 34: classes 0..3 have exactly 34 points and are dominating, class 4 has 33; "
                         "value precedence on 1,2,3")
    ap.add_argument("--maxmult", type=int, default=0, help="every vertex has at most m neighbours of each colour")
    ap.add_argument("--units", default="", help="extra unit literals, comma separated")
    a = ap.parse_args()
    cnf = CNF()
    base(cnf, amo=not a.no_amo)
    order = lex_order(circles=[int(c) for c in a.lexcircles.split(",") if c])
    if a.pin_edge:
        e = ADJ[0][0]
        cnf.add([X(0, 0)])
        cnf.add([X(e, 1)])
        vorder = [0, e] + [v for v in range(NV) if v not in (0, e)]
        value_precedence(cnf, [1, 2, 3, 4], vorder)
    if a.big0:
        totalizer_atleast(cnf, [X(v, 0) for v in range(NV)], a.big0)
    if a.cap:
        for c in range(K):
            totalizer_atmost(cnf, [X(v, c) for v in range(NV)], a.cap)
    if a.cap_others:
        for c in range(1, K):
            totalizer_atmost(cnf, [X(v, c) for v in range(NV)], a.cap_others)
    if a.lex0:
        s = {v: X(v, 0) for v in range(NV)}
        if a.lexstab:
            cnf.add([X(0, 0)])
            perms = stabiliser()[1:]
        else:
            perms = [p for p in automorphisms() if any(p[v] != v for v in range(NV))]
        lex_leader_sets(cnf, s, perms, order, a.lex0)
    if a.dom0:
        for v in range(NV):
            cnf.add([X(v, 0)] + [X(w, 0) for w in ADJ[v]])
    if a.vp:
        cols = [int(ch) for ch in a.vp]
        value_precedence(cnf, cols, order)
    if a.maxmult:
        for v in range(NV):
            for c in range(K):
                seq_atmost(cnf, [X(w, c) for w in ADJ[v]], a.maxmult)
    if a.domcap:
        for c in range(1, K):
            o = unary_count(cnf, [X(v, c) for v in range(NV)], a.domcap)
            for v in range(NV):
                cnf.add([-o[a.domcap - 1], X(v, c)] + [X(w, c) for w in ADJ[v]])
    if a.rigid34:
        for c in range(4):
            totalizer_atleast(cnf, [X(v, c) for v in range(NV)], 34)
            totalizer_atmost(cnf, [X(v, c) for v in range(NV)], 34)
            for v in range(NV):
                cnf.add([X(v, c)] + [X(w, c) for w in ADJ[v]])
        totalizer_atmost(cnf, [X(v, 4) for v in range(NV)], 33)
        value_precedence(cnf, [1, 2, 3], order)
    if a.nf:
        cnt = [unary_count(cnf, [X(v, c) for v in range(NV)], 43) for c in range(K)]
        for c in range(K - 1):
            for j in range(43):
                cnf.add([-cnt[c + 1][j], cnt[c][j]])      # |C_{c+1}| >= j+1  ->  |C_c| >= j+1
        cnf.add([cnt[0][33]])                              # |C0| >= 34
        for v in range(NV):
            for j in range(1, K):
                for i in range(j):
                    cnf.add([-X(v, j)] + [X(w, i) for w in ADJ[v]])
    for t in a.units.split(","):
        if t:
            cnf.add([int(t)])
    with open(a.out, "w") as fh:
        fh.write(cnf.text())
    print(f"{a.out}: {cnf.top} variables, {len(cnf.clauses)} clauses")


if __name__ == "__main__":
    main()
