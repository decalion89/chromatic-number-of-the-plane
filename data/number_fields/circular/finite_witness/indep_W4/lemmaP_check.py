#!/usr/bin/env python3
"""Numerical sanity check of Lemma P (a)-(d), Corollary P and the 'Use' paragraph on small examples.

Part 1 (brute force): random U (nonzero integer vectors in Z^2 or Z^3 with U and -U disjoint), a small
point set H, ALL proper 4-colourings with c(x0) = 0 enumerated by backtracking; those without a tight
cycle are kept; on random closed walks (random walk + BFS return path, so vertices get revisited) we test
 (a) per integer, (b) additivity and inverse, (c) backtrack insertion/deletion and square swaps
 (also the excluded cases t = s and t = -s), (d) the strict bounds.
 Negative controls: the same tests on colourings WITH a tight cycle, and on a U with an antipodal pair.
Part 2 (Use paragraph): U = {e_1..e_d} plus extra vectors, Z-basis r_j of the relation lattice, colourings
from characters x -> <a,x> mod 1 with all ||<a,u>|| > 1/4 and Markov-chain perturbations of them (kept only
if proper and tight-cycle free); chains 0 -> rho by +-r_j, every chain step bubble-sorted inside G;
checks per invariance at every move, per(W_rho) = sum_j a_j per(W_{r_j}), the ranges, and for character
colourings per(W_rho) = sum_u n_u f(u) with f(u) in (1/4, 3/4).
"""
import itertools
import random
import sys
from collections import deque
from fractions import Fraction as Fr

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 12345)
CNT = {}


def bump(k, v=1):
    CNT[k] = CNT.get(k, 0) + v


def vadd(x, s, k=1):
    return tuple(a + k * b for a, b in zip(x, s))


def vneg(s):
    return tuple(-a for a in s)


class Inst:
    def __init__(self, U, H):
        self.U = [tuple(u) for u in U]
        self.Uset = set(self.U)
        self.negU = set(vneg(u) for u in self.U)
        self.steps = list(self.Uset | self.negU)
        self.H = set(map(tuple, H))
        self.adj = {x: [vadd(x, s) for s in self.steps if vadd(x, s) in self.H] for x in self.H}

    def per(self, c, W):
        """W = vertex list x_0..x_L with x_L = x_0.  Returns (per as Fraction, P, N, Lambda)."""
        lam = 0
        P = N = 0
        for a, b in zip(W, W[1:]):
            d = (c[b] - c[a]) % 4
            assert d != 0
            lam += d
            s = tuple(q - p for p, q in zip(a, b))
            if s in self.negU:
                N += 1
            if s in self.Uset:
                P += 1
        return Fr(lam, 4) - N, P, N, lam


def tight_cycle_free(inst, c):
    """Kahn's algorithm on the digraph of tight arcs."""
    indeg = {x: 0 for x in inst.H}
    out = {x: [] for x in inst.H}
    for x in inst.H:
        for y in inst.adj[x]:
            if (c[y] - c[x]) % 4 == 1:
                out[x].append(y)
                indeg[y] += 1
    q = [x for x in inst.H if indeg[x] == 0]
    seen = 0
    while q:
        x = q.pop()
        seen += 1
        for y in out[x]:
            indeg[y] -= 1
            if indeg[y] == 0:
                q.append(y)
    return seen == len(inst.H)


def all_colourings(inst, x0, limit=None):
    order = sorted(inst.H)
    order.remove(x0)
    order = [x0] + order
    c = {}
    res = []

    def rec(i):
        if limit is not None and len(res) >= limit:
            return
        if i == len(order):
            res.append(dict(c))
            return
        x = order[i]
        for k in ([0] if i == 0 else range(4)):
            if all(c.get(y) != k for y in inst.adj[x]):
                c[x] = k
                rec(i + 1)
                del c[x]
    rec(0)
    return res


def bfs_path(inst, a, b):
    prev = {a: None}
    q = deque([a])
    while q:
        x = q.popleft()
        if x == b:
            break
        for y in inst.adj[x]:
            if y not in prev:
                prev[y] = x
                q.append(y)
    if b not in prev:
        return None
    p = [b]
    while prev[p[-1]] is not None:
        p.append(prev[p[-1]])
    return p[::-1]


def random_closed_walk(inst, x0, maxlen):
    W = [x0]
    for _ in range(rng.randint(1, maxlen)):
        nb = inst.adj[W[-1]]
        if not nb:
            break
        W.append(rng.choice(nb))
    back = bfs_path(inst, W[-1], x0)
    W = W + back[1:]
    if len(W) == 1:  # isolated start: no closed walk of positive length
        return None
    return W


def moves(inst, W):
    """All walks obtained from W by one move (c): backtrack insert/delete, square swap (incl. t = +-s)."""
    res = []
    L = len(W) - 1
    for i in range(L + 1):           # insert backtrack at position i
        x = W[i]
        for y in inst.adj[x]:
            res.append(("ins", W[:i + 1] + [y, x] + W[i + 1:]))
    for i in range(L - 1):           # delete backtrack W[i], W[i+1], W[i+2] = W[i]
        if W[i] == W[i + 2]:
            res.append(("del", W[:i + 1] + W[i + 3:]))
    for i in range(L - 1):           # swap the steps s = W[i+1]-W[i], t = W[i+2]-W[i+1]
        x, y, z = W[i], W[i + 1], W[i + 2]
        s = tuple(q - p for p, q in zip(x, y))
        t = tuple(q - p for p, q in zip(y, z))
        w = vadd(x, t)
        if w in inst.H:
            kind = "swap" if (s != t and s != vneg(t)) else "swap_degenerate"
            res.append((kind, W[:i + 1] + [w] + W[i + 2:]))
    return res


def test_walks(inst, c, x0, nwalks, tag, expect_ok=True):
    viol = {"a": 0, "b": 0, "c": 0, "c_deg": 0, "d": 0}
    for _ in range(nwalks):
        W = random_closed_walk(inst, x0, 8)
        if W is None:
            continue
        p, P, N, lam = inst.per(c, W)
        L = len(W) - 1
        bump(tag + " walks")
        if p.denominator != 1:
            viol["a"] += 1
        if not (Fr(P - 3 * N, 4) < p < Fr(3 * P - N, 4)) or P + N != L:
            viol["d"] += 1
        W2 = random_closed_walk(inst, x0, 8)
        if W2 is not None:
            p2 = inst.per(c, W2)[0]
            pc = inst.per(c, W + W2[1:])[0]
            pi = inst.per(c, W[::-1])[0]
            if pc != p + p2 or pi != -p:
                viol["b"] += 1
        for kind, W3 in moves(inst, W):
            q = inst.per(c, W3)[0]
            bump(tag + " moves " + kind)
            if q != p:
                if kind == "swap_degenerate":
                    viol["c_deg"] += 1
                else:
                    viol["c"] += 1
    for k, v in viol.items():
        bump(tag + " violations " + k, v)
    return viol


def part1(ninst):
    for it in range(ninst):
        d = rng.choice([2, 3])
        while True:
            k = rng.randint(2, 4)
            U = set()
            while len(U) < k:
                u = tuple(rng.randint(-1, 1) for _ in range(d))
                if any(u) and vneg(u) not in U:
                    U.add(u)
            box = list(itertools.product(*[range(3)] * d))
            H = rng.sample(box, rng.randint(6, 9))
            inst = Inst(sorted(U), H)
            if sum(len(v) for v in inst.adj.values()) >= 2 * len(H):
                break
        x0 = rng.choice(sorted(inst.H))
        cols = all_colourings(inst, x0)
        good = [c for c in cols if tight_cycle_free(inst, c)]
        bad = [c for c in cols if not tight_cycle_free(inst, c)]
        bump("P1 instances")
        bump("P1 proper colourings (c(x0)=0)", len(cols))
        bump("P1 tight-cycle-free colourings", len(good))
        for c in rng.sample(good, min(len(good), 6)):
            test_walks(inst, c, x0, 25, "P1[no tight cycle]")
        for c in rng.sample(bad, min(len(bad), 3)):
            test_walks(inst, c, x0, 25, "P1-control[has tight cycle]")
    # control: U with an antipodal pair
    U = [(1, 0), (-1, 0), (0, 1)]
    inst = Inst(U, list(itertools.product(range(3), range(3))))
    for c in all_colourings(inst, (0, 0), limit=2000):
        if tight_cycle_free(inst, c):
            test_walks(inst, c, (0, 0), 20, "P1-control[U meets -U]")
            break


# ---------------------------------------------------------------------------
def frac(x):
    return x - (x.numerator // x.denominator)


def part2(ninst):
    for it in range(ninst):
        d = rng.choice([2, 3])
        E = [tuple(1 if t == i else 0 for t in range(d)) for i in range(d)]
        extra = []
        while len(extra) < rng.randint(2, 3):
            u = tuple(rng.randint(-2, 2) for _ in range(d))
            Uall = E + extra
            if sum(1 for a in u if a) >= 2 and u not in Uall and vneg(u) not in Uall:
                extra.append(u)
        U = E + extra
        m = len(U)
        # Z-basis of the relation lattice: r_i = e_{d+i} - sum_t (u_{d+i})_t e_t  (coefficient space Z^m)
        basis = []
        for i, u in enumerate(extra):
            r = [0] * m
            r[d + i] = 1
            for t in range(d):
                r[t] = -u[t]
            basis.append(tuple(r))
        # a character with all ||<a,u>|| > 1/4 (rational a, exact arithmetic)
        a = None
        for _ in range(20000):
            cand = tuple(Fr(rng.randint(0, 59), 60) for _ in range(d))
            vals = [frac(sum(x * y for x, y in zip(cand, u))) for u in U]
            if all(Fr(1, 4) < v < Fr(3, 4) for v in vals):
                a = cand
                break
        if a is None:
            bump("P2 instances skipped (no character found)")
            continue
        fU = [frac(sum(x * y for x, y in zip(a, u))) for u in U]

        def walk_steps(rho, sign=1):
            st = []
            for j, n in enumerate(rho):
                st += [(j, 1 if n > 0 else -1)] * abs(n)
            if sign < 0:
                st = [(j, -s) for (j, s) in reversed(st)]
            return st

        def points(st):
            p = (0,) * d
            pts = [p]
            for j, s in st:
                p = vadd(p, U[j], s)
                pts.append(p)
            return pts

        # chains 0 -> rho for random rho in the span of the basis
        S = []
        for _ in range(8):
            coeffs = [rng.randint(-2, 2) for _ in basis]
            if any(coeffs):
                S.append(coeffs)
        triples = []      # (rho, j, s, rho') with W_rho W_{r_j}^s bubble-sorted into W_rho'
        Hset = {(0,) * d}
        for coeffs in S:
            cur = [0] * len(basis)
            seq = []
            for j, aj in enumerate(coeffs):
                seq += [(j, 1 if aj > 0 else -1)] * abs(aj)
            rng.shuffle(seq)
            for j, s in seq:
                rho = tuple(sum(cur[i] * basis[i][u] for i in range(len(basis))) for u in range(m))
                cur[j] += s
                rho2 = tuple(sum(cur[i] * basis[i][u] for i in range(len(basis))) for u in range(m))
                triples.append((rho, j, s, rho2))
        # build H: all vertices met, plus the fourth corners of the bubble-sort swaps; record walks
        sorted_walks = []
        for (rho, j, s, rho2) in triples:
            st = walk_steps(rho) + walk_steps(basis[j], s)
            pts = points(st)
            Hset.update(pts)
            seqs = [list(pts)]
            changed = True
            while changed:
                changed = False
                i = 0
                while i < len(st) - 1:
                    x, y = st[i], st[i + 1]
                    if x[0] == y[0] and x[1] == -y[1]:
                        del st[i:i + 2]
                        del pts[i + 1:i + 3]
                        seqs.append(list(pts))
                        changed = True
                        i = max(i - 1, 0)
                        continue
                    if x[0] > y[0]:
                        new = vadd(pts[i], U[y[0]], y[1])
                        Hset.add(new)
                        st[i], st[i + 1] = y, x
                        pts[i + 1] = new
                        seqs.append(list(pts))
                        changed = True
                    i += 1
            assert st == walk_steps(rho2)
            sorted_walks.append(seqs)
            for r in (rho, rho2, basis[j]):
                Hset.update(points(walk_steps(r)))
        inst = Inst(U, sorted(Hset))
        # colourings: the character colouring and Markov-chain perturbations of it
        c0 = {x: int(4 * frac(sum(p * q for p, q in zip(a, x)))) for x in inst.H}
        assert all(c0[x] != c0[y] for x in inst.H for y in inst.adj[x])
        assert tight_cycle_free(inst, c0)
        cols = [("char", c0)]
        c = dict(c0)
        Hl = sorted(inst.H)
        accepted = 0
        for _ in range(1500):
            x = rng.choice(Hl)
            k = rng.randrange(4)
            if k == c[x] or any(c[y] == k for y in inst.adj[x]):
                continue
            old = c[x]
            c[x] = k
            if tight_cycle_free(inst, c):
                accepted += 1
            else:
                c[x] = old
        cols.append(("mcmc", dict(c)))
        bump("P2 instances")
        bump("P2 |H| total", len(inst.H))
        bump("P2 mcmc accepted recolourings", accepted)
        for kind, c in cols:
            per = lambda rho: inst.per(c, points(walk_steps(rho)))[0]
            p = [per(r) for r in basis]
            # every intermediate walk of every bubble sort has the same per as the concatenation
            for seqs in sorted_walks:
                vals = set(inst.per(c, W)[0] for W in seqs)
                bump("P2 bubble-sort walks checked", len(seqs))
                if len(vals) != 1:
                    bump("P2 VIOLATION per changed during a bubble sort")
            for (rho, j, s, rho2) in triples:
                bump("P2 triples")
                if per(rho2) != per(rho) + s * p[j]:
                    bump("P2 VIOLATION triple")
            for coeffs in S:
                rho = tuple(sum(coeffs[i] * basis[i][u] for i in range(len(basis))) for u in range(m))
                y = per(rho)
                P = sum(n for n in rho if n > 0)
                N = -sum(n for n in rho if n < 0)
                bump("P2 relations rho in S")
                if y != sum(ci * pi for ci, pi in zip(coeffs, p)):
                    bump("P2 VIOLATION homomorphism")
                if not (Fr(P - 3 * N, 4) < y < Fr(3 * P - N, 4)):
                    bump("P2 VIOLATION range")
                if kind == "char":
                    chi = sum(n * f for n, f in zip(rho, fU))
                    if chi != y or chi.denominator != 1:
                        bump("P2 VIOLATION character value")
                    else:
                        bump("P2 per(W_rho) == sum n_u f(u) (character colouring)")
        # real-character claim on many random nonzero relations (no walks needed)
        for _ in range(200):
            coeffs = [rng.randint(-3, 3) for _ in basis]
            if not any(coeffs):
                continue
            rho = [sum(coeffs[i] * basis[i][u] for i in range(len(basis))) for u in range(m)]
            P = sum(n for n in rho if n > 0)
            N = -sum(n for n in rho if n < 0)
            v = sum(n * f for n, f in zip(rho, fU))
            bump("P2 character claim tested")
            if v.denominator != 1 or not (Fr(P - 3 * N, 4) < v < Fr(3 * P - N, 4)):
                bump("P2 VIOLATION character claim")


if __name__ == "__main__":
    part1(int(sys.argv[2]) if len(sys.argv) > 2 else 40)
    part2(int(sys.argv[3]) if len(sys.argv) > 3 else 12)
    for k in sorted(CNT):
        print("%-60s %d" % (k, CNT[k]))
