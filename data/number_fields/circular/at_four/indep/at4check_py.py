#!/usr/bin/env python3
"""Independent check (no code shared with at4check.c) of Lemma F13, Corollary F14 and Corollary F13'
of the draft q3_at4.md on small finite abelian groups, with exact arithmetic (int, Fraction).

For EVERY symmetric generating S of Gamma = Z/m_1 x ... x Z/m_d (no reduction by automorphisms) with
Cay(Gamma,S) 4-colourable, all proper 4-colourings with c(0) = 0 are enumerated, and per S we print
  kappa, #colourings, #with acyclic tight digraph, #without tight 4-cycle, #without tight square,
  #with A = B everywhere,
and assert:
  F13 : A != B at (x,s,t)  =>  the square at (x,s,t) is a tight cycle (one orientation);
        no tight square => A = B everywhere, a(-s) = 4 - a(s), s -> a(s)/4 is a character on S,
        and the windings of the closed walks (ord(s) x s) and (s, t, -(s+t)) do not depend on x;
  F14 : kappa <= 1/4 => no colouring has acyclic tight digraph; kappa > 1/4 => some has;
        an acyclic colouring has 1 < a(s) < 3 for all s;
  F13': some colouring has no tight square  <=>  kappa >= 1/4;
  hom : Cay -> K_{P/Q} (P = 4Q-1, Q = floor((n+1)/4)) <=> kappa > 1/4.
Tight 4-cycles and acyclicity are computed from powers of the 0/1 adjacency matrix of the tight digraph.
"""
import sys
import itertools
from fractions import Fraction
import numpy as np


def group(ms):
    els = list(itertools.product(*[range(m) for m in ms]))
    add = lambda a, b: tuple((x + y) % m for x, y, m in zip(a, b, ms))
    neg = lambda a: tuple((-x) % m for x, m in zip(a, ms))
    return els, add, neg


def frac_part(q):
    return q - (q.numerator // q.denominator)


def main(ms, maxprint=None):
    els, add, neg = group(ms)
    n = len(els)
    idx = {e: i for i, e in enumerate(els)}
    zero = tuple(0 for _ in ms)
    nonzero = [e for e in els if e != zero]
    classes = []
    for e in nonzero:
        c = frozenset([e, neg(e)])
        if c not in classes:
            classes.append(c)
    chars = els  # character j: x -> sum j_i x_i / m_i
    def chi(j, x):
        return frac_part(sum(Fraction(a * b, m) for a, b, m in zip(j, x, ms)))
    def norm(q):
        q = frac_part(q)
        return min(q, 1 - q)
    results = {}
    nbad = 0
    for r in range(1, len(classes) + 1):
        for cl in itertools.combinations(classes, r):
            S = sorted(set().union(*cl))
            # generation
            seen = {zero}
            frontier = [zero]
            while frontier:
                x = frontier.pop()
                for s in S:
                    y = add(x, s)
                    if y not in seen:
                        seen.add(y)
                        frontier.append(y)
            if len(seen) != n:
                continue
            kappa = max(min(norm(chi(j, s)) for s in S) for j in chars)
            # colourings in lexicographic element order (c(0)=0)
            order = els
            nbrs = {x: [add(x, s) for s in S] for x in els}
            cols = []
            col = {}
            def rec(i):
                if i == n:
                    cols.append(dict(col))
                    return
                x = order[i]
                for c in (range(1) if i == 0 else range(4)):
                    if all(col.get(y) != c for y in nbrs[x]):
                        col[x] = c
                        rec(i + 1)
                        del col[x]
            rec(0)
            if not cols:
                results[tuple(S)] = (kappa, 0)
                continue
            st = dict(col=0, acyc=0, noT4=0, noSq=0, ABeq=0)
            bad = []
            for c in cols:
                st['col'] += 1
                ell = {(x, s): (c[add(x, s)] - c[x]) % 4 for x in els for s in S}
                Amat = np.zeros((n, n), dtype=np.int64)
                for (x, s), v in ell.items():
                    if v == 1:
                        Amat[idx[x], idx[add(x, s)]] = 1
                A4 = np.linalg.matrix_power(Amat, 4)
                has_t4 = int(np.trace(A4)) > 0
                # acyclic iff nilpotent: boolean powers
                P = Amat.copy()
                for _ in range(n):
                    P = np.minimum(P @ Amat, 1)
                acyclic = not P.any()
                # tight squares and A = B
                has_sq = False
                abneq = False
                for x in els:
                    for s in S:
                        for t in S:
                            xs, xt = add(x, s), add(x, t)
                            xst = add(xs, t)
                            if ell[(x, s)] == 1 and ell[(xs, t)] == 1 and ell[(xst, neg(s))] == 1 and ell[(xt, neg(t))] == 1:
                                has_sq = True
                            Av = ell[(x, s)] + ell[(xs, t)]
                            Bv = ell[(x, t)] + ell[(xt, s)]
                            if Av != Bv:
                                abneq = True
                                if (Av, Bv) == (2, 6):
                                    ok = ell[(x, s)] == 1 and ell[(xs, t)] == 1 and ell[(xst, neg(s))] == 1 and ell[(xt, neg(t))] == 1
                                elif (Av, Bv) == (6, 2):
                                    ok = ell[(x, t)] == 1 and ell[(xt, s)] == 1 and ell[(xst, neg(t))] == 1 and ell[(xs, neg(s))] == 1
                                else:
                                    ok = False
                                if not ok:
                                    bad.append(('F13 square', c, x, s, t))
                if not has_t4:
                    st['noT4'] += 1
                if not has_sq:
                    st['noSq'] += 1
                if not abneq:
                    st['ABeq'] += 1
                if acyclic:
                    st['acyc'] += 1
                if has_sq and not has_t4:
                    bad.append(('square but no 4-cycle', c))
                a = {s: Fraction(sum(ell[(x, s)] for x in els), n) for s in S}
                if not has_sq:
                    if abneq:
                        bad.append(('no tight square but A != B', c))
                    if any(a[s] + a[neg(s)] != 4 for s in S):
                        bad.append(('a(-s) != 4 - a(s)', c))
                    if not any(all(frac_part(a[s] / 4) == chi(j, s) for s in S) for j in chars):
                        bad.append(('not a character', c))
                    # windings
                    for s in S:
                        o = 1
                        y = s
                        while y != zero:
                            y = add(y, s)
                            o += 1
                        vals = set()
                        for x in els:
                            tot, y = 0, x
                            for _ in range(o):
                                tot += ell[(y, s)]
                                y = add(y, s)
                            vals.add(tot)
                        if len(vals) != 1 or next(iter(vals)) % 4 or next(iter(vals)) != o * a[s]:
                            bad.append(('winding n.s', c, s))
                    for s in S:
                        for t in S:
                            u = add(s, t)
                            if u == zero or u not in S:
                                continue
                            vals = set()
                            for x in els:
                                y1 = add(x, s)
                                y2 = add(y1, t)
                                vals.add(ell[(x, s)] + ell[(y1, t)] + ell[(y2, neg(u))])
                            if len(vals) != 1 or next(iter(vals)) % 4 or next(iter(vals)) != a[s] + a[t] + a[neg(u)]:
                                bad.append(('winding s+t-u', c, s, t))
                if acyclic and not all(1 < a[s] < 3 for s in S):
                    bad.append(('acyclic but a(s) in {1,3}', c))
            if kappa <= Fraction(1, 4) and st['acyc']:
                bad.append(('(ii) fails',))
            if kappa > Fraction(1, 4) and not st['acyc']:
                bad.append(('(iii) fails',))
            if (st['noSq'] > 0) != (kappa >= Fraction(1, 4)):
                bad.append(("F13' fails",))
            # independent homomorphism test to K_{P/Q}
            Q = max(1, (n + 1) // 4)
            P = 4 * Q - 1
            hc = {}
            def hrec(i):
                if i == n:
                    return True
                x = order[i]
                for cc in (range(1) if i == 0 else range(P)):
                    if all(((cc - hc[y]) % P) in range(Q, P - Q + 1) for y in nbrs[x] if y in hc):
                        hc[x] = cc
                        if hrec(i + 1):
                            return True
                        del hc[x]
                return False
            hom = hrec(0)
            if hom != (kappa > Fraction(1, 4)):
                bad.append(('hom test', hom))
            if bad:
                nbad += 1
            results[tuple(S)] = (kappa, st['col'], st['acyc'], st['noT4'], st['noSq'], st['ABeq'], hom, len(bad))
            print('S=' + ''.join('(' + ','.join(map(str, s)) + ')' for s in S),
                  f"kappa={kappa} col={st['col']} acyc={st['acyc']} noT4={st['noT4']} noSq={st['noSq']} "
                  f"ABeq={st['ABeq']} hom{P}/{Q}={int(hom)} {'BAD ' + repr(bad[:3]) if bad else 'ok'}", flush=True)
    print(f'# summary Z/{ms}: generating S={len(results)}, 4-colourable={sum(1 for v in results.values() if v[1])}, BAD={nbad}')


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]])
