"""Exhaustive checker for a small finite witness of chi_c = p/q over F = Q(sqrt d) (no SAT solver, no proof file).

usage: python3 check_small.py WITNESS.json.gz [OUT.cnf]

The witness file contains
  d, denominator D, p, q     the field Q(sqrt d) (d > 1 squarefree), the common denominator, the value p/q
  points                     integer 4-tuples [a, b, c, e]: the point ((a + b r)/D, (c + e r)/D) of F^2, r = sqrt d
  edges                      pairs [i, j]
  colouring                  a (p, q)-colouring V -> Z/p
  cycles                     lists [v_0, ..., v_{m-1}] of vertices (directed cycles)
  critical_colourings        for every vertex v, a map V -> Z/p with v marked -1 (optional)
A (p, q)-colouring is a map c with q <= (c(y) - c(x)) mod p <= p - q on every edge xy, that is, a homomorphism to
K_{p/q}; the arc x -> y is tight when c(y) - c(x) = q mod p.  The checker verifies, with Python integers only:
  (1) the points are distinct, every edge joins two points at distance exactly 1, and the edges are ALL the pairs at
      distance 1: (Da)^2 + d (Db)^2 + (Dc)^2 + d (De)^2 = D^2 and Da*Db + Dc*De = 0 for the difference;
  (2) the colouring is a (p, q)-colouring, so chi_c(H) <= p/q;
  (3) it enumerates all p^n maps V -> Z/p and finds that every (p, q)-colouring has a tight cycle (a directed cycle
      of its tight digraph, found by depth-first search), indeed one of the listed cycles; by Lemma 20 of
      papers/three-colours (Guichard), chi_c(H) >= p/q;
  (4) independently of (3): H has no homomorphism to K_{p'/q'}, where p'/q' is the largest fraction below p/q with
      p' <= n (backtracking search); as chi_c(H) is a fraction with numerator at most n (Vince, Bondy-Hell), this
      gives chi_c(H) >= p/q again;
  (5) for every vertex v, critical_colourings[v] is a (p, q)-colouring of H - v whose tight digraph has no directed
      cycle, so chi_c(H - v) < p/q (perturb the colours along a topological order): H is vertex-critical.
With OUT.cnf it also writes the formula "a (p, q)-colouring in which every listed cycle has a non-tight arc"
(variables x(v,k) = p v + k + 1, one variable per arc of a listed cycle, the layout of check_witness.py) for a SAT
solver and a proof checker.  Every check is explicit (not an assert, which python -O would skip).
"""
import sys, json, gzip, itertools
from fractions import Fraction


def _req(ok, *msg):
    """An explicit check (not assert, so that python -O cannot skip it)."""
    if not ok:
        print('REJECTED:', *msg, file=sys.stderr)
        sys.exit(1)


def load(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        return json.load(f)


def unit(P1, P2, D, d):
    a, b, c, e = (P2[k] - P1[k] for k in range(4))
    return a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0


def proper(c, E, p, q):
    return all(q <= (c[j] - c[i]) % p <= p - q for i, j in E if c[i] >= 0 and c[j] >= 0)


def tight_cycle(c, n, E, p, q):
    """A directed cycle of the tight digraph of c (vertices with c = -1 ignored), or None."""
    out = {v: [] for v in range(n)}
    for i, j in E:
        if c[i] < 0 or c[j] < 0:
            continue
        if (c[j] - c[i]) % p == q:
            out[i].append(j)
        if (c[i] - c[j]) % p == q:
            out[j].append(i)
    state = [0] * n                       # 0 new, 1 on the stack, 2 done
    for s in range(n):
        if state[s]:
            continue
        stack = [(s, iter(out[s]))]; path = [s]; state[s] = 1
        while stack:
            v, it = stack[-1]
            w = next(it, None)
            if w is None:
                stack.pop(); path.pop(); state[v] = 2
            elif state[w] == 1:
                return path[path.index(w):]
            elif state[w] == 0:
                state[w] = 1; stack.append((w, iter(out[w]))); path.append(w)
    return None


def homomorphism(n, adj, p, q):
    """A homomorphism to K_{p/q} found by backtracking (colour 0 for vertex 0), or None."""
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    c = [-1] * n
    def ok(v, k):
        return all(c[w] < 0 or q <= (k - c[w]) % p <= p - q for w in adj[v])
    def go(t):
        if t == n:
            return True
        v = order[t]
        for k in ([0] if t == 0 else range(p)):
            if ok(v, k):
                c[v] = k
                if go(t + 1):
                    return True
                c[v] = -1
        return False
    return list(c) if go(0) else None


def cnf_text(n, E, cycles, p, q):
    x = lambda v, k: p * v + k + 1
    cl = [[x(v, k) for k in range(p)] for v in range(n)]
    bad = [t for t in range(p) if not q <= t <= p - q]
    for i, j in E:
        for k in range(p):
            for t in bad:
                cl.append([-x(i, k), -x(j, (k + t) % p)])
    tv = {}; nv = p * n
    for cyc in cycles:
        for k in range(len(cyc)):
            a, b = cyc[k], cyc[(k + 1) % len(cyc)]
            if (a, b) not in tv:
                nv += 1; tv[(a, b)] = nv
                for s in range(p):
                    cl.append([-x(a, s), -x(b, (s + q) % p), nv])
    for cyc in cycles:
        cl.append([-tv[(cyc[k], cyc[(k + 1) % len(cyc)])] for k in range(len(cyc))])
    return '\n'.join([f'p cnf {nv} {len(cl)}'] + [' '.join(map(str, c)) + ' 0' for c in cl]) + '\n'


def main():
    W = load(sys.argv[1])
    d, D, p, q = W['d'], W['denominator'], W['p'], W['q']
    P, E, col, cycles = W['points'], [tuple(e) for e in W['edges']], W['colouring'], W['cycles']
    n = len(P)
    _req(isinstance(d, int) and d > 1 and all(d % (k * k) for k in range(2, int(d ** 0.5) + 1)), 'd not squarefree')
    _req(0 < 2 * q < p and p ** n <= 10 ** 7, 'bad (p, q) or too many maps to enumerate')
    _req(all(len(t) == 4 and all(isinstance(z, int) for z in t) for t in P) and len(set(map(tuple, P))) == n,
         'bad or repeated points')
    Es = {(min(i, j), max(i, j)) for i, j in E}
    _req(len(Es) == len(E) and all(0 <= i < n and 0 <= j < n and i != j for i, j in E), 'bad or repeated edges')
    _req(all(unit(P[i], P[j], D, d) for i, j in E), 'an edge is not at distance 1')
    units = {(i, j) for i in range(n) for j in range(i + 1, n) if unit(P[i], P[j], D, d)}
    _req(units == Es, 'the edges are not all the unit pairs')
    print(f'(1) Q(sqrt{d}): {n} points, {len(E)} edges, every edge at distance exactly 1, and no other unit pair '
          f'(induced)')
    _req(len(col) == n and all(isinstance(k, int) and 0 <= k < p for k in col) and proper(col, E, p, q),
         'the colouring is not a (p, q)-colouring')
    print(f'(2) the colouring is a ({p},{q})-colouring: chi_c(H) <= {p}/{q}')
    for cyc in cycles:
        _req(len(cyc) >= 3 and len(set(cyc)) == len(cyc), 'a cycle is not simple', cyc)
        _req(all((min(a, b), max(a, b)) in Es for a, b in zip(cyc, cyc[1:] + cyc[:1])), 'a cycle uses a non-edge', cyc)
    count = 0
    for c in itertools.product(range(p), repeat=n):
        if not proper(c, E, p, q):
            continue
        count += 1
        _req(tight_cycle(c, n, E, p, q) is not None, 'a colouring without a tight cycle', c)
        _req(any(all((c[b] - c[a]) % p == q for a, b in zip(cyc, cyc[1:] + cyc[:1])) for cyc in cycles),
             'no listed cycle is tight', c)
    print(f'(3) all {p ** n} maps V -> Z/{p}: {count} ({p},{q})-colourings, each with a tight cycle (one of the '
          f'{len(cycles)} listed cycles): chi_c(H) = {p}/{q} by Lemma 20')
    lower = max(Fraction(a, b) for a in range(2, n + 1) for b in range(1, a) if Fraction(a, b) < Fraction(p, q))
    adj = {v: [] for v in range(n)}
    for i, j in E:
        adj[i].append(j); adj[j].append(i)
    _req(homomorphism(n, adj, lower.numerator, lower.denominator) is None, 'a homomorphism to K', lower)
    _req(homomorphism(n, adj, p, q) is not None, 'no homomorphism to K', p, q)
    print(f'(4) no homomorphism to K_{lower.numerator}/{lower.denominator}, the largest fraction below {p}/{q} with '
          f'numerator <= {n}: chi_c(H) = {p}/{q} again')
    if 'critical_colourings' in W:
        crit = W['critical_colourings']
        _req(len(crit) == n, 'one certificate per vertex is needed')
        for v in range(n):
            c = crit[v]
            _req(len(c) == n and c[v] == -1 and
                 all(isinstance(c[u], int) and 0 <= c[u] < p for u in range(n) if u != v),
                 'bad certificate', v)
            _req(proper(c, E, p, q), 'not a colouring of H - v', v)
            _req(tight_cycle(c, n, E, p, q) is None, 'tight cycle in H - v', v)
        print(f'(5) for every vertex v, a ({p},{q})-colouring of H - v with an acyclic tight digraph: H is '
              f'vertex-critical')
    if len(sys.argv) > 2:
        open(sys.argv[2], 'w').write(cnf_text(n, E, cycles, p, q))
        print('wrote', sys.argv[2])


if __name__ == '__main__':
    main()
