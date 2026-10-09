"""Stand-alone checker for the finite witnesses of chi_c = 4 over F = Q(sqrt3, sqrt11) and F = Q(sqrt2, sqrt3).

usage: python3 check_witness4.py WITNESS.json.gz [CNF] [PROOF.drat] [path/to/drat-trim]

The witness file contains
  field        "Q(sqrt3, sqrt11)" or "Q(sqrt2, sqrt3)", that is F = Q(sa, sb) with (a, b) = (3, 11) or (2, 3)
  basis        "(1, sqrt3, sqrt11, sqrt33)" or "(1, sqrt2, sqrt3, sqrt6)", that is (1, sa, sb, sab), sab = sqrt(ab)
  denominator  D (an integer)
  points       integer 8-tuples [a0, a1, a2, a3, b0, b1, b2, b3]: the point
               ((a0 + a1 sa + a2 sb + a3 sab)/D, (b0 + b1 sa + b2 sb + b3 sab)/D) of F^2
               (for Q(sqrt3, sqrt11): sa = sqrt3, sb = sqrt11, sab = sqrt33)
  edges        pairs [i, j]
  colouring    a map V -> Z/4
  cycles       lists [v_0, ..., v_{m-1}] of vertices, m divisible by 4
  fixed_vertex a vertex v_0 whose colour the formula fixes to 0
  generators   (optional) integer 8-tuples over D: a set U of unit vectors, one of each pair +-u
The checker verifies, with Python integers only:
  (1) the points are distinct and every edge joins two points at Euclidean distance exactly 1: for the difference
      (x0, ..., x3, y0, ..., y3), (x0 + x1 sa + x2 sb + x3 sab)^2 + (y0 + ...)^2 = D^2, computed in the basis
      (1, sa, sb, sab) with sa sb = sab, sa sab = a sb, sb sab = b sa (four integer equations; the basis is
      independent over Q, as a, b and ab are not rational squares and F has degree 4); it also reports whether
      `edges` is the set of ALL pairs at distance 1 (the induced unit-distance graph);
  (1b) if `generators` is given: every generator is a unit vector, no two are equal or opposite, and `edges` is
      exactly the set of pairs of points that differ by an element of U u -U, so that H is the subgraph of the
      Cayley graph Cay(F^2, U u -U) induced on the points; it also reports whether every point is joined to the
      fixed vertex by such steps inside H (then, if the fixed vertex is the origin, H lies in Cay(ZU, U u -U));
  (2) the colouring is a proper 4-colouring (so chi(H) <= 4 and chi_c(H) <= 4);
  (3) every listed cycle is a closed walk v_0 -> v_1 -> ... -> v_{m-1} -> v_0 along edges of H, with distinct
      vertices and m divisible by 4;
  (4) the CNF (if given) is exactly the formula built here from (edges, cycles, fixed_vertex):
        variables x(v,k) = 4v + k + 1 ("v has colour k")
        clauses  OR_k x(v,k)                          (every vertex gets a colour)
                 -x(v,k) v -x(v,l)  for k < l         (at most one colour)
                 -x(i,k) v -x(j,k)  for every edge ij  (adjacent vertices get different colours)
                 x(v_0,0)                             (the fixed vertex)
                 OR_t -x(v_t, k + t mod 4)            for every listed cycle and every k in Z/4
      the last clauses say that the cycle is not tight: a tight cycle is one with c(v_{t+1}) - c(v_t) = 1 mod 4
      for every t (indices mod m), that is c(v_t) = c(v_0) + t;
  (5) drat-trim (if given) verifies the DRAT proof of unsatisfiability of the CNF.
Logic: a proper 4-colouring c without a tight listed cycle would satisfy the CNF (x(v,k) iff c(v) - c(v_0) = k):
rotating the colours keeps a proper colouring proper and keeps every colour difference, hence every tight arc.  So
(5) shows that every 4-colouring of H has a tight cycle, and Lemma 20 of papers/three-colours (Guichard) gives
chi_c(H) >= 4; with (2), chi_c(H) = chi(H) = 4.  Every check is explicit (not an assert, which python -O would skip),
and (5) needs a line 's VERIFIED' from drat-trim.
"""
import sys, json, gzip, hashlib, os, subprocess, tempfile


def _req(ok, *msg):
    """An explicit check (not assert, so that python -O cannot skip it)."""
    if not ok:
        print('REJECTED:', *msg, file=sys.stderr)
        sys.exit(1)


def load(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        return json.load(f)


FIELDS = {'Q(sqrt3, sqrt11)': (3, 11, '(1, sqrt3, sqrt11, sqrt33)'),
          'Q(sqrt2, sqrt3)': (2, 3, '(1, sqrt2, sqrt3, sqrt6)')}


def square(x0, x1, x2, x3, a=3, b=11):
    """(x0 + x1 sa + x2 sb + x3 sab)^2 in the basis (1, sa, sb, sab), sa = sqrt(a), sb = sqrt(b), sab = sqrt(ab)."""
    return (x0 * x0 + a * x1 * x1 + b * x2 * x2 + a * b * x3 * x3, 2 * x0 * x1 + 2 * b * x2 * x3,
            2 * x0 * x2 + 2 * a * x1 * x3, 2 * x0 * x3 + 2 * x1 * x2)


def unit(p, q, D, a=3, b=11):
    d = [q[k] - p[k] for k in range(8)]
    s, t = square(*d[:4], a, b), square(*d[4:], a, b)
    return (s[0] + t[0], s[1] + t[1], s[2] + t[2], s[3] + t[3]) == (D * D, 0, 0, 0)


def cnf_text(n, E, cycles, fixed):
    x = lambda v, k: 4 * v + k + 1
    cl = []
    for v in range(n):
        cl.append([x(v, k) for k in range(4)])
        for k in range(4):
            for l in range(k + 1, 4):
                cl.append([-x(v, k), -x(v, l)])
    for i, j in E:
        for k in range(4):
            cl.append([-x(i, k), -x(j, k)])
    cl.append([x(fixed, 0)])
    for cyc in cycles:
        for k in range(4):
            cl.append([-x(cyc[t], (k + t) % 4) for t in range(len(cyc))])
    return '\n'.join([f'p cnf {4 * n} {len(cl)}'] + [' '.join(map(str, c)) + ' 0' for c in cl]) + '\n'


def main():
    W = load(sys.argv[1])
    D, P, E, col, cycles, fixed = (W['denominator'], W['points'], W['edges'], W['colouring'], W['cycles'],
                                   W['fixed_vertex'])
    n = len(P)
    field = W.get('field')
    _req(field in FIELDS, 'unknown field', field)
    a, b, basis = FIELDS[field]
    _req(W.get('basis') == basis, 'the basis should be', basis)
    _req(isinstance(D, int) and D > 0, 'bad denominator')
    _req(all(len(p) == 8 and all(isinstance(t, int) for t in p) for p in P), 'bad points')
    _req(isinstance(fixed, int) and 0 <= fixed < n, 'bad fixed_vertex')
    _req(len(set(map(tuple, P))) == n, 'repeated point')
    Es = set()
    for i, j in E:
        _req(isinstance(i, int) and isinstance(j, int) and 0 <= i < n and 0 <= j < n and i != j, 'bad edge', i, j)
        _req(unit(P[i], P[j], D, a, b), 'not a unit-distance edge', i, j)
        Es.add((min(i, j), max(i, j)))
    _req(len(Es) == len(E), 'repeated edge')
    print(f'(1) {field}: {n} points, {len(E)} edges, every edge at distance exactly 1')
    # all unit pairs, exactly; the filter compares the rational parts of the squared lengths, which must be D^2
    allpairs = missing = 0
    for i in range(n):
        pi = P[i]
        for j in range(i + 1, n):
            pj = P[j]
            d = [pj[k] - pi[k] for k in range(8)]
            r = (d[0] * d[0] + a * d[1] * d[1] + b * d[2] * d[2] + a * b * d[3] * d[3] + d[4] * d[4] + a * d[5] * d[5]
                 + b * d[6] * d[6] + a * b * d[7] * d[7])
            if r != D * D:
                continue
            if unit(pi, pj, D, a, b):
                allpairs += 1
                if (i, j) not in Es:
                    missing += 1
    print(f'    unit-distance pairs among the points: {allpairs}; not listed as edges: {missing}'
          f' ({"induced" if missing == 0 else "not induced"})')
    gens = W.get('generators')
    if gens is not None:
        _req(isinstance(gens, list) and len(gens) > 0 and all(
            isinstance(g, list) and len(g) == 8 and all(isinstance(t, int) for t in g) for g in gens), 'bad generators')
        for g in gens:
            _req(unit([0] * 8, g, D, a, b), 'a generator is not a unit vector', g)
        S = set()
        for g in gens:
            for h in (tuple(g), tuple(-t for t in g)):
                _req(h not in S, 'two generators are equal or opposite', g)
                S.add(h)
        idx = {tuple(p): k for k, p in enumerate(P)}
        cay = set()
        for k, p in enumerate(P):
            for h in S:
                j = idx.get(tuple(p[t] + h[t] for t in range(8)))
                if j is not None:
                    cay.add((min(k, j), max(k, j)))
        _req(cay == Es, 'the edges are not the pairs of points that differ by an element of U u -U')
        adj = [[] for _ in range(n)]
        for i, j in Es:
            adj[i].append(j)
            adj[j].append(i)
        seen, stack = {fixed}, [fixed]
        while stack:
            for v in adj[stack.pop()]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        print(f'(1b) {len(gens)} generators U: unit vectors, pairwise distinct up to sign; the edges are exactly the pairs '
              f'of points that differ by an element of U u -U; {len(seen)} of the {n} points are joined to the fixed '
              f'vertex (the {"origin" if P[fixed] == [0] * 8 else "point " + str(P[fixed])}) inside H')
    _req(len(col) == n and all(isinstance(c, int) and 0 <= c < 4 for c in col), 'bad colouring')
    for i, j in E:
        _req(col[i] != col[j], 'the colouring is not proper', i, j)
    print('(2) the colouring is a proper 4-colouring: chi(H) <= 4, chi_c(H) <= 4')
    for cyc in cycles:
        _req(len(cyc) >= 4 and len(cyc) % 4 == 0 and len(set(cyc)) == len(cyc), 'bad cycle', cyc)
        for k in range(len(cyc)):
            a, b = cyc[k], cyc[(k + 1) % len(cyc)]
            _req(isinstance(a, int) and 0 <= a < n and (min(a, b), max(a, b)) in Es, 'cycle uses a non-edge', cyc)
    print(f'(3) {len(cycles)} cycles, all closed walks along edges (lengths {min(map(len, cycles))}..'
          f'{max(map(len, cycles))}, divisible by 4)')
    txt = None
    if len(sys.argv) > 2:
        txt = cnf_text(n, E, cycles, fixed)
        op = gzip.open if sys.argv[2].endswith('.gz') else open
        with op(sys.argv[2], 'rt') as f:
            given = f.read()
        _req(given == txt, 'CNF differs from the formula built from the witness')
        print('(4) the CNF is the formula built from the edges, cycles and fixed vertex; sha256',
              hashlib.sha256(txt.encode()).hexdigest())
    if len(sys.argv) > 4:
        cnf, tmp = sys.argv[2], None
        if cnf.endswith('.gz'):                  # drat-trim reads plain text: a temporary copy, removed afterwards
            fd, tmp = tempfile.mkstemp(suffix='.cnf')
            with os.fdopen(fd, 'w') as f:
                f.write(txt)
            cnf = tmp
        try:
            r = subprocess.run([sys.argv[4], cnf, sys.argv[3]], capture_output=True, text=True)
        finally:
            if tmp:
                os.remove(tmp)
        ok = any(line.strip() == 's VERIFIED' for line in r.stdout.replace('\r', '\n').split('\n'))
        print('(5) drat-trim:', 's VERIFIED' if ok else r.stdout[-500:])
        _req(ok, 'drat-trim did not print the line "s VERIFIED"')
        print('=> every 4-colouring of H has a tight cycle, so chi_c(H) = chi(H) = 4')


if __name__ == '__main__':
    main()
