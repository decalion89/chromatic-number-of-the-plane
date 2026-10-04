"""Stand-alone checker for the finite witness of chi_c = 7/2 over F = Q(sqrt11).

usage: python3 check_witness.py WITNESS.json.gz [CNF] [PROOF.drat] [path/to/drat-trim]

The witness file contains
  denominator  D (an integer, prime to 7)
  points       integer 4-tuples [a, b, c, e]: the point ((a + b r)/D, (c + e r)/D) of F^2, r = sqrt11
  edges        pairs [i, j]
  colouring    a map V -> Z/7
  cycles       lists [v_0, ..., v_{m-1}] of vertices
The checker verifies, with Python integers only:
  (1) the points are distinct and every edge joins two points at Euclidean distance exactly 1:
      (Da)^2 + 11 (Db)^2 + (Dc)^2 + 11 (De)^2 = D^2 and Da*Db + Dc*De = 0 for the difference (Da, Db, Dc, De);
      it also reports whether `edges` is the set of ALL pairs at distance 1 (the induced unit-distance graph);
  (2) the colouring is a (7,2)-colouring: c(y) - c(x) in {2,3,4,5} mod 7 on every edge (so chi_c(H) <= 7/2);
  (3) every listed cycle is a closed walk v_0 -> v_1 -> ... -> v_{m-1} -> v_0 along edges of H, with distinct vertices;
  (4) the CNF (if given) is exactly the formula built here from (edges, cycles):
        variables x(v,k) = 7v + k + 1 ("v has colour k"), and one variable t(a,b) per arc of a listed cycle;
        clauses  OR_k x(v,k)                                   (every vertex gets a colour)
                 -x(i,k) v -x(j,k+d)  for d in {0,1,6}, every edge ij and k   (difference not in {0,1,6})
                 -x(a,k) v -x(b,k+2) v t(a,b)                   (t(a,b) is true when the arc a->b is tight)
                 OR over the arcs of C of -t(a,b)               (C is not tight), for every listed cycle C;
  (5) drat-trim (if given) verifies the DRAT proof of unsatisfiability of the CNF.
Logic: a (7,2)-colouring c without a tight listed cycle would satisfy the CNF (x(v,k) iff c(v) = k,
t(a,b) iff c(b) - c(a) = 2), so (5) shows that every (7,2)-colouring of H has a tight cycle (a directed cycle of the
tight digraph), and Lemma 20 of papers/three-colours (Guichard) gives chi_c(H) >= 7/2; with (2), chi_c(H) = 7/2.
"""
import sys, json, gzip, hashlib, subprocess, itertools


def load(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        return json.load(f)


def cnf_text(n, E, cycles):
    x = lambda v, c: 7 * v + c + 1
    cl = []
    for v in range(n):
        cl.append([x(v, c) for c in range(7)])
    for i, j in E:
        for c in range(7):
            for d in (0, 1, 6):
                cl.append([-x(i, c), -x(j, (c + d) % 7)])
    tv = {}; nv = 7 * n
    for cyc in cycles:
        for k in range(len(cyc)):
            a, b = cyc[k], cyc[(k + 1) % len(cyc)]
            if (a, b) not in tv:
                nv += 1; tv[(a, b)] = nv
                for c in range(7):
                    cl.append([-x(a, c), -x(b, (c + 2) % 7), nv])
    for cyc in cycles:
        cl.append([-tv[(cyc[k], cyc[(k + 1) % len(cyc)])] for k in range(len(cyc))])
    return '\n'.join([f'p cnf {nv} {len(cl)}'] + [' '.join(map(str, c)) + ' 0' for c in cl]) + '\n'


def unit(p, q, D):
    a, b, c, e = (q[k] - p[k] for k in range(4))
    return a * a + 11 * b * b + c * c + 11 * e * e == D * D and a * b + c * e == 0


def main():
    W = load(sys.argv[1])
    D, P, E, col, cycles = W['denominator'], W['points'], W['edges'], W['colouring'], W['cycles']
    n = len(P)
    assert D % 7 != 0 and all(len(p) == 4 and all(isinstance(t, int) for t in p) for p in P)
    assert len(set(map(tuple, P))) == n, 'repeated point'
    Es = set()
    for i, j in E:
        assert 0 <= i < n and 0 <= j < n and i != j
        assert unit(P[i], P[j], D), ('not a unit-distance edge', i, j)
        Es.add((min(i, j), max(i, j)))
    assert len(Es) == len(E), 'repeated edge'
    print(f'(1) {n} points, {len(E)} edges, every edge at distance exactly 1')
    # all unit pairs (exact; quadratic loop with a cheap integer filter on the first coordinate norm)
    allpairs = 0; missing = 0
    for i in range(n):
        pi = P[i]
        for j in range(i + 1, n):
            pj = P[j]
            a = pj[0] - pi[0]; c = pj[2] - pi[2]
            if a * a + c * c > D * D:
                continue
            if unit(pi, pj, D):
                allpairs += 1
                if (i, j) not in Es:
                    missing += 1
    print(f'    unit-distance pairs among the points: {allpairs}; not listed as edges: {missing}'
          f' ({"induced" if missing == 0 else "not induced"})')
    assert len(col) == n and all(0 <= c < 7 for c in col)
    for i, j in E:
        assert (col[j] - col[i]) % 7 in (2, 3, 4, 5), ('bad colouring', i, j)
    print('(2) the colouring is a (7,2)-colouring: chi_c(H) <= 7/2')
    for cyc in cycles:
        assert len(cyc) >= 3 and len(set(cyc)) == len(cyc)
        for k in range(len(cyc)):
            a, b = cyc[k], cyc[(k + 1) % len(cyc)]
            assert (min(a, b), max(a, b)) in Es, ('cycle uses a non-edge', cyc)
    print(f'(3) {len(cycles)} cycles, all closed walks along edges (lengths {min(map(len, cycles))}..{max(map(len, cycles))})')
    if len(sys.argv) > 2:
        txt = cnf_text(n, E, cycles)
        op = gzip.open if sys.argv[2].endswith('.gz') else open
        with op(sys.argv[2], 'rt') as f:
            given = f.read()
        assert given == txt, 'CNF differs from the formula built from the witness'
        print('(4) the CNF is the formula built from the edges and cycles; sha256',
              hashlib.sha256(txt.encode()).hexdigest())
    if len(sys.argv) > 4:
        cnf = sys.argv[2]
        if cnf.endswith('.gz'):
            cnf = 'witness_check_tmp.cnf'
            open(cnf, 'w').write(txt)
        r = subprocess.run([sys.argv[4], cnf, sys.argv[3]], capture_output=True, text=True)
        ok = 's VERIFIED' in r.stdout
        print('(5) drat-trim:', 's VERIFIED' if ok else r.stdout[-500:])
        assert ok
        print('=> every (7,2)-colouring of H has a tight cycle, so chi_c(H) = 7/2')


if __name__ == '__main__':
    main()
