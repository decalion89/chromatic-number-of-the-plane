#!/usr/bin/env python3
"""Second referee encoding (E2): two-sided tight indicators, then a random renaming of the variables with random
polarity flips, a random clause order and a random literal order (seeded), so that neither the numbering nor the
clause list resembles the repository's encodings or E1.

Usage:  python3 referee_encode2.py W.json.gz OUT.cnf [SEED]

Unrenamed semantics:
  x(v,k)  vertex v has colour k (k in Z/7);   T(a,b)  the arc a->b of a listed cycle is tight.
  ALO and AMO per vertex; for every edge {u,w} and k: not(x(u,k) and x(w,k+delta)), delta in {0,1,6};
  T(a,b) <-> tight:  (-x(a,k) | -x(b,k+2) | T)  and  (-T | -x(a,k) | x(b,k+2))  for every k;
  every listed cycle: OR over its arcs of -T;  unit x(fixed_vertex, 0).
Soundness of 'UNSAT => every (7,2)-colouring has a tight listed cycle': a (7,2)-colouring c with a non-tight arc on
every listed cycle, rotated so that c(fixed_vertex) = 0 (rotation keeps all colour differences), satisfies every
clause with x(v,k) = [c(v) = k], T = [arc tight]; the renaming is a bijection on assignments.
The edge list is taken from the file; referee_check.py has verified that it equals the set of all unit pairs.
"""
import gzip, json, sys, random, hashlib


def main():
    with gzip.open(sys.argv[1], 'rt') as f:
        W = json.load(f)
    out = sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261004
    n = len(W['points'])
    E = [tuple(e) for e in W['edges']]
    cyc = W['cycles']
    fv = W.get('fixed_vertex')
    var = {}

    def V(key):
        if key not in var:
            var[key] = len(var) + 1
        return var[key]

    x = lambda v, k: V(('x', v, k % 7))
    cls = []
    for v in range(n):
        cls.append([x(v, k) for k in range(7)])
        cls += [[-x(v, k), -x(v, l)] for k in range(7) for l in range(k + 1, 7)]
    for (u, w) in E:
        for k in range(7):
            for dl in (0, 1, 6):
                cls.append([-x(u, k), -x(w, k + dl)])
    arcs = []
    for C in cyc:
        m = len(C)
        for i in range(m):
            a = (C[i], C[(i + 1) % m])
            if ('T',) + a not in var:
                arcs.append(a)
                t = V(('T',) + a)
                for k in range(7):
                    cls.append([-x(a[0], k), -x(a[1], k + 2), t])
                    cls.append([-t, -x(a[0], k), x(a[1], k + 2)])
    for C in cyc:
        m = len(C)
        cls.append([-var[('T', C[i], C[(i + 1) % m])] for i in range(m)])
    if fv is not None:
        cls.append([x(fv, 0)])
    nv = len(var)
    rng = random.Random(seed)
    perm = list(range(1, nv + 1))
    rng.shuffle(perm)
    flip = [rng.choice((1, -1)) for _ in range(nv + 1)]
    ren = lambda lit: (1 if lit > 0 else -1) * flip[abs(lit)] * perm[abs(lit) - 1]
    newcls = []
    for c in cls:
        c2 = [ren(l) for l in c]
        rng.shuffle(c2)
        newcls.append(c2)
    rng.shuffle(newcls)
    txt = f"p cnf {nv} {len(newcls)}\n" + ''.join(' '.join(map(str, c)) + ' 0\n' for c in newcls)
    with open(out, 'w') as f:
        f.write(txt)
    print(f"E2 {out}: {nv} variables ({7*n} colour, {len(arcs)} arc), {len(newcls)} clauses, seed {seed}, "
          f"sha256 {hashlib.sha256(txt.encode()).hexdigest()}")


if __name__ == '__main__':
    main()
