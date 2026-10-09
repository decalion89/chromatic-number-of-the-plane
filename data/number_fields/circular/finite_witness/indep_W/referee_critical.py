#!/usr/bin/env python3
"""Referee check of the vertex-criticality certificates (written from scratch).

Usage:  python3 referee_critical.py W.json.gz W_critical.json.gz

For every vertex v, the certificate c_v : V -> Z/7 must satisfy c_v(v) = -1, c_v(x) in 0..6 for x != v,
be a (7,2)-colouring of H - v, and have an acyclic tight digraph (arcs x -> y with c_v(y) - c_v(x) = 2 mod 7).

From such a c_v we build EXPLICITLY a homomorphism of H - v to K_{7N/(2N+1)}, N = |V(H - v)|:
    f(x) = N * c_v(x) + pos(x)   (mod 7N),
pos = position in a topological order of the tight digraph (0..N-1), and check it edge by edge:
  for every edge xy of H - v, (f(y) - f(x)) mod 7N lies in [2N+1, 5N-1].
Hence chi_c(H - v) <= 7N/(2N+1) < 7/2, without any appeal to real perturbations or limits.
(Proof that it works: with delta = c_v(y)-c_v(x) mod 7 in {2,3,4,5} and e = pos(y)-pos(x), |e| <= N-1,
 f(y)-f(x) = N*delta + e; delta in {3,4} gives [2N+1, 5N-1]; delta = 2 is a tight arc x->y, so e >= 1 and
 N*2+e in [2N+1, 3N-1]; delta = 5 is a tight arc y->x, so e <= -1 and 5N+e in [4N+1, 5N-1].)
"""
import gzip, json, sys
from collections import defaultdict, deque


if not __debug__:
    sys.exit('run without -O: this script relies on assert statements')


def main():
    with gzip.open(sys.argv[1], 'rt') as f:
        W = json.load(f)
    with gzip.open(sys.argv[2], 'rt') as f:
        C = json.load(f)
    n = len(W['points'])
    edges = [tuple(e) for e in W['edges']]
    certs = C['critical_colourings']
    assert isinstance(certs, list) and len(certs) == n, "need one certificate per vertex"
    worst_ratio = None
    ntight_total = 0
    for v in range(n):
        c = certs[v]
        assert isinstance(c, list) and len(c) == n
        assert c[v] == -1, f"vertex {v}: c_v(v) != -1"
        assert all(isinstance(t, int) and not isinstance(t, bool) and 0 <= t < 7 for i, t in enumerate(c) if i != v), \
            f"vertex {v}: bad colour values"
        Ev = [(x, y) for (x, y) in edges if x != v and y != v]
        out = defaultdict(list)
        indeg = [0] * n
        for (x, y) in Ev:
            dlt = (c[y] - c[x]) % 7
            assert dlt in (2, 3, 4, 5), f"vertex {v}: edge {x}-{y} has difference {dlt}"
            if dlt == 2:
                out[x].append(y); indeg[y] += 1
            elif dlt == 5:
                out[y].append(x); indeg[x] += 1
        ntight_total += sum(len(l) for l in out.values())
        # Kahn's algorithm on V - v
        verts = [x for x in range(n) if x != v]
        q = deque(x for x in verts if indeg[x] == 0)
        pos = {}
        while q:
            x = q.popleft()
            pos[x] = len(pos)
            for y in out[x]:
                indeg[y] -= 1
                if indeg[y] == 0:
                    q.append(y)
        assert len(pos) == len(verts), f"vertex {v}: tight digraph has a directed cycle"
        N = len(verts)
        p, qq = 7 * N, 2 * N + 1
        f = {x: (N * c[x] + pos[x]) % p for x in verts}
        for (x, y) in Ev:
            dd = (f[y] - f[x]) % p
            assert qq <= dd <= p - qq, f"vertex {v}: explicit ({p},{qq})-colouring fails on {x}-{y}"
        worst_ratio = (p, qq)
    print(f"{sys.argv[1].split('/')[-1]}: {n} certificates OK: each c_v is a (7,2)-colouring of H-v with an acyclic tight "
          f"digraph; explicit (7N,2N+1)-colourings with N={n-1} verified, so chi_c(H-v) <= {worst_ratio[0]}/{worst_ratio[1]} "
          f"= {worst_ratio[0]/worst_ratio[1]:.6f} < 3.5 for every v; average tight arcs per certificate "
          f"{ntight_total/n:.1f}")


if __name__ == '__main__':
    main()
