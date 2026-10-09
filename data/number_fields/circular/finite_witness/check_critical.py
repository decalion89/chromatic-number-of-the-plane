"""Check the vertex-criticality certificates of a finite witness for chi_c = 7/2.

usage: python3 check_critical.py WITNESS.json[.gz] CRITICAL.json[.gz]
For every vertex v of H, CRITICAL gives a map c_v: V -> Z/7 with c_v(v) = -1.  Checked: c_v is a (7,2)-colouring of
H - v (c_v(y) - c_v(x) in {2, 3, 4, 5} mod 7 on every edge not at v), and the tight digraph of c_v on H - v (arcs
x -> y with c_v(y) - c_v(x) = 2 mod 7) has no directed cycle (Kahn's algorithm).  Then chi_c(H - v) < 7/2 for every v
(Guichard: a (7,2)-colouring with an acyclic tight digraph can be perturbed into a (p,q)-colouring with p/q < 7/2), so
no proper induced subgraph of H has circular chromatic number 7/2."""
import gzip, json, sys
from collections import deque


def _req(ok, *msg):
    """An explicit check (not assert, so that python -O cannot skip it)."""
    if not ok:
        print('REJECTED:', *msg, file=sys.stderr)
        sys.exit(1)


def load(p):
    return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))

W = load(sys.argv[1]); C = load(sys.argv[2])
n = len(W['points']); E = [tuple(e) for e in W['edges']]
certs = C['critical_colourings']
_req(len(certs) == n, 'one certificate per vertex')
for v, col in enumerate(certs):
    _req(len(col) == n and col[v] == -1 and all(isinstance(c, int) and 0 <= c < 7 for u, c in enumerate(col) if u != v),
         'certificate', v, 'is not a map V - v -> Z/7 with v marked -1')
    out = [[] for _ in range(n)]; indeg = [0] * n
    for i, j in E:
        if v in (i, j):
            continue
        dd = (col[j] - col[i]) % 7
        _req(dd in (2, 3, 4, 5), 'not a (7,2)-colouring of H - v:', v, 'edge', i, j)
        if dd == 2:
            out[i].append(j); indeg[j] += 1
        elif dd == 5:
            out[j].append(i); indeg[i] += 1
    q = deque(u for u in range(n) if u != v and indeg[u] == 0); seen = 0
    while q:
        u = q.popleft(); seen += 1
        for w in out[u]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    _req(seen == n - 1, 'the tight digraph of H - v has a directed cycle:', v)
print(f'{n} vertices: for every v, a (7,2)-colouring of H - v with an acyclic tight digraph, so chi_c(H - v) < 7/2;'
      ' H is vertex-critical for chi_c = 7/2')
