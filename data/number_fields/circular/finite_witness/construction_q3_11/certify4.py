"""certify4.py WITNESS.json OUT.cnf: the formula of 'every 4-colouring of H has a tight cycle among the listed ones'.

WITNESS.json: {"D", "units", "points", "cycles"} (grow8 output; points are integer 8-vectors over (1, s3, s11, s33)/D).
Edges are recomputed exactly: all pairs of points whose difference is a unit vector of Q(sqrt3, sqrt11)^2 (not only
the listed units).  Encoding without auxiliary variables: x(v, c) for v a point and c in Z/4; each point has a colour,
at most one; adjacent points differ; vertex 0 (the origin) has colour 0 (colour rotation keeps tightness); every
listed cycle v_0 -> ... -> v_{m-1} -> v_0 is checked to be a closed walk of edges and gets, for each c in Z/4, the
clause 'not (x(v_0, c) and x(v_1, c + 1) and ... and x(v_{m-1}, c + m - 1))', i.e. it is not tight starting at c
(a tight cycle has m = 0 mod 4; for m != 0 mod 4 no clause is needed).  UNSAT means every 4-colouring of the graph
has one of the listed cycles tight, so chi_c >= 4 (with a 4-colouring, = 4)."""
import sys, json, itertools
import numpy as np

W = json.load(open(sys.argv[1])); out = sys.argv[2]
P = np.array(W["points"], dtype=np.int64); n = len(P); D = W["D"]
idx = {tuple(p): i for i, p in enumerate(W["points"])}


def is_unit(v):
    a0, a1, a2, a3, b0, b1, b2, b3 = (int(x) for x in v)
    def sq(x0, x1, x2, x3):
        return (x0 * x0 + 3 * x1 * x1 + 11 * x2 * x2 + 33 * x3 * x3, 2 * x0 * x1 + 22 * x2 * x3,
                2 * x0 * x2 + 6 * x1 * x3, 2 * x0 * x3 + 2 * x1 * x2)
    s, t = sq(a0, a1, a2, a3), sq(b0, b1, b2, b3)
    return (s[0] + t[0], s[1] + t[1], s[2] + t[2], s[3] + t[3]) == (D * D, 0, 0, 0)


# all unit pairs: vectorised prefilter on the rational part, then exact check
edges = set()
B = 700
for i0 in range(0, n, B):
    A = P[i0:i0 + B]
    for j0 in range(i0, n, B):
        Bm = P[j0:j0 + B]
        Dl = A[:, None, :] - Bm[None, :, :]                       # differences
        x = Dl[..., :4]; y = Dl[..., 4:]
        r0 = (x[..., 0] ** 2 + 3 * x[..., 1] ** 2 + 11 * x[..., 2] ** 2 + 33 * x[..., 3] ** 2 +
              y[..., 0] ** 2 + 3 * y[..., 1] ** 2 + 11 * y[..., 2] ** 2 + 33 * y[..., 3] ** 2)
        cand = np.argwhere(r0 == D * D)
        for a, b in cand:
            i, j = i0 + int(a), j0 + int(b)
            if i < j and is_unit(P[i] - P[j]):
                edges.add((i, j))
adj = [set() for _ in range(n)]
for i, j in edges:
    adj[i].add(j); adj[j].add(i)
print(f"{n} points, {len(edges)} unit pairs", flush=True)

v = lambda i, c: 4 * i + c + 1
cls = []
for i in range(n):
    cls.append([v(i, c) for c in range(4)])
    for c1, c2 in itertools.combinations(range(4), 2):
        cls.append([-v(i, c1), -v(i, c2)])
for i, j in sorted(edges):
    for c in range(4):
        cls.append([-v(i, c), -v(j, c)])
zero = idx[(0,) * 8]
cls.append([v(zero, 0)])
ncyc = 0
for cyc in W["cycles"]:
    m = len(cyc)
    for a, b in zip(cyc, cyc[1:] + cyc[:1]):
        assert b in adj[a], "a listed cycle uses a non-edge"
    assert len(set(cyc)) == m
    if m % 4:
        continue
    for c in range(4):
        cls.append([-v(cyc[k], (c + k) % 4) for k in range(m)])
    ncyc += 1
with open(out, "w") as fh:
    fh.write(f"p cnf {4 * n} {len(cls)}\n")
    for c in cls:
        fh.write(" ".join(map(str, c)) + " 0\n")
print(f"wrote {out}: {4 * n} variables, {len(cls)} clauses ({ncyc} cycles of length 0 mod 4)", flush=True)
