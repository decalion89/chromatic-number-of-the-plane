"""The Minty-type characterisation, checked before it is trusted.

Goddyn, Tarsi and Zhang: chi_c(G) is the minimum over orientations D of the
maximum over cycles C of |C| / min(|C+|, |C-|), where C+ and C- are the
edges of C agreeing and disagreeing with a traversal.  An orientation with a
directed cycle gives min = 0 and is excluded automatically.

If it holds it changes what the target looks like.  chi_c > 9/2 becomes: in
EVERY acyclic orientation some cycle has fewer than 2|C|/9 edges the
minority way -- a 9-cycle with 2, a 5-cycle with 1, a 14-cycle with 3.  That
is a statement about the cycles of the graph, with no geometry in it at all.

Brute force over all orientations and all cycles, on graphs small enough to
do both, against values already measured by the solver.
"""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.geometry import DEGREY_FIELD as K, Point
from hn.graph import build_graph


def cycles(n, E):
    """Every cycle as a list of vertices, each once up to rotation/reversal."""
    adj = {v: set() for v in range(n)}
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    out, seen = [], set()
    for start in range(n):
        stack = [(start, [start], {start})]
        while stack:
            v, path, vis = stack.pop()
            for u in adj[v]:
                if u == start and len(path) >= 3:
                    key = frozenset((min(x, y), max(x, y)) for x, y in
                                    zip(path, path[1:] + [path[0]]))
                    if key not in seen:
                        seen.add(key)
                        out.append(list(path))
                elif u not in vis and u > start:
                    stack.append((u, path + [u], vis | {u}))
    return out


def chi_c_minty(n, E):
    E = sorted(set((min(a, b), max(a, b)) for a, b in E))
    C = cycles(n, E)
    best = None
    for bits in range(1 << len(E)):
        # bit set means the edge points from a to b
        d = {e: bool(bits >> i & 1) for i, e in enumerate(E)}
        worst = Fr(0)
        for cyc in C:
            fwd = 0
            L = len(cyc)
            for x, y in zip(cyc, cyc[1:] + [cyc[0]]):
                e = (min(x, y), max(x, y))
                if d[e] == (e[0] == x):
                    fwd += 1
            lo = min(fwd, L - fwd)
            if lo == 0:
                worst = None
                break
            worst = max(worst, Fr(L, lo))
        if worst is not None and (best is None or worst < best):
            best = worst
    return best


tests = []
for m in (3, 5, 7):
    tests.append((f"C_{m}", m, [(i, (i + 1) % m) for i in range(m)]))
tests.append(("K_4", 4, [(i, j) for i in range(4) for j in range(i + 1, 4)]))
half, r3, r11 = K.rational(Fr(1, 2)), K.sqrt(3), K.sqrt(11)
A = Point(K.zero(), K.zero())
D = Point(r3, K.zero())
B = Point(r3 * half, half)
Cp = Point(r3 * half, -half)
c, s = K.rational(Fr(5, 6)), r11 * K.rational(Fr(1, 6))
sp = [A, B, Cp, D] + [Point(p.x * c - p.y * s, p.x * s + p.y * c)
                      for p in (B, Cp, D)]
g = build_graph(sp)
tests.append(("Moser spindle", g.n,
              sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))))
known = {"C_3": Fr(3), "C_5": Fr(5, 2), "C_7": Fr(7, 3), "K_4": Fr(4),
         "Moser spindle": Fr(7, 2)}
print(f"{'graph':16s} {'edges':>6s} {'Minty':>8s} {'solver':>8s}  agree")
for name, n, E in tests:
    v = chi_c_minty(n, E)
    k = known[name]
    print(f"{name:16s} {len(E):6d} {str(v):>8s} {str(k):>8s}  "
          f"{'YES' if v == k else 'NO'}")
