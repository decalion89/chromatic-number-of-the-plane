"""Kempe swaps of coset colourings, exactly.

In c0 = psi + k on M, the Kempe component of a point x in colours alpha,
beta = alpha + delta is x + L  u  x + p0 + L, with P = {units u : psi(u) = delta}
and L = <p - p' : p, p' in P> inside ker psi.  Swapping one component is a
proper colouring of the whole unit-distance graph on M.  It joins a, a + 2e
(alpha = c0(a), delta = 2 psi(e)) iff 2e - p0 is NOT in L, and it splits
a, a + 5e (any delta) iff 5e is NOT in L.  So everything hinges on the index
[ker psi : L]; index 1 means every such swap is trivial.
"""
import sys, json, itertools
from fractions import Fraction as Fr
exec(open(__import__("os").path.join("/home/user/darwin-50/research/hadwiger-nelson/scripts", "gate.py")).read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
name = sys.argv[1]
d = json.load(open(f"{ROOT}/data/{name}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P))
B = echelon(E); C = [coords(B, v) for v in E]; r = len(B)
Cs = C + [[-x for x in c] for c in C]                   # all units, both signs
def hnf_index(vecs, r):
    """|Z^r / span(vecs)| (0 if not full rank), by integer echelon + product of pivots."""
    rows = [list(v) for v in vecs if any(v)]
    basis = []
    for v in rows:
        v = v[:]
        while any(v):
            p = next(i for i, t in enumerate(v) if t)
            row = next((b for b in basis if next(i for i, t in enumerate(b) if t) == p), None)
            if row is None:
                basis.append(v if v[p] > 0 else [-t for t in v]); break
            a, b = row[p], v[p]
            x0, x1, y0, y1, aa, bb = 1, 0, 0, 1, a, b
            while bb:
                q = aa // bb; aa, bb = bb, aa - q * bb; x0, x1, y0, y1 = x1, x0 - q * x1, y1, y0 - q * y1
            new = [x0 * s + y0 * w for s, w in zip(row, v)]
            v = [(a // aa) * w - (b // aa) * s for s, w in zip(row, v)]
            row[:] = new if new[p] > 0 else [-t for t in new]
    if len(basis) < r: return 0, basis
    idx = 1
    for b in basis: idx *= b[next(i for i, t in enumerate(b) if t)]
    return abs(idx), basis
def member(v, basis):
    v = list(v)
    for b in sorted(basis, key=lambda b: next(i for i, t in enumerate(b) if t)):
        p = next(i for i, t in enumerate(b) if t)
        if v[p] % b[p]: return False
        k = v[p] // b[p]; v = [x - k * y for x, y in zip(v, b)]
    return not any(v)
adm = [psi for psi in itertools.product(range(5), repeat=r) if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in C)]
print(f"{name}: {len(C)} directions, rank {r}, {len(adm)} admissible psi", flush=True)
from collections import Counter
hist = Counter(); refute2 = set(); refute5 = set()
for psi in adm:
    val = lambda c: sum(a * b for a, b in zip(psi, c)) % 5
    for delta in range(1, 5):
        Pd = [c for c in Cs if val(c) == delta]
        p0 = Pd[0]
        idx, basis = hnf_index([[x - y for x, y in zip(p, p0)] for p in Pd[1:]], r)
        rel = idx // 5 if idx else 0
        hist[rel] += 1
        if rel == 1: continue
        for k, c in enumerate(C):
            for sgn in (1, -1):
                e = [sgn * x for x in c]
                if val(e) * 2 % 5 == delta and not member([2 * x - y for x, y in zip(e, p0)], basis): refute2.add(k)
                if not member([5 * x for x in e], basis): refute5.add(k)
print(f"  [ker psi : L] over all (psi, delta): {dict(hist)}  (0 = L not of full rank)")
print(f"  Kempe swaps join a 2e pair in {len(refute2)} directions, split a 5e pair in {len(refute5)} directions")
