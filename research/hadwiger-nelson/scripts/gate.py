"""The gate, exactly: coset colourings over the module the edge vectors generate.

The ambient test (phi on Z^dim) can disagree with the true one (phi on the
module M the edge vectors span) when M is not saturated.  So: integer echelon
basis of M, each edge vector in that basis, then has_homomorphism at n = 2..5,
and periodic_screen at n = 6..12.
"""
import sys, json, time
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.homcol import has_homomorphism, periodic_screen

def edge_vectors(g):
    out = set()
    for a, b in g.edges():
        d = (g.vertices[b].x - g.vertices[a].x, g.vertices[b].y - g.vertices[a].y)
        out.add(tuple(d[0].c) + tuple(d[1].c))
    den = 1
    for v in out:
        for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
    vs = sorted({tuple(int(Fr(q) * den) for q in v) for v in out})
    half = {}
    for v in vs:
        k = max(v, tuple(-x for x in v)); half[k] = True
    return list(half)

def echelon(vecs):
    B = []   # rows with distinct pivots, pivot entry > 0
    for v in vecs:
        v = list(v)
        for r in B:
            p = next(i for i, x in enumerate(r) if x)
            if v[p] == 0: continue
            # gcd-combine v into r on column p
            a, b = r[p], v[p]
            # extended gcd
            x0, x1, y0, y1, aa, bb = 1, 0, 0, 1, a, b
            while bb:
                q = aa // bb; aa, bb = bb, aa - q * bb
                x0, x1 = x1, x0 - q * x1; y0, y1 = y1, y0 - q * y1
            gg = aa  # = x0*a + y0*b
            newr = [x0 * ri + y0 * vi for ri, vi in zip(r, v)]
            v = [(a // gg) * vi - (b // gg) * ri for ri, vi in zip(r, v)]
            r[:] = newr
            if newr[p] < 0: r[:] = [-x for x in r]
        if any(v):
            if v[next(i for i, x in enumerate(v) if x)] < 0: v = [-x for x in v]
            B.append(v)
            B.sort(key=lambda r: next(i for i, x in enumerate(r) if x))
        # keep pivots distinct: re-reduce until stable
        changed = True
        while changed:
            changed = False
            for i in range(len(B)):
                for j in range(len(B)):
                    if i != j and B[i] and B[j]:
                        pi = next(t for t, x in enumerate(B[i]) if x)
                        pj = next(t for t, x in enumerate(B[j]) if x)
                        if pi == pj:
                            w = B.pop(j); changed = True
                            B2 = echelon([w]) if False else None
                            # merge w into B[i] via one gcd step
                            r = B[i if i < j else i - 1]; a, b = r[pi], w[pi]
                            x0, x1, y0, y1, aa, bb = 1, 0, 0, 1, a, b
                            while bb:
                                q = aa // bb; aa, bb = bb, aa - q * bb
                                x0, x1 = x1, x0 - q * x1; y0, y1 = y1, y0 - q * y1
                            gg = aa
                            newr = [x0 * ri + y0 * wi for ri, wi in zip(r, w)]
                            rest = [(a // gg) * wi - (b // gg) * ri for ri, wi in zip(r, w)]
                            r[:] = newr if newr[pi] > 0 else [-x for x in newr]
                            if any(rest):
                                if rest[next(t for t, x in enumerate(rest) if x)] < 0: rest = [-x for x in rest]
                                B.append(rest)
                            B.sort(key=lambda r: next(t for t, x in enumerate(r) if x))
                            break
                    if changed: break
                if changed: break
    return B

def coords(B, v):
    v = list(v); c = []
    for r in B:
        p = next(i for i, x in enumerate(r) if x)
        assert v[p] % r[p] == 0, "not in the module"
        k = v[p] // r[p]; c.append(k)
        v = [vi - k * ri for vi, ri in zip(v, r)]
    assert not any(v)
    return c

def gate(g, label):
    t0 = time.time()
    E = edge_vectors(g)
    B = echelon(E)
    C = [coords(B, v) for v in E]
    print(f"  {label}: n={g.n} m={g.m}, {len(E)} directions (up to sign), module rank {len(B)}", flush=True)
    for n in (2, 3, 4, 5):
        phi, why = has_homomorphism(C, n)
        print(f"    mod {n}: {'coset colouring exists' if phi else 'BLOCKS (' + why + ')'}", flush=True)
    for n in range(6, 13):
        r = periodic_screen(C, n, cap=5, rounds=60)
        print(f"    Z/{n}: {'periodic 5-colouring FOUND' if r.get('colourable') else 'no periodic 5-colouring (' + ('exhausted' if r.get('exhausted') else 'rounds out') + ')'}   [{time.time()-t0:.0f}s]", flush=True)

if __name__ == "__main__":
    exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/screen_lambda.py").read().split("for name, pts in")[0])
    gate(build_graph(VK), "E-I K, unit edges")
    ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
    d = json.load(open(f"{ROOT}/data/five_247_c.json"))
    F2 = Field(tuple(d["field_generators"]))
    P2 = [Point(F2.element([Fr(a, b) for a, b in x]), F2.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P2)
    gate(g, "803 alone")
    v0 = max(range(g.n), key=lambda v: len(g.adj[v])); c = P2[v0]
    ca = F2.rational(Fr(49, 50)); sa = F2.rational(Fr(3, 50)) * F2.sqrt(11)
    g2 = build_graph(list(dict.fromkeys(P2 + [rot_about(c, p) for p in P2])))
    gate(g2, "803 u lambda(803)")
