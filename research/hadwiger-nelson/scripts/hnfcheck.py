"""Test blocking on the lattice the directions generate, not on Z^d.

A coset colouring is a homomorphism phi from the additive group the edge
vectors generate -- call it M -- to Z/n, nonzero on every vector.  Searching
instead over phi : Z^d -> Z/n only finds the ones that extend to the ambient
lattice, and Z/n is not injective, so some do not.  The gap is real and it
points the dangerous way, towards claiming blocking that is not there:

    d = 1, D = {2}, n = 2.  Over M = 2Z the map phi(2) = 1 escapes.  Over Z
    every psi has psi(2) = 0, so the test reports blocking.

Dividing out the global content kills that example, which is why the pipeline
has done it from the start, but content 1 does not imply M = Z^d.  The fix is
to put the vectors in a Z-basis of M, by Hermite reduction, so that M IS Z^r
by construction and the two questions coincide.

Everything this work has claimed to block is re-checked here.
"""
import sys, time
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import has_homomorphism

t0 = time.time()


def _ext_gcd(a, b):
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def hermite(rows, d):
    """A Z-basis of the lattice the rows generate, in echelon form."""
    basis = {}
    for r0 in rows:
        w = list(r0)
        c = 0
        while c < d:
            if w[c] == 0:
                c += 1
                continue
            if c not in basis:
                basis[c] = [-x for x in w] if w[c] < 0 else w
                break
            b = basis[c]
            g, x, y = _ext_gcd(b[c], w[c])
            nb = [x * bb + y * ww for bb, ww in zip(b, w)]
            nw = [(b[c] // g) * ww - (w[c] // g) * bb
                  for bb, ww in zip(b, w)]
            basis[c] = nb
            w = nw
            c += 1
    piv = sorted(basis)
    return [basis[c] for c in piv], piv


def in_basis(vecs, d):
    """Rewrite each vector in a Z-basis of the lattice they generate."""
    B, piv = hermite(vecs, d)
    out = []
    for v in vecs:
        w, co = list(v), []
        for b, c in zip(B, piv):
            q, rem = divmod(w[c], b[c])
            if rem:
                raise AssertionError("not in the lattice -- Hermite is wrong")
            co.append(q)
            if q:
                w = [x - q * y for x, y in zip(w, b)]
        if any(w):
            raise AssertionError("residue left over")
        out.append(tuple(co))
    return sorted(set(out)), len(B)


def verdict(vecs, name):
    d = len(vecs[0])
    raw = [n for n in (2, 3, 4, 5) if has_homomorphism(vecs, n)[0] is None]
    red, r = in_basis(vecs, d)
    lat = [n for n in (2, 3, 4, 5) if has_homomorphism(red, n)[0] is None]
    flag = "  <-- DIFFERENT" if raw != lat else ""
    print(f"  {name}: {len(vecs)} vectors, dim {d} -> rank {r}; "
          f"over Z^d blocks {raw}, over the lattice blocks {lat}{flag}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return raw, lat


print("re-checking every blocking claim on the true module", flush=True)
from hn.degrey import build_G
from hn.graph import build_graph
from hn.homcol import edge_vectors
G = build_graph(build_G(as_graph=False))
verdict(edge_vectors(G), "de Grey's G")

from hn.spindle import moser_spindle
try:
    verdict(edge_vectors(moser_spindle()), "Moser spindle")
except Exception as e:
    print(f"  (moser: {e})", flush=True)
