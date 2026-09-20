"""How far is the direction lattice from the ambient one?

If M has finite index in Z^d and that index is prime to n, every
homomorphism M -> Z/n extends: from 0 -> M -> Z^d -> Q -> 0 the sequence
Hom(Z^d, Z/n) -> Hom(M, Z/n) -> Ext^1(Q, Z/n) is exact, and multiplication by
|Q| is invertible on Z/n while killing Q, so the Ext term vanishes and the
restriction is onto.  So the two tests AGREE exactly when gcd(index, n) = 1,
and every old verdict at n stands unless n divides the index.

No SAT needed for that -- the index is the product of the Hermite pivots.
"""
import sys, time
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import lattice_basis, denominator_29_directions, edge_vectors

t0 = time.time()


def report(vecs, name):
    d = len(vecs[0])
    B = lattice_basis(vecs)
    piv = [next(i for i, x in enumerate(b) if x) for b in B]
    if len(B) < d:
        print(f"{name}: {len(vecs)} vectors, dim {d}, rank {len(B)} < d -- "
              f"the span is a proper subspace, so the ambient lattice was "
              f"never the right one  [{time.time()-t0:.0f}s]", flush=True)
        idx = None
    else:
        idx = 1
        for b, c in zip(B, piv):
            idx *= abs(b[c])
        bad = [n for n in (2, 3, 4, 5) if gcd(idx, n) > 1]
        print(f"{name}: {len(vecs)} vectors, dim {d}, full rank, index "
              f"{idx}; verdicts at risk for n = {bad or 'none'}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    return idx


report(denominator_29_directions(3), "denominator-29 set (box 3)")
from hn.degrey import build_G
from hn.graph import build_graph
report(edge_vectors(build_graph(build_G(as_graph=False))), "de Grey's G")
