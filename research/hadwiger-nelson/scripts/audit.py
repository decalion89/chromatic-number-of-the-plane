"""Audit every blocking claim against the lattice the directions generate.

The search for phi has been running over Z^d, the ambient lattice of the
coordinates, and the question is about M, the group the edge vectors actually
generate.  M sits inside Z^d, possibly properly, and Z/n is not injective, so
a phi on M need not extend -- the test can report blocking that is not there.
The gap is not hypothetical: the G u rho(G) direction set blocks at 5 over
Z^32 and does not over its own lattice.

So every claim gets re-run with the directions rewritten in a Hermite basis of
M, which makes M = Z^r by construction and the two questions the same.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import (has_homomorphism, denominator_29_directions,
                       edge_vectors)
exec(open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/hnfcheck.py").read()
     .split('print("re-checking')[0].split("t0 = time.time()")[1])

t0 = time.time()


def audit(vecs, name):
    d = len(vecs[0])
    raw = [n for n in (2, 3, 4, 5) if has_homomorphism(vecs, n)[0] is None]
    red, r = in_basis(vecs, d)
    lat = [n for n in (2, 3, 4, 5) if has_homomorphism(red, n)[0] is None]
    mark = "  *** DIFFERENT ***" if raw != lat else "  (agrees)"
    print(f"{name}: {len(vecs)} vectors, dim {d}, lattice rank {r}\n"
          f"    over Z^d {raw}   over the lattice {lat}{mark}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return raw, lat


audit(denominator_29_directions(3), "denominator-29 set (box 3)")
from hn.degrey import build_G
from hn.graph import build_graph
audit(edge_vectors(build_graph(build_G(as_graph=False))), "de Grey's G")
