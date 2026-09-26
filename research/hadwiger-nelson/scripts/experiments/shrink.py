"""What makes a unit-distance graph TIGHT at its chromatic number?

Sa and Y have chi_c = 4 exactly, where the Moser spindle is loose at 7/2.
Tightness is what G lacks at five, and it is the whole of the remaining gap,
so the mechanism behind it is worth seeing rather than guessing at.  Shrink
Sa to a minimal subgraph that still refuses K(35/9) = 3.8889 -- the largest
ratio below 4 with denominator at most 9 -- and whatever is left is the
structure doing the work.

Removal is by delta-debugging rather than one vertex at a time: try dropping
a whole block, keep the drop if the graph still refuses, otherwise put it
back and halve the block.  That turns a few hundred solver calls into a few
dozen.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
WHICH = sys.argv[1] if len(sys.argv) > 1 else "Sa"
R = Fr(*map(int, (sys.argv[2] if len(sys.argv) > 2 else "35/9").split("/")))
t0 = time.time()
P = {"Sa": build_Sa, "Y": build_Y,
     "G": lambda k: build_G(k, as_graph=False)}[WHICH](K)
b = IntBasis.covering(P)
rows = b.rows(P)
adj = defaultdict(set)
for a, c in fast_edges_complete(b, rows):
    adj[a].add(c)
    adj[c].add(a)
p, q = R.numerator, R.denominator
print(f"{WHICH}: {len(P)} points, "
      f"{sum(len(v) for v in adj.values())//2} edges; refusing "
      f"K({p}/{q}) = {float(R):.4f}?  [{time.time()-t0:.0f}s]", flush=True)


def refuses(keep):
    """True when the induced subgraph has NO homomorphism to K(p/q)."""
    ks = sorted(keep)
    idx = {v: i for i, v in enumerate(ks)}
    cl = [[1 + i * p + j for j in range(p)] for i in range(len(ks))]
    for v in ks:
        for u in adj[v]:
            if u in idx and u > v:
                for j in range(p):
                    for d in range(-(q - 1), q):
                        cl.append([-(1 + idx[v] * p + j),
                                   -(1 + idx[u] * p + (j + d) % p)])
    cl += [[-(1 + j)] for j in range(1, p)] + [[1]]
    s = Solver(name="cd15", bootstrap_with=cl)
    ok = s.solve()
    s.delete()
    return not ok


keep = set(range(len(P)))
assert refuses(keep), "carrier does not refuse this ratio"
print(f"yes -- now shrinking  [{time.time()-t0:.0f}s]", flush=True)
# peel first: at p/q a neighbour forbids 2q-1 positions, so a vertex with
# fewer than p/(2q-1) neighbours can always be coloured last and cannot
# belong to a minimal refusing subgraph.
need = -(-p // (2 * q - 1))
while True:
    drop = {v for v in keep if len(adj[v] & keep) < need}
    if not drop:
        break
    keep -= drop
print(f"peeled to {len(keep)} (min degree {need})  [{time.time()-t0:.0f}s]",
      flush=True)
block = max(1, len(keep) // 2)
calls = 0
while block >= 1:
    order = sorted(keep)
    random.shuffle(order)
    moved = False
    for i in range(0, len(order), block):
        chunk = set(order[i:i + block]) & keep
        if not chunk:
            continue
        calls += 1
        if refuses(keep - chunk):
            keep -= chunk
            moved = True
            while True:
                d = {v for v in keep if len(adj[v] & keep) < need}
                if not d:
                    break
                keep -= d
            print(f"   dropped {len(chunk)} -> {len(keep)} left "
                  f"({calls} calls)  [{time.time()-t0:.0f}s]", flush=True)
    if not moved or block == 1:
        if block == 1:
            break
        block = max(1, block // 2)
    else:
        block = max(1, min(block, len(keep) // 2))
ks = sorted(keep)
m = sum(1 for v in ks for u in adj[v] if u in keep) // 2
degs = sorted(len(adj[v] & keep) for v in ks)
print(f"\nminimal refusing subgraph: {len(ks)} points, {m} edges, "
      f"degrees {degs[0]}..{degs[-1]}, mean {2*m/len(ks):.2f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
pickle.dump([P[v] for v in ks],
            open(SC + f"tight_{WHICH}_{p}_{q}.pkl", "wb"))
print("DONE", flush=True)
