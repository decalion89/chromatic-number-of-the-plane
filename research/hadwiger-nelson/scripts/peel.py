"""Squeeze the tight core further, and look at what is left.

The UNSAT core converged at 91 points, but a core is only as small as the
proof the solver happened to find.  Greedy deletion -- drop a vertex, keep
the drop if the rest still refuses -- takes it to vertex-minimality, and at
91 points each refusal test is quick enough to afford one per vertex.

Then describe it: chromatic number, how many Moser spindles it carries, the
distances that appear, and whether the points sit in recognisable rings.
The object is whatever makes a unit-distance graph tight at its chromatic
number, and that property is the one missing at five.
"""
import sys, time, pickle, math
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
R = Fr(35, 9)
t0 = time.time()
P = pickle.load(open(SC + "tight_Sa_35_9.pkl", "rb"))
p, q = R.numerator, R.denominator


def build(pts):
    b = IntBasis.covering(pts)
    r = b.rows(pts)
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    return E


def refuses(pts):
    E = build(pts)
    n = len(pts)
    cl = [[1 + v * p + j for j in range(p)] for v in range(n)]
    for a, c in E:
        for j in range(p):
            for d in range(-(q - 1), q):
                cl.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
    cl += [[-(1 + j)] for j in range(1, p)] + [[1]]
    s = Solver(name="cd15", bootstrap_with=cl)
    ok = s.solve()
    s.delete()
    return not ok


def chrom(pts, k):
    E = build(pts)
    n = len(pts)
    cl = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cl.append([-(1 + a * k + col), -(1 + c * k + col)])
    s = Solver(name="cd15", bootstrap_with=cl)
    ok = s.solve()
    s.delete()
    return ok


print(f"start: {len(P)} points, {len(build(P))} edges"
      f"  [{time.time()-t0:.0f}s]", flush=True)
cur = list(P)
changed = True
while changed:
    changed = False
    for i in range(len(cur) - 1, -1, -1):
        trial = cur[:i] + cur[i + 1:]
        if len(trial) < 10:
            continue
        if refuses(trial):
            cur = trial
            changed = True
    print(f"   pass done: {len(cur)} points  [{time.time()-t0:.0f}s]",
          flush=True)
E = build(cur)
n = len(cur)
pickle.dump(cur, open(SC + "tightmin.pkl", "wb"))
xy = [(float(t.x), float(t.y)) for t in cur]
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
degs = sorted(len(adj[v]) for v in range(n))
tri = sum(len(adj[a] & adj[c]) for a, c in E) // 3
k3 = chrom(cur, 3)
k4 = chrom(cur, 4)
diam = max(math.hypot(a[0] - c[0], a[1] - c[1]) for a in xy for c in xy)
d2 = Counter(round((a[0] - c[0]) ** 2 + (a[1] - c[1]) ** 2, 6)
             for i, a in enumerate(xy) for c in xy[i + 1:])
print(f"\nvertex-minimal tight core: {n} points, {len(E)} edges"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"   degrees {degs[0]}..{degs[-1]} (mean {2*len(E)/n:.2f}), "
      f"{tri} triangles, diameter {diam:.3f}", flush=True)
print(f"   3-colourable: {k3}   4-colourable: {k4}   "
      f"(chi_c = 4 exactly, since it refuses {R})", flush=True)
print(f"   commonest squared distances: "
      f"{[(f'{k:.3f}', v) for k, v in d2.most_common(8)]}", flush=True)
print("DONE", flush=True)
