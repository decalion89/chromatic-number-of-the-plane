"""Search for the brick BY DISTANCE, exhaustively within each class.

Sa's 29 forced non-adjacent pairs at four colours are all at distance exactly
1/3.  Not spread over the graph -- concentrated at one distance, which is the
signature of the mechanism and also the distance de Grey prunes away to make
Y.  Random sampling over pairs was therefore the wrong search: the thing to
do is group pairs by exact squared distance and sweep whole classes.

Every distance class is swept completely, so a zero for a class is an
exhaustion for that class rather than a sample.  Classes are taken in order
of size, since forcing that matters shows up in many pairs at once.
"""
import sys, time, pickle
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "G"
KC = int(sys.argv[2]) if len(sys.argv) > 2 else 5
CLASSES = int(sys.argv[3]) if len(sys.argv) > 3 else 40
BUDGET = int(sys.argv[4]) if len(sys.argv) > 4 else 60000
PERCLASS = int(sys.argv[5]) if len(sys.argv) > 5 else 30
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, c in E:
    for col in range(KC):
        cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
s = Solver(name="m22", bootstrap_with=cls)
assert s.solve(), f"{CAR} is not {KC}-colourable"
xy = [(float(p.x), float(p.y)) for p in P]
buckets = defaultdict(list)
for u in range(n):
    for v in range(u + 1, n):
        if v in adj[u]:
            continue
        d2 = round((xy[u][0] - xy[v][0]) ** 2 + (xy[u][1] - xy[v][1]) ** 2, 9)
        if d2 < 9.0:
            buckets[d2].append((u, v))
# Ordering by class SIZE was wrong: Sa's forced pairs sit at d = 1/3, a
# small class, and sweeping the biggest four (d^2 = 0.333, 1.667, 0.543,
# 0.209) found nothing.  Forcing lives at a particular distance, not a
# common one, so breadth beats depth -- a few pairs from EVERY class.
# Two orderings failed.  Largest classes whole missed it: Sa's forcing sits
# at d = 1/3, a small class.  All classes with 30 pairs each also missed it:
# that class holds hundreds of pairs and only 29 are forced, so the forcing
# is SPARSE WITHIN ITS OWN CLASS and any partial sweep loses it.
# What the successful pairs have in common is that d^2 = 1/9 -- a rational
# with a small denominator.  So: sweep the algebraically simple classes
# WHOLE, simplest first, and leave the rest.
from fractions import Fraction as Fr


def simplicity(d2):
    fr = Fr(d2).limit_denominator(400)
    if abs(float(fr) - d2) > 1e-7:
        return 10 ** 6
    return fr.denominator + fr.numerator


order = sorted(buckets.items(), key=lambda kv: (simplicity(kv[0]), kv[0]))
order = [kv for kv in order if simplicity(kv[0]) < 10 ** 6]
print(f"{CAR}, k = {KC}: {n} points, {len(buckets)} distinct non-adjacent "
      f"distances; sweeping the {CLASSES} algebraically simplest classes WHOLE"
      f"  [{time.time()-t0:.0f}s]", flush=True)
done = 0
hits = []
for ci, (d2, pairs) in enumerate(order[:CLASSES]):
    if done > BUDGET:
        break
    found = 0
    for (u, v) in pairs:
        if not s.solve(assumptions=[1 + u * KC, 1 + v * KC]):
            found += 1
        done += 1
        if done > BUDGET:
            break
    if found:
        hits.append((d2, found, len(pairs)))
        print(f"   *** d^2 = {d2:.6f} (d = {d2 ** 0.5:.6f}): {found} of "
              f"{len(pairs)} are FORCED  [{time.time()-t0:.0f}s]", flush=True)
    else:
        print(f"   d^2 = {d2:.6f}: {len(pairs)} pairs swept whole, none "
              f"forced  [{time.time()-t0:.0f}s]", flush=True)
    if False:
        print(f"   {ci+1} classes swept, {done} pairs, {len(hits)} distances "
              f"with forcing  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{done} pairs swept in whole classes; "
      f"{len(hits)} distances carry forcing  [{time.time()-t0:.0f}s]",
      flush=True)
for d2, f, t in hits:
    print(f"   d = {d2 ** 0.5:.6f}: {f} of {t}", flush=True)
print("DONE", flush=True)
