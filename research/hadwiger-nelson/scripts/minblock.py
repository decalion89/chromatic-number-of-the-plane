"""The smallest set of unit directions that blocks every coset 5-colouring.

Blocking is exactly a covering problem in the dual.  A direction d kills the
maps phi with phi(d) = 0, which is the hyperplane d-perp, so a step set blocks
iff its hyperplanes cover all 3906 points of PG(5,5).  Each hyperplane holds
781 of them, so five is the counting floor; a projective line of six would do
it in principle, but unit vectors all share one norm and no line over F_5 has
all six points in one square class -- so seven is the real floor here.

Greedy covering finds a small one and a local exchange pass tries to shrink it
further, which turns the abstract obstruction into an explicit short list of
steps that a construction can actually carry.
"""
import sys, json, itertools, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import has_homomorphism

S = [tuple(v) for v in json.load(open(
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/blocking.json"))]
print(f"{len(S)} candidate directions (denominator <= 29)")

D = 6
# Every projective point of the dual, one representative per class.
pts = []
for v in itertools.product(range(5), repeat=D):
    if any(v) and next(x for x in v if x) == 1:
        pts.append(v)
print(f"{len(pts)} points of PG(5,5) to cover")

cover = {}
for d in S:
    r = tuple(x % 5 for x in d)
    cover[d] = frozenset(i for i, p in enumerate(pts)
                         if sum(a * b for a, b in zip(r, p)) % 5 == 0)
sizes = {len(c) for c in cover.values()}
print(f"each direction covers {sorted(sizes)} points\n")

full = frozenset(range(len(pts)))


def greedy(order):
    got, chosen = set(), []
    for d in order:
        if cover[d] - got:
            chosen.append(d)
            got |= cover[d]
            if len(got) == len(pts):
                return chosen
    return None


best = None
rng = random.Random(7)
for trial in range(400):
    order = sorted(S, key=lambda d: -len(cover[d])) if trial == 0 else S[:]
    if trial:
        rng.shuffle(order)
    ch = greedy(order)
    if ch and (best is None or len(ch) < len(best)):
        # drop any direction the rest already covers
        changed = True
        while changed:
            changed = False
            for d in list(ch):
                if frozenset().union(*[cover[e] for e in ch if e != d]) == full:
                    ch.remove(d)
                    changed = True
        best = ch
        print(f"  trial {trial}: {len(ch)} directions", flush=True)

print(f"\nminimum found: {len(best)} directions")
phi, why = has_homomorphism(list(best), 5)
print(f"independent check with the SAT screen: "
      f"{'STILL 5-COLOURABLE (bug)' if phi else 'NO homomorphism -- blocks'}")
for d in best:
    print("   ", d)
json.dump([list(d) for d in best], open(
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/minblock.json", "w"))
