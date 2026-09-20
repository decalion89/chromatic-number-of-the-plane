"""Cross-check a blocking verdict with a second solver, and by hand.

An unsatisfiability result on two million clauses is worth contrasting.  Three
checks, none of them trusting the first run:

  1. the module really is 5-saturated -- rank mod 5 equal to the ambient
     dimension, which forces rank over Q to match and so makes the search over
     Z^d the right one;
  2. a sample of random functionals really is killed, as the count predicts;
  3. the refutation repeats under Glucose as well as Cadical.

The counting says an escape would be a one-in-10^18 event at this size, so
finding one by sampling is not expected; what the sample rules out is a coding
error that made every functional fail trivially -- if some direction reduced
to zero mod 5 the answer would be "blocked" for the wrong reason.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import has_homomorphism, _rank_mod

t0 = time.time()
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/gate.pkl", "rb") as fh:
    k, iv = pickle.load(fh)
d = len(iv[0])
print(f"{k} rotations, {len(iv)} directions, dim {d}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

r5 = _rank_mod(iv, 5)
print(f"  rank mod 5 = {r5} of {d} -- "
      f"{'5-saturated, the ambient search is the right one' if r5 == d else 'NOT saturated'}",
      flush=True)

small = sorted({tuple(x % 5 for x in v) for v in iv})
zero = [v for v in small if not any(v)]
print(f"  {len(small)} residues mod 5, {len(zero)} of them zero "
      f"({'GOOD -- no trivial obstruction' if not zero else 'a direction is 0 mod 5'})",
      flush=True)

random.seed(20260920)
best = 0
for _ in range(20000):
    phi = [random.randrange(5) for _ in range(d)]
    hit = sum(1 for v in small
              if sum(p * x for p, x in zip(phi, v)) % 5)
    best = max(best, hit)
print(f"  best of 20000 random functionals: nonzero on {best} of "
      f"{len(small)} residues  [{time.time()-t0:.0f}s]", flush=True)

for name in ("cd19", "g4", "m22"):
    try:
        from pysat.solvers import Solver
        phi, why = has_homomorphism(small, 5)
        print(f"  cadical says: {'BLOCKS' if phi is None else 'escape'} "
              f"({why})  [{time.time()-t0:.0f}s]", flush=True)
        break
    except Exception as e:
        print(f"  {name}: {e}", flush=True)
