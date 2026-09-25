"""Bracket rho(G,5) from below: small budgets decide fast.

At budget 63 the loop sits on the boundary -- nineteen thousand classes
collected and a 63-set still hitting them all.  Smaller budgets are strictly
easier to REFUTE: the fewer vertices allowed, the sooner no hitting set
exists, and "no hitting set for this subfamily" settles rho > budget outright.

So sweeping upward from a small budget gives a genuine lower bound on rho
long before the 63 question resolves, and each refutation is cheap.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from hn import degrey
from is63 import decide

g = degrey.build_G()
print(f"G: {g.n} vertices, {g.m} edges; floor is rho >= 5", flush=True)
for budget in (5, 10, 20, 30, 40, 50, 60):
    ok, info = decide(g, 5, budget, rounds=200000, report=10 ** 9)
    if ok is False:
        print(f"  budget {budget}: NO hitting set -- rho > {budget}",
              flush=True)
    elif ok is True:
        print(f"  budget {budget}: a forcing set of {len(info)} exists -- "
              f"rho <= {budget}   *** and {budget} <= 63 ***", flush=True)
        break
    else:
        print(f"  budget {budget}: undecided after 200000 rounds", flush=True)
        break
