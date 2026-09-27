"""Exact box subdivision for the circular gate in low rank.

F = { phi in T^r : frac(C_u . phi) in [1/k, 1 - 1/k] for every unit u }.  Boxes of the unit cube (dyadic,
exact integer arithmetic: a box is [a_i, a_i + 1] * 2^-L per coordinate) are
  * pruned when some unit maps the whole box into the open bad arc (-1/k, 1/k) + Z,
  * certified feasible when every unit maps the whole box into [1/k, 1 - 1/k] + Z (then F has
    interior: a circular colouring exists),
  * split otherwise (longest side first).
If the stack empties with no certified box and no undecided leaf, F is empty on the given units:
an exact proof that the module admits no circular k-colouring along those directions.
Leaves below a depth limit are reported as undecided (F may still be empty, or contain
lower-dimensional pieces such as coset points).

usage: python3 boxgate.py <module.json> [k] [max depth] [max boxes]"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, os, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
src = open(os.path.join(HN_DIR, "scripts", "circgate.py")).read()
pre = src[:src.index("# optional: also require the multiples")]
ARGS = list(sys.argv)
sys.argv = ["circgate.py", ARGS[1], "1"]
exec(pre)
KC = int(ARGS[2]) if len(ARGS) > 2 else 5
MAXD = int(ARGS[3]) if len(ARGS) > 3 else 40
MAXB = int(ARGS[4]) if len(ARGS) > 4 else 20_000_000
reps, seenr = [], set()
for i in range(len(CU)):
    k_ = tuple(CU[i]); nk_ = tuple(-CU[i])
    if nk_ in seenr or k_ in seenr: continue
    seenr.add(k_); reps.append(i)
C = [[int(x) for x in CU[i]] for i in reps]; m = len(C)
Cp = [[max(c, 0) for c in row] for row in C]; Cn = [[min(c, 0) for c in row] for row in C]
order = sorted(range(m), key=lambda j: sum(abs(c) for c in C[j]))      # short units first: they prune early
print(f"{os.path.basename(ARGS[1])}: {m} directions, rank {r}, k = {KC}", flush=True)
# level-by-level (breadth first), vectorised over the boxes of a level with numpy int64 arithmetic
Cm = np.array(C, dtype=np.int64); CPs = np.maximum(Cm, 0).sum(axis=1); CNs = np.minimum(Cm, 0).sum(axis=1)
alive = np.zeros((1, r), dtype=np.int64)       # corners at level L (box = [a, a+1] / 2^L)
pruned = certified = undecided = 0; boxes = 0
for L in range(0, MAXD + 1):
    S = 1 << L
    if len(alive) == 0: break
    base = alive @ Cm.T                                        # (boxes x units)
    lo = (base + CNs) * KC; hi = (base + CPs) * KC
    n = np.floor_divide(lo + S, KC * S)
    bad = (lo > KC * S * n - S) & (hi < KC * S * n + S)
    n2 = np.floor_divide(lo - S, KC * S)
    good = (lo >= KC * S * n2 + S) & (hi <= KC * S * n2 + (KC - 1) * S)
    is_pruned = bad.any(axis=1); is_good = good.all(axis=1) & ~is_pruned
    boxes += len(alive); pruned += int(is_pruned.sum())
    if is_good.any():
        certified += 1
        print(f"  CERTIFIED feasible box at level {L}: corner {tuple(alive[np.argmax(is_good)])} / 2^{L}", flush=True)
        break
    keep = alive[~is_pruned]
    print(f"  level {L:2d}: {len(alive):>10} boxes, {int(is_pruned.sum()):>10} pruned, {len(keep):>9} alive   [{time.time()-t0:.0f}s]", flush=True)
    if L == MAXD: undecided = len(keep); break
    if len(keep) * (1 << r) > MAXB: print("  box budget exhausted", flush=True); undecided = len(keep); break
    kids = [2 * keep + np.array([(bits >> i) & 1 for i in range(r)], dtype=np.int64) for bits in range(1 << r)]
    alive = np.concatenate(kids)
    if np.abs(alive).max() * np.abs(Cm).sum(axis=1).max() * KC > 2 ** 62: print("  int64 range reached"); undecided = len(keep); break
stack = []
print(f"boxes {boxes}; pruned {pruned}; certified {certified}; undecided leaves (depth {MAXD}) {undecided}; stack left {len(stack)}   [{time.time()-t0:.0f}s]")
if certified == 0 and undecided == 0 and not stack: print("==> EXACT: no circular colouring along these units (F is empty)")
