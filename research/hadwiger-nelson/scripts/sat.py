import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/fastblock.py").read()
     .split('print(f"necklace')[0])
from hn.homcol import saturated_at, _rank_q, _rank_mod, has_homomorphism
t1 = time.time()
rq = _rank_q(ints)
print(f"necklace: {len(ints)} directions, rank_Q {rq}  "
      f"[{time.time()-t1:.1f}s]", flush=True)
for p in (2, 3, 5):
    print(f"  rank mod {p} = {_rank_mod(ints, p)}  [{time.time()-t1:.1f}s]",
          flush=True)
for n in (2, 3, 4, 5):
    sat = saturated_at(ints, n)
    small = sorted({tuple(x % n for x in v) for v in ints})
    zero = any(not any(v) for v in small)
    amb = True if zero else has_homomorphism(small, n)[0] is None
    print(f"  n = {n}: saturated {sat}, ambient test says "
          f"{'blocks' if amb else 'escapes'}"
          f"{'' if sat else '  <-- needs the module test'}  "
          f"[{time.time()-t1:.1f}s]", flush=True)
