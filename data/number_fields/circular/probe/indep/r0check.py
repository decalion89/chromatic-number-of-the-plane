import sys, time
from fractions import Fraction as Fr
from sn import lift_levels, classify, make_refs
K = int(sys.argv[1])
for rs in sys.argv[2:]:
    r = Fr(rs)
    t = time.time()
    lev = lift_levels(K, r, keep_all=True)
    out = []
    for k in range(1, K + 1):
        refs = make_refs(k, r)
        junk = sum(1 for P in lev[k] if classify(P, k, r, refs) is None)
        out.append(f"k={k}: {len(lev[k])} pieces, junk {junk}")
    print(f"r = {r} = {float(r):.9f}: " + "; ".join(out) + f"  ({time.time()-t:.0f}s)")
    sys.stdout.flush()
