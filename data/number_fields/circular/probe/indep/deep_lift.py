"""Lift S_{5^k}^(r) level by level up to K (keeping one level), print piece/junk counts per level."""
import sys, time
from fractions import Fraction as Fr
from hgeom import funcs
from sn import square, classify, make_refs

def run(K, r):
    lo, hi = r, 1 - r
    cur = [square(r)]
    for k in range(1, K + 1):
        t = time.time()
        Np = 5 ** (k - 1)
        fs = [f for j in (k, -k) for f in funcs(j)]
        new = []
        for P in cur:
            for a in range(5):
                for b in range(5):
                    pieces = [P.translate((Fr(Np * a), Fr(Np * b)))]
                    for f in fs:
                        nxt = []
                        for R in pieces:
                            nxt.extend(R.strip_split(f, lo, hi))
                        pieces = nxt
                        if not pieces:
                            break
                    new.extend(pieces)
        cur = new
        refs = make_refs(k, r)
        junk = [P for P in cur if classify(P, k, r, refs) is None]
        print(f"  r={float(r):.9f} k={k}: {len(cur)} pieces, junk {len(junk)} ({sum(1 for P in junk if len(P.V)==1)} points)  [{time.time()-t:.1f}s]")
        sys.stdout.flush()
    return cur

K = int(sys.argv[1])
for rs in sys.argv[2:]:
    run(K, Fr(rs))
