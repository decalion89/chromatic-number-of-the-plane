"""gatescan.py d,d,... DMAX PMAX: for each d, the denominators D <= DMAX (multiples of 6, prime factors <= PMAX)
with 40 to 400 directions and every gate open (gates2.report), sorted by the number of directions."""
import sys
import numpy as np
from units_fast import units_fast
from gates2 import report
ds = [int(x) for x in sys.argv[1].split(",")]; DMAX, PMAX = int(sys.argv[2]), int(sys.argv[3])
def smooth(n, pmax):
    for p in range(2, pmax + 1):
        while n % p == 0:
            n //= p
    return n == 1
for d in ds:
    ok = []
    for D in range(6, DMAX + 1, 6):
        if not smooth(D, PMAX):
            continue
        U = np.asarray(units_fast(d, D), dtype=np.int64).reshape(-1, 4)
        if len(U) < 40 or len(U) > 400:
            continue
        r = report(d, D)
        if any(r.get(g) == 0 for g in ("2", "3", "5")):
            continue
        ok.append((len(U), D))
    ok.sort(reverse=True)
    print(f"d={d}: {len(ok)} denominators with every gate open; most directions: {ok[:12]}", flush=True)
