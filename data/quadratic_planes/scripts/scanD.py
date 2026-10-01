import sys
from units_fast import units_fast
def smooth(n, P=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)):
    for p in P:
        while n % p == 0:
            n //= p
    return n == 1
ds = [int(x) for x in sys.argv[1].split(",")]
Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
for d in ds:
    res = []
    for D in range(2, Dmax + 1):
        if not smooth(D):
            continue
        n = len(units_fast(d, D))
        res.append((n, D))
    good = sorted([(n, D) for n, D in res if 80 <= n <= 260], key=lambda x: (x[1]))
    best = sorted(res)[-5:]
    print(f"d={d}: D with 80..260 units: {good[:14]}  largest: {best}", flush=True)
