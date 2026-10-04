"""Like test_new.py, but V also contains vectors that are p-adically large at the split places above 2 and 3
(p close to sqrt(-1) there), which the proof needs (the torus is non-compact at a split place)."""
import sys, random
from fractions import Fraction as Fr
from nf import Field, unit_from, lconj
from limit_test import feasible

def sqrt_mod(a, p, k):
    """x with x^2 = a mod p^k (p odd, or p = 2 with a = 1 mod 8)"""
    m = p ** k
    for x in range(m):
        if (x * x - a) % m == 0:
            return x
    raise ValueError

name = sys.argv[1]; tries = int(sys.argv[2]); nrand = int(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed)
DENS = [1, 1, 2, 5] if len(sys.argv) > 5 and sys.argv[5] == 'int3' else [1, 1, 2, 3, 5]
if name == "3_5":
    Fd = Field([[-3, 0, 1], [-5, 0, 1]])
    s3, s5 = Fd.gen(0), Fd.gen(1); s15 = Fd.mul(s3, s5)
    # at 3: i = sqrt(-5)/sqrt5 = r sqrt5 / 5 with r^2 = -5 (3-adic); at 2: i = sqrt(-15)/sqrt15 = r sqrt15/15, r^2 = -15
    special = []
    for k in (2, 3, 4):
        r = sqrt_mod(-5, 3, k); special.append(Fd.scal(Fr(r, 5), s5))
    for k in (4, 5, 6):
        r = sqrt_mod(-15, 2, k); special.append(Fd.scal(Fr(r, 15), s15))
elif name == "c7_7":
    Fd = Field([[-1, -2, 1, 1], [-7, 0, 1]])
    s7 = Fd.gen(1)
    special = []
    for k in (4, 5, 6):
        r = sqrt_mod(-7, 2, k); special.append(Fd.scal(Fr(r, 7), s7))   # i = sqrt(-7)/sqrt7 at 2... up to sign
elif name == "c7_2":
    Fd = Field([[-1, -2, 1, 1], [-2, 0, 1]])
    special = []
elif name == "2_7":
    Fd = Field([[-2, 0, 1], [-7, 0, 1]])
    s2, s7 = Fd.gen(0), Fd.gen(1)
    special = []
    for k in (2, 3, 4):
        r = sqrt_mod(-2, 3, k); special.append(Fd.scal(Fr(r, 2), s2))     # i = sqrt(-2)/sqrt2 at 3
    for k in (4, 5, 6):
        r = sqrt_mod(-7, 2, k); special.append(Fd.scal(Fr(r, 7), s7))     # i = sqrt(-7)/sqrt7 at 2
elif name == "3_7":
    Fd = Field([[-3, 0, 1], [-7, 0, 1]])
    s7 = Fd.gen(1)
    special = []
    for k in (4, 5, 6):
        r = sqrt_mod(-7, 2, k); special.append(Fd.scal(Fr(r, 7), s7))
else:
    raise SystemExit
basis = Fd.basis
def rand_elt():
    a = {}
    for b in basis:
        if rng.random() < 0.6:
            c = Fr(rng.randint(-4, 4), rng.choice(DENS))
            if c: a[b] = c
    return a
pool = []
while len(pool) < 120:
    try:
        pool.append(unit_from(Fd, rand_elt()))
    except Exception:
        pass
spec = [unit_from(Fd, p) for p in special]
found = 0
for t in range(tries):
    V = [(Fd.one(), {})] + spec + rng.sample(pool, nrand)
    ok, pat = feasible(Fd, V, want=True)
    print(f"try {t}: |V| = {len(V)}: {'feasible ' + pat if ok else 'INFEASIBLE'}", flush=True)
    if not ok:
        found += 1
        import json
        enc = lambda a: {",".join(map(str, k)): str(v) for k, v in a.items()}
        json.dump({"field": name, "polys": [[str(c) for c in pl] for pl in Fd.polys],
                   "V": [[enc(X), enc(Y)] for (X, Y) in V]}, open(f"V_{name}_{seed}_{t}.json", "w"))
        if found >= 2: break
print("done", name, "infeasible:", found)
