"""Search small finite quotients, where the coordinate-wise ones were too big.

The homomorphism-to-Z_5 test ruled out one family; the coordinate-wise
quotients (Z_m)^r would rule out far more, but at rank eight the useful moduli
are 5^8 = 390625 and 7^8 = 5.7 million cosets, and the small ones all have a
generator vanishing.  So the test never actually ran for the two
configurations that matter.

Any finite abelian quotient will do, and it need not be (Z_m)^r.  Pick a small
group A -- Z_n, or Z_a x Z_b, a few dozen elements -- and a homomorphism
phi : Z^r -> A given by one image per basis vector.  If no generator maps to
zero, 5-colour Cay(A, phi(U)); it has |A| vertices, so the call is instant.
A colouring there lifts to the WHOLE infinite group, which is a proof.

There are |A|^r homomorphisms, far too many to enumerate at rank eight, so
they are sampled.  A hit is a proof; a miss is a miss, and is reported as one.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from quotient import units, coords
from pysat.solvers import Solver

k = 5
t0 = time.time()
rng = random.Random(20240921)


def colourable(A, imgs):
    """Is Cay(A, imgs) 5-colourable?  A is a list of moduli."""
    size = 1
    for a in A:
        size *= a
    pw, acc = [], 1
    for a in A:
        pw.append(acc)
        acc *= a
    def code(v):
        return sum((x % a) * p for x, a, p in zip(v, A, pw))
    E = set()
    for v in itertools.product(*[range(a) for a in A]):
        i = code(v)
        for g in imgs:
            j = code([x + y for x, y in zip(v, g)])
            if i != j:
                E.add((min(i, j), max(i, j)))
    cls = [[1 + v * k + c for c in range(k)] for v in range(size)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok, size, len(E)


SHAPES = [(7,), (8,), (9,), (11,), (13,), (16,), (5, 5), (4, 6), (3, 9),
          (5, 7), (6, 6), (4, 4, 3), (3, 3, 5)]
for rots in ((3, 4), (3, 4, 7)):
    U = units(rots, 1)
    r, ints = coords(U)
    print(f"\nrotations {rots}: {len(U)} unit vectors, rank {r}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    found = False
    for A in SHAPES:
        tries = 0
        for _ in range(400):
            phi = [[rng.randrange(a) for a in A] for _ in range(r)]
            imgs = []
            dead = False
            for co in ints:
                g = [sum(c * phi[t][s] for t, c in enumerate(co)) % A[s]
                     for s in range(len(A))]
                if not any(g):
                    dead = True
                    break
                imgs.append(tuple(g))
            if dead:
                continue
            tries += 1
            ok, size, e = colourable(A, sorted(set(imgs)))
            if ok:
                print(f"   A={A}: a homomorphism works -- {size} cosets, "
                      f"{e} edges, 5-colourable, so the WHOLE GROUP is"
                      f"  [{time.time()-t0:.0f}s]", flush=True)
                found = True
                break
        if found:
            break
        print(f"   A={A}: {tries} usable homomorphisms of 400 sampled, none "
              f"gave a 5-colourable quotient  [{time.time()-t0:.0f}s]",
              flush=True)
    if not found:
        print(f"   no small quotient found for {rots} -- a miss, not a proof",
              flush=True)
print("\nDONE", flush=True)
