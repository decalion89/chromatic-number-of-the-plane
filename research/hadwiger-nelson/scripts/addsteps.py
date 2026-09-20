"""Close under rotations for the chromatic number, translate for the blocking.

The closure of the Moser spindle under the six rotations of
K = Q(m, sqrt-3, sqrt-11) reaches 597 points and chi = 4, and blocks at 2, 3
and 4 -- but not at 5, with 594 directions, where the necklace blocked at 5
with 138.  So it is not a matter of how many directions but which.

The field can block: 5 is inert in Q(m) and in Q(sqrt33), so its residue
degree in F is 6, far above the bound of 3.  What is missing is steps of the
right kind, and those come from Hilbert 90: every unit step is a/abar.

Adding one as a TRANSLATION -- X u (X + w) -- keeps chi >= 4, since X is
still there, and brings w's direction into the set.  So translate by the steps
that kill the escaping homomorphisms, one at a time, until nothing escapes.
"""
import sys, time, itertools
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.homcol import has_homomorphism

t0 = time.time()


def kinv(x):
    """Inverse in K, by solving the 12x12 multiplication system."""
    basis = []
    for i in range(4):
        for j in range(3):
            e = [list(LZ) for _ in range(4)]
            e[i][j] = Fr(1)
            basis.append(tuple(tuple(t) for t in e))
    cols = [kmul(x, b) for b in basis]

    def flat(v):
        out = []
        for part in v:
            out.extend(part)
        return out

    A = [[flat(cols[j])[i] for j in range(12)]
         + [Fr(1 if i == 0 else 0)] for i in range(12)]
    for c in range(12):
        p = next(t for t in range(c, 12) if A[t][c])
        A[c], A[p] = A[p], A[c]
        sc = Fr(1) / A[c][c]
        A[c] = [v * sc for v in A[c]]
        for t in range(12):
            if t != c and A[t][c]:
                f = A[t][c]
                A[t] = [u - f * v for u, v in zip(A[t], A[c])]
    co = [A[i][12] for i in range(12)]
    out = []
    for i in range(4):
        out.append(tuple(co[3 * i:3 * i + 3]))
    return tuple(out)


seeds = [K1, W, S, LMK := kmul(K1, K1)]
seeds = [K1, W, S, (LM, LZ, LZ, LZ), kadd(K1, W), kadd(K1, S),
         kadd((LM, LZ, LZ, LZ), K1), ksub((LM, LZ, LZ, LZ), K1),
         kadd(W, S), kadd(kmul(krat(2), W), K1)]
steps = set()
for a in seeds:
    for b in seeds:
        for z in (kadd(a, b), ksub(a, b), kmul(a, b), kadd(a, kmul(krat(2), b))):
            if z == kz():
                continue
            try:
                u = kmul(z, kinv(kconj(z)))
            except (StopIteration, ZeroDivisionError):
                continue
            if knorm2(u) == K1:
                v = u
                for _ in range(6):
                    steps.add(v)
                    v = kmul(v, Z6)
steps = sorted(steps)
print(f"{len(steps)} unit steps of K by Hilbert 90  [{time.time()-t0:.0f}s]",
      flush=True)


def dirs_of(vecs):
    dn = 1
    for v in vecs:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    out = set()
    for v in vecs:
        w = tuple(int(q * dn) for q in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w) if g > 1 else w)
    return sorted(out)


def flat(x):
    out = []
    for part in x:
        out.extend(part)
    return tuple(out)


rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
X = close(seed, [r for _, r in rots])
print(f"closure: {len(X)} points  [{time.time()-t0:.0f}s]", flush=True)

base = set()
for p in X:
    for q in X:
        if p != q and knorm2(ksub(q, p)) == K1:
            base.add(flat(ksub(q, p)))
cur = list(base)
chosen = []
for rnd in range(40):
    iv = dirs_of(cur)
    blocks = [n for n in (2, 3, 4, 5)
              if has_homomorphism(iv, n)[0] is None]
    print(f"  round {rnd}: {len(chosen)} translations, {len(iv)} directions, "
          f"blocks at {blocks}  [{time.time()-t0:.0f}s]", flush=True)
    if 5 in blocks and len(blocks) == 4:
        print("  *** blocks at 2, 3, 4 and 5 ***", flush=True)
        break
    phi, _ = has_homomorphism(iv, 5)
    if phi is None:
        phi, _ = has_homomorphism(iv, 4)
    best, bk = None, -1
    for w in steps:
        if w in chosen:
            continue
        fw = flat(w)
        d = dirs_of([fw])[0]
        k = 1 if sum(x * y for x, y in zip(phi, d)) % 5 == 0 else 0
        if k > bk:
            best, bk = w, k
            if k:
                break
    if best is None:
        print("  no step left", flush=True)
        break
    chosen.append(best)
    cur.append(flat(best))
