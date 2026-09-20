"""Check the ramified claim instead of asserting it.

`cyclotomic_can_block` returns "cannot block" whenever 5 divides n, on the
grounds that 5 ramifies there.  But the orbit classification enumerated the
UNRAMIFIED types, so that return value was an assertion, not a computation --
only Q(zeta_5) had actually been decided.  It matters, because de Grey's field
contains sqrt5 and so has 5 ramified.

Decided here by the same exhaustion over A = O/5 that settled the unramified
types: build N = {x : x.sigma(x) = 1} and ask whether every F_5-hyperplane
meets it.
"""
import sys, itertools, time
from fractions import Fraction
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from hn.cyclotomic import CycloField

P = 5


def decide(name, K):
    t0 = time.time()
    d = K.degree

    def red(x):
        out = []
        for t in x:
            if t.denominator % P == 0:
                return None
            out.append(int(t.numerator * pow(t.denominator, -1, P) % P))
        return tuple(out)

    def e(i):
        return tuple(Fraction(1 if j == i else 0) for j in range(d))

    C = [[red(K.mul(e(i), e(j))) for j in range(d)] for i in range(d)]
    S = [red(K.conj(e(i))) for i in range(d)]
    if any(r is None for row in C for r in row) or any(r is None for r in S):
        print(f"{name}: basis not integral at 5", flush=True)
        return None
    one = tuple([1] + [0] * (d - 1))

    def mul(a, b):
        out = [0] * d
        for i in range(d):
            if a[i]:
                ai = a[i]
                for j in range(d):
                    if b[j]:
                        f = ai * b[j]
                        row = C[i][j]
                        for k in range(d):
                            if row[k]:
                                out[k] += f * row[k]
        return tuple(x % P for x in out)

    def conj(a):
        out = [0] * d
        for i in range(d):
            if a[i]:
                row = S[i]
                for k in range(d):
                    out[k] += a[i] * row[k]
        return tuple(x % P for x in out)

    N = [a for a in itertools.product(range(P), repeat=d)
         if mul(a, conj(a)) == one]
    seen, hyps = set(), []
    for c in itertools.product(range(P), repeat=d):
        if not any(c) or c in seen:
            continue
        for k in range(1, P):
            seen.add(tuple((k * x) % P for x in c))
        hyps.append(c)
    missed = 0
    for c in hyps:
        if not any(sum(x * y for x, y in zip(c, u)) % P == 0 for u in N):
            missed += 1
    ok = missed == 0
    print(f"{name}: degree {d}, |N| = {len(N)}, {len(hyps)} hyperplanes -> "
          + ("CAN BLOCK" if ok else f"cannot block ({missed} miss)")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
    return ok


for n in (5, 15, 20, 25):
    decide(f"Q(zeta_{n}), 5 ramified", CycloField(n))
