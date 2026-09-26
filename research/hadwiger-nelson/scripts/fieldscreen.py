"""Screen a multiquadratic CM field K = Q(sqrt-d_1, ..., sqrt-d_n) by its non-split places.

A place v of K+ that does not split in K colours the unit-distance graph on all of K:
  * ramified in K/K+  ->  chi(K) <= 3 (the norm-one residues are +-1);
  * unramified, residue field F_q of K+_v  ->  chi(K) <= chi(Cay(F_{q^2}, N1)), the finite plane.
(hn/adelic.py, notes/rigidity.md sections 8-9.)  Local Galois groups from Kummer theory: the
decomposition group D_p is the image of Hom(Q_p^*/Q_p^*2, +-1) acting on the square classes of the
-d_i, inertia I_p the image of the characters trivial on the unramified unit class; complex
conjugation c = (1, ..., 1).  v non-split iff c in D_p, ramified iff c in I_p.

usage: python3 fieldscreen.py d1,d2,... [bound]"""
import sys, itertools
from sympy import primerange, factorint

def square_class(a, p):
    """Coordinates of a in Q_p^*/Q_p^*2 over F_2: odd p -> (v mod 2, unit non-residue);
    p = 2 -> (v mod 2, unit = -1 class, unit = 5 class) with units mod 8 = (-1)^x 5^y."""
    v = 0
    while a % p == 0:
        a //= p; v += 1
    if p != 2:
        return (v % 2, 0 if pow(a % p, (p - 1) // 2, p) == 1 else 1)
    u = a % 8
    x = 1 if u in (3, 7) else 0            # 3 = -5, 7 = -1
    y = 1 if u in (3, 5) else 0
    return (v % 2, x, y)

def local_groups(ds, p):
    S = [square_class(-d, p) for d in ds]
    m = len(S[0])
    unram = [0, 1] if p != 2 else [0, 2]   # index of the unramified unit class in the coordinates
    # a character tau of the square-class group is a vector in F_2^m; sigma_i = <tau, s_i>
    D, I = set(), set()
    for tau in itertools.product((0, 1), repeat=m):
        sig = tuple(sum(t * s for t, s in zip(tau, si)) % 2 for si in S)
        D.add(sig)
        # inertia: tau kills the unramified unit class (odd p: coordinate 1; p = 2: coordinate 2)
        if (p != 2 and tau[1] == 0) or (p == 2 and tau[2] == 0):
            I.add(sig)
    return D, I

def screen(ds, bound=400):
    n = len(ds); c = tuple([1] * n)
    out = []
    for p in primerange(2, bound):
        D, I = local_groups(ds, p)
        if c not in D:
            continue
        e = len(I); f = len(D) // e
        if c in I:
            out.append((p, "ramified", None))
        else:
            out.append((p, "unramified", p ** (f // 2)))
    return out

if __name__ == "__main__":
    ds = [int(x) for x in sys.argv[1].split(",")]
    if any(d <= 0 for d in ds):
        sys.exit("each d_i must be a positive integer, for the field Q(sqrt-d_1, ..., sqrt-d_n); write a "
                 "real generator sqrt(m) as sqrt-(3m) next to sqrt-3, for example")
    bound = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    res = screen(ds, bound)
    print(f"K = Q({', '.join('sqrt-%d' % d for d in ds)}): non-split places of K+ below {bound}:")
    for p, kind, q in res:
        print(f"  p = {p}: {kind}" + (f", residue field of K+_v: F_{q}" if q else "  -> chi(K) <= 3"))
