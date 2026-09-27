"""Every arithmetic type, decided once.

A = O/5 depends only on how 5 splits.  Write K for the edge field, F for its
maximal real subfield, sigma for complex conjugation.  When 5 is unramified,

    A = prod_i F_{5^{f_i}},    sum f_i = [K:Q],

and sigma permutes the factors with order 2.  Each of its orbits is either a
single factor of EVEN degree f, on which sigma must be the unique involution
Frob^{f/2}, or a PAIR of factors of equal degree swapped.  So the type is a
partition of the degree into such orbits -- finitely many per degree.

N = {x : x.sigma(x) = 1} is determined by the type: a fixed factor of degree f
contributes the norm-one subgroup of F_{5^f} over F_{5^{f/2}}, of order
5^{f/2} + 1, and a swapped pair contributes {(a, a^{-1})}, of order 5^f - 1.

Blocking over K needs every hyperplane of A to meet N.  Deciding every type of
a given degree therefore decides EVERY field of that degree at once.
"""
import sys, itertools, time

P = 5


def irreducible(f):
    """A monic irreducible of degree f over F_5, as a coefficient list."""
    if f == 1:
        return [0, 1]
    for tail in itertools.product(range(P), repeat=f):
        poly = list(tail) + [1]
        roots_free = True
        # trial-divide by all monic polys of degree 1..f//2
        for g in range(1, f // 2 + 1):
            for t2 in itertools.product(range(P), repeat=g):
                div = list(t2) + [1]
                r = poly[:]
                for i in range(len(r) - 1, g - 1, -1):
                    c = r[i]
                    if c:
                        for j in range(g + 1):
                            r[i - g + j] = (r[i - g + j] - c * div[j]) % P
                if not any(r[:g]):
                    roots_free = False
                    break
            if not roots_free:
                break
        if roots_free:
            return poly
    raise RuntimeError(f)


def fq(f):
    """Multiplication and Frobenius on F_{5^f} in the basis 1, t, .., t^{f-1}."""
    poly = irreducible(f)

    def mul(a, b):
        r = [0] * (2 * f - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        r[i + j] = (r[i + j] + x * y) % P
        for i in range(2 * f - 2, f - 1, -1):
            c = r[i]
            if c:
                for j in range(f + 1):
                    r[i - f + j] = (r[i - f + j] - c * poly[j]) % P
        return tuple(r[:f])

    one = tuple([1] + [0] * (f - 1))

    def power(a, e):
        out, base = one, a
        while e:
            if e & 1:
                out = mul(out, base)
            base = mul(base, base)
            e >>= 1
        return out

    return mul, one, power


def decide(orbits, show=True):
    """orbits: list of ('fix', f) or ('swap', f)."""
    t0 = time.time()
    parts, d = [], 0
    for kind, f in orbits:
        mul, one, power = fq(f)
        parts.append((kind, f, mul, one, power))
        d += f if kind == "fix" else 2 * f
    if d > 8:
        return None

    def split(vec):
        out, i = [], 0
        for kind, f, _, _, _ in parts:
            if kind == "fix":
                out.append((tuple(vec[i:i + f]),)); i += f
            else:
                out.append((tuple(vec[i:i + f]), tuple(vec[i + f:i + 2 * f])))
                i += 2 * f
        return out

    def is_norm_one(vec):
        for (kind, f, mul, one, power), comp in zip(parts, split(vec)):
            if kind == "fix":
                a = comp[0]
                if mul(a, power(a, P ** (f // 2))) != one:
                    return False
            else:
                a, b = comp
                if mul(a, b) != one:
                    return False
        return True

    N = [v for v in itertools.product(range(P), repeat=d) if is_norm_one(v)]
    seen, hyps = set(), []
    for c in itertools.product(range(P), repeat=d):
        if not any(c) or c in seen:
            continue
        for k in range(1, P):
            seen.add(tuple((k * x) % P for x in c))
        hyps.append(c)
    missed = sum(1 for c in hyps
                 if not any(sum(x * y for x, y in zip(c, u)) % P == 0
                            for u in N))
    label = " + ".join(f"{'fixed' if k == 'fix' else 'swapped pair'} f={f}"
                       for k, f in orbits)
    if show:
        print(f"  degree {d}: {label:46s} |N| = {len(N):5d}  "
              + ("CAN BLOCK" if missed == 0 else
                 f"cannot block ({missed} hyperplanes miss)")
              + f"   [{time.time()-t0:.0f}s]", flush=True)
    return missed == 0


def types(d):
    """All sigma-orbit multisets summing to d, with fixed orbits of even f."""
    out = []

    def rec(rem, cur, minsize):
        if rem == 0:
            out.append(list(cur))
            return
        for f in range(1, rem + 1):
            if f % 2 == 0 and f <= rem and ("fix", f) >= minsize:
                rec(rem - f, cur + [("fix", f)], ("fix", f))
            if 2 * f <= rem and ("swap", f) >= minsize:
                rec(rem - 2 * f, cur + [("swap", f)], ("swap", f))
    rec(d, [], ("", 0))
    return out


for d in (2, 4, 6):
    print(f"-- degree {d}: every unramified arithmetic type --", flush=True)
    any_ok = False
    for t in types(d):
        if decide(t):
            any_ok = True
    print(f"  => degree {d} " + ("HAS a blocking type" if any_ok else
                                 "CANNOT BLOCK, in any field, ever"),
          flush=True)
