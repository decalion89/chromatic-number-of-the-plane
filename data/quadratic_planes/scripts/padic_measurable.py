"""padic_measurable.py PMAX: lower bounds for the measurable chromatic number of the p-adic plane (note section 5).

For p = 3 (mod 4) the unit circle S of Q_p^2 (x^2 + y^2 = 1) lies in Z_p^2, so the compact group Z_p^2 is a union
of edges' components: x ~ x + s with s in S. If a measurable set A in Z_p^2 contains no pair x, x + s, its
autocorrelation f(x) = mu(A cap (A - x)) is positive definite with f(0) = mu(A) and f = 0 on S; integrating f
against the Haar measure sigma of S gives mu(A) <= -m / (1 - m), where m is the least value of the Fourier
transform of sigma at a nontrivial character of Z_p^2 (Delsarte, Hoffman). A character of exact level 1 gives
lambda(xi) / (p + 1), an eigenvalue of the unit-distance graph of F_p^2 over its degree; a character of level
l >= 2 gives at most 2 / (p + 1) in absolute value (average over the rotations = 1 mod p^ceil(l/2) first; it
leaves the two arcs where s is parallel to xi mod p^(l - ceil(l/2)), of measure at most 2 / (p + 1)). So every
measurable colouring of Q_p^2 (each coset of Z_p^2 is coloured) has a class of measure >= 1/k in Z_p^2, and
    chi_m(Q_p^2) >= 1 + 1/|m|,   m = min(lambda_min / (p + 1), -2 / (p + 1)).
lambda depends only on the norm of xi (rotations), so one xi per norm class is enough. Everything is computed in
interval arithmetic (mpmath.iv); the script prints the rigorous bound ceil(lower end of 1 + 1/|m|) and checks
Weil's bound |lambda| <= 2 sqrt(p), which gives chi_m >= 1 + (p + 1)/(2 sqrt p) for every such p."""
import sys
from math import isqrt
from mpmath import iv, mpf, floor as mfloor

iv.prec = 120


def isprime(n):
    return n > 1 and all(n % q for q in range(2, isqrt(n) + 1))


def bound(p):
    C = [(x, y) for x in range(p) for y in range(p) if (x * x + y * y) % p == 1]
    assert len(C) == p + 1
    reps = {}
    for x in range(p):
        for y in range(p):
            nrm = (x * x + y * y) % p
            if nrm and nrm not in reps:
                reps[nrm] = (x, y)
    cosv = [iv.cos(2 * iv.pi * t / p) for t in range(p)]
    lam_min = None
    for a, b in reps.values():
        lam = sum((cosv[(a * x + b * y) % p] for x, y in C), iv.mpf(0))
        assert lam.b <= 2 * iv.sqrt(p) and lam.a >= -2 * iv.sqrt(p)        # Weil
        lam_min = lam if lam_min is None or lam.a < lam_min.a else lam_min
    # |m| <= max(|lambda_min| / (p + 1), 2 / (p + 1)); use the upper end of |lambda_min| (rigorous)
    absm = max(-lam_min.a / (p + 1), mpf(2) / (p + 1))
    low = 1 + 1 / iv.mpf(absm)               # chi_m >= this (an interval; take its lower end)
    k = int(mfloor(low.a)) + 1 if low.a != mfloor(low.a) else int(low.a)
    return lam_min, low, k


if __name__ == "__main__":
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    print("p    | least eigenvalue of F_p^2   | 1 + 1/|m| (lower end) | chi_m(Q_p^2) >=")
    for p in range(3, PMAX, 4):
        if not isprime(p):
            continue
        lam_min, low, k = bound(p)
        print(f"{p:<4} | {float(lam_min.a):>+.10f}             | {float(low.a):.6f}              | {k}")
