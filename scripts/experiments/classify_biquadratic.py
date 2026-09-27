"""Local upper bounds for the planes over L = Q(sqrt a, sqrt b), from the places of L that do not split
in K = L(i) (Proposition A of notes/local_colourings.md):
  * K_w / L_v ramified: the norm-one residues are +-1, so chi(L^2) <= 3; above 2, where +1 = -1, the
    target is a perfect matching and chi(L^2) <= 2;
  * unramified, residue field F_q: chi(L^2) <= chi(G_q), G_q = Cay(F_{q^2}, N_1).
Known values (SAT; tests/test_split_places.py and tests/test_biquadratic_bounds.py):
  chi(G_q) = 4, 4, 3, 4, 5, 5 for q = 2, 4, 3, 7, 11, 19.
Places come from scripts/fieldscreen.py (Kummer theory) with d = (1, a, b), since K = Q(sqrt-1, sqrt-a, sqrt-b).

usage: classify_biquadratic.py [pmax] [placebound]    (pairs of primes a < b < pmax)
"""
import os, sys, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fieldscreen import screen
from sympy import primerange

CHI_G = {2: 4, 4: 4, 3: 3, 7: 4, 11: 5, 19: 5}


def local_bound(a, b, placebound=200):
    """(best bound, place giving it) over the non-split places below placebound, or (None, None)"""
    best, why = None, None
    for p, kind, q in screen([1, a, b], placebound):
        if kind == "ramified":
            val = 2 if p == 2 else 3
        else:
            val = CHI_G.get(q)
            if val is None:
                continue
        if best is None or val < best:
            best, why = val, (p, kind, q)
    return best, why


if __name__ == "__main__":
    pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    pb = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    primes = list(primerange(2, pmax))
    rows = []
    for i, a in enumerate(primes):
        for b in primes[i + 1:]:
            bd, why = local_bound(a, b, pb)
            rows.append((a, b, 3 in (a, b), bd, why))
    cnt = collections.Counter((r[2], r[3]) for r in rows)
    print(f"pairs of primes below {pmax}; places below {pb}")
    for has3 in (True, False):
        print(("with" if has3 else "without") + " sqrt3:",
              {str(k[1]): v for k, v in sorted(cnt.items(), key=str) if k[0] == has3})
    for a, b, has3, bd, why in rows:
        print(f"  Q(sqrt{a}, sqrt{b}): {bd if bd else '-'}" + (f"  via p = {why[0]} ({why[1]}, F_{why[2]})" if why else ""))
