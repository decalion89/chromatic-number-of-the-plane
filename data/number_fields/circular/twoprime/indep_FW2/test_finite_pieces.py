"""For examples with 2 < chi_c = p/q < 4 (kappa computed exactly), grow finite pieces of Cay(Gamma, S)
until every homomorphism of the piece to K_{p/q} has a tight cycle (SAT: no acyclic-tight homomorphism),
i.e. until the piece has chi_c = p/q (Lemma 22 predicts that this happens).  For small pieces cross-check
with an independent test: chi_c(H) < p/q iff H -> K_{a/b} for the largest a/b < p/q with a <= |V(H)|."""
import sys, time, json
from fractions import Fraction as Fr
from kappa_exact import Group, symmetrize, half, kappa
from sat_tight import piece, hom_exists, acyclic_hom_exists, farey_pred

EXAMPLES = [
    # (name, r, m, S_half (with torsion coordinate), growth: list of boxes)
    ("Z, D={3,4,9,12}", 1, 1, [(3, 0), (4, 0), (9, 0), (12, 0)], [[(0, n)] for n in range(6, 80)]),
    ("Z, D={1,2,5}", 1, 1, [(1, 0), (2, 0), (5, 0)], [[(0, n)] for n in range(2, 60)]),
    ("Z, D={5,3}+... (16/5 via Z x Z/2)", 1, 2, [(5, 1), (3, 0), (2, 1)], [[(0, n)] for n in range(2, 60)]),
    ("Z, D={1,4}", 1, 1, [(1, 0), (4, 0)], [[(0, n)] for n in range(2, 30)]),
    ("Z, D={3,4}", 1, 1, [(3, 0), (4, 0)], [[(0, n)] for n in range(2, 40)]),
    ("Z, D={2,5,9}", 1, 1, [(2, 0), (5, 0), (9, 0)], [[(0, n)] for n in range(2, 80)]),
    ("Z, D={1,3,8}", 1, 1, [(1, 0), (3, 0), (8, 0)], [[(0, n)] for n in range(2, 80)]),
    ("Z^2, S={(1,0),(0,1),(1,2)}", 2, 1, [(1, 0, 0), (0, 1, 0), (1, 2, 0)], [[(0, n), (0, n)] for n in range(1, 12)]),
    ("Z^2, S={(2,-1),(1,1),(0,2),(1,0)}", 2, 1, [(2, -1, 0), (1, 1, 0), (0, 2, 0), (1, 0, 0)],
     [[(0, n), (0, n)] for n in range(1, 12)]),
    ("Z^2, S={(3,2),(0,1),(2,2)}", 2, 1, [(3, 2, 0), (0, 1, 0), (2, 2, 0)], [[(0, n), (0, n)] for n in range(1, 12)]),
    ("Z^2, S={(1,0),(0,1),(2,3)}", 2, 1, [(1, 0, 0), (0, 1, 0), (2, 3, 0)], [[(0, n), (0, n)] for n in range(1, 12)]),
    ("Z x Z/4, S={(4,1),(2,2),(1,2)}", 1, 4, [(4, 1), (2, 2), (1, 2)], [[(0, n)] for n in range(1, 40)]),
    ("Z x Z/3, S={(4,0),(3,0),(0,1)}", 1, 3, [(4, 0), (3, 0), (0, 1)], [[(0, n)] for n in range(1, 40)]),
    ("Z x Z/3, S={(1,1),(2,0)}", 1, 3, [(1, 1), (2, 0)], [[(0, n)] for n in range(1, 40)]),
]


def run(names=None, cross_check_max=36):
    results = []
    for (name, r, m, Sh, boxes) in EXAMPLES:
        if names and not any(nm in name for nm in names):
            continue
        G = Group(r, m)
        S = symmetrize(G, Sh)
        kap, opt = kappa(G, S)
        chi = 1 / kap
        p, q = chi.numerator, chi.denominator
        rec = dict(name=name, kappa=str(kap), chi_c=str(chi))
        print(f"== {name}: kappa={kap}, chi_c={chi}", flush=True)
        if not (2 < chi < 4):
            print("   (not in (2,4), skipped)")
            continue
        t0 = time.time()
        found = None
        for box in boxes:
            V, E = piece(r, m, Sh, box)
            N = len(V)
            col = hom_exists(N, E, p, q)
            assert col is not None, "piece must map to K_{p/q}"
            acyc = acyclic_hom_exists(N, E, p, q)
            status = "acyclic-tight hom EXISTS (chi_c(piece) < p/q)" if acyc else "every hom has a tight cycle (chi_c(piece) = p/q)"
            cc = ""
            if N <= cross_check_max:
                fp = farey_pred(p, q, N)
                if fp is not None:
                    f, a, b = fp
                    lower = hom_exists(N, E, a, b)
                    agree = (lower is None) == (acyc is None)
                    cc = f" | cross-check K_{a}/{b}: {'hom' if lower else 'no hom'} -> {'agree' if agree else 'DISAGREE'}"
                    assert agree, (name, box)
            print(f"   box={box} |V|={N} |E|={len(E)}: {status}{cc}", flush=True)
            if acyc is None:
                found = dict(box=box, V=N, E=len(E))
                break
        rec['first_piece_with_chi_c'] = found
        rec['time_s'] = round(time.time() - t0, 1)
        results.append(rec)
        print("   ->", rec, flush=True)
    return results


if __name__ == '__main__':
    names = sys.argv[1:] or None
    res = run(names)
    with open('finite_pieces_results.json', 'a') as f:
        f.write(json.dumps(res) + "\n")
