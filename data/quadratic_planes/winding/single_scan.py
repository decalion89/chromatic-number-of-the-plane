"""single_scan.py d1,d2,... TL NPROC: floating-point guide for Theorem W over single denominators.
For each d, every D <= 6000 with 6 | D and at least 24 unit vectors (up to sign) is tested with the relation-space
MILP of kapparel.py (time limit TL seconds, NPROC at a time), fewest vectors first. One line per test; at the first
INFEASIBLE denominator it prints 'FOUND d D' and moves to the next d."""
import sys, time
from multiprocessing import Pool
from fractions import Fraction as Fr
from kappaD import units, halve
from kapparel import relation_basis, solve_rel


def test(job):
    D, V, TL = job
    t0 = time.time()
    ker = relation_basis(V)
    name, _ = solve_rel(len(V), ker, Fr(1, 3), TL)
    return D, len(V), len(ker), name, time.time() - t0


if __name__ == "__main__":
    ds = [int(x) for x in sys.argv[1].split(',')]
    TL, NP = int(sys.argv[2]), int(sys.argv[3])
    for d in ds:
        t0 = time.time()
        cand = sorted((len(V), D, V) for D in range(6, 6001, 6) for V in [halve(units(d, D))] if len(V) >= 24)
        print(f"d={d}: {len(cand)} denominators with >= 24 vectors, sizes {cand[0][0]}..{cand[-1][0]} [{time.time()-t0:.0f}s]",
              flush=True)
        with Pool(NP) as pool:
            for D, n, r, name, t in pool.imap_unordered(test, [(D, V, TL) for n, D, V in cand]):
                print(f"d={d} D={D} |U|={n} relations={r}: {name} [{t:.0f}s]", flush=True)
                if name == "INFEASIBLE":
                    print(f"FOUND {d} {D}", flush=True)
                    pool.terminate()
                    break
        print(f"d={d} scan ended [{time.time()-t0:.0f}s]", flush=True)
