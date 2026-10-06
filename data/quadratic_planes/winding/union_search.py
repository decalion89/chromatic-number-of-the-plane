"""union_search.py d N TL NPROC: floating-point guide for Theorem W over many small unions.
Ranks the denominators D <= 6000 (6 | D) by the number of unit vectors (up to sign), keeps the N richest, and tests
every single one, every pair and every triple of the eight richest, smallest union first, with the relation-space
MILP of kapparel.py (time limit TL seconds each, NPROC at a time). One line per test. At the first INFEASIBLE test
it prints 'FOUND d D1,D2,...' and stops; union_cert.py then writes and checks the exact certificate."""
import sys, math, time, itertools
from multiprocessing import Pool
from fractions import Fraction as Fr
from kappaD import units, halve
from kapparel import relation_basis, solve_rel


def union(Ds, cache):
    L = 1
    for D in Ds: L = L * D // math.gcd(L, D)
    seen, V = set(), []
    for D in Ds:
        for u in cache[D]:
            w = tuple(x * (L // D) for x in u)
            key = max(w, tuple(-x for x in w))
            if key not in seen:
                seen.add(key); V.append(list(w))
    return L, V


def test(job):
    Ds, L, V, TL = job
    t0 = time.time()
    ker = relation_basis(V)
    name, _ = solve_rel(len(V), ker, Fr(1, 3), TL)
    return Ds, L, len(V), len(ker), name, time.time() - t0


if __name__ == "__main__":
    d, N, TL, NP = map(int, sys.argv[1:5])
    t0 = time.time()
    cache = {D: halve(units(d, D)) for D in range(6, 6001, 6)}
    ranked = sorted(((len(v), D) for D, v in cache.items() if len(v) >= 12), reverse=True)[:N]
    top = [D for n, D in ranked]
    print(f"d={d} the {len(top)} richest D: " + ", ".join(f"{D}:{n}" for n, D in ranked) + f"  [{time.time()-t0:.0f}s]",
          flush=True)
    cands = [(D,) for D in top] + list(itertools.combinations(top, 2)) + list(itertools.combinations(top[:8], 3))
    jobs, sizes = [], set()
    for Ds in cands:
        L, V = union(Ds, cache)
        key = frozenset(tuple(Fr(x, L) for x in v) for v in V)
        if key in sizes: continue          # the same set again (one D divides another)
        sizes.add(key); jobs.append((len(V), Ds, L, V))
    jobs.sort(key=lambda x: x[0])
    print(f"d={d}: {len(jobs)} distinct unions, sizes {jobs[0][0]}..{jobs[-1][0]}", flush=True)
    with Pool(NP) as pool:
        for Ds, L, n, r, name, t in pool.imap_unordered(test, [(Ds, L, V, TL) for n, Ds, L, V in jobs]):
            print(f"d={d} union {list(Ds)} L={L} |U|={n} relations={r}: {name} [{t:.0f}s]", flush=True)
            if name == "INFEASIBLE":
                print(f"FOUND {d} {','.join(map(str, Ds))}", flush=True)
                pool.terminate()
                break
    print(f"d={d} search ended [{time.time()-t0:.0f}s]", flush=True)
