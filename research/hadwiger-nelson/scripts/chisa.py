import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver
t0 = time.time()
for name, P in (("Sa", build_Sa(K)), ("Y", build_Y(K))):
    b = IntBasis.covering(P)
    r = b.rows(P)
    E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
    n = len(P)
    cands = sorted({Fr(p, q) for q in range(1, 10)
                    for p in range(2, 5 * q + 1)
                    if Fr(p, q).denominator == q and 3 <= Fr(p, q) <= 5})
    print(f"{name}: {n} pts, {len(E)} edges  [{time.time()-t0:.0f}s]",
          flush=True)
    for rr in cands:
        p, q = rr.numerator, rr.denominator
        cl = [[1 + v * p + j for j in range(p)] for v in range(n)]
        for a, c in E:
            for j in range(p):
                for d in range(-(q - 1), q):
                    cl.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
        cl += [[-(1 + j)] for j in range(1, p)] + [[1]]
        s = Solver(name="cd15", bootstrap_with=cl)
        ok = s.solve()
        s.delete()
        if ok:
            print(f"   chi_c({name}) = {rr} = {float(rr):.4f}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            break
        print(f"   not {rr} = {float(rr):.4f}  [{time.time()-t0:.0f}s]",
              flush=True)
print("DONE", flush=True)
