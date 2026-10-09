"""periodic4.py -- periodic K-colourings of the unit-distance graph over Q(sqrt d) with directions U_D:
is Cay(M/mM, U) K-colourable, where M is the Z-span of U_D (a lattice in Z^4)?  (periodicq.py at K colours)
A yes for some m gives a K-colouring of every graph built from U_D: growth with these directions cannot succeed.
usage: KCOL=4 KISSAT=path periodic4.py d D TL m1 m2 ...   (kissat, time limit TL seconds per m; a colouring found
is checked and saved in OUTDIR)"""
import sys, os, itertools, time, subprocess
import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from units_fast import units_fast
KISSAT = os.environ.get("KISSAT", "kissat")
OUT = os.environ.get("OUTDIR", ".")
K = int(os.environ.get("KCOL", "4"))
d, D, TL = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ms = [int(x) for x in sys.argv[4:]]
U = [tuple(int(x) for x in u) for u in np.asarray(units_fast(d, D)).reshape(-1, 4)]
H = hermite_normal_form(Matrix(U).T)
B = Matrix.hstack(*[H[:, j] for j in range(H.shape[1]) if any(H[:, j])])
r = B.shape[1]
Binv = B.inv()
C = [tuple(int(x) for x in Binv * Matrix(u)) for u in U]
assert all(B * Matrix(c) == Matrix(u) for c, u in zip(C, U))
print(f"d={d} D={D} K={K}: {len(U)} directions, rank {r}", flush=True)
for m in ms:
    t = time.time()
    G = sorted(set(tuple(x % m for x in c) for c in C))
    if tuple([0] * r) in G:
        print(f"  m={m}: a direction is 0 mod m: no periodic colouring", flush=True)
        continue
    G = [g for g in G if g <= tuple((-x) % m for x in g)]          # one of each pair +-g
    N = m ** r
    P = np.array(list(itertools.product(range(m), repeat=r)), dtype=np.int64)
    w = m ** np.arange(r - 1, -1, -1)
    path = f"/dev/shm/per4_{d}_{D}_{m}.cnf"
    ncl = N + 1 + len(G) * N * K
    with open(path, "w") as f:
        f.write(f"p cnf {N * K} {ncl}\n")
        f.write("".join(" ".join(str(K * i + c + 1) for c in range(K)) + " 0\n" for i in range(N)))
        f.write("1 0\n")
        for g in G:
            J = ((P + np.array(g)) % m) @ w
            I = np.arange(N)
            lines = []
            for c in range(K):
                lines.append(np.char.add(np.char.add(np.char.add((-(K * I + c + 1)).astype(str), " "), (-(K * J + c + 1)).astype(str)), " 0"))
            f.write("\n".join(np.concatenate(lines)) + "\n")
    out = subprocess.run([KISSAT, f"--time={TL}", path], capture_output=True, text=True).stdout
    st = [l for l in out.splitlines() if l.startswith("s ")]
    res = st[0] if st else "s UNKNOWN (time-out)"
    print(f"  m={m}: {N} classes, {len(G)} generator pairs: {res}  [{time.time() - t:.0f}s]", flush=True)
    if res == "s SATISFIABLE":
        lits = [int(x) for l in out.splitlines() if l.startswith("v ") for x in l[2:].split()]
        col = np.full(N, -1)
        for x in lits:
            if x > 0 and col[(x - 1) // K] < 0:
                col[(x - 1) // K] = (x - 1) % K
        # check the colouring independently
        bad = sum(int(np.sum(col == col[((P + np.array(g)) % m) @ w])) for g in G)
        print(f"    colouring check: {bad} monochromatic edges", flush=True)
        np.save(f"{OUT}/periodic4_col_{d}_{D}_{m}.npy", col)
    os.remove(path)
