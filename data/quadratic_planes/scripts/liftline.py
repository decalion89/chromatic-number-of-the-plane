"""liftline.py p k IN.npy|all OUT [K] [budget]: is the preimage at level k+1 of a set S of level-k points of the
p-adic plane (p = 3 mod 4) K-colourable?  Line encoding of the level-(k+1) graph over level k:
a point of level k+1 is z~ + p^k x with z in S (integer lift z~ in [0, p^k)^2) and x in F_p^2 (the next digit).
For z, w in S with w - z = t a unit vector mod p^k, put d = w~ - z~ (an integer vector, N(d) = 1 mod p^k),
u = d mod p and e = (N(d) - 1) / p^k mod p.  Since N(d + p^k v) = N(d) + 2 p^k d.v mod p^(k+1) (k >= 1),
z~ + p^k x and w~ + p^k y are adjacent iff u.(y - x) = -e/2 mod p: the line u.x = s of block z is completely
joined to the line u.y = s - e/2 of block w, and there are no other edges (points of one block differ by a
multiple of p^k, which is not a unit vector).  CNF: point variables P[z,x,c] (at least one colour each), line
variables L[z,u,s,c] with P[z,x,c] -> L[z,u,u.x,c] for every label u at z, and not L[z,u,s,c] or not
L[w,u,s-e/2,c] for every edge and every s, c.  A proper colouring gives a model (L = 'some point of the line
has colour c'), so UNSAT means the preimage has no proper K-colouring.  One edge of level k + 1 is pinned
(colours 0, 1).  Writes OUT.cnf, runs kissat (DRAT to OUT.drat) and drat-trim if UNSAT.
IN.npy holds codes x * p^k + y of level-k points; 'all' takes all of (Z/p^k)^2."""
import sys, os, subprocess, time
import numpy as np
from collections import defaultdict
p, k, IN, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
K = int(sys.argv[5]) if len(sys.argv) > 5 else 4; budget = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # kissat/ and drat-trim/ live here
t0 = time.time(); m = p ** k
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
S = np.arange(m * m, dtype=np.int64) if IN == "all" else np.unique(np.load(IN).astype(np.int64))
n = len(S); sx, sy = S // m, S % m
pos = -np.ones(m * m, dtype=np.int64); pos[S] = np.arange(n)
EA, EB, EU, EC, ESG = [], [], [], [], []
inv2 = pow(2, -1, p)
for i in range(0, len(T), 64):
    tb = T[i:i + 64]
    wx = (sx[:, None] + tb[None, :, 0]) % m; wy = (sy[:, None] + tb[None, :, 1]) % m
    j = pos[wx * m + wy]
    a = np.broadcast_to(np.arange(n)[:, None], j.shape)
    keep = j > a
    a, j, wx, wy = a[keep], j[keep], wx[keep], wy[keep]
    dx, dy = wx - sx[a], wy - sy[a]
    N = dx * dx + dy * dy
    assert np.all((N - 1) % m == 0)
    e = ((N - 1) // m) % p
    EA.append(a); EB.append(j); cu, nu = (dx % p) * p + dy % p, ((-dx) % p) * p + (-dy) % p
    EU.append(np.minimum(cu, nu)); ESG.append(np.where(cu <= nu, 1, -1)); EC.append((-e * inv2) % p)
EA, EB, EU, EC, ESG = (np.concatenate(v) for v in (EA, EB, EU, EC, ESG))
EC = (EC * ESG) % p      # u = sign * canonical label (u and -u give the same lines): canonical.(y - x) = sign * c
print(f"p={p}: {n} points at level {k}, {len(EA)} edges among them; preimage at level {k + 1}: {n * p * p} points"
      f"  [{time.time() - t0:.1f}s]", flush=True)
# variables
P = lambda z, x, c: (z * p * p + x) * K + c + 1          # x = x0 * p + x1 in [0, p^2)
nP = n * p * p * K
labels = defaultdict(set)
for a, b, u in zip(EA.tolist(), EB.tolist(), EU.tolist()):
    labels[a].add(u); labels[b].add(u)
lid = {}; nxt = nP + 1
for z in range(n):
    for u in sorted(labels[z]):
        lid[z, u] = nxt; nxt += p * K
L = lambda z, u, s, c: lid[z, u] + s * K + c
X0 = np.arange(p * p) // p; X1 = np.arange(p * p) % p
nclauses = 0
with open(out + ".cnf", "w") as f:
    f.write(" " * 64 + "\n")
    buf = []
    for z in range(n):
        for x in range(p * p):
            buf.append(" ".join(str(P(z, x, c)) for c in range(K)) + " 0\n")
        for u in sorted(labels[z]):
            ux, uy = u // p, u % p
            svals = ((ux * X0 + uy * X1) % p).tolist()
            for x in range(p * p):
                for c in range(K):
                    buf.append(f"-{P(z, x, c)} {L(z, u, svals[x], c)} 0\n")
        if len(buf) > 200000:
            f.write("".join(buf)); nclauses += len(buf); buf = []
    for a, b, u, cc in zip(EA.tolist(), EB.tolist(), EU.tolist(), EC.tolist()):
        for s in range(p):
            s2 = (s + cc) % p
            for c in range(K):
                buf.append(f"-{L(a, u, s, c)} -{L(b, u, s2, c)} 0\n")
        if len(buf) > 200000:
            f.write("".join(buf)); nclauses += len(buf); buf = []
    # pin one edge of level k+1: block a point x=0, block b a point y on the matched line
    a, b, u, cc = int(EA[0]), int(EB[0]), int(EU[0]), int(EC[0])
    ux, uy = u // p, u % p
    y = next(y for y in range(p * p) if (ux * (y // p) + uy * (y % p)) % p == cc % p)
    buf += [f"{P(a, 0, 0)} 0\n", f"{P(b, y, 1)} 0\n"]
    f.write("".join(buf)); nclauses += len(buf)
nvars = nxt - 1
with open(out + ".cnf", "r+") as f:
    f.write(f"p cnf {nvars} {nclauses}")
print(f"  CNF {out}.cnf: {nvars} variables, {nclauses} clauses  [{time.time() - t0:.1f}s]", flush=True)
NOPROOF = os.environ.get("NOPROOF") == "1"          # decide first without a proof file (proofs of big runs fill the disk)
r = subprocess.run([f"{SC}/kissat/build/kissat", f"--time={budget}", out + ".cnf"] + ([] if NOPROOF else [out + ".drat"]), capture_output=True, text=True)
res = [l for l in r.stdout.splitlines() if l.startswith("s ")]
print(f"  kissat: {res}  [{time.time() - t0:.0f}s]", flush=True)
if res == ["s UNSATISFIABLE"] and not NOPROOF:
    o = subprocess.run([f"{SC}/drat-trim/drat-trim", out + ".cnf", out + ".drat", "-t", "50000"], capture_output=True, text=True).stdout
    print(f"  drat-trim: {[l for l in o.splitlines() if l.startswith('s ')]}  [{time.time() - t0:.0f}s]", flush=True)
elif res == ["s SATISFIABLE"]:
    if not NOPROOF:
        os.unlink(out + ".drat")
    # keep the colouring: point (z, x) gets the first colour c with P[z, x, c] true
    val = set(int(t) for l in r.stdout.splitlines() if l.startswith("v ") for t in l.split()[1:] if int(t) > 0)
    col = np.array([[next((c for c in range(K) if P(z, x, c) in val), -1) for x in range(p * p)] for z in range(n)])
    codes = np.array([(int(sx[z]) + m * (x // p)) * (m * p) + int(sy[z]) + m * (x % p) for z in range(n) for x in range(p * p)])
    np.savez(out + "_colouring.npz", codes=codes, colours=col.ravel())
    print(f"  saved the colouring of {len(codes)} level-{k + 1} points to {out}_colouring.npz (check it with an independent script)", flush=True)
