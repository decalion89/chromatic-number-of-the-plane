"""certify.py -- exact checks + kissat/DRAT/drat-trim certificate of non-4-colourability.

usage: python3 certify.py PTS.npy NAME [kissat_time] [drat_time]
Writes NAME.json (vertices: exact 12-integer vectors over 7 in the power basis of zeta21, float coordinates
for zeta21 -> exp(2 pi i/21); edges with direction index into U), NAME.edges, NAME.cnf, NAME.kissat.log,
NAME.drat-trim.log. The proof goes to /dev/shm and is deleted after the check.
"""
import sys, os, json, time, subprocess, itertools, shutil
import numpy as np
from flat import *
from satutil import write_dimacs

SC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the scratchpad directory
KISSAT, DRAT = f"{SC}/kissat/build/kissat", f"{SC}/drat-trim/drat-trim"
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("pts"); ap.add_argument("name")
ap.add_argument("--kt", type=int, default=3000); ap.add_argument("--dt", type=int, default=20000)
ap.add_argument("--tri", default=None, help="a,b,c: triangle to fix (default: first found)")
ap.add_argument("--reuse", default=None, help="existing proof: skip kissat, CNF must equal NAME.cnf, NAME.kissat.log must say UNSAT")
ap.add_argument("--core", action="store_true", help="also let drat-trim write the clause core (NAME.core.cnf)")
args = ap.parse_args()
P = np.load(args.pts).astype(np.int64)
name = args.name
KT, DT = args.kt, args.dt
D, CD, U = directions()
Uset = {u: j for j, u in enumerate(U)}
n = len(P)
rep = []
# 1. vertices distinct
assert len(set(map(tuple, P.tolist()))) == n
rep.append(f"vertices: {n}, all distinct (exact): True")
# 2. edges: every pair whose difference is in U (exact lookup), then re-check each edge by exact arithmetic
E, J = build_edges(P, U)
ok_dir = all(tuple((P[b] - P[a]).tolist()) in Uset for a, b in E)
ok_norm = all(mul(list((P[b] - P[a]).tolist()), list(conj(tuple((P[b] - P[a]).tolist())))) == [49] + [0] * 11 for a, b in E)
rep.append(f"edges: {len(E)}; every edge vector is one of the 126 directions: {ok_dir}; "
           f"diff*conj(diff) = 1 exactly for every edge: {ok_norm}")
assert ok_dir and ok_norm
# brute-force numeric check: all pairs at distance 1 (float) are exactly the edges
Z = np.array([cval(p) for p in P])
cnt, mind = 0, 9.0
Es = set(map(tuple, np.sort(E, axis=1).tolist()))
extra = 0
for i in range(0, n, 500):
    d = np.abs(Z[i:i + 500, None] - Z[None, :])
    for a_, b_ in zip(*np.nonzero(np.abs(d - 1) < 1e-9)):
        a, b = i + a_, b_
        if a < b:
            cnt += 1
            if (a, b) not in Es:
                extra += 1
    dd = d[d > 0]
    mind = min(mind, dd.min())
rep.append(f"float check: unit-distance pairs (|d-1|<1e-9) = {cnt}, not among the edges: {extra}; "
           f"min distance between distinct vertices {mind:.6g}")
# degrees
deg = np.bincount(E.ravel(), minlength=n)
rep.append(f"degree min/mean/max: {deg.min()}/{deg.mean():.2f}/{deg.max()}")
# classes of the vertices in M/O = F7 omega + F7 conj(omega): count vertices in Z[zeta21]
inO = int(np.all(P % 7 == 0, axis=1).sum())
rep.append(f"vertices in Z[zeta21]: {inO}")
dir_use = np.bincount(J, minlength=126)
rep.append(f"directions used: mu42 {int((dir_use[[j for j,u in enumerate(U) if all(c % 7 == 0 for c in u)]] > 0).sum())}/42, "
           f"total {int((dir_use > 0).sum())}/126; edges by family: mu42 {int(sum(dir_use[j] for j,u in enumerate(U) if j < 84 and j % 2 == 0))}, "
           f"omega*mu42 {int(sum(dir_use[j] for j in range(1, 84, 2)))}, conj(omega)*mu42 {int(dir_use[84:].sum())}")
# 3. CNF: triangle fixed
tri = tuple(int(x) for x in args.tri.split(",")) if args.tri else find_triangle(n, E.tolist())
Es_ = set(map(tuple, np.sort(E, axis=1).tolist()))
assert all((min(a, b), max(a, b)) in Es_ for a, b in ((tri[0], tri[1]), (tri[0], tri[2]), (tri[1], tri[2]))), "not a triangle"
cls = cnf_clauses(n, E.tolist(), tri)
if args.reuse:
    old_cls = set()
    for line in open(f"{name}.cnf"):
        if line.startswith(("c", "p")):
            if line.startswith("p"):
                nv_old = int(line.split()[2])
            continue
        old_cls.add(tuple(sorted(int(x) for x in line.split()[:-1])))
    same = nv_old == 4 * n and old_cls == set(tuple(sorted(c)) for c in cls) and len(old_cls) == len(cls)
    rep.append(f"existing {name}.cnf equals the CNF regenerated from the points: {same}")
    assert same
else:
  write_dimacs(f"{name}.cnf", 4 * n, cls,
               [f"{name}: 4-colourability of a unit-distance graph with {n} vertices and {len(E)} edges",
              "variable 4v+c+1 (v = 0..n-1, c = 0..3): vertex v has colour c",
              f"clauses: one at-least-one clause per vertex, one clause per edge and colour, "
              f"units fixing the triangle {tri} to colours 0,1,2"])
rep.append(f"CNF {name}.cnf: {4 * n} variables, {len(cls)} clauses; triangle {tri} fixed")
for r in rep:
    print(r, flush=True)
# 4. kissat + DRAT
free = shutil.disk_usage("/dev/shm").free / 2 ** 30
if args.reuse:
    proof = args.reuse
    rc = 20 if any(l.startswith("s UNSATISFIABLE") for l in open(f"{name}.kissat.log")) else -1
    kt = float("nan")
    print(f"reusing kissat run {name}.kissat.log and proof {proof}", flush=True)
else:
    proof = f"/dev/shm/{name}.drat"
    print(f"free in /dev/shm: {free:.1f} GiB; running kissat (time {KT}s) with proof {proof}", flush=True)
    t = time.time()
    with open(f"{name}.kissat.log", "w") as f:
        rc = subprocess.run(["nice", "-n", "19", KISSAT, f"--time={KT}", f"{name}.cnf", proof], stdout=f, stderr=subprocess.STDOUT).returncode
    kt = time.time() - t
st = [l for l in open(f"{name}.kissat.log") if l.startswith("s ")]
psize = os.path.getsize(proof) / 2 ** 20 if os.path.exists(proof) else 0
print(f"kissat: {st[0].strip() if st else 'no answer'} (exit {rc}) in {kt:.1f}s wall; proof {psize:.1f} MiB", flush=True)
res = {"kissat": st[0].strip() if st else None, "kissat_exit": rc, "kissat_wall_s": round(kt, 1), "proof_MiB": round(psize, 1)}
if rc == 20:
    t = time.time()
    with open(f"{name}.drat-trim.log", "w") as f:
        rc2 = subprocess.run(["nice", "-n", "19", DRAT, f"{name}.cnf", proof, "-t", str(DT)] + (["-c", f"{name}.core.cnf"] if args.core else []),
                             stdout=f, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t
    out = open(f"{name}.drat-trim.log").read()
    verified = "s VERIFIED" in out
    print(f"drat-trim: {'VERIFIED' if verified else 'NOT verified'} (exit {rc2}) in {dt:.1f}s", flush=True)
    res.update({"drat_trim": "VERIFIED" if verified else "NOT VERIFIED", "drat_trim_wall_s": round(dt, 1)})
if os.path.exists(proof):
    os.unlink(proof)
# 5. JSON
verts = []
for p, z in zip(P.tolist(), Z):
    verts.append({"exact7": p, "x": repr(float(z.real)), "y": repr(float(z.imag))})
out = {"description": "unit-distance graph in the plane; vertex = (c_0 + c_1 z + ... + c_11 z^11)/7 with exact7 = [c_0..c_11], "
                      "z = zeta21, embedded by z -> exp(2 pi i/21); edges join vertices whose difference is one of the 126 unit "
                      "vectors U = mu42 u omega*mu42 u conj(omega)*mu42 (listed, scaled by 7, in 'U').",
       "n": n, "m": int(len(E)), "checks": rep, "certificate": res, "triangle_fixed": list(tri),
       "U": [list(u) for u in U], "vertices": verts, "edges": [[int(a), int(b), int(j)] for (a, b), j in zip(E, J)]}
json.dump(out, open(f"{name}.json", "w"))
with open(f"{name}.edges", "w") as f:
    f.write(f"# {n} vertices, {len(E)} edges: a b (0-based vertex indices into {name}.json)\n")
    for a, b in E:
        f.write(f"{a} {b}\n")
print(f"wrote {name}.json, {name}.edges, {name}.cnf, {name}.kissat.log" + (", " + name + ".drat-trim.log" if rc == 20 else ""))
