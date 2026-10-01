"""certify_q.py -- certify a vertex-critical non-3-colourable unit-distance graph over Q(sqrt R).

usage: python3 certify_q.py R CRIT.json OUTDIR
Checks, with its own code: the points are distinct elements (a + b sqrt R)/D, (c + d sqrt R)/D; every edge has
squared length exactly 1 (in Z[sqrt R], over D^2); the edge list equals all unit pairs among the given directions;
no triangle. Then: kissat + DRAT + drat-trim on two encodings of "3-colourable with one edge fixed"; a proper
4-colouring; for every vertex v a proper 3-colouring of G - v (vertex-criticality). Writes OUTDIR/certificate.json
and the logs."""
import sys, os, json, subprocess, itertools
R = int(sys.argv[1]); g = json.load(open(sys.argv[2])); out = sys.argv[3]
os.makedirs(out, exist_ok=True)
SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KISSAT, DRAT = f"{SC}/kissat/build/kissat", f"{SC}/drat-trim/drat-trim"
D, pts, E = g["D"], [tuple(p) for p in g["points"]], [tuple(e) for e in g["edges"]]
n = len(pts)
assert len(set(pts)) == n


def sqlen(p, q):
    """|q - p|^2 * D^2 as an element x + y sqrt R of Z[sqrt R]: (dx)^2 + (dy)^2 with dx = a + b s, dy = c + d s"""
    a, b, c, d = (qq - pp for pp, qq in zip(p, q))
    return (a * a + R * b * b + c * c + R * d * d, 2 * a * b + 2 * c * d)


unit = (D * D, 0)
assert all(sqlen(pts[a], pts[b]) == unit for a, b in E), "an edge is not of length 1"
allu = {(a, b) for a in range(n) for b in range(a + 1, n) if sqlen(pts[a], pts[b]) == unit}
extra = allu - {(min(a, b), max(a, b)) for a, b in E}
adj = [set() for _ in range(n)]
for a, b in E:
    adj[a].add(b); adj[b].add(a)
tri = sum(1 for a, b in E for c in adj[a] & adj[b])
print(f"{n} points, {len(E)} edges of length exactly 1; unit pairs not listed as edges: {len(extra)}; triangles: {tri}")


def write(path, k, enc, fixed_edge=True):
    var = (lambda v, c: k * v + c + 1) if enc == 0 else (lambda v, c: c * n + v + 1)
    cl = [[var(v, c) for c in range(k)] for v in range(n)]
    cl += [[-var(a, c), -var(b, c)] for a, b in E for c in range(k)]
    if fixed_edge:
        a0, b0 = E[0]
        cl += [[var(a0, 0)], [var(b0, 1)]]
    with open(path, "w") as f:
        f.write(f"p cnf {k * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    return var


res = {}
for enc in (0, 1):
    cnf, proof = f"{out}/g3_enc{enc}.cnf", f"/dev/shm/g3_{R}_{os.getpid()}_enc{enc}.drat"
    write(cnf, 3, enc)
    k = subprocess.run([KISSAT, cnf, proof], capture_output=True, text=True).stdout
    open(f"{out}/g3_enc{enc}.kissat.log", "w").write(k)
    d = subprocess.run([DRAT, cnf, proof, "-t", "20000"], capture_output=True, text=True).stdout
    open(f"{out}/g3_enc{enc}.drat-trim.log", "w").write(d)
    os.remove(proof)
    res[enc] = ("s UNSATISFIABLE" in k, "s VERIFIED" in d)
    print(f"encoding {enc}: kissat UNSAT {res[enc][0]}, drat-trim VERIFIED {res[enc][1]}")


def solve(k, drop=None):
    """a proper k-colouring of G (or of G - drop), by kissat"""
    keep = [v for v in range(n) if v != drop]
    pos = {v: i for i, v in enumerate(keep)}
    Es = [(pos[a], pos[b]) for a, b in E if drop not in (a, b)]
    cl = [[k * i + c + 1 for c in range(k)] for i in range(len(keep))]
    cl += [[-(k * a + c + 1), -(k * b + c + 1)] for a, b in Es for c in range(k)]
    p = f"/dev/shm/certq_col_{R}_{os.getpid()}.cnf"
    with open(p, "w") as f:
        f.write(f"p cnf {k * len(keep)} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    o = subprocess.run([KISSAT, p], capture_output=True, text=True).stdout
    if "s SATISFIABLE" not in o:
        return None
    lits = set(int(x) for l in o.splitlines() if l.startswith("v ") for x in l.split()[1:] if int(x) > 0)
    col = [-1] * n
    for v, i in pos.items():
        col[v] = next(c for c in range(k) if k * i + c + 1 in lits)
    assert all(col[a] != col[b] for a, b in E if drop not in (a, b))
    return col


c4 = solve(4)
print("proper 4-colouring:", c4 is not None)
crit = {}
for v in range(n):
    c = solve(3, drop=v)
    assert c is not None, f"G - {v} is not 3-colourable: not vertex-critical"
    crit[str(v)] = c
print(f"vertex-critical: a proper 3-colouring of G - v for all {n} vertices")
json.dump({"field": f"Q(sqrt{R})", "D": D, "points": [list(p) for p in pts], "edges": [list(e) for e in E],
           "four_colouring": c4, "critical_3_colourings": crit,
           "checks": {"edges_unit_exact": True, "unlisted_unit_pairs": len(extra), "triangles": tri,
                      "not_3_colourable": {f"encoding_{e}": {"kissat_unsat": r[0], "drat_trim_verified": r[1]}
                                           for e, r in res.items()}}},
          open(f"{out}/certificate.json", "w"))
print("CERTIFIED" if all(r[0] and r[1] for r in res.values()) and c4 else "NOT CERTIFIED")
