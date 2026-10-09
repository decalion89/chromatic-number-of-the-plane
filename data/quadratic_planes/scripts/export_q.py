"""export_q.py d CERTDIR OUTDIR: write OUTDIR/q{d}.json, q{d}.cnf and q{d}.logs/ from a certify_q.py directory."""
import json, os, shutil, sys
d, cdir, out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
c = json.load(open(os.path.join(cdir, "certificate.json")))
assert c["field"] == f"Q(sqrt{d})"
ch = c["checks"]
assert ch["edges_unit_exact"] and ch["unlisted_unit_pairs"] == 0 and ch["triangles"] == 0
for enc in ("encoding_0", "encoding_1"):
    assert ch["not_3_colourable"][enc] == {"kissat_unsat": True, "drat_trim_verified": True}
n = len(c["points"])


def digits(col, drop=None):
    return "".join("-" if v == drop else str(int(col[v])) for v in range(n))


crit = c["critical_3_colourings"]
g = {
    "d": d,
    "D": c["D"],
    "description": f"unit-distance graph in the plane over Q(sqrt{d}): point [a, b, c, e] is ((a + b sqrt{d})/D, "
                   f"(c + e sqrt{d})/D); edges are all pairs at distance 1",
    "points": c["points"],
    "edges": c["edges"],
    "fixed_edge": c["edges"][0],
    "four_colouring": digits(c["four_colouring"]),
    "critical_3_colourings": {str(v): digits(crit[str(v)], drop=v) for v in range(n)},
    "certificate": {"formula": f"q{d}.cnf: 3-colourable with the fixed edge coloured 0, 1 (variable 3v + c + 1)",
                    "logs": f"q{d}.logs/",
                    "result": "kissat UNSATISFIABLE, drat-trim VERIFIED (two encodings)"},
}
os.makedirs(os.path.join(out, f"q{d}.logs"), exist_ok=True)
json.dump(g, open(os.path.join(out, f"q{d}.json"), "w"))
shutil.copy(os.path.join(cdir, "g3_enc0.cnf"), os.path.join(out, f"q{d}.cnf"))
for e in (0, 1):
    for kind in ("kissat", "drat-trim"):
        shutil.copy(os.path.join(cdir, f"g3_enc{e}.{kind}.log"), os.path.join(out, f"q{d}.logs", f"encoding{e}.{kind}.log"))
print(f"q{d}: {n} points, {len(c['edges'])} edges, D = {c['D']}")
