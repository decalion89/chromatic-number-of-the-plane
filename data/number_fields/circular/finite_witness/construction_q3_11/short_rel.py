"""short_rel.py B OUT.json: all relations rho (in coordinates over the certificate's Z-basis r_j of the relation
lattice of the 27 vectors) with sum rho_u^2 <= B, one per +- pair, via PARI qfminim on the Gram matrix."""
import json, subprocess, sys, gzip, os
B = int(sys.argv[1]); out = sys.argv[2]
C = json.load(gzip.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "at_four", "cert311_open_4.json.gz"), "rt"))
R = C["relations"]
n, m = len(R), len(R[0])
G = [[sum(R[i][u] * R[j][u] for u in range(m)) for j in range(n)] for i in range(n)]
gp = "G = [" + ";".join(",".join(map(str, row)) for row in G) + "];\n"
gp += f"V = qfminim(G, {B}, , 2)[3];\n"
gp += 'f = "coef.txt"; for (k = 1, #V, write1(f, Vec(V[,k]))); quit;\n'
open("q.gp", "w").write(gp)
import os
if os.path.exists("coef.txt"):
    os.remove(os.path.abspath("coef.txt"))
subprocess.run(["gp", "-q", "-s", "2000000000", "q.gp"], check=True)
txt = open("coef.txt").read().replace("][", "]\n[")
rows = [json.loads(line) for line in txt.split("\n") if line.strip()]
res = []
for a in rows:
    rho = [sum(a[j] * R[j][u] for j in range(n)) for u in range(m)]
    res.append({"a": a, "rho": rho, "l1": sum(abs(x) for x in rho)})
json.dump(res, open(out, "w"))
from collections import Counter
print(len(res), "relations with l2^2 <=", B, "; l1 histogram:", sorted(Counter(r["l1"] for r in res).items()))
