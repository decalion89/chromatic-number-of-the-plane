"""Write a saved vector set V_*.json as a GP file: X and Y as polynomials in a, b (the two generators)."""
import json, sys
d = json.load(open(sys.argv[1]))
polys = d["polys"]
assert len(polys) == 2
def poly(coeffs, var):
    return "+".join(f"({c})*{var}^{k}" for k, c in enumerate(coeffs))
def elt(dct):
    terms = []
    for key, c in dct.items():
        e = [int(x) for x in key.split(",")]
        terms.append(f"({c})*a^{e[0]}*b^{e[1]}")
    return "+".join(terms) if terms else "0"
out = [f"fa = {poly(polys[0], 'a')};", f"fb = {poly(polys[1], 'b')};", f"da = {len(polys[0]) - 1}; db = {len(polys[1]) - 1};", "V = ["]
body = ",".join(f"[{elt(X)}, {elt(Y)}]" for X, Y in d["V"])
lines = out[:3] + ["V = [" + body + "];"]
open(sys.argv[2], "w").write("\n".join(lines) + "\n")
