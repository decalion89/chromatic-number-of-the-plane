"""union_cert.py d D1,D2,... [OUTDIR]: write the union of U_D1, U_D2, ... at the common denominator L (one vector per
+- pair) as a certify_w2 input {"d", "D": L, "units"} (OUTDIR/in_d_L.json), certify it (OUTDIR/cert_d_L.json, then
gzipped) and check it with check_w.py. A single denominator is the union of one. OUTDIR defaults to the current
directory."""
import sys, math, json, subprocess, os
from kappaD import units, halve
W = os.path.dirname(os.path.abspath(__file__))
d = int(sys.argv[1]); Ds = [int(x) for x in sys.argv[2].split(',')]
OUT = sys.argv[3] if len(sys.argv) > 3 else os.getcwd()
L = 1
for D in Ds: L = L * D // math.gcd(L, D)
seen, V = set(), []
for D in Ds:
    for u in halve(units(d, D)):
        w = tuple(x * (L // D) for x in u)
        key = max(w, tuple(-x for x in w))
        if key not in seen:
            seen.add(key); V.append(list(w))
for a, b, c, e in V:
    assert a * a + d * b * b + c * c + d * e * e == L * L and a * b + c * e == 0
inp = f"{OUT}/in_{d}_{L}.json"; cert = f"{OUT}/cert_{d}_{L}.json"
json.dump({"d": d, "D": L, "units": V}, open(inp, "w"))
print(f"d={d} L={L} |U|={len(V)} -> certifying", flush=True)
r = subprocess.run(["python3", f"{W}/certify_w2.py", inp, cert], capture_output=True, text=True, timeout=30 * 3600)
if os.path.exists(cert):
    c = subprocess.run(["python3", f"{W}/check_w.py", cert], capture_output=True, text=True)
    print(c.stdout.strip() or c.stderr.strip()[-300:], flush=True)
    subprocess.run(["gzip", "-f", "-9", cert])
else:
    print("certificate failed:", (r.stdout + r.stderr)[-400:], flush=True)
