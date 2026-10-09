"""scan_all.py DMIN DMAX OUT [workers] [CERTDIR]: for every squarefree d = 11 (mod 12) in [DMIN, DMAX], rank the
denominators D <= 2400 with 6 | D by the number of unit vectors with denominator D (keeping 24 <= |U_D/+-| <= 110),
run the relation-space Lemma W test on the 12 richest (120 s each), and record the first INFEASIBLE one.  With CERTDIR,
build the exact certificate for that (d, D) with certify_w2 (full unit set) and run the independent checker on it;
the line then ends with the checker's verdict.  One output line per d; resumable (skips d already in OUT)."""
import sys, os, time
from fractions import Fraction as Fr
from multiprocessing import Pool
from kappaD import units, halve
from kapparel import relation_basis, solve_rel

def squarefree(n):
    k = 2
    while k * k <= n:
        if n % (k * k) == 0: return False
        k += 1
    return True

CERTDIR = sys.argv[5] if len(sys.argv) > 5 else None

def one(d):
    t0 = time.time()
    ranked = sorted(((len(halve(units(d, D))), D) for D in range(6, 2401, 6)), reverse=True)
    ranked = [(n, D) for n, D in ranked if 24 <= n <= 110]
    tried = []
    for n, D in ranked[:12]:
        V = halve(units(d, D))
        name, _ = solve_rel(len(V), relation_basis(V), Fr(1, 3), 120)
        tried.append(f"{D}:{name[:4]}")
        if name == "INFEASIBLE":
            verdict = ""
            if CERTDIR:
                import json, subprocess
                inp = f"{CERTDIR}/in_{d}_{D}.json"; cert = f"{CERTDIR}/cert_{d}_{D}.json"
                json.dump({"d": d, "D": D, "units": [list(u) for u in V]}, open(inp, "w"))
                here = os.path.dirname(os.path.abspath(__file__))
                r = subprocess.run(["python3", f"{here}/certify_w2.py", inp, cert], capture_output=True, text=True, timeout=3 * 3600)
                if os.path.exists(cert):
                    c = subprocess.run(["python3", f"{here}/check_w.py", cert], capture_output=True, text=True)
                    verdict = " | " + (c.stdout.strip() or c.stderr.strip()[-200:])
                    subprocess.run(["gzip", "-f", "-9", cert])
                else:
                    verdict = " | certificate failed: " + (r.stdout + r.stderr).strip()[-200:]
            return f"d={d} INFEASIBLE at D={D} (|U|={n}) tried {' '.join(tried)} [{time.time()-t0:.0f}s]{verdict}"
    return f"d={d} no infeasible D found; tried {' '.join(tried)} [{time.time()-t0:.0f}s]"

if __name__ == "__main__":
    lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    CERTDIR = sys.argv[5] if len(sys.argv) > 5 else None
    done = set()
    if os.path.exists(out):
        for line in open(out):
            if line.startswith("d="): done.add(int(line.split()[0][2:]))
    ds = [d for d in range(lo, hi + 1) if d % 12 == 11 and squarefree(d) and d not in done]
    with Pool(workers) as p:
        for line in p.imap_unordered(one, ds):
            with open(out, "a") as f: f.write(line + "\n")
            print(line, flush=True)
