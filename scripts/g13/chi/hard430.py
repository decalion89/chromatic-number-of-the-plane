"""Certify the hard sub-leaf F34_split430_split210_split90_split35_leaf1 of G13 share 33 (FAILED at 16:37: at that
depth every colour-0 variable is assigned, so the re-splits on vars_c0_lex had one cube each). Same resolve() as
share_driver.py, with the plan-D style name (CASE_leafI, which verify_plan_D.py tries first) and the cuts taken on
the colour-1 then colour-2 variables in lex order (work/vars_c12_lex.txt); the cover of every cube file is checked
by cuber2.py and again by verify_plan_D.py. usage: hard430.py JOBS"""
import hashlib, os, sys
D = os.path.dirname(os.path.abspath(__file__))
os.chdir(D)
sys.path.insert(0, D)
import share_driver as sd
want = dict(reversed(l.split()) for l in open("cases/SHA256SUMS") if l.strip())
assert hashlib.sha256(open("cases/F34.cnf", "rb").read()).hexdigest() == want["F34.cnf"], "F34.cnf differs"
leaf = "cases/F34_split430_split210_split90_split35_leaf1.cnf"
assert hashlib.sha256(open(leaf, "rb").read()).hexdigest().startswith("150ee2639f1f516f"), "leaf formula differs"
sd.VARS["F34"] = "work/vars_c12_lex.txt"
jobs = int(sys.argv[1])
sd.say(f"hard leaf F34_split430_split210_split90_split35_leaf1: start, {jobs} jobs, cuts on colours 1 and 2")
sd.resolve("F34", "cases/F34_split430_split210_split90_split35.cnf", [1], 0, jobs, style="leaf",
           limits=[120] + [1200] * 10)
sd.finish("DONE", "F34_split430_split210_split90_split35_leaf1 certified (re-splits included)")
