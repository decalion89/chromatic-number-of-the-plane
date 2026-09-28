"""Plan C for 'G_13 has no proper 5-colouring': split by the largest colour class.

(Plan D, run_plan_D.sh, uses the same lemma but refutes the case formulas F34/F35/F36 and E37 directly by cube
and conquer on 'v in C0' along lex_order(); it needs no orbit lists and is the recommended route.  Plan C, the
explicit enumeration below, is the fallback: its R(I) checks are cheap, but there are >= 6 500 orbits of
maximal 34-sets.)

Notation.  G = G_13, Aut(G) the 4 732 maps z -> l z + a, l conj(z) + a (N(l) = 1).  For a 5-colouring chi, M(chi) is
the size of its largest class, and M* the maximum of M(chi) over all 5-colourings (if there is one).

Lemma (the case split).  Suppose G has a 5-colouring, and let s = M*.  Then 34 <= s <= alpha(G) and there are a
5-colouring chi and an automorphism g such that, with C0 = g(largest class of chi):
  (1) C0 is an independent dominating (= maximal independent) set with exactly s points;
  (2) the indicator of C0, read along lex_order() on its first L positions, is >= that of every image p(C0),
      p in Aut(G)  (lex-leader: take g maximising this prefix over the orbit);
  (3) the other four classes 4-colour G - C0, each has at most s points, and each class with s points is
      dominating in G;
  (4) after g, the colours 1..4 can be renamed by their first appearance along lex_order() (value precedence
      1 < 2 < 3 < 4; unused colours last); this does not touch C0, so (1)-(3) still hold;
  (5) if s = 34, the class sizes are exactly 34, 34, 34, 34, 33 (they are at most 34 and sum to 169), every
      34-class is dominating by (3), and C0 can be any 34-class: then the 33-class is named 4 and the other
      three 34-classes are renamed 1..3 by first appearance (g13cnf.py --rigid34).
Proof. s >= 34 since 5 * 33 < 169.  Take chi with M(chi) = s and C0 a largest class.  If some vertex outside a
class C with |C| = s had no neighbour in C, moving it into C would give a 5-colouring with a class of s + 1
points, against the choice of s; so every class with s points (C0 included) is dominating.  Automorphisms map
colourings to colourings and preserve all of this, so g can be chosen for (2); the renaming of (4) and (5)
comes after g and changes neither C0 nor the class sizes.  QED

Certificates (every one a kissat DRAT proof checked by drat-trim):
  A(37)   E(37): no independent dominating set with >= 37 points satisfies (2)   [alpha <= 36, so s <= 36]
  B(s)    E(s) with the list L_s blocked: every set satisfying (1) and (2) is one of the lex-leader images of a
          listed orbit representative (enum_cert.py; s = 34, 35, 36)
  Each E formula is split (enum_cert.py --rosette c / --norosette): either the set contains a point p with its
  whole circle p + {N = c}, c in {6, 7, 9, 11} (translate p to 0: Part A_c, no lex-leader, every image of a
  listed set containing {0} + C_c blocked; formulas of ~1 200 variables, E(37)'s four refuted in < 1 s), or it
  contains no such point (Part B: a property invariant under Aut, so lex-leader stays sound).
  R(I)    for every listed representative I of size s: G - I has no 4-colouring with classes of at most s
          points (rest4.py; this implies (3) is impossible)
If a 5-colouring existed, the lemma, B(s) and an automorphism would give a listed I with R(I) satisfiable.

The case s = 34 could instead be refuted in one formula (g13cnf.py --rigid34 --lex0 25: classes 0..3 have 34
points and are dominating, class 4 has 33), avoiding the thousands of orbits of maximal 34-sets; kissat did not
finish it in 300 s, so it needs cube and conquer too (cnc_case.sh with work/vars_c0_lex.txt).

usage: plan_C.py lists [SAMPLES [SEED]]  -- sample maximal independent sets (work/misample), merge the orbit
                                            representatives into lists/L{34,35,36}.txt
       plan_C.py formulas OUTDIR         -- write E(t)_A{6,7,9,11}.cnf and E(t)_B.cnf for t = 37, 36, 35, 34
Then certify_R.py for the lists, certify.py for the Part A formulas, cnc_case.sh for the Part B formulas.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def lists(samples=20000, seed=1):
    """sample maximal independent sets of 34, 35, 36 points (work/misample) and merge their canonical forms into
    lists/L{t}.txt (seeded from work/reps{t}.txt, the lists of the test run).  Repeat with other seeds until the
    counts stop growing."""
    sys.path.insert(0, os.path.join(HERE, "work"))
    from orbits import canon
    os.makedirs(os.path.join(HERE, "lists"), exist_ok=True)
    for t in (36, 35, 34):
        out = subprocess.run([os.path.join(HERE, "work", "misample"), os.path.join(HERE, "work", "g13.adj"),
                              str(t), str(samples), str(seed * 1000 + t)], capture_output=True, text=True).stdout
        reps = {}
        for fn in (os.path.join(HERE, "lists", f"L{t}.txt"), os.path.join(HERE, "work", f"reps{t}.txt")):
            if os.path.exists(fn):
                for line in open(fn):
                    S = list(map(int, line.split()))
                    if S:
                        reps.setdefault(canon(S), S)
        before = len(reps)
        for line in out.splitlines():
            S = list(map(int, line.split()))
            if S:
                reps.setdefault(canon(S), S)
        with open(os.path.join(HERE, "lists", f"L{t}.txt"), "w") as fh:
            for c in sorted(reps):
                fh.write(" ".join(map(str, c)) + "\n")
        print(f"t = {t}: {len(reps)} orbits ({len(reps) - before} new from {samples} samples)")


def formulas(outdir):
    """E formulas, each split into Part A_c (0 and its whole circle N = c in S, c = 6, 7, 9, 11; tiny) and
    Part B (no point with its whole circle; the hard part, for cube and conquer).  The R(I) formulas are made
    on the fly by certify_R.py (7 000+ of them would not fit on the disk)."""
    from enum_cert import build
    os.makedirs(outdir, exist_ok=True)
    for t, tmax in ((37, 0), (36, 36), (35, 35), (34, 34)):
        sets = []
        if t <= 36:
            sets = [list(map(int, l.split())) for l in open(os.path.join(HERE, "lists", f"L{t}.txt")) if l.strip()]
        for c in (6, 7, 9, 11):
            # Part A_c: no lex-leader (L = 0); block every image of each listed set that contains 0 and C_c
            cnf, nb = build(t, tmax, 0, [], rosette=c)
            nb = add_rosette_blocks(cnf, sets, c)
            open(os.path.join(outdir, f"E{t}_A{c}.cnf"), "w").write(cnf.text())
        cnf, nb = build(t, tmax, 25, sets, norosette=True)
        open(os.path.join(outdir, f"E{t}_B.cnf"), "w").write(cnf.text())
        print(f"t = {t}: {len(sets)} listed orbits; E{t}_A{{6,7,9,11}}.cnf, E{t}_B.cnf ({nb} blocking clauses)")


def add_rosette_blocks(cnf, sets, c):
    """block, in Part A_c, every image of a listed set that contains 0 and the circle N = c"""
    import numpy as np
    from enum_cert import AUT
    from g13 import INDEX, NV, circle
    core = [INDEX[(0, 0)]] + [INDEX[q] for q in circle(c)]
    nb = 0
    for S in sets:
        ind = np.zeros((len(AUT), NV), dtype=bool)
        ind[np.repeat(np.arange(len(AUT)), len(S)), AUT[:, S].ravel()] = True
        rows = np.nonzero(ind[:, core].all(axis=1))[0]
        for img in {frozenset(np.nonzero(ind[g])[0].tolist()) for g in rows}:
            cnf.add([-(1 + v) for v in sorted(img)])
            nb += 1
    return nb


if __name__ == "__main__":
    if sys.argv[1] == "lists":
        lists(int(sys.argv[2]) if len(sys.argv) > 2 else 20000, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    elif sys.argv[1] == "formulas":
        formulas(sys.argv[2])
