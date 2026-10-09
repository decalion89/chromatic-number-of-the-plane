# Review of plan D: "G_13 has no proper 5-colouring"

This is an adversarial review, written inside the project (not independent of it), of `run_plan_D.sh` and the code it calls. I read `g13.py`, `g13cnf.py`, `enum_cert.py`, `plan_C.py`, `cuber2.py`, `certify.py`, `cnc_case.sh`, `resplit.sh` and `summarize.py`. I also compared the plan with `scripts/g17_alpha.py` and `scripts/g17_part_b.py` in the repository. Every claim below was checked against the code and, where possible, tested (section 5).

## Verdict

| Point | Verdict |
|---|---|
| 1. Lemma / exhaustiveness of the case split | **Sound.** One documentation gap: the written lemma does not justify value precedence or the F34 normal form. |
| 2. Encodings (F36, F35, F34, E37) | **Sound.** Every constraint says what its case needs, and brute-force tests pass. |
| 3. E37 (rosette cases A6/A7/A9/A11 + rest B) | **Sound.** Together they cover every independent set of at least 37 points. |
| 4. Cube split and certification | **Each leaf is certified soundly, but the final bookkeeping is not a valid check.** `summarize.py` can exit 0 while leaves, or whole cases, were never certified. It also fails in two situations where the proof is complete. `certify.py`'s resume logic can skip leaves that were never certified. Fixes are given below. |

The mathematics is sound. The weak link is the last step, which should turn "all logs look fine" into "every leaf of every cover is certified". Until that is fixed, "chi(G_13) = 6 is certified" depends on the operator's process, not on a machine check.

---

## 1. The lemma: sound

**Statement used.** Suppose a 5-colouring exists. Let M* be the maximum of the largest class size over all 5-colourings. Then 34 ≤ M*, because 5·33 = 165 < 169. There are four cases.

- **M* ≥ 37.** Some class has at least 37 points, so alpha ≥ 37. Every independent set extends to a maximal one, so a maximal independent set of at least 37 points exists. This is exactly what E37 refutes (section 3). So M* ≤ 36.
- **M* = s ∈ {35, 36}.** Take chi with M(chi) = s and name a largest class colour 0. Every class C with s points is dominating. Otherwise a vertex outside C with no neighbour in C could be moved into C. Its old class stays independent, so this would give a proper 5-colouring with a class of s+1 points, contradicting the maximality of M*. No class is modified, so the other sizes are untouched: dominance comes for free from the global maximality of M*. This handles the "making a class maximal changes the other sizes" trap. The normalisation then goes as follows:
  1. Apply an automorphism g that maximises the 25-prefix of g(C0) over the orbit. Automorphisms preserve properness, class sizes and dominance, so (1) and (3) of the lemma still hold.
  2. Rename colours 1..4 by first appearance along `lex_order()`, with unused colours last. This does not touch C0, so the lex-leader property still holds. The caps and "s-class ⇒ dominating" are symmetric in colours 1..4. The result satisfies `--vp 1234`.

  The order matters: automorphism first, colour renaming second. The code needs exactly that order: vp involves colours 1..4 only and lex0 involves colour 0 only. The classic conflict, where the automorphism used for lex-leader breaks value precedence, cannot arise. Precedence is re-established afterwards without touching C0. Colour 0 is correctly excluded from vp (`--vp 1234`, not `01234`).
- **M* = 34.** Every 5-colouring has all classes of at most 34 points summing to 169 = 5·34 − 1. The sizes are therefore exactly 34, 34, 34, 34, 33, and each 34-class is dominating by the argument above. Call any 34-class C0 and apply the lex-leader automorphism. Colour the 33-class 4, and rename the other three 34-classes 1..3 by first appearance. This is exactly `--rigid34`. Class 4 is correctly *not* assumed dominating, since moving a vertex into it gives 34, which is no contradiction.

**Lex-leader on a prefix, over the whole group.** `lex_leader_sets` compares S with p⁻¹(S) = {v : p(v) ∈ S} on `order[:L]`, for every non-identity p among the 4 732 maps.
- These maps are all automorphisms. I checked that each maps the 1 183 edges onto edges, that all are distinct bijections, and (on a 2 000-pair sample) that they are closed under composition.
- Let g maximise the prefix over the orbit. For every p, p⁻¹g(C0) is in the orbit, so its prefix is at most that of g(C0).
- Positions fixed by p are skipped. They compare a literal with itself, so this is correct.
- The prefix makes the constraint weaker than full lex-leader, never unsound.

**Documentation gap (not a soundness gap).** The lemma in `plan_C.py` lists (1)-(3) only. It justifies neither value precedence nor the F34 normal form (sizes 34, 34, 34, 34, 33, which 34-class becomes C0, vp on 1..3 only). `g13cnf.py`'s docstring points to that lemma for soundness. Add two lines to the lemma:

- (4) After g, rename colours 1..4 (1..3 when s = 34) by first appearance along `lex_order()`. This fixes C0 and keeps (1)-(3).
- (5) If s = 34, the sizes are forced to be 34, 34, 34, 34, 33, and every 34-class is dominating by (3).

## 2. The encodings: sound

What each case formula contains (read in `g13cnf.main()`):

- **F36 / F35.** They contain:
  - the base (ALO, AMO, and one edge clause per colour);
  - class 0 ≥ s (`totalizer_atleast`) and every class ≤ s (`totalizer_atmost` = at least 169−s negations);
  - class 0 dominating;
  - lex-leader of the class-0 indicator over all 4 731 non-identity maps on the first 25 positions of `lex_order()`;
  - vp for 1<2<3<4 along the full `lex_order()` (169 vertices);
  - `--domcap s`: `unary_count` with outputs o_j ↔ (count ≥ j), plus the clause ¬o_s ∨ x(v,c) ∨ ⋁_{w~v} x(w,c) for c = 1..4.
- **F34.** Classes 0..3 are each ≥ 34, ≤ 34 and dominating. Class 4 is ≤ 33, which is exactly 33 given ALO/AMO. There is vp on 1,2,3 and lex0 with L = 25.
- **E37_A_c.** It says: independent, dominating, ≥ 37, with units s₀ and s_q for q ∈ C_c, and no lex (L = 0).
- **E37_B.** It says: independent, dominating, ≥ 37, lex-leader (L = 25, full group), and for each c ∈ {6,7,9,11} and each point p, the clause ¬s_p ∨ ⋁_{q∈C_c} ¬s_{p+q}.

Tests, each passing with 0 failures (details in section 5):

- **Totalizer at-least / at-most.** Exhaustive for n ≤ 7 and all t, including t = 0 and t > n, and with mixed-sign inputs. The cap on node outputs is handled correctly.
- **`unary_count`.** Both directions, o_j ↔ (count ≥ j), exhaustive for n ≤ 7 and all t.
- **`lex_leader_sets`.** 4 020 random cases (random perms, orders, prefixes) against the definition. On the **real** group and order with L = 25: 300 random subsets, where the clauses accept S **iff** its 25-prefix is maximal in its orbit. The lex-max image is always accepted.
- **`value_precedence`.** Exhaustive over all 5-colourings of orders of length ≤ 5, for [1,2,3,4] and [1,2,3].
- **Whole-pipeline analogues.** I loaded the real `g13.py` / `g13cnf.py` in memory with q = 5 and q = 7 (only Q and N0 changed; no file edited).
  - For G_7, M* = 14, and `--big0 14 --cap 14 --dom0 --lex0 L --vp 1234 --domcap 14` is SAT for L = 1, 10, 25, 49. It is UNSAT for 15.
  - For s = 10..14, "big0 s + cap s" and "big0 s + cap s + lex0 + vp" agree.
  - 87 random 5-colourings with a 14-class were scrambled by a random automorphism and colour permutation. After the lemma's normalisation (largest class → 0, lex-max g, rename 1..4), all 87 are accepted by the formula's clauses.
  - I also tested `--rigid34` with 34/33 textually replaced by 10/9 (n = 49 = 4·10 + 9). 27 colourings with sizes (10,10,10,10,9) and dominating 10-classes were normalised with a randomly chosen 10-class as C0. All 81 normalised forms were accepted for L = 10, 25, 49.
- **Real G_13 formulas.** `F36`, `F35` and `F34` were regenerated with the plan's exact options. I dropped only the 4 732 edge clauses of colours 1..4 and fixed C0 to the lex-max image of a listed 36-, 35- or 34-point maximal independent set. All three are SAT: the non-colouring parts accept a genuine normal-form C0.
- **Determinism.** The formulas regenerate byte-for-byte identically to the test copies:
  - F36 = `work/F36.cnf`, F35 = `work/F35.cnf`;
  - F34 = `work/v3_s34_L25.cnf`, whose sampled leaves are in `work/certQ_F34.log`;
  - E37_* = `work/E37*.cnf` = `work/cases_test/E37_*.cnf`;
  - the E37_A hashes equal the VERIFIED hashes in `work/partA.log`.
- **Variables and order.** `work/vars_c0_lex.txt` = [x(v,0) for v in `lex_order()`] and `work/vars_s_lex.txt` = [s_v for v in `lex_order()`]. `lex_order()` is a permutation of the 169 vertices.

## 3. E37: sound

Take any independent set of at least 37 points and extend it to a maximal one, S. The split works as follows:

- **S has a rosette:** some p ∈ S with p + C_c ⊆ S for some c ∈ {6,7,9,11}. Then the translate S − p satisfies E37_A_c. Translations preserve independence, dominance and size, and lex is correctly dropped (L = 0).
- **S has no rosette:** this property is Aut-invariant. l·C_c = C_c because N is multiplicative with N(l) = 1, and conj(C_c) = C_c. So the lex-max image of S satisfies E37_B.

The split is tautological: B's clauses forbid exactly "p ∈ S and p + C_c ⊆ S". So exhaustiveness does not depend on which circles are listed. I computed the circles with no unit distance inside: {1, 6, 7, 9, 11}. C_1 needs no case, because p ∈ S excludes its neighbours p + C_1. Circles with a unit distance cannot lie inside S.

This matches G_17's part A and part B (`g17_alpha.py`): the same "centre + whole independent circle" split. G_17's A_c counts ≥ 58 − 19 in the region, which is equivalent to the units used here. G_17's B uses 0 ∈ S plus the stabiliser lex-leader. G_13's B uses full-group lex-leader on a prefix, plus dominance. Both are sound.

Test: every listed maximal set was normalised as above and fed to the E formula at t = |S|. That covers all 15 orbits of L36, all 170 of L35 and 120 sampled from L34. Rosette sets were translated to A_c and accepted; the others were lex-maxed and accepted by B. B rejected every rosette set. The E37_A formulas that plan D will write have exactly the hashes already VERIFIED in `work/partA.log`. `lists/` is not used by plan D (`enum_cert.py` is called without `--sets`).

## 4. Cube split and certification: per-leaf sound, global check not sound

**What is right**
- **`cuber2.check_cover`** accepts exactly the leaf sets of a complete binary tree that branches on one variable both ways at every node. It rejects missing leaves, duplicates and prefix pairs. I tested adversarial inputs, and removing one leaf from any of the four test icnf files makes it fail. Closed leaves are included and are certified too.
- **`resplit.sh`** writes the sub-formula as CASE.cnf plus the leaf's units. This is byte-identical to certify's leaf formula. The sub-tree is cover-checked, so parent plus sub-tree still covers every assignment.
- **`certify.py`** runs kissat `--unsat` with a binary DRAT proof, then drat-trim on the same file.
  - It logs VERIFIED only on `^s VERIFIED`. This works because `text=True` turns drat-trim's `\r` into newlines, and `NOT VERIFIED` does not match.
  - The header count is correct.
  - kissat's default NORMAL parsing rejects a truncated formula ("trailing zero missing" / "clause missing"). drat-trim reads formula and proof from separate files, so a lemma can never be taken as an input clause.
  - drat-trim ignoring unit deletions is sound.
  - SAT leaves are logged and their model is saved, and `check_colouring.py` checks it independently.

**Gaps (bookkeeping), with concrete demonstrations** (synthetic logs, section 5)

- **G1. `summarize.py` does not check completeness, and the plan never passes it the main logs.** Step 5 runs `summarize.py logs/E37_A.log cases/*_leaf*.certlog`.
  - `cases/E37_B.certlog`, `F36.certlog`, `F35.certlog` and `F34.certlog` are never read, and no `.icnf` is read.
  - Demo (scenario A): F36 has 4 leaves; leaf 3 has no line at all, and E37_B/F35/F34 have no logs. `summarize.py` prints "4 VERIFIED … 2 VERIFIED" and **exits 0**.

  In a clean single run, `set -e` makes a missing line unlikely. But the gaps below can produce one, and the final check would not see it.
- **G2. False failures make exit 0 unreachable in normal situations, which invites manual overrides.**
  - With **no** re-split at all, the glob stays literal. `summarize.py` then crashes with FileNotFoundError and exits 1 (scenario C).
  - After a **second-level** re-split, as step 5 advises, the intermediate log `cases/F34_leafN.certlog` matches the glob and keeps its timeout line. `summarize.py` then exits 1 forever (scenario B).

  The "superseded" rule exists only as a comment.
- **G3. `certify.py` resumes by name only.** `done` is the set of names with a VERIFIED line; the logged sha256 is never compared. If `CASE.icnf` is regenerated while `CASE.certlog` is kept, the old VERIFIED names skip the new leaves. That happens with another depth, another `CUBER_ARGS` such as `--minpos` (advertised in `cnc_case.sh`), another var file, or a regenerated case formula.
  - Demo: I re-cubed a tiny case with a different variable order and kept the log. `certify.py` then ran **0** leaves, and the new tree was never certified. The plan's summary would not notice (G1).
- **G4. The cover is checked only when it is written, and the file is not written atomically.** `cnc_case.sh` reuses any existing `CASE.icnf` (`[ -f … ] ||`) without re-checking it.
  - If cuber2 is killed, or the disk fills during the write, a truncated icnf remains. The documented "rerun, it is resumable" then certifies only its leaves.
  - The check is also an `assert`, which disappears under `PYTHONOPTIMIZE` / `-O`.
- **G5. Nothing ties the certified formulas to the reviewed code at the end.** `cases/SHA256SUMS` is written but never checked. The generation block is guarded only by `cases/E37_B.cnf`, so an interrupted generation is never redone.
- **G6 (minor hardening).**
  - Proofs are deleted, so the logs are the only evidence.
  - No kissat or drat-trim version or binary hash is logged; kissat here is 4.0.4.
  - Two concurrent `certify.py` runs on the same case would share leaf file names in `proofs/`.

**Smallest fix.** Replace step 5 with a tree-aware verifier. It is cheap: about 1 s per case formula to regenerate, plus hashing and cover checks. The core below was tested on a tiny certified case: complete run → exit 0; a leaf with no line → FAIL; a certified re-split → exit 0; stale log after re-cubing → FAIL.

```python
# verify_plan_D.py (run in g13/ after run_plan_D.sh): exit 0 only if every leaf of every cover is certified
import hashlib, os, re, subprocess, sys, tempfile
sys.setrecursionlimit(100000)
from cuber2 import check_cover

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

def verified(log):                       # name -> set of sha256 with drat-trim VERIFIED
    ok = {}
    for l in (open(log) if os.path.exists(log) else []):
        m = re.match(r"(\S+): sha256 ([0-9a-f]{64}); .*; drat-trim VERIFIED", l)
        if m: ok.setdefault(m.group(1), set()).add(m.group(2))
    return ok

def cubes_of(icnf):
    out = []
    for line in open(icnf):
        if line.startswith("a "): out.append(line.split()[1:-1])
        elif line.startswith("c closed"): out.append(line.split()[2:-1])
    return out

def check_case(cnf):
    stem = cnf[:-4]; tag = os.path.basename(stem)
    head, body = open(cnf, "rb").read().split(b"\n", 1)
    nv, nc = map(int, head.split()[2:4])
    cubes = cubes_of(stem + ".icnf")
    if not check_cover([tuple(map(int, c)) for c in cubes]): sys.exit(f"FAIL {stem}.icnf is not a cover")
    ok, pre, n = verified(stem + ".certlog"), {}, 0
    for i, cube in enumerate(cubes):
        if len(cube) not in pre:
            pre[len(cube)] = hashlib.sha256(f"p cnf {nv} {nc + len(cube)}\n".encode() + body)
        h = pre[len(cube)].copy(); h.update("".join(f"{l} 0\n" for l in cube).encode())
        name, H = f"{tag}_leaf{i}", h.hexdigest()          # exactly the leaf file certify.py hashes
        if H in ok.get(name, ()): n += 1; continue
        sub = f"{stem}_leaf{i}.cnf"                       # resplit.sh: byte-identical to the leaf formula
        if not (os.path.exists(sub) and sha(sub) == H): sys.exit(f"FAIL {name}: not VERIFIED and not re-split")
        n += check_case(sub)
    return n

GEN = {f"E37_A{c}": f"enum_cert.py {{}} 37 --L 0 --rosette {c}" for c in (6, 7, 9, 11)}
GEN.update({"E37_B": "enum_cert.py {} 37 --norosette",
            "F36": "g13cnf.py {} --big0 36 --cap 36 --dom0 --lex0 25 --vp 1234 --domcap 36",
            "F35": "g13cnf.py {} --big0 35 --cap 35 --dom0 --lex0 25 --vp 1234 --domcap 35",
            "F34": "g13cnf.py {} --rigid34 --lex0 25"})
with tempfile.TemporaryDirectory() as d:                  # the certified formulas are what the code writes
    for k, cmd in GEN.items():
        out = os.path.join(d, k + ".cnf")
        subprocess.run(["python3"] + cmd.format(out).split(), check=True, capture_output=True)
        if sha(out) != sha(f"cases/{k}.cnf"): sys.exit(f"FAIL cases/{k}.cnf differs from the code's output")
okA = {h for s in verified("logs/E37_A.log").values() for h in s}
for c in (6, 7, 9, 11):
    if sha(f"cases/E37_A{c}.cnf") not in okA: sys.exit(f"FAIL E37_A{c} not VERIFIED")
for case in ("E37_B", "F36", "F35", "F34"):
    print(case, check_case(f"cases/{case}.cnf"), "leaves VERIFIED (re-split sub-leaves included)")
print("PLAN D FULLY CERTIFIED")
```

Also make these three small changes:

1. **`certify.py`:** skip a leaf only if a VERIFIED line has the same name **and** the sha256 of the leaf formula about to be written. Compute the hash before skipping.
2. **`cuber2.py`:** write the icnf to `OUT.tmp`, then `os.replace`. Replace `assert check_cover(...)` with an explicit `sys.exit` on failure.
3. **`run_plan_D.sh`:** guard generation on `cases/SHA256SUMS` and run `sha256sum -c`. Log `kissat --version` and `sha256sum` of both binaries. Optionally use a per-process temporary workdir in `certify.py`.

## 5. What I ran (all under `nice -n 10`; each a few seconds, except the list check in section 3, which took 85 s)

- **Graph and group sanity:** 169 vertices, 14 units, 1 183 edges, 14-regular, GEN of order 14. The 4 732 maps are edge-preserving distinct bijections, closed under composition on a sample. `lex_order()` is a permutation, and the var files match it.
- **Encoding brute force:** totalizer at-least/at-most, `unary_count`, `lex_leader_sets` (random, and exact on the real group with L = 25), `value_precedence`. 0 failures.
- **Analogues G_5 and G_7:** the real modules loaded in memory with Q/N0 changed.
  - F_s SAT exactly at s = M*; lex0 plus vp agree with the plain formula for s = 10..14.
  - 87/87 normalised random colourings accepted.
  - The rigid34 logic with 10/9 accepted 81/81 normalised colourings.
- **Real G_13 formulas:**
  - regeneration is byte-identical to the test copies;
  - relaxed F36/F35/F34 are SAT with a genuine normal-form C0;
  - E formulas accept the normalised form of all L36 and L35 orbits and 120 L34 orbits, and B rejects rosette sets.
- **Certification tools on a tiny UNSAT case** (pigeonhole 5→4): real `cuber2.py` and `certify.py`, with a separate `--workdir` and log. The verifier sketch passes the complete run and fails on a missing leaf and on stale logs. It passes a certified re-split, and `certify.py` skips all re-cubed leaves by name (G3).
- **Synthetic logs through the real `summarize.py`:**
  - scenario A (missing leaf, three cases absent): exit 0;
  - scenario B (nested re-split, complete): exit 1;
  - scenario C (no re-split, complete): FileNotFoundError, exit 1.

**Housekeeping.** Temporary files went to `scratchpad/review_tmp/` and were removed. Importing modules from `g13/` refreshed `g13/__pycache__/g13cnf.cpython-311.pyc`, which was stale (older than the 00:34 edit of `g13cnf.py`). I deleted the `cuber2` `.pyc` my import created. No other file in `g13/` was touched except this report.

## Open questions

- **Q1.** The plan will certify every leaf again rather than reuse the test-run certificates. F34 is byte-identical to `work/v3_s34_L25.cnf` and F36/F35 to `work/F36.cnf` / `work/F35.cnf`, so reuse would be possible. That is fine; just do not mix test logs into `cases/` without the hash-matching verifier.
- **Q2.** alpha(G_13) ≤ 36 rests entirely on E37_B being refuted. If E37_B turns out SAT, the model gives a 37-point independent set. The plan would then need a colouring case for s ≥ 37. That would not make it unsound, only incomplete.
