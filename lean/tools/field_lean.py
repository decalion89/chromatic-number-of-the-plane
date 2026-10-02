"""field_lean.py: write lean/Sqrt{d}.lean, the Lean proof that the plane over Q(sqrt d) has chromatic number 4, from
the published data. The upper bound is a reduction at 7 (d a nonzero square mod 7, or d = 7d' with 7 not dividing d')
or at 2 (d = 3 mod 8).

    python3 lean/tools/field_lean.py              # writes Sqrt{d}.lean for every d in FIELDS
    python3 lean/tools/field_lean.py 11 191       # writes Sqrt11.lean and Sqrt191.lean
    python3 lean/tools/field_lean.py --check      # exits 1 if a file differs from what the data give, or if the
                                                  # 4-colouring of F_7^2 in QuadraticPlanes.lean is not the one
                                                  # of data/quadratic_planes/finite_planes.json

The lower bound uses the graph of data/quadratic_planes/q{d}.json (points [a, b, c, e] with denominator D, edges,
fixed edge), its formula q{d}.cnf and the LRAT proof q{d}.lrat (kissat, then drat-trim -L). The general parts are
in QuadraticPlanes.lean."""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data", "quadratic_planes")
# The fields with a Lean proof: those whose LRAT proof the kernel checks in a few minutes and a few GB.
FIELDS = [11, 119, 131, 179, 191, 251, 431, 455, 935]


def chunks(items, per):
    return [items[k:k + per] for k in range(0, len(items), per)]


def source(d):
    g = json.load(open(os.path.join(DATA, f"q{d}.json")))
    assert g["d"] == d
    D = g["D"]
    P = [tuple(p) for p in g["points"]]; E = [tuple(e) for e in g["edges"]]; n = len(P)
    u0, v0 = g["fixed_edge"]
    assert (u0, v0) == tuple(E[0])
    for i, j in E:
        a, b, c, e = (P[i][k] - P[j][k] for k in range(4))
        assert a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0
    if d % 7 == 0:
        assert (d // 7) % 7 != 0
        lemma = "QuadraticPlanes.colorable_four_ramified"
        upper_doc = f"`{d} = 7 · {d // 7}`, so Moorhouse's reduction at 7 applies (the ramified case)."
        upper = f"colorable_four_ramified (d' := {d // 7}) (by norm_num) (by norm_num)"
    elif any((s * s - d) % 7 == 0 for s in range(1, 4)):
        s = [s for s in range(1, 4) if (s * s - d) % 7 == 0][0]
        lemma = "QuadraticPlanes.colorable_four"
        upper_doc = f"`{d} ≡ {s}² (mod 7)`, so Moorhouse's reduction at 7 applies."
        upper = f"colorable_four (s := {s}) (by norm_num) (by norm_num)"
    else:
        assert d % 8 == 3, f"{d} is not 0 or a nonzero square modulo 7, nor 3 modulo 8"
        lemma = "QuadraticPlanes.colorable_four_two"
        upper_doc = f"`{d} ≡ 3 (mod 8)`, so the reduction at 2 applies (Fischer, Theorem 10)."
        upper = "colorable_four_two (by norm_num)"
    with open(os.path.join(DATA, f"q{d}.cnf")) as fh:
        header = fh.readline().split()
    nvars, nclauses = int(header[2]), int(header[3])
    assert nvars == 3 * n
    pts = ",\n    ".join(", ".join(f"({a}, {b}, {c}, {e})" for a, b, c, e in row) for row in chunks(P, 4))
    edges = ",\n    ".join(", ".join(f"({i}, {j})" for i, j in row) for row in chunks(E, 10))
    args = "\n    ".join(" ".join(f"(V {i})" for i in row) for row in chunks(list(range(nvars)), 12))
    return f'''import QuadraticPlanes

/-!
# The plane over `ℚ(√{d})` has chromatic number 4

**Theorem** (`Sqrt{d}.chromaticNumber_eq_four`). The unit-distance graph on `ℚ(√{d})²` has chromatic number 4
(`notes/quadratic_planes.md`). This file is written from the data by `tools/field_lean.py`.

*Lower bound* (`not_colorable_three`). The graph of `data/quadratic_planes/q{d}.json` has {n} vertices `[a, b, c, e]`,
standing for `((a + b√{d})/{D}, (c + e√{d})/{D})`, and {len(E)} edges. The kernel checks that each edge has length 1
(`checkEdges_E`, `QuadraticPlanes.adj_pt`), and Mathlib's `lrat_proof` checks the LRAT proof
`data/quadratic_planes/q{d}.lrat` that the colouring formula `q{d}.cnf` ({nvars} variables, {nclauses} clauses) is
unsatisfiable. A 3-colouring would satisfy it (`QuadraticPlanes.clause_facts`).

*Upper bound* (`{lemma}`): {upper_doc}
-/

namespace Sqrt{d}

open LocalColouring QuadraticPlanes

set_option maxHeartbeats 4000000 in
/-- The {n} points of `data/quadratic_planes/q{d}.json`, as `[a, b, c, e]` with denominator {D}. -/
def points : List (ℤ × ℤ × ℤ × ℤ) := [
    {pts}]

/-- The point with index `i`. -/
def P (i : Fin {n}) : ℤ × ℤ × ℤ × ℤ := points.getD i.val (0, 0, 0, 0)

set_option maxHeartbeats 4000000 in
/-- The {len(E)} edges of `data/quadratic_planes/q{d}.json`: all the unit distances among the points. -/
def E : List (Fin {n} × Fin {n}) := [
    {edges}]

lemma checkEdges_E : checkEdges {d} {D} P E = true := by decide +kernel

lemma noLoops_E : noLoops E = true := by decide +kernel

-- `q{d}.cnf` is unsatisfiable: for all propositions `x₀, …, x_{nvars - 1}`, some clause is false. The statement is
-- a disjunction, over the clauses, of the negations of the clauses.
lrat_proof refuted
  (include_str "../data/quadratic_planes/q{d}.cnf")
  (include_str "../data/quadratic_planes/q{d}.lrat")

set_option maxHeartbeats 10000000 in
/-- The graph of `q{d}.json` is not 3-colourable. -/
theorem graph_not_colorable : ¬ (edgeGraph E).Colorable 3 := by
  rintro ⟨C⟩
  obtain ⟨V, hVert, hEdge, hVu, hVv⟩ := clause_facts noLoops_E (u₀ := {u0}) (v₀ := {v0}) (List.mem_cons_self ..) C
  have H := refuted
    {args}
  casesm* _ ∨ _
  all_goals first
    | exact hVert _ ‹_› (by norm_num) (by norm_num)
    | exact hEdge _ _ ‹_› (by decide +kernel)
    | exact ‹¬V {3 * u0}› hVu
    | exact ‹¬V {3 * v0 + 1}› hVv

/-- **Lower bound.** The unit-distance graph of `ℚ(√{d})²` is not 3-colourable. -/
theorem not_colorable_three : ¬ (unitDistGraph (L {d})).Colorable 3 := fun h =>
  graph_not_colorable (h.of_hom (edgeGraph.hom (fun i => pt {d} {D} (P i)) (pt_adj (by norm_num) checkEdges_E)))

/-- **Theorem.** The unit-distance graph of `ℚ(√{d})²` has chromatic number 4. -/
theorem chromaticNumber_eq_four : (unitDistGraph (L {d})).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of ({upper}) not_colorable_three

end Sqrt{d}
'''


def f7_table_matches():
    """The table f7Colour of QuadraticPlanes.lean is the 4-colouring of F_7^2 in finite_planes.json (vertex 7x + y)."""
    planes = json.load(open(os.path.join(DATA, "finite_planes.json")))
    f7 = planes["planes"]["7"]["colouring"]
    src = open(os.path.join(HERE, "..", "QuadraticPlanes.lean"), encoding="utf-8").read()
    start = src.index("def f7Colour : Fin 49 → Fin 4 := ![") + len("def f7Colour : Fin 49 → Fin 4 := ![")
    table = [int(x) for x in src[start:src.index("]", start)].replace(",", " ").split()]
    return table == f7


if __name__ == "__main__":
    args = sys.argv[1:]
    check = "--check" in args
    ds = [int(a) for a in args if a != "--check"] or FIELDS
    bad = 0
    for d in ds:
        src = source(d)
        out = os.path.join(HERE, "..", f"Sqrt{d}.lean")
        if check:
            same = os.path.exists(out) and open(out, encoding="utf-8").read() == src
            print(f"Sqrt{d}.lean is up to date" if same else f"Sqrt{d}.lean differs from the data")
            bad += not same
        else:
            open(out, "w", encoding="utf-8").write(src)
            print(f"wrote Sqrt{d}.lean")
    if check:
        same = f7_table_matches()
        print("the colouring of F_7^2 in QuadraticPlanes.lean is that of finite_planes.json" if same
              else "the colouring of F_7^2 in QuadraticPlanes.lean differs from finite_planes.json")
        bad += not same
    sys.exit(1 if bad else 0)
