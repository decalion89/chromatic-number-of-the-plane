"""chi(G_13) = 6 (notes/g13_chi.md): the code in scripts/g13/chi and the checks around it, without kissat or drat-trim.

- The code: the copies of g13.py and enum_cert.py in scripts/g13/chi are those of scripts/g13, which prove
  alpha(G_13) = 36; g13cnf.py and certify.py differ from theirs only in their docstrings, and cuber2.py only in how it
  propagates the unit clauses of a leaf formula. The cube file of E37_B is certificates/g13_alpha_part_b_cubes.icnf.
  The maps of the lex-leader clauses are automorphisms, and lex_order() starts as the note says.
- The eight case formulas, written by the commands of share_driver.py, have the SHA-256 of cases/SHA256SUMS, and the
  four formulas E37_A_c those of certificates/g13_alpha_part_a_checks.txt.
- F36, F35 and F34, read from their text. Every clause on the colour variables alone is of a known kind, in the right
  number, and every other clause has the auxiliary variables of one count, of the chains or of value precedence. The
  counts and the chains pass the functions of scripts/g13/g13_audit.py, which audits the formulas of alpha(G_13).
- F36, F35 and F34 against intended models. A colouring is built from a known maximal independent set: the set,
  moved to its lex-leader image, is class 0; the other vertices, along lex_order(), go greedily to the class with the
  fewest neighbours among those not yet full; the colours are then renamed by first appearance. Every variable gets
  its intended value, the auxiliary ones in the order in which g13cnf.py numbers them, and the clauses that fail are
  exactly the edge clauses of the monochromatic edges and the domination clauses of the vertices that a class which
  must dominate does not dominate. For colourings far from a normal form, the clauses that fail are exactly those that
  the meaning of each kind of clause predicts.
- The lex-leader clauses: for each orbit of the 15 known 36-point sets and of eight maximal sets of 35 and 34 points,
  the lex-leader image satisfies every lex-leader clause of F36, F35 or F34, and the lex-least image does not.
- F36 with one wrong clause of each kind fails one of these checks.
- The four cube files are covers, with the counts of notes/g13_chi.md.
- check_colouring.py accepts the stored 6-colouring and rejects an improper 5-colouring.
- verify_plan_D.py, on a tiny plan made up here (stand-ins for the formula writers, and trees with both kinds of
  re-split), accepts the certified plan and rejects each of nine broken copies. pack_certificates.py packs it
  deterministically, and scripts/verify_g13_chi.py checks the archive, prints its numbers, fills placeholders, and
  rejects a changed archive and the archive of a broken plan.
- Once certificates/g13_chi_certlogs.tar.gz is there: it matches its SHA256SUMS file, has only logs of certify.py
  (none with a satisfiable leaf) and cube files of re-splits, and its logs of E37 are the certificates of
  alpha(G_13) = 36.
"""
import ast
import collections
import gzip
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
G13 = os.path.join(ROOT, "scripts", "g13")
CHI = os.path.join(G13, "chi")
sys.path.insert(0, os.path.join(ROOT, "scripts", "experiments"))
from largest_sets_whole_circles import KNOWN_G13  # noqa: E402

# maximal independent sets of 35 and 34 points, in different orbits, from a sample of maximal independent sets
# (plan C, whose lists hold 170 orbits of 35 points and 2 862 of 34); the tests check what they need of them
MAXIMAL_35 = [
    [0, 1, 2, 3, 4, 5, 7, 9, 11, 19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 33, 37, 47, 49, 59, 61, 63, 87, 92, 95, 126,
     150, 152, 154, 164, 166],
    [0, 1, 2, 3, 5, 11, 17, 25, 34, 45, 49, 57, 59, 60, 61, 63, 66, 67, 69, 77, 78, 81, 83, 84, 88, 89, 92, 93, 95, 99,
     103, 138, 147, 155, 164],
    [0, 1, 2, 4, 6, 20, 24, 27, 29, 30, 31, 38, 48, 58, 60, 62, 63, 65, 69, 79, 87, 88, 93, 96, 104, 107, 121, 126,
     127, 131, 150, 151, 153, 155, 165],
    [0, 1, 3, 6, 8, 11, 15, 18, 22, 25, 26, 27, 29, 33, 37, 58, 60, 72, 81, 83, 87, 89, 91, 92, 98, 106, 108, 114, 116,
     137, 148, 152, 160, 163, 166]]
MAXIMAL_34 = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 22, 25, 36, 37, 69, 80, 81, 83, 84, 91, 92, 98, 99, 100, 101, 102, 103, 108, 127, 128,
     146, 147, 148, 166, 167],
    [0, 1, 2, 3, 5, 11, 17, 21, 23, 27, 31, 33, 37, 47, 49, 57, 59, 63, 66, 67, 69, 85, 89, 92, 95, 124, 125, 128, 144,
     150, 152, 154, 160, 164],
    [0, 1, 2, 4, 5, 11, 21, 23, 27, 30, 33, 37, 44, 47, 49, 59, 61, 63, 66, 69, 92, 108, 124, 125, 126, 128, 131, 135,
     140, 147, 150, 152, 154, 164],
    [0, 1, 3, 5, 11, 15, 21, 23, 30, 32, 33, 35, 38, 47, 61, 79, 82, 90, 91, 94, 96, 97, 101, 105, 106, 108, 111, 116,
     125, 126, 128, 136, 140, 168]]


def load(name, folder=CHI):
    """folder/name.py as a module of its own name: the copies in scripts/g13/chi have the names of the modules of
    scripts/g13, which the other tests import, so the modules they import on the way are taken out again"""
    saved_path, saved = list(sys.path), set(sys.modules)
    spec = importlib.util.spec_from_file_location(f"g13chi_{name}", os.path.join(folder, name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, folder)
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path[:] = saved_path
        for m in set(sys.modules) - saved:
            if m in ("g13", "g13cnf", "cuber2", "enum_cert"):
                del sys.modules[m]
    return module


g13 = load("g13")
gc = load("g13cnf")
NV, K, X = g13.NV, g13.K, gc.X
EDGES, ADJ = g13.EDGES, g13.ADJ
ORDER = gc.lex_order()
MAPS = np.array(g13.automorphisms(), dtype=np.int64)
PERMS = [list(p) for p in MAPS.tolist() if any(p[v] != v for v in range(NV))]    # as g13cnf.py --lex0 takes them
SUMS = dict(reversed(line.split()) for line in open(os.path.join(CHI, "cases", "SHA256SUMS")) if line.strip())


def sha256(text):
    return hashlib.sha256(text.encode() if isinstance(text, str) else text).hexdigest()


# ------------------------------------------------------------------------------------------------------ the code

def without_docstring(path):
    tree = ast.parse(open(path).read())
    tree.body = tree.body[1:] if isinstance(tree.body[0], ast.Expr) and isinstance(tree.body[0].value, ast.Constant) \
        else tree.body
    return ast.dump(tree)


def test_the_copies_are_the_code_of_alpha():
    for name in ("g13.py", "enum_cert.py"):
        assert open(os.path.join(CHI, name), "rb").read() == open(os.path.join(G13, name), "rb").read()
    for name in ("g13cnf.py", "certify.py"):
        assert without_docstring(os.path.join(CHI, name)) == without_docstring(os.path.join(G13, name))
    # cuber2.py: the same cover check; the tree builder passes the unit clauses of a formula as assumptions
    chi, alpha = (ast.parse(open(os.path.join(d, "cuber2.py")).read()) for d in (CHI, G13))
    function = lambda tree, name: ast.dump(next(f for f in tree.body if getattr(f, "name", None) == name))
    assert function(chi, "check_cover") == function(alpha, "check_cover")
    assert function(chi, "build") != function(alpha, "build") and "units + path" in open(os.path.join(CHI, "cuber2.py")).read()
    assert (open(os.path.join(CHI, "cases", "E37_B.icnf"), "rb").read()
            == open(os.path.join(ROOT, "certificates", "g13_alpha_part_b_cubes.icnf"), "rb").read())


def test_the_graph_the_maps_and_the_order():
    """the maps of the lex-leader clauses are automorphisms, and the prefix of lex_order() is 0, the circle N = 2 and
    ten points of the circle N = 3"""
    adjacency = np.zeros((NV, NV), dtype=bool)
    for a, b in EDGES:
        adjacency[a, b] = adjacency[b, a] = True
    e = np.array(EDGES)
    assert len(EDGES) == 1183 and len(MAPS) == 4732 and len({tuple(p) for p in PERMS}) == 4731
    assert adjacency[MAPS[:, e[:, 0]], MAPS[:, e[:, 1]]].all()
    assert sorted(ORDER) == list(range(NV))
    assert [g13.norm(g13.POINTS[v]) for v in ORDER[:25]] == [0] + [2] * 14 + [3] * 10


# -------------------------------------------------------------------------------------------------- the formulas

@pytest.fixture(scope="module")
def formulas(tmp_path_factory):
    """the eight case formulas, written by the commands of share_driver.py (GEN) with the code in scripts/g13/chi:
    name -> text"""
    folder = tmp_path_factory.mktemp("g13_chi_formulas")
    out = {}
    for name, command in load("share_driver").GEN.items():
        path = str(folder / f"{name}.cnf")
        args = command.format(path).split()
        subprocess.run([sys.executable, os.path.join(CHI, args[0])] + args[1:], check=True, capture_output=True)
        with open(path) as fh:
            out[name] = fh.read()
        os.remove(path)
    return out


def test_the_case_formulas_are_the_certified_ones(formulas):
    assert sorted(formulas) == sorted(name[:-4] for name in SUMS) and len(SUMS) == 8
    for name, text in formulas.items():
        assert sha256(text) == SUMS[name + ".cnf"], name
    alpha = {}
    for line in open(os.path.join(ROOT, "certificates", "g13_alpha_part_a_checks.txt")):
        m = re.match(r"(E37_A\d+)\.cnf: sha256 ([0-9a-f]{64}); ", line)
        if m:
            alpha[m.group(1)] = m.group(2)
    assert alpha == {f"E37_A{c}": SUMS[f"E37_A{c}.cnf"] for c in (6, 7, 9, 11)}
    heads = {name: tuple(map(int, text.split("\n", 1)[0].split()[2:4])) for name, text in formulas.items()}
    assert (heads["F36"], heads["F35"], heads["F34"], heads["E37_B"]) == (
        (125416, 484970), (125386, 483890), (124192, 451390), (114055, 351802))


class Formula:
    """the clauses of a DIMACS text, as one array of literals and the ends of the clauses"""

    def __init__(self, text):
        head, body = text.split("\n", 1)
        self.nv, nc = map(int, head.split()[2:4])
        self.lits = np.array(body.split(), dtype=np.int64)
        self.ends = np.flatnonzero(self.lits == 0)
        self.starts = np.concatenate(([0], self.ends[:-1] + 1))
        assert len(self.ends) == nc and len(self.lits) == self.ends[-1] + 1

    def satisfied(self, val):
        """for each clause, whether it holds when variable i has the value val[i]"""
        val = np.array(val, dtype=bool)
        assert len(val) == self.nv + 1
        truth = val[np.abs(self.lits)] ^ (self.lits < 0)
        truth[self.ends] = False
        return np.logical_or.reduceat(truth, self.starts)

    def clause(self, k):
        return frozenset(self.lits[self.starts[k]:self.ends[k]].tolist())

    def violated(self, val):
        return {self.clause(k) for k in np.flatnonzero(~self.satisfied(val))}


@pytest.fixture(scope="module")
def parsed(formulas):
    return {name: Formula(formulas[name]) for name in ("F36", "F35", "F34")}


class Intended:
    """the intended value of every variable of a formula of g13cnf.py, for a colouring col (any map of the vertices
    to the colours 0..4): x(v, c) = [col[v] = c], then the auxiliary variables, in the order in which g13cnf.py
    numbers them; it also records the root of each count, the variables of each chain, and those of value
    precedence"""

    def __init__(self, col):
        self.col = col
        self.val = [False] + [col[v] == c for v in range(NV) for c in range(K)]
        self.units, self.counts, self.chain, self.precede, self.ranges = [], [], [], [], {}

    def value(self, lit):
        return self.val[lit] if lit > 0 else not self.val[-lit]

    def new(self, value):
        self.val.append(bool(value))
        return len(self.val) - 1

    def counter(self, lits, t, unit=True):
        """the tree of totalizer_atleast (and of unary_count, unit=False): halves, the left one first; each node gets
        its outputs after its children, o_1, ..., o_m with m = min(t, the outputs of its children), o_j = [at least j
        of its leaves are true]; totalizer_atleast adds the unit clause o_t of the root. Returns the root's outputs."""
        def build(items):
            if len(items) == 1:
                return [items[0]], int(self.value(items[0]))
            (left, a), (right, b) = build(items[:len(items) // 2]), build(items[len(items) // 2:])
            return [self.new(a + b >= j) for j in range(1, min(t, len(left) + len(right)) + 1)], a + b
        start = len(self.val)
        root = build(list(lits))[0]
        self.counts.append((start, len(self.val), list(lits), t, unit))
        if unit:
            self.units.append(root[t - 1])
        return root

    def chains(self):
        """lex_leader_sets on class 0 and the first 25 positions of lex_order(): for each map, the pairs (v, p(v)) of
        the positions it moves, and e_k = [the first k + 1 pairs agree] for every pair but the last"""
        start = len(self.val)
        for p in PERMS:
            pairs = [(v, p[v]) for v in ORDER[:25] if p[v] != v]
            agree, es = True, []
            for v, w in pairs[:-1]:
                agree = agree and self.val[X(v, 0)] == self.val[X(w, 0)]
                es.append(self.new(agree))
            self.chain.append((pairs, es))
        self.ranges["chains"] = (start, len(self.val))

    def precedence(self, colours):
        """value_precedence along lex_order(): for each pair of consecutive colours a, b and each position i,
        y_i = [a is used among the first i + 1 positions]"""
        start = len(self.val)
        for a, b in zip(colours, colours[1:]):
            used, ys = False, []
            for v in ORDER:
                used = used or self.val[X(v, a)]
                ys.append(self.new(used))
            self.precede.append((b, ys))
        self.ranges["precedence"] = (start, len(self.val))


def intended(name, col):
    """the intended assignment of F36, F35 or F34 for the colouring col, in the order of g13cnf.main()"""
    m = Intended(col)
    cls = lambda c: [X(v, c) for v in range(NV)]
    if name in ("F36", "F35"):          # --big0 s --cap s --dom0 --lex0 25 --vp 1234 --domcap s
        s = int(name[1:])
        m.size0 = m.counter(cls(0), s)                          # totalizer_atleast
        for c in range(K):
            m.counter([-x for x in cls(c)], NV - s)             # totalizer_atmost: at least NV - s negations
        m.chains()
        m.precedence([1, 2, 3, 4])
        m.domcap = [m.counter(cls(c), s, unit=False)[s - 1] for c in range(1, K)]     # unary_count
    else:                               # --rigid34 --lex0 25
        m.chains()
        for c in range(4):
            m.counter(cls(c), 34)
            m.counter([-x for x in cls(c)], NV - 34)
        m.counter([-x for x in cls(4)], NV - 33)
        m.precedence([1, 2, 3])
    return m


def predicted(name, m):
    """the clauses of formula name that fail in the intended assignment m, from their meaning alone:
    - the edge clause of each monochromatic edge;
    - the domination clause of each vertex that a class which must dominate does not dominate: class 0 (in F34
      classes 0 to 3), and in F36 and F35 a class of colours 1 to 4 with at least s points, whose clauses carry the
      output o_s of its counter;
    - the unit clause of each count whose bound the colouring breaks;
    - for each map whose image of class 0 is lex-greater on the prefix: at the first position where they differ,
      the comparison clause and the two clauses that define the chain variable there;
    - for value precedence of colours a < b: the clause of each position where b is used before a is."""
    col = m.col
    out = {frozenset((-X(a, col[a]), -X(b, col[a]))) for a, b in EDGES if col[a] == col[b]}
    size = collections.Counter(col)
    guard = {c: [] for c in ((0, 1, 2, 3) if name == "F34" else (0,))}
    if name != "F34":
        guard.update({c: [-o] for c, o in zip(range(1, K), m.domcap) if size[c] >= int(name[1:])})
    for c, extra in guard.items():
        for v in range(NV):
            if col[v] != c and all(col[w] != c for w in ADJ[v]):
                out.add(frozenset(extra + [X(v, c)] + [X(w, c) for w in ADJ[v]]))
    out |= {frozenset([o]) for o in m.units if not m.val[o]}
    for pairs, es in m.chain:
        k = next((k for k, (v, w) in enumerate(pairs) if (col[v] == 0) != (col[w] == 0)), None)
        if k is not None and col[pairs[k][1]] == 0:
            (v, w), before = pairs[k], [-es[k - 1]] if k else []
            out.add(frozenset(before + [X(v, 0), -X(w, 0)]))
            if k < len(es):
                out |= {frozenset(before + [X(v, 0), es[k]]), frozenset(before + [-X(w, 0), es[k]])}
    for b, ys in m.precede:
        out |= {frozenset(([ys[i - 1]] if i else []) + [-X(v, b)]) for i, v in enumerate(ORDER)
                if col[v] == b and not (i and m.val[ys[i - 1]])}
    return out


def images(S):
    """the indicator vectors of the images of S under the 4 732 maps, and their first 25 entries along lex_order()
    as numbers (larger is lex-greater)"""
    ind = np.zeros((len(MAPS), NV), dtype=bool)
    ind[np.repeat(np.arange(len(MAPS)), len(S)), MAPS[:, S].ravel()] = True
    return ind, ind[:, ORDER[:25]].astype(np.int64) @ (1 << np.arange(24, -1, -1))


def lex_leader(S):
    ind, key = images(S)
    return set(np.flatnonzero(ind[int(np.argmax(key))]).tolist())


def lex_least(S):
    ind, key = images(S)
    return set(np.flatnonzero(ind[int(np.argmin(key))]).tolist())


def near_colouring(C0, sizes, renamed):
    """class 0 is C0; the other vertices, along lex_order(), go to the class 1, 2, ... not yet full (it holds
    sizes[c - 1] points) with the fewest neighbours already in it; then the colours `renamed` are renamed by their
    first appearance along lex_order()"""
    col = [0 if v in C0 else None for v in range(NV)]
    room = dict(enumerate(sizes, 1))
    for v in ORDER:
        if col[v] is None:
            c = min((c for c in room if room[c]), key=lambda c: (sum(col[w] == c for w in ADJ[v]), c))
            col[v], room[c] = c, room[c] - 1
    return by_first_appearance(col, renamed)


def by_first_appearance(col, renamed):
    first = {c: next((i for i, v in enumerate(ORDER) if col[v] == c), NV) for c in renamed}
    new = dict(zip(sorted(renamed, key=lambda c: (first[c], c)), renamed))
    return [new.get(c, c) for c in col]


SHAPES = {"F36": (KNOWN_G13, (36, 35, 35, 27), (1, 2, 3, 4)), "F35": (MAXIMAL_35, (35, 34, 33, 32), (1, 2, 3, 4)),
          "F34": (MAXIMAL_34, (34, 34, 34, 33), (1, 2, 3))}


@pytest.mark.parametrize("name", ["F36", "F35", "F34"])
def test_every_clause_holds_in_the_intended_models(parsed, name):
    """two near-colourings in normal form per formula: every clause holds but the edge clauses of the monochromatic
    edges and the domination clauses of the vertices that a class which must dominate does not dominate. In F36
    class 1 has 36 points, and in F35 35, so that its domination clauses, guarded by its counter, apply."""
    f = parsed[name]
    sets, sizes, renamed = SHAPES[name]
    for S in sets[:2]:
        col = near_colouring(lex_leader(S), sizes, renamed)
        assert sorted(collections.Counter(col).values(), reverse=True) == sorted((len(S),) + sizes, reverse=True)
        model = intended(name, col)
        assert len(model.val) == f.nv + 1                         # every variable, and no other, gets a value
        bad = f.violated(model.val)
        assert bad == predicted(name, model)
        assert 0 < sum(col[a] == col[b] for a, b in EDGES) < 80   # a near-colouring: a few dozen monochromatic edges


def colourings(name):
    """colourings far from a normal form: at random, with a class 0 of the density 0.1, 0.5 or 0.9, or with the
    class 0 of a normal form; a normal form whose class 0 is lex-least in its orbit, whose colours 1 and 2 are
    swapped, or whose class 0 lacks one point; and one colour on every vertex"""
    sets, sizes, renamed = SHAPES[name]
    rng = np.random.default_rng(13)
    good = near_colouring(lex_leader(sets[2]), sizes, renamed)
    p = max((v for v in range(NV) if good[v] == 0 and v not in ORDER[:25]), key=ORDER.index)   # still lex-leader
    tries = {"random": rng.integers(0, K, NV).tolist(),
             "random but class 0": [0 if c == 0 else int(rng.integers(1, K)) for c in good],
             "lex-least class 0": near_colouring(lex_least(sets[2]), sizes, renamed),
             "colours 1 and 2 swapped": [{1: 2, 2: 1}.get(c, c) for c in good],
             "class 0 one point short": by_first_appearance(
                 [min(renamed, key=good.count) if v == p else c for v, c in enumerate(good)], renamed)}
    for d in (0.1, 0.5, 0.9):                   # class 0 sparse or dense, the others at random
        tries[f"class 0 of density {d}"] = [0 if r < d else int(rng.integers(1, K)) for r in rng.random(NV)]
    for c in range(K):                          # every count at its least or its greatest
        tries[f"every vertex of colour {c}"] = [c] * NV
    return tries


def check_predictions(name, f, tries):
    """for each colouring, the clauses of f that fail in its intended assignment are those that predicted() names;
    returns how many of them are lex-leader clauses, value precedence clauses and unit clauses of counts"""
    kinds = collections.Counter()
    for what, col in tries.items():
        model = intended(name, col)
        bad = f.violated(model.val)
        assert bad == predicted(name, model), what
        lo, hi = model.ranges["chains"]
        kinds[what, "chain"] = sum(any(lo <= abs(x) < hi for x in c) or (len(c) == 2 and all(abs(x) <= NV * K and
                                   (abs(x) - 1) % K == 0 for x in c) and min(c) < 0 < max(c)) for c in bad)
        lo, hi = model.ranges["precedence"]
        kinds[what, "precedence"] = sum(any(lo <= abs(x) < hi for x in c) or c in [frozenset([-X(0, b)])
                                        for b in range(2, K)] for c in bad)
        kinds[what, "unit"] = sum(len(c) == 1 and max(map(abs, c)) > NV * K for c in bad)
    return kinds


@pytest.mark.parametrize("name", ["F36", "F35", "F34"])
def test_the_failing_clauses_are_predicted_for_any_colouring(parsed, name):
    """for the colourings of colourings(), the clauses that fail are exactly those that predicted() names from the
    meaning of each kind of clause; and they are of the kinds that each change of a normal form should break"""
    kinds = check_predictions(name, parsed[name], colourings(name))
    assert kinds["lex-least class 0", "chain"] and kinds["colours 1 and 2 swapped", "precedence"]
    assert not kinds["lex-least class 0", "unit"] and not kinds["colours 1 and 2 swapped", "chain"]
    # the count of class 0 fails; in F34 also the count of the class that grows to 35
    assert kinds["class 0 one point short", "unit"] == (2 if name == "F34" else 1)
    assert not kinds["class 0 one point short", "chain"]
    assert all(kinds["random", k] for k in ("chain", "precedence", "unit"))


def is_chain_clause(f, lo, hi):
    """for each clause of f: a clause with a chain variable, or a first comparison x(v, 0) | -x(w, 0)"""
    out = np.zeros(len(f.ends), dtype=bool)
    chain = (np.abs(f.lits) >= lo) & (np.abs(f.lits) < hi)
    out[np.unique(np.searchsorted(f.ends, np.flatnonzero(chain)))] = True
    for k in np.flatnonzero(f.ends - f.starts == 2):
        a, b = f.lits[f.starts[k]:f.ends[k]]
        out[k] |= a * b < 0 and max(abs(a), abs(b)) <= NV * K and (abs(a) - 1) % K == 0 and (abs(b) - 1) % K == 0
    return out


@pytest.mark.parametrize("name", ["F36", "F35", "F34"])
def test_every_orbit_has_an_image_that_satisfies_the_lex_leader_clauses(parsed, name):
    """for each set, class 0 its lex-leader image satisfies every lex-leader clause, and its lex-least image does not;
    the sets are maximal independent, of the size of class 0, and in different orbits"""
    f, (sets, _, _) = parsed[name], SHAPES[name]
    size = int(name[1:])
    adjacency = [set(a) for a in ADJ]
    canon = set()
    for S in sets:
        assert len(set(S)) == size and not any(adjacency[v] & set(S) for v in S)
        assert all(v in S or adjacency[v] & set(S) for v in range(NV))
        canon.add(min(tuple(sorted(np.flatnonzero(row).tolist())) for row in images(S)[0]))
    assert len(canon) == len(sets)
    lo, hi = intended(name, [0] * NV).ranges["chains"]
    chains = is_chain_clause(f, lo, hi)
    assert chains.sum() == sum(3 * len([v for v in ORDER[:25] if p[v] != v]) - 2 for p in PERMS)
    for S in sets:
        for T, want in ((lex_leader(S), True), (lex_least(S), False)):
            model = intended(name, [0 if v in T else 1 for v in range(NV)])
            assert f.satisfied(model.val)[chains].all() == want


@pytest.fixture(scope="module")
def au():
    """scripts/g13/g13_audit.py, the audit of the formulas of alpha(G_13) <= 36"""
    return load("g13_audit", G13)


def check_counts(name, f, au):
    """every count of at least t of f is, for au.audit_totalizer, a standard totalizer over the 169 variables of its
    class with the unit clause o_t at its root; a count of at most t, one of at least 169 - t over the negations"""
    absolute = np.abs(f.lits)
    for start, end, lits, t, unit in intended(name, [0] * NV).counts:
        if not unit:
            continue                            # the counts of --domcap have clauses both ways: not a totalizer
        leaf = {x: v + 1 for v, x in enumerate(lits)}
        leaf.update({-x: -(v + 1) for v, x in enumerate(lits)})
        ks = np.unique(np.searchsorted(f.ends, np.flatnonzero((absolute >= start) & (absolute < end))))
        clauses = [tuple(leaf.get(x, x) for x in f.lits[f.starts[k]:f.ends[k]].tolist()) for k in ks]
        assert au.audit_totalizer(clauses, list(range(1, NV + 1)), t) == NV - 1, (name, t)


def chain_literals(name, f):
    """the literals of the lex-leader clauses of f, with x(v, 0) read as the vertex variable v + 1 and the chain
    variables numbered from 170"""
    lo, hi = intended(name, [0] * NV).ranges["chains"]
    lits = np.concatenate([f.lits[f.starts[k]:f.ends[k] + 1] for k in np.flatnonzero(is_chain_clause(f, lo, hi))])
    a = np.abs(lits)
    assert ((a - 1) % K == 0)[(a > 0) & (a <= NV * K)].all()
    return np.where(a > NV * K, a - lo + NV + 1, (a - 1) // K + 1) * np.sign(lits)


def check_chains(lits, au):
    """au.read_chains finds the chains; they compare the positions of ORDER[:25], each with its image under a map of
    the group that fixes the positions it skips (au.chain_maps), one chain for each map other than the identity"""
    clauses = [tuple(c[:-1]) for c in np.split(lits, np.flatnonzero(lits == 0) + 1)[:-1]]
    firsts = [(max(c) - 1, -min(c) - 1) for c in clauses if max(map(abs, c)) <= NV]
    read = au.read_chains([c for c in clauses if max(map(abs, c)) > NV], firsts, NV)
    prefix = au.common_order(read)
    maps = au.chain_maps(read, prefix)
    assert prefix == ORDER[:25] and len(read) == len(set(maps)) == 4731
    assert all(au.MAPS[i] != tuple(range(NV)) for i in maps)
    assert sum(map(len, read)) == sum(len([v for v in ORDER[:25] if p[v] != v]) for p in PERMS)


def test_the_counts_and_the_chains_pass_the_audit_of_alpha(parsed, au):
    """the counts and the lex-leader clauses of F36, F35 and F34, read from the text of each formula by the functions
    of scripts/g13/g13_audit.py, which audits the formulas of alpha(G_13) <= 36 (check_counts, check_chains); the
    chains are the same in the three formulas"""
    for name, f in parsed.items():
        check_counts(name, f, au)
    chains = {name: chain_literals(name, f) for name, f in parsed.items()}
    assert (chains["F36"] == chains["F35"]).all() and (chains["F35"] == chains["F34"]).all()
    check_chains(chains["F34"], au)


def check_place(name, f):
    """each clause of f is either on the colour variables x(v, c) alone, and then one of: at least one colour for a
    vertex, at most one, an edge clause, the domination clause of a class that must dominate (class 0; in F34 classes
    0 to 3), the first comparison of a chain, or value precedence at the first position; each kind in its number. Or
    its other variables all belong to one count, to the chains, or to value precedence."""
    edges = {frozenset(e) for e in EDGES}
    model = intended(name, [0] * NV)
    kinds = collections.Counter()
    vertex, colour = lambda x: (abs(x) - 1) // K, lambda x: (abs(x) - 1) % K
    dominating = (0, 1, 2, 3) if name == "F34" else (0,)
    pairs = (1, 2, 3) if name == "F34" else (1, 2, 3, 4)
    for k in np.flatnonzero(np.maximum.reduceat(np.abs(f.lits), f.starts) <= NV * K):
        c = f.lits[f.starts[k]:f.ends[k]].tolist()
        pos, neg = sorted(x for x in c if x > 0), sorted(-x for x in c if x < 0)
        if not neg and len(pos) == K and {vertex(x) for x in pos} == {vertex(pos[0])}:
            kinds["at least one colour"] += 1
        elif not pos and len(neg) == 2 and vertex(neg[0]) == vertex(neg[1]) and colour(neg[0]) != colour(neg[1]):
            kinds["at most one colour"] += 1
        elif not pos and len(neg) == 2 and colour(neg[0]) == colour(neg[1]) and frozenset(map(vertex, neg)) in edges:
            kinds["edge"] += 1
        elif (not neg and len(pos) == 15 and len({colour(x) for x in pos}) == 1 and colour(pos[0]) in dominating
              and any(sorted(map(vertex, pos)) == sorted([v] + ADJ[v]) for v in map(vertex, pos))):
            kinds["domination"] += 1
        elif len(pos) == len(neg) == 1 and colour(pos[0]) == colour(neg[0]) == 0:
            kinds["first comparison"] += 1
        elif not pos and len(neg) == 1 and neg[0] in [X(ORDER[0], b) for b in pairs[1:]]:
            kinds["value precedence at the first position"] += 1
        else:
            raise AssertionError(f"{name}: a clause of no known kind: {c}")
    assert kinds == {"at least one colour": NV, "at most one colour": NV * 10, "edge": len(EDGES) * K,
                     "domination": NV * len(dominating), "first comparison": 4731,
                     "value precedence at the first position": len(pairs) - 1}, (name, kinds)
    runs = sorted([(start, end) for start, end, *_ in model.counts] + [model.ranges["chains"],
                                                                       model.ranges["precedence"]])
    assert runs[0][0] == NV * K + 1 and runs[-1][1] == f.nv + 1
    assert all(a[1] == b[0] for a, b in zip(runs, runs[1:]))                         # every auxiliary variable
    run = np.searchsorted([start for start, _ in runs], np.abs(f.lits), side="right") - 1
    run[np.abs(f.lits) <= NV * K] = -1
    high = np.maximum.reduceat(run, f.starts)
    low = np.minimum.reduceat(np.where(run < 0, len(runs), run), f.starts)
    assert ((high == -1) | (low == high)).all()                                      # one run per clause


def test_every_clause_has_its_place(parsed):
    for name, f in parsed.items():
        check_place(name, f)


def replaced(f, k, clause):
    """a copy of f with its clause k replaced (None: removed), or with the clause appended (k = None)"""
    g = Formula.__new__(Formula)
    new = np.array(list(clause) + [0], dtype=np.int64) if clause is not None else np.zeros(0, dtype=np.int64)
    g.lits = (np.concatenate([f.lits, new]) if k is None
              else np.concatenate([f.lits[:f.starts[k]], new, f.lits[f.ends[k] + 1:]]))
    g.nv, g.ends = f.nv, np.flatnonzero(g.lits == 0)
    g.starts = np.concatenate(([0], g.ends[:-1] + 1))
    return g


def first(f, test):
    return next(k for k in range(len(f.ends)) if test(f.lits[f.starts[k]:f.ends[k]].tolist()))


def wrong_edge(f, m):
    k = first(f, lambda c: len(c) == 2 and c[0] < 0 and c[1] < 0 and (-c[0] - 1) % K == (-c[1] - 1) % K)
    a, colour = divmod(-f.lits[f.starts[k]] - 1, K)
    w = next(v for v in range(NV) if v != a and v not in ADJ[a])
    return replaced(f, k, [-X(a, colour), -X(w, colour)])


def counter_clause_negated(f, m):
    start, end = m.counts[0][:2]                            # the count of class 0
    k = first(f, lambda c: len(c) == 3 and any(start <= abs(x) < end for x in c))
    c = f.lits[f.starts[k]:f.ends[k]].tolist()
    return replaced(f, k, [-c[0]] + c[1:])


def root_one_lower(f, m):
    k = first(f, lambda c: c == [m.size0[35]])
    return replaced(f, k, [m.size0[34]])


def chain_clause_negated(f, m):
    lo, hi = m.ranges["chains"]
    k = first(f, lambda c: len(c) == 3 and sum(lo <= abs(x) < hi for x in c) == 1 and min(c) < 0 < max(c)
              and max(abs(x) for x in c if not lo <= abs(x) < hi) <= NV * K and
              sum(x > 0 for x in c if abs(x) <= NV * K) == 1)
    c = f.lits[f.starts[k]:f.ends[k]].tolist()
    return replaced(f, k, [-x if 0 < x <= NV * K else x for x in c])


def precedence_clause_negated(f, m):
    lo, hi = m.ranges["precedence"]
    k = first(f, lambda c: len(c) == 2 and any(lo <= x < hi for x in c) and any(-NV * K <= x < 0 for x in c))
    return replaced(f, k, [abs(x) for x in f.lits[f.starts[k]:f.ends[k]].tolist()])


def guard_negated(f, m):
    k = first(f, lambda c: len(c) == 16 and -m.domcap[0] in c)
    return replaced(f, k, [-x if x == -m.domcap[0] else x for x in f.lits[f.starts[k]:f.ends[k]].tolist()])


def stray_unit(f, m):
    return replaced(f, None, [X(0, 1)])


WRONG = [(wrong_edge, "place"), (stray_unit, "place"), (counter_clause_negated, "counts"),
         (root_one_lower, "counts"), (chain_clause_negated, "chains"), (precedence_clause_negated, "models"),
         (guard_negated, "models")]


@pytest.mark.parametrize("change, check", WRONG, ids=[w.__name__ for w, _ in WRONG])
def test_one_wrong_clause_is_caught(parsed, au, change, check):
    """F36 with one wrong clause of each kind fails one of the checks above: an edge clause on a pair that is not an
    edge, a stray unit clause, a clause of a count with a literal negated, the unit clause of a count one lower, a
    lex-leader clause with a literal negated, a clause of value precedence with a literal negated, and a domination
    clause of --domcap with its guard negated"""
    f, m = parsed["F36"], intended("F36", [0] * NV)
    wrong = change(f, m)
    assert len(wrong.ends) in (len(f.ends), len(f.ends) + 1) and not np.array_equal(wrong.lits, f.lits)
    runs = {"place": lambda g: check_place("F36", g), "counts": lambda g: check_counts("F36", g, au),
            "chains": lambda g: check_chains(chain_literals("F36", g), au),
            "models": lambda g: check_predictions("F36", g, {c: [c] * NV for c in range(K)})}
    if check != "chains":                           # the check passes on F36 itself (the chains: the test above)
        runs[check](f)
    with pytest.raises((AssertionError, au.AuditError)):
        runs[check](wrong)


# ------------------------------------------------------------------------------------------------ the cube files

def read_cubes(path):
    out, closed = [], 0
    for line in open(path):
        if line.startswith("a "):
            out.append(tuple(map(int, line.split()[1:-1])))
        elif line.startswith("c closed "):
            out.append(tuple(map(int, line.split()[2:-1])))
            closed += 1
    return out, closed


@pytest.mark.parametrize("name, leaves, closed, variables", [("E37_B", 4822, 564, "vars_s_lex.txt"),
                                                              ("F36", 4823, 565, "vars_c0_lex.txt"),
                                                              ("F35", 4823, 564, "vars_c0_lex.txt"),
                                                              ("F34", 4823, 564, "vars_c0_lex.txt")])
def test_the_cube_files_are_covers(name, leaves, closed, variables):
    """each cube file of a case is the set of leaves of a binary tree (check_cover), with at most 14 decisions, on the
    variables of its list (x(v, 0) or s_v along lex_order()) in the order of the list"""
    cuber2 = load("cuber2")
    cubes, n = read_cubes(os.path.join(CHI, "cases", name + ".icnf"))
    assert (len(cubes), n) == (leaves, closed) and cuber2.check_cover(cubes)
    rank = {int(x): k for k, x in enumerate(open(os.path.join(CHI, "work", variables)).read().split(","))}
    want = [X(v, 0) for v in ORDER] if name != "E37_B" else [1 + v for v in ORDER]
    assert sorted(rank, key=rank.get) == want
    assert all(len(c) <= 14 and all(rank[abs(a)] < rank[abs(b)] for a, b in zip(c, c[1:])) for c in cubes)
    assert not cuber2.check_cover(cubes[:-1]) and not cuber2.check_cover(cubes + cubes[-1:])
    if name == "F35":           # F34 and F35 propagate alike on these decisions: the same tree
        assert open(os.path.join(CHI, "cases", "F35.icnf")).read() == open(os.path.join(CHI, "cases", "F34.icnf")).read()


def test_check_colouring(tmp_path):
    """check_colouring.py, which does not import g13.py: the stored 6-colouring is proper; a 5-colouring in the form
    of a kissat model is not"""
    col = json.load(open(os.path.join(ROOT, "data", "small_plane_colourings.json")))["G_q"]["13"]
    (tmp_path / "six").write_text(" ".join(map(str, col)))
    r = subprocess.run([sys.executable, os.path.join(CHI, "check_colouring.py"), str(tmp_path / "six")],
                       capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == "6 colours, 1183 edges checked, 0 monochromatic"
    model = [X(v, c) if c == v % 5 else -X(v, c) for v in range(NV) for c in range(K)]
    (tmp_path / "five").write_text("s SATISFIABLE\nv " + " ".join(map(str, model)) + " 0\n")
    r = subprocess.run([sys.executable, os.path.join(CHI, "check_colouring.py"), str(tmp_path / "five")],
                       capture_output=True, text=True)
    assert r.returncode == 1 and "monochromatic edge" in r.stdout


# ------------------------------------------------------------------------------ a tiny plan, and the checkers

STUB = '''"""a stand-in for NAME.py, for tests/test_g13_chi.py: a small formula that depends on the arguments"""
import hashlib
import sys

out, args = sys.argv[1], sys.argv[2:]
b = hashlib.sha256(" ".join(["NAME"] + args).encode()).digest()
with open(out, "w") as fh:
    fh.write("p cnf 9 3\\n" + "".join(f"{b[k] % 9 + 1} -{b[k + 1] % 9 + 1} 0\\n" for k in (0, 2, 4)))
'''

# each cube file: its leaves, and those that have only a timeout line in the log. F35 is certified by a re-split
# CASE_splitI; F34 by the re-split F34_split1, after the re-split F34_leaf1, which is tried first, fails
SYNTHETIC = {"E37_B": (["a 1 0", "a -1 0"], ()),
             "F36": (["c closed 1 0", "a -1 2 0", "a -1 -2 0"], ()),
             "F35": (["a 3 0", "a -3 0"], (1,)),
             "F35_split1": (["a 4 0", "a -4 0"], ()),
             "F34": (["a 5 0", "a -5 0"], (1,)),
             "F34_leaf1": (["a 6 0", "a -6 0"], (1,)),
             "F34_split1": (["a 7 0", "a -7 0"], (0,)),
             "F34_split1_split0": (["a 8 0", "a -8 0"], ())}


def leaf_text(text, cube_line):
    """the formula of a leaf, as certify.py writes it: the formula with the literals of the leaf as unit clauses"""
    head, body = text.split("\n", 1)
    nv, nc = map(int, head.split()[2:4])
    lits = cube_line.split()[1:-1] if cube_line.startswith("a ") else cube_line.split()[2:-1]
    return f"p cnf {nv} {nc + len(lits)}\n" + body + "".join(f"{x} 0\n" for x in lits)


def log_line(name, text, verified=True):
    if verified:
        return f"{name}: sha256 {sha256(text)}; kissat UNSAT in 0.5 s (proof 0.1 MB); drat-trim VERIFIED, 3 of 7 " \
               f"lemmas in core in 0.2 s\n"
    return f"{name}: sha256 {sha256(text)}; kissat UNKNOWN in 1200.0 s (proof 90.0 MB); drat-trim not run in 0.0 s\n"


def synthetic_plan(root):
    """a tiny plan D: root/code holds the real verify_plan_D.py, share_driver.py, cuber2.py and pack_certificates.py,
    stand-ins for enum_cert.py and g13cnf.py, the base cube files and cases/SHA256SUMS; root/assembly is a copy with
    the logs, the cube files of the re-splits and logs/E37_A.log of a certified run"""
    code = root / "code"
    (code / "cases").mkdir(parents=True)
    for name in ("verify_plan_D.py", "share_driver.py", "cuber2.py", "pack_certificates.py"):
        shutil.copy(os.path.join(CHI, name), code / name)
    for name in ("enum_cert", "g13cnf"):
        (code / f"{name}.py").write_text(STUB.replace("NAME", name))
    texts = {}
    for name, command in load("share_driver").GEN.items():
        args = command.format(str(root / "f.cnf")).split()
        subprocess.run([sys.executable, str(code / args[0])] + args[1:], check=True)
        texts[name] = (root / "f.cnf").read_text()
    assert len(set(texts.values())) == 8
    (code / "cases" / "SHA256SUMS").write_text("".join(f"{sha256(t)}  {n}.cnf\n" for n, t in texts.items()))
    for stem in ("E37_B", "F36", "F35", "F34"):
        (code / "cases" / f"{stem}.icnf").write_text("".join(x + "\n" for x in SYNTHETIC[stem][0]))
    assembly = root / "assembly"
    shutil.copytree(code, assembly)
    (assembly / "logs").mkdir()
    (assembly / "logs" / "E37_A.log").write_text("".join(log_line(f"E37_A{c}.cnf", texts[f"E37_A{c}"])
                                                         for c in (11, 6, 9, 7)))
    for stem in sorted(SYNTHETIC, key=len):
        m = re.fullmatch(r"(.*)_(?:leaf|split)(\d+)", stem)
        if m:
            texts[stem] = leaf_text(texts[m.group(1)], SYNTHETIC[m.group(1)][0][int(m.group(2))])
            (assembly / "cases" / f"{stem}.icnf").write_text("".join(x + "\n" for x in SYNTHETIC[stem][0]))
        leaves, timeouts = SYNTHETIC[stem]
        (assembly / "cases" / f"{stem}.certlog").write_text("".join(
            log_line(f"{stem}_leaf{i}", leaf_text(texts[stem], x), i not in timeouts) for i, x in enumerate(leaves)))
    return code, assembly


@pytest.fixture(scope="module")
def plan(tmp_path_factory):
    """the tiny plan, with `share_driver.py regen` run in the assembly"""
    root = tmp_path_factory.mktemp("g13_chi_plan")
    code, assembly = synthetic_plan(root)
    r = subprocess.run([sys.executable, "share_driver.py", "regen"], cwd=assembly, capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == "4 re-split formulas written", r.stdout + r.stderr
    return code, assembly


def verify_plan_d(folder):
    return subprocess.run([sys.executable, "verify_plan_D.py"], cwd=folder, capture_output=True, text=True)


def test_verify_plan_d_accepts_a_certified_plan(plan):
    r = verify_plan_d(plan[1])
    assert r.returncode == 0, r.stdout + r.stderr
    assert r.stdout.splitlines() == [f"{case} {n} leaves VERIFIED (re-split sub-leaves included)"
                                     for case, n in (("E37_B", 2), ("F36", 3), ("F35", 3), ("F34", 4))] + [
                                         "PLAN D FULLY CERTIFIED"]


def edit(path, old, new):
    text = path.read_text()
    assert old in text
    path.write_text(text.replace(old, new, 1))


def drop_line(path, start):
    lines = path.read_text().splitlines(keepends=True)
    assert sum(x.startswith(start) for x in lines) == 1
    path.write_text("".join(x for x in lines if not x.startswith(start)))


def stale(d):
    line = next(x for x in (d / "cases" / "F36.certlog").read_text().splitlines() if x.startswith("F36_leaf2:"))
    edit(d / "cases" / "F36.certlog", line, line.replace(line[18:26], "00000000"))


BROKEN = {
    "a leaf deep down without a line": (
        lambda d: drop_line(d / "cases" / "F34_split1_split0.certlog", "F34_split1_split0_leaf1:"),
        "FAIL F34: F34_leaf1: not VERIFIED and not re-split into a certified case"),
    "a line for another formula": (stale, "FAIL F36: F36_leaf2: not VERIFIED"),
    "drat-trim NOT VERIFIED": (
        lambda d: edit(d / "cases" / "F36.certlog", "drat-trim VERIFIED", "drat-trim NOT VERIFIED"),
        "FAIL F36: F36_leaf0: not VERIFIED"),
    "kissat SAT": (lambda d: edit(d / "cases" / "E37_B.certlog", "kissat UNSAT in 0.5 s (proof 0.1 MB); drat-trim "
                                  "VERIFIED, 3 of 7 lemmas in core", "kissat SAT in 0.5 s (proof 0.1 MB); drat-trim not run"),
                   "FAIL E37_B: E37_B_leaf0: not VERIFIED"),
    "a base cube file that is not a cover": (lambda d: drop_line(d / "cases" / "F36.icnf", "a -1 -2 0"),
                                             "FAIL F36: cases/F36.icnf is not a cover"),
    "a re-split that is not a cover": (lambda d: drop_line(d / "cases" / "F35_split1.icnf", "a -4 0"),
                                       "FAIL F35: F35_leaf1: not VERIFIED"),
    "a re-split formula that is not the leaf formula": (
        lambda d: edit(d / "cases" / "F35_split1.cnf", "p cnf 9 4\n", "p cnf 9 5\n1 0\n"),
        "FAIL F35: F35_leaf1: not VERIFIED"),
    "a formula E37_A without a line": (lambda d: drop_line(d / "logs" / "E37_A.log", "E37_A9.cnf:"),
                                       "FAIL E37_A9 not VERIFIED"),
    "a case formula that the code does not write": (lambda d: edit(d / "cases" / "F34.cnf", "p cnf 9 3\n", "p cnf 9 4\n1 0\n"),
                                                    "FAIL cases/F34.cnf differs from the code's output")}


@pytest.mark.parametrize("breakage", list(BROKEN), ids=list(BROKEN))
def test_verify_plan_d_rejects_a_broken_plan(plan, tmp_path, breakage):
    change, message = BROKEN[breakage]
    copy = tmp_path / "copy"
    shutil.copytree(plan[1], copy)
    change(copy)
    r = verify_plan_d(copy)
    assert r.returncode != 0 and "PLAN D FULLY CERTIFIED" not in r.stdout
    assert (r.stdout + r.stderr).strip().splitlines()[-1] == message or message in r.stderr, r.stdout + r.stderr


def pack(code, assembly, out):
    return subprocess.run([sys.executable, str(code / "pack_certificates.py"), str(assembly), "--out-dir", str(out)],
                          capture_output=True, text=True)


def verify_chi(archive, code, workdir, *extra):
    return subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "verify_g13_chi.py"), str(archive), "--chi",
                           str(code), "--workdir", str(workdir), "--force"] + [str(x) for x in extra],
                          capture_output=True, text=True)


def test_pack_and_check_the_certificates(plan, tmp_path):
    code, assembly = plan
    for out in ("one", "two"):
        r = pack(code, assembly, tmp_path / out)
        assert r.returncode == 0, r.stdout + r.stderr
    archive = tmp_path / "one" / "g13_chi_certlogs.tar.gz"
    assert archive.read_bytes() == (tmp_path / "two" / "g13_chi_certlogs.tar.gz").read_bytes()     # deterministic
    sums = (tmp_path / "one" / "g13_chi_SHA256SUMS.txt").read_text().splitlines()
    assert sums[0] == f"{sha256(archive.read_bytes())}  g13_chi_certlogs.tar.gz"
    names = sorted([f"cases/{s}.certlog" for s in SYNTHETIC] + ["logs/E37_A.log"]
                   + [f"cases/{s}.icnf" for s in SYNTHETIC if "_leaf" in s or "_split" in s])
    assert [x.split("  ")[1] for x in sums[1:]] == names
    with tarfile.open(archive) as tar:
        infos = tar.getmembers()
        assert [i.name for i in infos] == names
        assert all((i.mode, i.mtime, i.uid, i.gid, i.uname, i.gname, i.isfile()) == (0o644, 1790467200, 0, 0, "", "", True)
                   for i in infos)
        assert all(sha256(tar.extractfile(i).read()) == x.split("  ")[0] for i, x in zip(infos, sums[1:]))

    note = tmp_path / "note.md"
    note.write_text("{{F34_CERTIFIED}}, {{TOTAL_CERTIFIED}}, {{F34_DEPTH}}, {{DATE}}\n{{TREE_TABLE}}\n")
    r = verify_chi(archive, code, tmp_path, "--stats", "--fill", note)
    assert r.returncode == 0, r.stdout + r.stderr
    out = r.stdout
    assert "PLAN D FULLY CERTIFIED" in out and out.strip().splitlines()[-1].startswith("CONFIRMED: ")
    for line in ["F34: 2 leaves (0 closed by unit propagation); 1 VERIFIED, 1 split again (1 as CASE_splitI, 0 as "
                 "CASE_leafI); certified in all: 4; the deepest re-split is at level 2",
                 "level 2: 1 re-split formulas, 2 leaves (0 closed); 2 VERIFIED, 0 split again",
                 "all cases: 16 leaves and formulas certified, each by one VERIFIED line (4 formulas E37_A_c and 12 "
                 "leaves); 3 re-split formulas",
                 "cube files of re-splits: 4, of which 1 are not used by the certified trees",
                 "{{TIMEOUTS_1200}} 4", "{{UNUSED_RESPLITS}} 1", "{{KISSAT_HOURS}} 0.0", "{{ARCHIVE_FILES}} 13"]:
        assert line in out, line
    assert note.read_text() == "4, 16, 2, {{DATE}}\n| `F35` | 1 | 1 | 2 | 2 | 0 |\n| `F34` | 1 | 1 | 2 | 1 | 1 |\n" \
                               "| `F34` | 2 | 1 | 2 | 2 | 0 |\n"
    assert "placeholders written; left: DATE" in out

    changed = tmp_path / "changed"
    changed.mkdir()
    data = bytearray(archive.read_bytes())
    data[len(data) // 2] ^= 1
    (changed / "g13_chi_certlogs.tar.gz").write_bytes(bytes(data))
    shutil.copy(tmp_path / "one" / "g13_chi_SHA256SUMS.txt", changed)
    r = verify_chi(changed / "g13_chi_certlogs.tar.gz", code, tmp_path)
    assert r.returncode == 1 and "not the" in r.stdout and r.stdout.strip().splitlines()[-1].startswith("NOT CONFIRMED")

    broken = tmp_path / "broken"
    shutil.copytree(assembly, broken)
    drop_line(broken / "cases" / "F35_split1.certlog", "F35_split1_leaf0:")
    assert pack(code, broken, tmp_path / "three").returncode == 0
    r = verify_chi(tmp_path / "three" / "g13_chi_certlogs.tar.gz", code, tmp_path, "--stats")
    assert r.returncode == 1 and "FAIL F35: F35_leaf1: not VERIFIED" in r.stdout
    assert "F35: NOT fully certified" in r.stdout and r.stdout.strip().splitlines()[-1].startswith("NOT CONFIRMED")


# ------------------------------------------------------------------------------------- the certificates themselves

ARCHIVE = os.path.join(ROOT, "certificates", "g13_chi_certlogs.tar.gz")
ARCHIVE_SUMS = os.path.join(ROOT, "certificates", "g13_chi_SHA256SUMS.txt")
CERTIFY_LINE = re.compile(r"(\S+): sha256 [0-9a-f]{64}; kissat (UNSAT|SAT|UNKNOWN) in [\d.]+ s \(proof [\d.]+ MB\); "
                          r"drat-trim (VERIFIED|NOT VERIFIED|not run)(?:, \d+ of \d+ lemmas in core)? in [\d.]+ s$")


def check_archive(archive, sums_path):
    """the archive against its SHA256SUMS file: the SHA-256 of the archive and of each file; only logs of certify.py
    and cube files of re-splits, packed as pack_certificates.py packs them; no satisfiable leaf; and the logs of the
    formulas E37 are the certificates of alpha(G_13) = 36"""
    sums = [line.split("  ") for line in open(sums_path).read().splitlines()]
    assert sums[0] == [sha256(open(archive, "rb").read()), os.path.basename(archive)]
    member = re.compile(r"cases/(?:E37_B|F36|F35|F34)(?:_(?:leaf|split)\d+)*\.certlog|"
                        r"cases/(?:E37_B|F36|F35|F34)(?:_(?:leaf|split)\d+)+\.icnf|logs/E37_A\.log")
    texts = {}
    with tarfile.open(archive) as tar:
        infos = tar.getmembers()
        assert [i.name for i in infos] == [name for _, name in sums[1:]] == sorted(i.name for i in infos)
        for info, (digest, name) in zip(infos, sums[1:]):
            assert member.fullmatch(name) and info.isfile() and (info.mode, info.uid, info.gid) == (0o644, 0, 0)
            data = tar.extractfile(info).read()
            assert sha256(data) == digest, name
            if not name.endswith(".icnf"):
                texts[name] = data.decode().splitlines()
    for name, lines in texts.items():
        for line in lines:
            m = CERTIFY_LINE.match(line)
            assert m and m.group(2) != "SAT", (name, line)
    alpha = {}
    for part, opener in (("a", open), ("b", gzip.open)):
        path = os.path.join(ROOT, "certificates", f"g13_alpha_part_{part}_checks.txt" + (".gz" if part == "b" else ""))
        with opener(path, "rt") as fh:
            alpha[part] = sorted(line.rstrip("\n") for line in fh if not line.startswith("#"))
    assert sorted(texts["logs/E37_A.log"]) == alpha["a"] and sorted(texts["cases/E37_B.certlog"]) == alpha["b"]
    return len(infos)


@pytest.mark.skipif(not os.path.exists(ARCHIVE), reason="certificates/g13_chi_certlogs.tar.gz is not there yet")
def test_the_certificate_archive():
    assert check_archive(ARCHIVE, ARCHIVE_SUMS) > 4
