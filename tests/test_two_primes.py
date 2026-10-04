"""Theorems E and F (notes/circular_planes.md §6, papers/three-colours Section 10): the computer-assisted steps and
the finite facts at 7, rerun from the stored programs and certificates (data/number_fields/circular/twoprime/ and the
referees' subfolders indep_E/, indep_F/, indep_S10/, indep_FW/ and indep_FW2/)."""
import gzip, hashlib, json, lzma, os, re, shutil, subprocess, sys

import pytest

P = os.path.join(os.path.dirname(__file__), "..", "data", "number_fields", "circular", "twoprime")
FW = os.path.join(os.path.dirname(__file__), "..", "data", "number_fields", "circular", "finite_witness")
DRAT = shutil.which("drat-trim") or os.environ.get("DRAT_TRIM")


def run(cwd, *args, timeout=900, check=True):
    return subprocess.run([sys.executable, "-B", *args], cwd=cwd, capture_output=True, text=True, check=check,
                          timeout=timeout)


def stage(tmp_path, names, trees=()):
    """Copy programs and certificates into tmp_path and decompress the stored trees there."""
    for n in names:
        shutil.copy(os.path.join(P, n), tmp_path / n)
    for t in trees:
        with gzip.open(os.path.join(P, t + ".gz"), "rb") as src, open(tmp_path / t, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return tmp_path


def window65(tmp_path):
    return stage(tmp_path, ["verify_tree.py", "verify_prop7.py", "prop7_certificates.txt"], ["tree_r7_25.txt"])


def test_window_65_certificates(tmp_path):
    """The window G(1,1) modulo 65 (Proposition E1 of the note): completeness tree and the 13 index vectors."""
    d = window65(tmp_path)
    out = run(d, "verify_tree.py").stdout
    assert "tree verified: 4225 cells, 15620 branch nodes, 4364 Farkas-closed branches, 13 leaves" in out
    out = run(d, "verify_prop7.py").stdout
    assert "5 main components" in out and "8 extra components with kappa = 2/7 exactly" in out


def corrupt(path, pattern, repl):
    s = path.read_text()
    t = re.sub(pattern, repl, s, count=1, flags=re.M)
    assert t != s
    path.write_text(t)


@pytest.mark.parametrize("what", ["multiplier", "dropped branch", "certificate"])
def test_window_65_rejects_corrupted_copies(tmp_path, what):
    d = window65(tmp_path)
    if what == "multiplier":
        corrupt(d / "tree_r7_25.txt", r"^(D -?\d+ \S+ \S+ )(\d+)/(\d+)", lambda m: f"{m[1]}{int(m[2]) + 1}/{m[3]}")
        prog = "verify_tree.py"
    elif what == "dropped branch":
        corrupt(d / "tree_r7_25.txt", r"^D [^\n]*\n", "")
        prog = "verify_tree.py"
    else:
        corrupt(d / "prop7_certificates.txt", r"( ; cert -?\d+ -?\d+ \S+ \S+ )(\d+)/(\d+)",
                lambda m: f"{m[1]}{int(m[2]) + 1}/{m[3]}")
        prog = "verify_prop7.py"
    r = run(d, prog, check=False)
    assert r.returncode != 0 and "verified" not in r.stdout


@pytest.mark.slow
def test_window_325_certificates(tmp_path):
    """The window G(2,1) modulo 325 (Proposition F1 of the note): 29 index vectors, 16 of them excluded at 1/4
    (about 45 seconds)."""
    d = stage(tmp_path, ["verify_window.py", "cert_K2M1.txt"], ["tree_K2M1.txt"])
    out = run(d, "verify_window.py", "2", "1").stdout
    assert "tree verified (105625 cells, 550474 branch nodes, 199920 Farkas-closed branches, 29 leaves)" in out
    assert "labels {'7': 8, 'extra': 16, 'Q': 4, 'C': 1}" in out and "EXTRA components = 1/4" in out


@pytest.mark.parametrize("prog", ["seven_local", "seven_patterns"])
def test_facts_at_seven(prog):
    """A'_7: eight elements, one mu_8-orbit, no coset of a subgroup; the digit patterns of the 7-adic characters."""
    with open(os.path.join(P, prog + ".txt")) as fh:
        assert run(P, prog + ".py").stdout == fh.read()


def test_facts_at_seven_referee():
    out = run(os.path.join(P, "indep_F"), "f49_and_digits.py").stdout
    with open(os.path.join(P, "indep_F", "f49_and_digits.txt")) as fh:
        assert out == fh.read()
    for line in ["(F1) single mu8-orbit: True ; 0 not in A'_7: True", "(F2) cosets of non-zero subgroups inside "
                 "A'_7: none", "the 8 shifts are distinct and exhaust Z/8: True", "2 consecutive digits suffice: True",
                 "7-type points consistent with the characters: True"]:
        assert line in out


@pytest.mark.parametrize("prog", ["b4_model_check", "b5_f3_check"])
def test_lemmas_at_seven_on_models(prog):
    """Lemma F6 on the finite models F_49[pi]/(pi^e), e <= 3, and Lemma F7 for f = 3."""
    d = os.path.join(P, "indep_F")
    with open(os.path.join(d, prog + ".txt")) as fh:
        assert run(d, prog + ".py").stdout == fh.read()


def test_window_65_referee_vertices():
    """The referee's vertex enumeration: 13 components at 7/25 and 2/7, 5 above 2/7."""
    d = os.path.join(P, "indep_E")
    counts = {}
    for a, b in [(7, 25), (2, 7), (201, 700), (3, 10), (1, 3)]:
        head = run(d, "ref_vertex.py", str(a), str(b)).stdout.splitlines()[0]
        counts[(a, b)] = int(head.rsplit("components:", 1)[1])
    assert counts == {(7, 25): 13, (2, 7): 13, (201, 700): 5, (3, 10): 5, (1, 3): 5}


def test_window_65_referee_clip_and_kappa(tmp_path):
    """The referee's clipping at 7/25 and the exact largest least margin of each of the 13 components."""
    for n in ["ref_clip.py", "ref_kappa.py", "ref_rot.py"]:
        shutil.copy(os.path.join(P, "indep_E", n), tmp_path / n)
    run(tmp_path, "ref_clip.py", "7", "25", "clip.txt")
    out = run(tmp_path, "ref_kappa.py", "clip.txt").stdout
    vals = re.findall(r"kappa primal = (\S+), dual = (\S+), equal: True", out)
    assert sorted(p for p, q in vals if p == q) == sorted(["2/7"] * 8 + ["1/3"] * 4 + ["1/2"])


def test_torsion_characters():
    """One prime: orders 41 and 76 beat 2/7; two primes: only types c and q (orders <= 400)."""
    d = os.path.join(P, "indep_E")
    one = run(d, "ref_torsion.py", "one", "400").stdout
    two = run(d, "ref_torsion.py", "two", "400").stdout
    assert "exact orders with an orbit of kappa > 2/7: [2, 3, 41, 76]" in one
    assert "M=41: max kappa 12/41" in one and "M=76: max kappa 11/38" in one
    assert "exact orders with an orbit of kappa > 2/7: [2, 3]" in two


def test_padic_points_killed_by_sigma():
    out = run(os.path.join(P, "indep_E"), "ref_ck.py").stdout
    assert out.count("min margin of c_k over G_N = 3/10") == 5
    assert "k=1: min margin of c_k over G_N = 3/10 ; max over m of min margin of c_k+Nm over G(k,1) = 41/250" in out


def test_window_325_seven_certificates(tmp_path):
    """The 8 index vectors of family 7 in the window G(2,1) have largest least margin exactly 2/7, so for theta > 2/7
    that window leaves only the types c and q (Theorem E again, and Theorem F from Proposition F1 alone)."""
    d = stage(tmp_path, ["verify_seven_K2M1.py", "cert_K2M1.txt", "cert_K2M1_seven.txt"])
    assert "largest least margin exactly 2/7" in run(d, "verify_seven_K2M1.py").stdout
    corrupt(d / "cert_K2M1_seven.txt", r"(\| kappa = 2/7 ; cert -?\d+ -?\d+ \S+ \S+ )(\d+)/(\d+)",
            lambda m: f"{m[1]}{int(m[2]) + 1}/{m[3]}")
    r = run(d, "verify_seven_K2M1.py", check=False)
    assert r.returncode != 0 and "verified" not in r.stdout


def test_section10_referee(tmp_path):
    """The referee of Section 10 of the paper: the facts at 7 and the type points modulo 325, the points c_k, torsion,
    an archimedean point of Q(i), the family-7 vectors of G(2,1) restricted to G(1,1), and its own clipping of the
    window G(1,1) at 2/7 (the same 13 index vectors as the certificates)."""
    d = os.path.join(P, "indep_S10")
    for n in os.listdir(d):
        shutil.copy(os.path.join(d, n), tmp_path / n)
    (tmp_path / "run").mkdir()
    for n in ["cert_K2M1.txt", "prop7_certificates.txt"]:
        shutil.copy(os.path.join(P, n), tmp_path / "run" / n)
    for prog in ["check_seven", "check_ck", "check_qi_char", "check_torsion"]:
        assert run(tmp_path, prog + ".py").stdout == (tmp_path / (prog + ".txt")).read_text()
    out = run(tmp_path, "check_remark1.py").stdout.splitlines()
    assert {x: out.count(x) for x in out} == {"MAIN -> MAIN": 5, "SEVEN -> EXTRA": 8,
                                               "EXTRA -> not among the 13 (1,1)-vectors": 16}
    out = run(tmp_path, "clip_window.py", "1", "1", "2/7", "run/prop7_certificates.txt").stdout
    assert "(K,M)=(1,1) N=65 r0=2/7: 13 index vectors with nonempty polygons" in out
    assert "same set of index vectors as the certificate: True" in out


def test_tight_walks(tmp_path):
    """Corollary on finite witnesses: the eight certificates of value 2/7 give positive integer relations among tight
    rotations (closed walks of length 140 or 42), type q has relations 3, 4, 5 in G_5, and the Moser spindle has
    chi = 4 and maps to K_{7/2}; the program reproduces its stored output."""
    d = stage(tmp_path, ["tight_walks.py", "cert_K2M1_seven.txt"])
    out = run(d, "tight_walks.py").stdout
    assert out == open(os.path.join(P, "tight_walks.txt")).read()
    assert out.count("multiplicities [63, 52, 25], length 140") == 4
    assert out.count("multiplicities [13, 14, 15], length 42") == 4
    assert out.count("multiplicities (3, 4, 5)") == 8
    assert "proper 3-colourings: 0; a homomorphism to K_7/2: (0, 2, 4, 6, 3, 5, 1)" in out


def test_finite_ball(tmp_path):
    """Remark after the corollary on finite witnesses: the ball of radius 2 in the Cayley graph of the 140 vectors of
    Theorem C is bipartite (19 600 edges), and the unit-distance graph it induces (216 more edges) has a 5-cycle, no
    triangle and a homomorphism to K_{5/2}; the program reproduces its stored output."""
    d = stage(tmp_path, ["finite_ball.py"])
    out = run(d, "finite_ball.py").stdout
    assert out == open(os.path.join(P, "finite_ball.txt")).read()
    assert "9941 points; unit-distance pairs 19816; Cayley edges 19600; others 216" in out
    assert "Cayley ball bipartite: True" in out and "no triangle" in out
    assert "circular chromatic number 5/2" in out


def test_lemma13_referee():
    """The referee of the hand proof of Lemma 13 (twoprime/indep_L13/): Lemma 13 exhaustively, the new number-theoretic
    sentences of the paper, the radius-2 ball and the 628-vertex union; each program reproduces its stored output."""
    d = os.path.join(P, "indep_L13")
    for name, passes in [("check_lemma13", 68), ("check_diff_claims", 26), ("check_ball", 13), ("check_628", 4)]:
        out = run(d, name + ".py").stdout
        assert out == open(os.path.join(d, name + ".out")).read()
        assert out.count("PASS") == passes and out.rstrip().endswith("FAILED: none")


def test_finite_witness_referee(tmp_path):
    """The referee of the corollary on finite witnesses: its checks of the certificates, of type q, of the tight
    relations of the 7-adic colourings, of the Moser spindle and of the tight-cycle and winding lemmas reproduce."""
    src = os.path.join(P, "indep_FW")
    d = tmp_path / "indep_FW"
    shutil.copytree(src, d)
    for n in ["cert_K2M1.txt", "cert_K2M1_seven.txt"]:
        shutil.copy(os.path.join(P, n), tmp_path / n)
    for prog in ["check_seven_certs", "check_type_q", "seven_adic_tight", "moser", "test_lemmas_AB"]:
        assert run(d, prog + ".py").stdout == (d / (prog + ".out")).read_text()


def test_finite_connection_sets_referee(tmp_path):
    """The referee of Lemma 22 (finite connection sets): the radii of Lemma 9 (sup r_theta = sqrt10/14 < 0.23 on
    (2/7, 1/3]) and the 7-adic characters on the 140 vectors of Theorem C, whose tight sets carry positive relations."""
    d = tmp_path / "indep_FW2"
    shutil.copytree(os.path.join(P, "indep_FW2"), d)
    out = run(d, "check_radii.py").stdout
    assert "max |y|^2 over X = 2" in out and "sup s*sqrt10/3 = 0.2258769757" in out
    out = run(d, "thmC_tight7.py").stdout
    assert out.count("all values in {2..5}/7: True;  tight-set sizes: [35]  positive relation found & verified exactly for all: True") == 2


def _stage_witness(tmp_path):
    for n in ["check_witness.py", "witness_q11sum.json.gz", "witness_q11sum.cnf.gz"]:
        shutil.copy(os.path.join(FW, n), tmp_path / n)
    return tmp_path


def test_finite_witness_q11(tmp_path):
    """The explicit finite witness for chi_c(Q(sqrt11)^2) = 7/2 (finite_witness/): both checkers, which share no code,
    accept the graph (A + A for the 76-vertex graph, every unit-distance pair an edge), the (7,2)-colouring and the
    180 cycles; the stored formula is the one check_witness.py builds, and the second encoding is the one that kissat,
    drat-trim and cake_lpr refuted (verification.txt)."""
    d = _stage_witness(tmp_path)
    out = run(d, "check_witness.py", "witness_q11sum.json.gz", "witness_q11sum.cnf.gz").stdout
    assert "(1) 2237 points, 11300 edges" in out and "not listed as edges: 0 (induced)" in out
    assert "(2) the colouring is a (7,2)-colouring" in out and "(3) 180 cycles" in out
    assert "(4) the CNF is the formula built from the edges and cycles; sha256 426f67eec64a8034" in out
    out = run(FW, "verify_independent.py", str(tmp_path / "H.cnf")).stdout
    assert "witness points = A + A: yes" in out and "unit pairs: 11300 ; witness edges: 11300 ; equal: True" in out
    assert "(7,2)-colouring: yes" in out and "cycles: 180 simple cycles" in out
    digest = hashlib.sha256((tmp_path / "H.cnf").read_bytes()).hexdigest()
    assert digest == "ee36797440e0c7b52f16592f2bfc18a4292853d65c1aefe744db1e0d432088ab"


@pytest.mark.slow
@pytest.mark.skipif(not DRAT, reason="drat-trim not found (PATH or DRAT_TRIM)")
def test_finite_witness_q11_proof(tmp_path):
    """drat-trim verifies the stored DRAT proof that every (7,2)-colouring of the witness has a tight listed cycle."""
    d = _stage_witness(tmp_path)
    with lzma.open(os.path.join(FW, "witness_q11sum.drat.xz"), "rb") as src, open(d / "witness_q11sum.drat", "wb") as dst:
        shutil.copyfileobj(src, dst)
    out = run(d, "check_witness.py", "witness_q11sum.json.gz", "witness_q11sum.cnf.gz", "witness_q11sum.drat", DRAT,
              timeout=3600).stdout
    assert "(5) drat-trim: s VERIFIED" in out and "chi_c(H) = 7/2" in out


# The grown and minimised witnesses (finite_witness/README.md): file stem, d, vertices, edges, listed cycles, sha256 of
# the formula of check_witness.py, sha256 of the formula of verify_independent.py.
GROWN = [
    ("witness_q11", 11, 170, 468, 879, "f36a7e6a79e5459d6c7f735be8143129af0b16651278b3f27e02d7c8a6f68517",
     "f8b9ac2dc0113f731406c1715521b4ba4c1c554a749182dd37856d4ab5d9186d"),
    ("witness_q191", 191, 293, 803, 489, "e62238610df1e6ec1899e619ecb94d8968d67dbc527dfc7dfa272866d4f7e9fd",
     "2c343503bcc20307c4bf16edcf1d98afbfd52a535efdaff441e2400d3644229c"),
    ("witness_q455", 455, 175, 434, 91, "0e06b8498615d6dea0bc54cbd39ad7ae979134369c1288587768550b61663839",
     "9994a1a50544107d07c5fa8677c9e9e1c6bf27c609ccc276c2ff2aff43afa022"),
]


@pytest.mark.parametrize("stem,d,nv,ne,nc,sha1,sha2", GROWN)
def test_finite_witness_grown(tmp_path, stem, d, nv, ne, nc, sha1, sha2):
    """The witnesses grown from the 4-chromatic graphs and minimised: both checkers, which share no code, accept the
    graph (every unit-distance pair an edge), the (7,2)-colouring and the cycles; the stored formula is the one
    check_witness.py builds, the second encoding is the one that kissat, drat-trim and cake_lpr refuted
    (verification.txt), and every H - v has a (7,2)-colouring with an acyclic tight digraph."""
    for n in ["check_witness.py", "check_critical.py", stem + ".json.gz", stem + ".cnf.gz", stem + "_critical.json.gz"]:
        shutil.copy(os.path.join(FW, n), tmp_path / n)
    out = run(tmp_path, "check_witness.py", stem + ".json.gz", stem + ".cnf.gz").stdout
    head = "" if d == 11 else f"Q(sqrt{d}): "
    assert f"(1) {head}{nv} points, {ne} edges, every edge at distance exactly 1" in out
    assert "not listed as edges: 0 (induced)" in out
    assert "(2) the colouring is a (7,2)-colouring" in out and f"(3) {nc} cycles" in out
    assert f"(4) the CNF is the formula built from the edges and cycles; sha256 {sha1}" in out
    out = run(FW, "verify_independent.py", str(tmp_path / "H.cnf"), os.path.join(FW, stem + ".json.gz")).stdout
    assert f"field Q(sqrt{d})" in out and "points of A in the witness:" in out
    assert f"unit pairs: {ne} ; witness edges: {ne} ; equal: True" in out and "(7,2)-colouring: yes" in out
    assert hashlib.sha256((tmp_path / "H.cnf").read_bytes()).hexdigest() == sha2
    out = run(tmp_path, "check_critical.py", stem + ".json.gz", stem + "_critical.json.gz").stdout
    assert f"{nv} vertices: for every v" in out and "vertex-critical" in out


def test_finite_witness_critical_rejects(tmp_path):
    """check_critical.py rejects a certificate whose colouring of H - v gives some edge a forbidden difference, and
    one whose tight digraph has a directed cycle (the stored colouring of H itself, with v uncoloured, has one when v
    is off every tight listed cycle)."""
    stem = GROWN[0][0]
    for n in ["check_critical.py", stem + ".json.gz"]:
        shutil.copy(os.path.join(FW, n), tmp_path / n)
    W = json.load(gzip.open(os.path.join(FW, stem + ".json.gz"), "rt"))
    C = json.load(gzip.open(os.path.join(FW, stem + "_critical.json.gz"), "rt"))
    i, j = W["edges"][0]
    v = next(u for u in range(len(W["points"])) if u not in (i, j))
    bad = json.loads(json.dumps(C))
    bad["critical_colourings"][v][j] = bad["critical_colourings"][v][i]
    with gzip.open(tmp_path / "bad1.json.gz", "wt") as f:
        json.dump(bad, f)
    r = run(tmp_path, "check_critical.py", stem + ".json.gz", "bad1.json.gz", check=False)
    assert r.returncode != 0 and "not a (7,2)-colouring of H - v" in r.stderr
    tight = [C2 for C2 in W["cycles"] if all((W["colouring"][C2[(k + 1) % len(C2)]] - W["colouring"][C2[k]]) % 7 == 2
                                               for k in range(len(C2)))]
    assert tight, "the stored colouring has a tight listed cycle"
    v = next(u for u in range(len(W["points"])) if u not in tight[0])
    bad = json.loads(json.dumps(C))
    bad["critical_colourings"][v] = [-1 if u == v else c for u, c in enumerate(W["colouring"])]
    with gzip.open(tmp_path / "bad2.json.gz", "wt") as f:
        json.dump(bad, f)
    r = run(tmp_path, "check_critical.py", stem + ".json.gz", "bad2.json.gz", check=False)
    assert r.returncode != 0 and "has a directed cycle" in r.stderr


@pytest.mark.slow
@pytest.mark.skipif(not DRAT, reason="drat-trim not found (PATH or DRAT_TRIM)")
@pytest.mark.parametrize("stem", [g[0] for g in GROWN])
def test_finite_witness_grown_proof(tmp_path, stem):
    """drat-trim verifies the stored DRAT proof that every (7,2)-colouring of the witness has a tight listed cycle."""
    for n in ["check_witness.py", stem + ".json.gz", stem + ".cnf.gz"]:
        shutil.copy(os.path.join(FW, n), tmp_path / n)
    with lzma.open(os.path.join(FW, stem + ".drat.xz"), "rb") as src, open(tmp_path / (stem + ".drat"), "wb") as dst:
        shutil.copyfileobj(src, dst)
    out = run(tmp_path, "check_witness.py", stem + ".json.gz", stem + ".cnf.gz", stem + ".drat", DRAT,
              timeout=3600).stdout
    assert "(5) drat-trim: s VERIFIED" in out and "chi_c(H) = 7/2" in out


def test_finite_witness_q7(tmp_path):
    """The nine-point witness for chi_c = 3 over Q(sqrt7) (finite_witness/README.md): check_small.py checks the points,
    all unit pairs, all 3^9 maps (each of the 84 proper 3-colourings has a tight listed cycle), that there is no
    homomorphism to K_{8/3}, and the criticality certificates; the formula it writes is the one that kissat, drat-trim
    and cake_lpr refuted (verification.txt).  The graph is the Wagner graph with one chord subdivided: the 8-cycle
    P_0 ... P_7, the chords P_0P_4, P_1P_5, P_2P_6 and the path P_3 S P_7 (S = vertex 8)."""
    out = run(FW, "check_small.py", "witness_q7.json.gz", str(tmp_path / "q7.cnf")).stdout
    assert "(1) Q(sqrt7): 9 points, 13 edges, every edge at distance exactly 1, and no other unit pair" in out
    assert "(3) all 19683 maps V -> Z/3: 84 (3,1)-colourings, each with a tight cycle" in out
    assert "(4) no homomorphism to K_8/3" in out and "H is vertex-critical" in out
    digest = hashlib.sha256((tmp_path / "q7.cnf").read_bytes()).hexdigest()
    assert digest == "258bcc7fc1dea349bc33aa8c896dda337e42e8b90074c5f1e32c09434d8ee6fb"
    W = json.load(gzip.open(os.path.join(FW, "witness_q7.json.gz"), "rt"))
    wagner = [(i, (i + 1) % 8) for i in range(8)] + [(0, 4), (1, 5), (2, 6), (3, 8), (7, 8)]
    assert sorted(tuple(sorted(e)) for e in W["edges"]) == sorted(tuple(sorted(e)) for e in wagner)


@pytest.mark.parametrize("what,message", [("point", "an edge is not at distance 1"),
                                          ("edge", "the edges are not all the unit pairs"),
                                          ("cycles", "no listed cycle is tight"),
                                          ("certificate", "tight cycle in H - v")])
def test_finite_witness_q7_rejects(tmp_path, what, message):
    """check_small.py rejects a moved point, a missing edge, a cycle list that misses some 3-colouring, and a
    criticality certificate whose tight digraph has a directed cycle."""
    W = json.load(gzip.open(os.path.join(FW, "witness_q7.json.gz"), "rt"))
    if what == "point":
        W["points"][8][2] += 1
    elif what == "edge":
        W["edges"].pop()
    elif what == "cycles":
        W["cycles"] = W["cycles"][:1]
    else:
        col = W["colouring"]
        tight = [c for c in W["cycles"] if all((col[c[(k + 1) % len(c)]] - col[c[k]]) % 3 == 1 for k in range(len(c)))]
        v = next(u for u in range(9) if u not in tight[0])
        W["critical_colourings"][v] = [-1 if u == v else c for u, c in enumerate(col)]
    with gzip.open(tmp_path / "bad.json.gz", "wt") as f:
        json.dump(W, f)
    shutil.copy(os.path.join(FW, "check_small.py"), tmp_path / "check_small.py")
    r = run(tmp_path, "check_small.py", "bad.json.gz", check=False)
    assert r.returncode != 0 and message in r.stderr


def test_finite_witness_q7_found_by_growth(tmp_path):
    """grow.py with (p, q) = (3, 1) from q7_seed.json, then minimise.py and critical.py, give the stored nine-point
    witness and its certificates, up to the translation by (1, 0) and the order of the vertices."""
    for n in ["grow.py", "minimise.py", "critical.py", "q7_seed.json"]:
        shutil.copy(os.path.join(FW, n), tmp_path / n)
    out = run(tmp_path, "grow.py", "q7_seed.json", "A", "G7", "2000", "200", "3", "1").stdout
    assert "UNSAT with 607 points" in out and "chi_c(H) = 3/1" in out
    run(tmp_path, "minimise.py", "G7.json", "W7.json")
    run(tmp_path, "critical.py", "W7.json", "C7.json")
    W = json.load(open(tmp_path / "W7.json"))
    C = json.load(open(tmp_path / "C7.json"))
    S = json.load(gzip.open(os.path.join(FW, "witness_q7.json.gz"), "rt"))
    pos = {tuple(p): i for i, p in enumerate(S["points"])}
    m = [pos[(p[0] + S["denominator"], p[1], p[2], p[3])] for p in W["points"]]    # found vertex -> stored vertex
    assert sorted(m) == list(range(9)) and (W["p"], W["q"]) == (3, 1)
    assert sorted(tuple(sorted((m[a], m[b]))) for a, b in W["edges"]) == sorted(tuple(sorted(e)) for e in S["edges"])
    assert all(S["colouring"][m[v]] == W["colouring"][v] for v in range(9))
    assert [[m[v] for v in c] for c in W["cycles"]] == S["cycles"]
    assert C["not_critical"] == []
    assert all(S["critical_colourings"][m[v]][m[u]] == C["critical_colourings"][v][u] for v in range(9) for u in range(9))


def test_small_triangle_free_six():
    """small_triangle_free.py on 6 vertices: the count of triangle-free graphs agrees with a direct enumeration."""
    out = run(FW, "small_triangle_free.py", "6").stdout
    assert ("triangle-free graphs on 6 labelled vertices: 5789; maximal: 211; not bipartite: 180; without a "
            "homomorphism to K_8/3: 0") in out


@pytest.mark.slow
def test_small_triangle_free_eight():
    """Nine vertices are the fewest for chi_c = 3 in a plane without unit triangles: every triangle-free graph with at
    most 8 vertices maps to K_{8/3}."""
    out = run(FW, "small_triangle_free.py").stdout
    assert ("triangle-free graphs on 8 labelled vertices: 4682270; maximal: 15247; not bipartite: 15120; without a "
            "homomorphism to K_8/3: 0") in out
