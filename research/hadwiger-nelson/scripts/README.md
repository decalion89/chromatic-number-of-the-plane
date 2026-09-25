# Scripts: the maintained tools

The 84 scripts here are the project's working tools. The 726 one-off experiments behind the
research log are in [`experiments/`](experiments/). Run every tool from `research/hadwiger-nelson/`:
each finds the `hn` package from its own location.

**Start here.**
- `verify_pair.py` rebuilds a claimed witness or gadget from its JSON file, re-derives every edge
  exactly and asks three SAT solvers (plus kissat and drat-trim when given).
- `grow_lean.py` is the colouring-guided growth used by the search for six.
- `fieldscreen.py` and `module_gate.py` decide, before any search, whether a field or a module can
  hold a 6-chromatic graph at all.
- `worker_setup.sh` installs kissat, drat-trim and tabu2 on a fresh machine.

## Verification — run these before believing any claim

| script | what it does |
|---|---|
| `verify_pair.py` | Independent check of a claimed obstruction, before a word is said about it. |
| `verify6.py` | Independent verification of any claimed 6-chromatic graph, before a word is said. |
| `verify_gadget.py` | Independent verification of a claimed distance-2 gadget, before a word is said. |
| `orbit_witness_test.py` | Is a graph, with its unit edges plus the edges at one or more Galois orbits of distances, 5-colourable? |
| `lattice_witness.py` | Build, verify (3 pysat solvers + kissat/DRAT) and shrink a lattice witness: points of Z[w]/sqrt(-3) (written as Eisenstein integers a + b w, real distance^2 = N/3), edges at the forbidden norms. |
| `verify_forced_same.py` | Verify a claimed forced pair c(A) = c(B) at /AB/ = 5 and its lambda-closure, from scratch. |
| `verify_twisted.py` | Build the twisted coset colouring EXPLICITLY on a saved graph and check it. |
| `verifyfield.py` | Verify the field that has both: de Grey's spindle and blocking. |
| `verifygate.py` | Cross-check a blocking verdict with a second solver, and by hand. |
| `verifyquot.py` | Check the periodic colouring on real points, not just on the quotient. |
| `circverify.py` | Exact check of a circular colouring c(x) = floor(5 * frac(phi(x))) on a module and on a grown graph. |
| `verify229.py` | Verify the 229-pair obstruction independently before it is claimed. |
| `verify4.py` | Verify the four-colour claim that the whole contrast rests on. |
| `verify49.py` | Verify the colouring, independently of the search that produced it. |
| `verify5.py` | Independent check of rho = 5 on Sa u rot(Sa). |
| `verify7.py` | Independent check: are those seven vertices really rainbow-forcing? |
| `verifySa.py` | Verify the colouring, and note what the local search was telling us. |
| `ei_rebuild.py` | Rebuild Exoo-Ismailescu's 6-chromatic {1,2}-graph from the paper and re-verify it. |

## Growth and search for obstructions

| script | what it does |
|---|---|
| `grow_lean.py` | Colouring-guided growth with a numpy candidate pool (grow_kw.py, about 10x less memory). |
| `grow_kw.py` | Colouring-guided growth with local search first (lean memory). |
| `grow_ls2.py` | Colouring-guided growth with local search first (lean memory). |
| `asym_grow.py` | Grow unit-distance graphs asymmetrically, because every generator here symmetrises and symmetry is what floors the core. |
| `make_pair_seed.py` | Make a MODE=apart / MODE=same growth seed from a graph JSON: the same points, units and colouring, with the target pair A, B set (grow_lean.py reads the keys A and B). usage: make_pair_seed.py graph.json A B out.json |
| `mk_kseed.py` | Exoo-Ismailescu's K inside our modules. |
| `backbone.py` | Backbone of the 5-colourings of a grown graph: is ANY pair already forced? |
| `search_forced.py` | Hunt for a forced monochromatic pair under k colours, then spindle it. |
| `search_disjunction.py` | Search for a forced disjunction, the weakest hypothesis the spindle needs. |
| `gadget.py` | Search the space of LEMMAS, not the space of graphs. |
| `six.py` | The same lemma, one floor up. |
| `measure_fk.py` | Measure f(4): the smallest unit-distance configuration with a pair forced monochromatic under four colours. |
| `is63.py` | Ask the target question directly: is rho(G,5) at most 63? |
| `rhotight.py` | The decision loop, fixed: shrink the hitting set before asking for an escape. |
| `bigfilter.py` | Forced-pair search that does not enumerate the pairs first. |
| `pairsat.py` | Every pair at one given squared distance, by the same bucketing the unit edge finder uses. |
| `blocked.py` | The blocked point: one SAT call per candidate, and equivalent to the answer. |

## Colouring by local search

| script | what it does |
|---|---|
| `tabucol.c` | tabucol: K-colouring by tabu search (Hertz–de Werra, with Galinier–Hao incremental tables). Compile with gcc -O2. |
| `tabu2.c` | tabu2: a proper K-colouring on hard edges that minimises conflicts on soft pairs without breaking a hard edge. Compile with gcc -O2. |
| `tabucol.py` | TabuCol, written properly this time, and calibrated before it is believed. |
| `tabu.py` | Local search, to tell a hard SAT from an UNSAT without waiting for either. |
| `tabuSaP.py` | Settle the grown seed at five colours with the calibrated local search. |
| `worker_setup.sh` | Set up a fresh machine for the growth searches: Python packages, kissat, drat-trim and tabu2. Idempotent. |

## Gates: periodic, circular and local colourings that rule a search out

| script | what it does |
|---|---|
| `gate.py` | The gate, exactly: coset colourings over the module the edge vectors generate. |
| `module_gate.py` | Split places matter for modules: reduce a graph's module at a place and colour the finite image. |
| `residuegate.py` | Residue gate: colourings of a module through one place above an unramified prime p. |
| `ramgate.py` | Residue gate at a prime p that may ramify: embed Q(sqrt g_1, ..., sqrt g_n, i) in L = Q_p(sqrt r, sqrt p) (r a non-residue unit), for every choice of signs, and reduce the unit vectors mod sqrt p.  Residues lie in F_{p^2} = F_p(sqrt r).  If every unit is integral at the place, any colouring of Cay(F_{p^2}, residues) colours every graph along these units. |
| `circgate.py` | Circular colourings through a real character: c(x) = floor(5 * frac(phi(x))), phi in Hom(M, R/Z). |
| `circrel.py` | The circular gate in relation space. |
| `circsat.py` | Rational circular colourings by CP-SAT: characters phi = y/Q, y in Z^r, with every unit in the window. |
| `circextend.py` | Does a circular colouring of a module M1 extend to a bigger module M' (e.g. its lambda-closure)? |
| `cellsample.py` | Sample the cell of F(M1) around a known point and test every sample's extensions to M'. |
| `kempegate.py` | Kempe swaps of coset colourings, exactly. |
| `prodgate.py` | The product gate: periodic colourings through NON-cyclic quotients. |
| `quotgate.py` | The full-quotient gate: 5-colour the Cayley graph of M/qM on the unit vectors. |
| `quotient.py` | Periodic colourings beyond homomorphisms to Z_5. |
| `idealquot.py` | Periodic colourings through an ideal of Z[omega]: Cay(M / alpha M, U). |
| `idealquot2.py` | Colourings periodic under ker(psi) cap alpha*M, and whether they can split a 5e pair. |
| `twistquot.py` | Twisted periodic colourings psi(x) + k(x mod alpha M) for an ideal alpha of Z[omega]. |
| `torusgate.py` | The torus gate: a PERIODIC colouring of the physical plane that is proper for the module's unit directions only.  f : R^2 -> [5], constant on the cells of an h-grid of the torus R^2 / (P Z)^2; for every unit vector v of the module and every cell C, every cell meeting C + v must get another colour (conservative, so any SAT answer is a genuine colouring).  Then c(x) = f(position of x) properly colours the whole unit-distance graph of the module -- and it has nothing to do with any coset colouring.  A 5-colouring here would refute rigidity outright. usage: torusgate.py <module.json> <P> <G (cells per side)> [K] |
| `torecore.py` | Unwrap the UNSAT core of an ideal torus. |
| `wallgate.py` | The wall gate: 1-dimensional colourings c(x) = G(phi(x), psi(x)). |
| `nearadm.py` | Blocked module: coset colourings fail only on the few unit vectors lying in 5M.  Look for 'nearly admissible' psi (zero exactly on those) and repair them with a small second coordinate chi : M -> Z/m that is nonzero there: any proper 5-colouring of Cay(Z/5 x Z/m, (psi, chi)(U)) would be a periodic colouring of the whole blocked module -- and would end the plain search. |
| `screen_lambda.py` | Does a rotation that is NOT integral at 5 escape the coset colourings? |
| `stiemke.py` | Exact certificates that NO twisted colouring exists at all on a module. |
| `farkas.py` | Exact certificates that NO twisted colouring keeps a given 2e pair alike. |
| `kappa_check.py` | kappa = omega * rho7 = (-11 + 5 sqrt(-3)) / 14 is congruent to 1 mod 5 in Z[omega, 1/7]. Check on the module: for every admissible psi and every unit u with kappa^k u also a unit of the module (k = +-1, and also rho7, omega alone), compare psi(kappa u) with psi(u). |
| `cosetfit.py` | How close is a grown graph's colouring to a coset colouring, and where are the defects? For every admissible psi, the best relabelling pi of the 5 colours maximises #{x : pi(c(x)) = psi(x)}; report the best fit, then where the mismatched vertices sit (distance from the pair). |
| `periodfind.py` | Does a grown colouring hide a period?  For many small lattice vectors t (differences and sums of unit vectors, small multiples of units) measure how often c(x + t) = pi(c(x)) for the best relabelling pi, over all x with both ends in the graph.  A colouring that is the restriction of a periodic colouring of the whole module shows rate 1 on its period lattice. |
| `quot2d.py` | The lattice N spanned by the translations that a grown colouring nearly respects, and the quotient M / N it would factor through: Smith form, the images of the unit vectors, and whether the quotient Cayley graph (made finite by a torus in the free part) is 5-colourable. |
| `rigid.py` | Rigidity at finite scale: colour a graph by tabu search (several seeds) and measure the best coset fit R = max_{psi, relabelling} #{x : pi(c(x)) = psi(x)} / n.  R = 1 for a coset colouring, about 0.2 for an unstructured one.  Units are taken from the module file given (e.g. five_rho7.json). |

## Fields and unit vectors

| script | what it does |
|---|---|
| `fieldscreen.py` | Screen a multiquadratic CM field K = Q(sqrt-d_1, ..., sqrt-d_n) by its non-split places. |
| `fieldtypes.py` | Every arithmetic type, decided once. |
| `field24.py` | Q(zeta_21, sqrt-11): the smallest field with triangles, blocking and a spindle. |
| `reduce11.py` | Colour unit-distance graphs in Q(sqrt-3, sqrt-11, sqrt-247) by reduction at the place above 11. |
| `moser2adic.py` | The 2-adic 4-colouring of the Moser field. |
| `msqrt.py` | Square roots in a multiquadratic field, by the obvious recursion. |
| `sqrtK.py` | Square roots in a multiquadratic field, done by characters not by search. |
| `allunits.py` | Every unit vector of an edge module, not just the ones the graph uses. |
| `lamclosure_units.py` | All the unit vectors of the lambda-closure M' = M1 + lambda M1 that come for free: U1 (units of M1, closed under the 60-degree turn), lambda U1, and the Exoo-Ismailescu units w_e = 5(1 - lambda) e (/w_e/ = 5 /1 - lambda/ = 1, and w_e = 5e - 5 lambda e lies in M'). Writes a module file (units + one dummy point) for circgate.py / circextend.py. |
| `rotunits.py` | Unit vectors generated by rotations of a field, written as a module file for the circular gates. |
| `unitsinqm.py` | Unit vectors of a graph's edge module: the rank of the module the edge vectors span, and which unit vectors it contains. |
| `tune.py` | The field of a 5-chromatic unit-distance graph, as a parameter. |

## Measurements on balls and growth curves

| script | what it does |
|---|---|
| `ball2.py` | The complete combinatorial 2-ball of five_rho7: O + u + v for all unit vectors u, v. Colour it by tabu search and measure the best coset fit R.  Does rigidity show up in a graph that contains EVERY point within two unit steps of a vertex? |
| `ballr.py` | Combinatorial r-balls of five_rho7 in exact lattice coordinates, and what their colourings do on the unit circles at the centre. |
| `freecurve.py` | free@5 against mean degree, averaged instead of sampled once. |
