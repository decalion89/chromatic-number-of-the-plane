# Hadwiger–Nelson: research log

> **About this file.** This is the project's chronological research log, kept as written:
> results, dead ends, corrections and retractions, in the order they happened. It was the
> project README until 25 September 2026. For a summary of what is established, see
> [`../README.md`](../README.md). Paths below are relative to the root of the
> repository, which held the project in `research/hadwiger-nelson/` until release 1.1.0.
>
> **Corrections.** Claims that later turned out wrong are kept where they were made. A
> correction made on the spot is marked **Corrected**, **Withdrawn** or **Retraction**, and
> the early ones are collected in "Corrections to my own claims, kept rather than edited
> away". These sections correct claims made further back: "A correction worth its own
> section: the gap is 10.8, not 700"; "Corrected: density was not the missing ingredient";
> "Two mistakes in the instrument, and what they were hiding"; "The null model was wrong
> twice, and the raw numbers say it better"; "Exoo–Ismailescu rebuilt, and the denominator
> that the project filtered out" (the multiquadratic claim of "Why every known
> construction stops at five"); "The circular gate: colourings through a real
> character" (the sampled cyclic gates); "Correction: `χ(ℚ(√3, √11)²) = 4` is Fischer's
> theorem (1994)" (the status entry below and "An open question closed"); and "Pre-release
> audit (26 September)" (the whole-field results, the local colourings, the note, and the
> descriptions of the data and certificates). The
> state of every result is given by the README and the notes, not by this log.

**How many colours does the plane need, so that no two points at distance exactly 1
share a colour?**

Posed around 1950. Still open. The answer, written χ(ℝ²), is known only to lie in
**{5, 6, 7}**.

| bound | value | who, when | how |
|---|---|---|---|
| lower | ≥ 4 | Nelson, 1950; L. and W. Moser, 1961 | the 7-vertex Moser spindle |
| lower | ≥ 5 | Aubrey de Grey, 2018 | a 1581-vertex unit-distance graph with no 4-colouring |
| lower | ≥ 5 | Jaan Parts, 2020 | the same, down to 509 vertices |
| upper | ≤ 7 | Isbell, 1950 | a hexagonal tiling of diameter just under 1 |

Nothing has moved since.

## What this is

Lowering the upper bound to 6 would mean constructing a 6-colouring of the entire
plane — an infinite object, not something a search can produce. **Raising the lower
bound to 6 is different: it needs a single finite unit-distance graph with no proper
5-colouring.** That is a finite object. Finite objects can be searched for, and found
ones can be checked by anybody.

So that is the target here, and the whole package is built around making any claim
falsifiable by a stranger:

- coordinates are **exact** — elements of ℚ(√3, √11), never floating point, so
  "distance exactly 1" is a decidable predicate rather than a tolerance;
- colourability is decided by a **SAT solver**;
- results ship as **certificates**: a colouring re-checked against edges re-derived
  in O(n²) from the published coordinates, or a **DRAT proof verified by drat-trim**
  against a formula rebuilt from those same coordinates.

A certificate never asks you to trust this code. It asks you to run `drat-trim`.

## Status

**The machinery works and is validated. The open problem is not solved.** Concretely:

- ✅ χ(ℝ²) ≥ 4 reproduced from scratch and **machine-checked end to end** —
  `certificates/moser_spindle_no3coloring.json`, verified by drat-trim.
- ✅ The de Grey spindling argument is **automated**, and rediscovers the classic
  step unaided: given only a unit rhombus, it finds that 3 colours force the two far
  tips to match, spindles on that pair, and lands on the Moser spindle.
- ✅ A **three-copy pigeonhole variant** of the argument, which needs only a forced
  *disjunction* rather than a forced pair, implemented and validated.
- ✅ χ(ℝ²) ≥ 5 **reproduced and machine-verified.** de Grey's 1581-vertex graph
  was rebuilt from the published 39-point set S and its recipe; every count
  matches (39 → 397 → 1581, with exactly one coincident point). kissat 4.0.4
  returns UNSATISFIABLE for 4-colourability under two configurations, and
  **drat-trim verifies the proof**: `s VERIFIED`, 13 140 458 lemmas, 2 016 499 in
  core, 130 857 426 resolution steps, 522 s.
  `certificates/degrey_1581_no4coloring.json`.

  Exactly what is verified, stated precisely: *no proper 4-colouring of this
  graph assigns colours 0, 1, 2 to one pinned triangle.* One step remains outside
  the proof — a triangle's three vertices are pairwise adjacent, so any proper
  4-colouring gives them three distinct colours, and colours are interchangeable,
  so permuting them to read 0, 1, 2 loses nothing. Standard, and still not
  something drat-trim checked. The plain unbroken formula removes even that step
  and was still solving.

  > **Later.** That run was stopped after 110 minutes, on 19 September, to free
  > a core, without a verdict. The unpinned formula is unsolved here.

- ❌ χ(ℝ²) ≥ 6 — the actual goal. Not found.

- ✅ **An open question from the literature, closed: `χ(ℚ(√3, √11)²) = 4`.**
  `ℚ(√3, √11)` is the smallest field whose plane holds a Moser spindle, and
  the question of its chromatic number had stood since 2010.
  - Moorhouse (2010) left the value undetermined.
  - Madore (2015) proved `4 ≤ χ ≤ 5`.
  - Exoo–Ismailescu (2018) and Polymath16 asked whether a 5-chromatic
    unit-distance graph embeds there.
  - Voronov (Polymath16, 2021) wrote that `χ = 4` "seems likely" and that
    nobody had proved it.

  It does not. At a place over 2 the field is inert in `ℚ(i, √3, √11)`, so
  every unit vector reduces to a nonzero element of `𝔽₄`, and that residue is
  a 4-colouring of the whole plane. The 2-adic idea is Speyer's: in Polymath16
  (2018) he used it to colour the Moser ring. What is new is the extension to
  the whole plane. See *An open question closed* below and
  `tests/test_q311.py`.

  > **Corrected later.** K. G. Fischer had proved `χ(ℚ(√3, √11)²) = 4` in 1994;
  > see "Correction: `χ(ℚ(√3, √11)²) = 4` is Fischer's theorem (1994)". The
  > proof here is a short alternative.

- ✅ A **19-vertex, 33-edge graph with no 3-colouring**, drat-trim verified and
  **vertex-critical** — `certificates/genuine_pair_19_no3coloring.json`. Not a
  record; the Moser spindle reaches χ ≥ 4 on seven. What is new is that its
  forced pair is forced **only jointly** — neither leg is forced on its own, so
  there is no single target to spindle and the classical argument has nothing to
  grip. Built rather than found, over ℚ(√3,√11)(√v) with v = (−66 + 30√33)/256;
  that extension is forced, since v has a negative conjugate and every
  multiquadratic field is totally real.

- ✅ The **first graph here with no coset 5-colouring** — 15 313 vertices over
  ℚ(ζ₇), built from every modulus-one step of denominator ≤ 29. Every other
  construction in this package, de Grey's included, admits a colouring that is
  a homomorphism of the edge module to ℤ/5; this one admits none. It is
  nevertheless **3-chromatic**, because ζ₆ ∉ ℚ(ζ₇) and so it has no triangle.

- ✅ **Why the method stops at five, in one number.** Squeezing a pivot's circle
  into two colours gives every auxiliary that sees both colours the *same* list
  of size k−2, and that subgraph — the confined set — is measured
  **2-degenerate** at every size from 72 to 420 points. A 2-degenerate graph is
  not always 2-choosable (an odd cycle is not), which is what makes k = 4 work;
  it is always 3-choosable, which is what kills k = 5. The mechanism falls on
  the wrong side of that gap by exactly one.

- ✅ **A module where coset colourings are provably rigid at every large
  scale.** `ρ₇ = (1+4√−3)/7` and `ω` combine into `κ = ωρ₇ ≡ 1 (mod 5)`, an
  irrational rotation that every coset colouring is blind to. On the 803-type
  graph together with its `ρ₇`-images, **3 840 of 3 840** exact Stiemke
  certificates show that every colour class positively spans. From that:
  - **no twisted colouring** exists;
  - every coarse colouring is a coset colouring;
  - **rigidity propagates:** a colouring that is a coset colouring on a
    half-space, or on a large enough ball, is that coset colouring everywhere.

  Exoo–Ismailescu's graph sits inside the module, turned by 150°. See *κ, the
  rotation every coset colouring is blind to* and the sections after it.
  *Update:* nowhere-locally-coset colourings **do** exist there. A circular
  colouring through a character of order 280 is one, so the strong rigidity
  conjecture is false. It still colours Exoo–Ismailescu's pair alike. See
  *The circular gate*.
- ✅ **The circular gate.** Colourings `⌊5·frac(φ(x))⌋` through a real
  character `φ` generalise coset colourings, and a MILP finds them.
  - They 5-colour the blocked module `803 ∪ λ(803)`: exact check on a
    46 496-point graph. That search was futile from the start.
  - A 6-chromatic unit-distance graph needs a module with no such colouring.
  - If every finite set of unit vectors has such a colouring, then `χ(ℝ²) = 5`.

## Pressure, and why every k = 5 search was dead on arrival

Write **pressure(p)** for the least number of colours N(p) can be squeezed into,
over all k-colourings. `hn/forced.py` computes it in at most k incremental
solves, and two theorems make it the right invariant.

**A core of size r at p requires pressure(p) ≥ k − r.** Take a colouring
attaining the minimum: at least k − pressure(p) colours are free for p while
c(T) offers at most r values, so if k − pressure(p) > r, recolour p to a free
colour outside c(T). Properness at p asks only that its colour avoid c(N(p)), so
the colouring stays proper and T was not a core.

**A k-vertex-critical graph has no core of any size, at any vertex.** Colour
G − p with k−1 colours and give p the kth: p is then the *only* vertex carrying
it, so c(p) lies outside c(T) for every T at once. de Grey's G is
**not** 5-vertex-critical — and that claim, which this package leaned on
throughout, is **retracted**:

> **G − 1420 is still 5-chromatic.** Vertex 1420 has degree 4; removing it and
> asking for a 4-colouring with a triangle pinned to 0, 1, 2 comes back UNSAT —
> CaDiCaL in 1581 s, **and Glucose independently in 1900 s on a formula rebuilt
> from scratch.** A proper subgraph on 1580 vertices is already 5-chromatic.

What was offered as evidence before was *separability*, which is a
**consequence** of criticality rather than a proof of it — a consequence read
backwards. The warning sign was there: others have published 5-chromatic
unit-distance graphs of several hundred vertices.

The theorem is unaffected — a k-vertex-critical graph does have no core of any
size. It simply says nothing about G. What survives, because it was measured
rather than inferred: pressure exactly 2 at all 1581 vertices, no forced pair
in 29 930 queries, no forced-different non-edge in 40 539 pairs, and 1201
vertices not rainbow-forcing. **Those results stand and are now unexplained** —
criticality was the explanation, and it is gone.

Measured pressure at k = 5: **exactly 2 at all 1581 vertices of G**, the two
hubs of degree 60 included, in two seconds — and still 2 on every union tried,
up to 5533 vertices. So the smallest core available at five colours is **three**.

### The core condition is itself a pressure measurement

T is a core of p exactly when **min |c(N(p) ∪ T)| = k**. If some colouring left a
colour unused on both, recolouring p to it is proper and puts c(p) outside c(T);
conversely the colouring witnessing a failure leaves c(p) itself unused on both.
So cores can be *built forwards* instead of shrunk: ask whether a colouring
leaves the first colour free on N(p) ∪ T — by symmetry that is the whole question
— and either it proves T is a core or it hands back a colouring, and any vertex
carrying that colour elsewhere kills it when added. One solve per vertex added.

It reproduces the 19-vertex construction's pair in six solves, and it **fails** on
de Grey's G at k = 5, which is exactly what the criticality corollary demands.

### Where blocking stops

Every leg's conflict graph on the copies has maximum degree 2 — a point of a
circle is one apart from at most two points of that circle — so α ≥ m/3, and the
counting certificate Σ α(G_L) < m needs r·m/3 < m, hence **r ≤ 2**. The
19-vertex configuration sits exactly on that boundary: alphas [2, 3] against 6
copies, blocking by one.

So the two bounds meet. Pressure 2 puts every core at five colours at three or
more; counting blocks only up to two.

**A core of three can still be blocked, by misalignment.** Rotations alone never
manage it — searched exhaustively over every N ≤ 30 and every triple of angles on
the Nth roots of unity, an escape always exists, while the same search finds 298
blocking *pairs*. Reflections break it, because each leg picks up its **own**
angle in the cross-orbit conflicts: 27 344 configurations block. The smallest is
N = 3, three legs at d² = 1/3, six copies over ℚ(√3) — alphas [2,2,2] against 6
copies, exactly the counting floor, so counting says nothing while SAT and a
brute force over all 3⁶ choices agree that it blocks.

**And three is the only size that blocks at all.** Four, five and six legs were
searched the same way — forty thousand random (N, angles, shifts) draws apiece —
and not one blocks, while three blocks on the fifth draw. The counting slack says
why: with each conflict graph of maximum degree 2 the independent sets available
total r·m/3, which at r = 3 is exactly the copy count, the knife edge where
overlap can still decide the question, and at r = 4 is a third more room than
there are copies to place.

Put beside the pressure bound this **pins the size exactly**. A core of r needs
pressure ≥ k − r, so pressure 2 at five colours forces r ≥ 3; blocking forces
r ≤ 3. Any configuration that could give χ(ℝ²) ≥ 6 through a pivot and its
isometries has a core of **exactly three**, with its legs' radii and angles
matching one of the 27 344 blocking patterns. That is a far smaller target than
"a forced pair somewhere".

What is missing is the core itself, not the block.

### The machine that makes pressure, taken apart

Exactly one configuration in this package reaches pressure 3, and it is worth
naming rather than measuring. Delete everything from Sa that can go while the
pivot's neighbourhood still refuses to be squeezed into two colours, and
**47 vertices** survive — `certificates/pressure3_witness_47.json`:

- the pivot;
- its circle of 30, which splits into **five hexagons** — the 60° orbits, as the
  bipartiteness argument requires — each with exactly two alternating
  2-colourings, so five independent orientation bits;
- **sixteen** further points, at d² = 1/3, (7 ± √33)/6 and (3 ± √33)/6, each
  adjacent to exactly two circle points — and in all sixteen cases those two lie
  in **different** hexagons, so each one reads the relative orientation of a pair
  of them.

  > **Corrected later** (pre-release audit, 26 September). The last pair is
  > (9 ± √33)/6; (3 − √33)/6 is negative. Exactly: 1/3 and (7 − √33)/6 five
  > times each, (7 + √33)/6, (9 − √33)/6 and (9 + √33)/6 twice each.

Squeeze the circle into two colours. Every one of the sixteen that sees two
differently-coloured circle points is barred from both, so it is confined to the
remaining **k − 2**. And the sixteen among themselves form a 16-vertex, 27-edge
unit-distance graph that is 3-chromatic and **not bipartite**: an odd cycle.

**The mechanism is list colouring.** Squeezing the circle into two colours hands
each auxiliary point a *list* — the colours its circle neighbours leave it. Two
differently-coloured neighbours give a list of k − 2; two of the same colour give
k − 1. The auxiliary graph never faces one palette; it faces a list assignment,
and the squeeze is refused exactly when that assignment admits no proper
colouring. Over all 32 orientations:

| | lists | refused |
|---|---|---|
| k = 4 | sizes 2 and 3 | **32 of 32** → pressure 3 |
| k = 5 | sizes 3 and 4 | **0 of 32** → pressure 2 |

The odd cycle is the special case where every list is the same pair, and it is
load-bearing where it applies — deleting three of the sixteen makes their graph
bipartite and the pressure drops from 3 to 2 at once. But it covers only half the
orientations. In the other sixteen the confined set *is* bipartite, and the
minimal refusal there is 22 vertices using **eleven** auxiliaries of which only
six are confined; the other five carry lists of size 3 and sit at d² = 1/3. A
frame with only "confined" and "free" in it cannot express that. Lists can.

### de Grey's forcing, rebuilt from one angle — and a 127-vertex gadget

Three hexagon offsets suffice where de Grey uses five. Put hexagons at
**0, θ/2 and θ** on the pivot's unit circle, with cos θ = 5/6, and take as
auxiliaries every **u + v** with u, v in different hexagons — which is every
auxiliary there is, since a point one away from u and v is the pivot reflected
across the chord uv, i.e. u + v. That is 127 points and 528 edges over
ℚ(√3, √11), and it measures **pressure 3 at four colours** — a core of one is
possible, the classical spindle regime — and pressure 2 at five.
`certificates/three_hexagon_pressure3.json`.

### The hexagons are generated by the Moser angle

Reading the offsets off the witness, the five hexagons sit at

> **0, ±θ/2, ±(60° − θ)** modulo 60°, where **cos θ = 5/6**

and cos θ = 5/6 is the Moser rotation itself: the chord it subtends on the unit
circle squares to 2 − 2cos θ = **1/3**, the split-prime quotient at 3 in
ℚ(√−11) that the spindle's angle turned out to be. Verified in the field rather
than in degrees — the cosines between hexagon 0 and the rest come out exactly
1, √33/6 twice and (5 + √33)/12 twice, with 2(√33/6)² − 1 = 5/6 the half-angle
identity and (5 + √33)/12 = cos 60° cos θ + sin 60° sin θ.

So the same angle runs through the whole construction: it is the rotation the
spindle uses, the radius its targets sit on, and the offset that spaces the
hexagons whose orientation bits carry the forcing.

Two structural facts constrain any redesign. Auxiliaries are **Minkowski sums** —
the point one away from circle points u and v is the pivot reflected across the
chord, which is just u + v — so the auxiliary set is H₁ + H₂ + …, and two of them
are one apart exactly when |Δu + Δv| = 1, an equation with a *discrete* solution
set in the hexagons' relative angle. Scanning that angle finds nothing, because
the good values are measure zero; solving it gives ten per pair, and the
half-angles only appear once there are three hexagons, since they come from
edges *between* different hexagon pairs. And the always-confined points — those
whose two circle neighbours are adjacent — are exactly the ones at √3 from the
pivot, where the chord-1 angle is 2 arcsin(1/(2√3)) = θ, which does not divide
360°, so their graph is a union of paths: **bipartite, never even an odd cycle.**
Orientation-dependence is forced, not a quirk.

**So the lift is an exact and standard question, and it is about choosability
rather than chromatic number:** at five colours the auxiliary graph must fail to
be list-colourable with lists of sizes 3 and 4. A 4-chromatic subgraph all of
whose vertices carry the *same* list of three is one sufficient way — which is
why the Moser spindle and the 19-vertex jointly forced construction above are the
pieces to try — but it is not the only way, and the non-uniform lists are where
de Grey's construction gets its extra reach.

Searched over that family: hexagons chosen from de Grey's own offsets, all
Minkowski-sum auxiliaries, every orientation. Three hexagons at 0, θ/2, θ
refuse **8 of 8** orientations at four colours — that is the 127-vertex gadget
above. At five colours, with up to five hexagons and 360 auxiliaries, **0 of 32**
orientations refuse. The lists of size 3 and 4 are simply too generous for a
depth-one design.

The scope of that negative is worth stating exactly: it covers auxiliaries
adjacent to two circle points and nothing deeper.

**The second level was built, and it does not help.** Every point one away from
two circle points is a Minkowski sum u + v, so the depth-one design is already
exhausted — no further confined point exists without enlarging the circle. Points
attached to *auxiliaries* instead have no list at all, only propagation, so only
a solver sees them: adding all 498 of them to the gadget, 625 vertices and 3324
edges, leaves the pressure at 2.

### What a core of three actually needs, and why the gadget cannot host one

Pressure 3 was never the target — the theorems pin the core size at **exactly
three**, and pressure 2 is the condition that makes three right, not an obstacle.
But the gadget is **4-colourable**, so at five colours the pivot can always take
a fifth colour nobody else has and no core exists at all, of any size. That is
the criticality corollary again, from the other side: **a core at k = 5 needs a
host that is 5-chromatic with the pivot non-critical.** The gadget supplies
pressure; it cannot supply that.

And one honest measurement about how loose the pressure bound is. The gadget is
4-chromatic — not 3-colourable — with pressure 3 at four colours, which *permits*
a core of one. Its smallest core over all 127 pivots is **seven**. Necessary is a
long way from sufficient.

### The wall, stated precisely

Two halves that have never been in the same graph:

| | pressure | can host a core at k = 5 |
|---|---|---|
| three-hexagon gadget | 3 at k = 4, 2 at k = 5 | **no** — it is 4-colourable, so the pivot always finds a fifth colour of its own |
| de Grey's G | flat 2 | **no** — measured inert at every pivot; the criticality explanation is retracted, so this is now an observation without a reason |

Grafting was tried both ways and neither cures the other. Landing the gadget's
pivot on G's hub leaves that hub critical, because G is 5-vertex-critical and so
G − hub is 4-colourable, and the gadget minus its pivot is 4-colourable too;
the counterexample builder runs out, as it must. Placing a translated copy of G
so it *avoids* the pivot does make the pivot non-critical in principle, but over
204 placements the two never share more than a vertex or two, and near-disjoint
copies colour independently — their colours can be permuted against each other,
so no fixed target set is ever forced. Adding the whole second level, 625
vertices, leaves the gadget 4-colourable still.

So the requirement is sharp and unmet: **a 5-chromatic host, with the pivot
non-critical, substantially overlapping a pressure gadget.** Each piece exists
separately. Nothing here puts them in one graph.

### The obvious resolution, and why it is dead

There is a real tension. Rigid graphs have few colourings and so small cores —
but the most rigid are the vertex-critical ones, and criticality kills cores
outright. Loose graphs escape criticality, but their cores are huge: on
`G ∪ f(G)` a counterexample can park the free colour anywhere, and hundreds of
rounds of counterexample-killing close nothing.

The obvious way out is to add **one** point to a k-critical graph. In W = G + v
every old vertex u still has W − u ⊇ G − u, which is (k−1)-colourable, so u
stays critical; only v has W − v = G. So v is the unique non-critical vertex,
the unique possible pivot, and G keeps all its rigidity.

**It is provably useless.**

> **Theorem.** Let G be k-vertex-critical and W = G + v. Then every core of v
> contains all of V(G) ∖ N(v) — so the minimal core is the whole graph but for
> the pivot's neighbours.
>
> *Proof.* Fix u ∉ N(v). G − u is (k−1)-colourable, and v's neighbours all lie
> in it, so W − u is too; colour it with 1…k−1 and give u the colour k. Then u
> is the *only* vertex carrying k and none of v's neighbours does, so recolouring
> v to k stays proper. In that colouring c(v) = k = c(u) and nothing else has
> it, so any target set avoiding u misses c(v). ∎

Verified on the Moser spindle, which is 4-critical: adding any exact point one
away from two of its vertices gives a pivot of degree 2 whose minimal core is
**exactly** the five vertices outside its neighbourhood, every one forced in.

So a usable core needs a **third** condition on top of the other two: no vertex
may be leavable alone in a colour. Nothing here satisfies all three.

### What the three conditions really ask for

They have a classical name between them. What makes a core hard to shrink is
that colour *classes* move: T has to meet every class that could be free at the
pivot, and in a loose graph the classes are large and mobile. **Pin the classes
and T shrinks to one representative each.**

> In a **uniquely k-colourable** graph — one whose k-colouring is unique up to
> permuting colours — a pivot of pressure q has a core of exactly **k − q**.

At five colours with the pressure 2 that every graph here measures, that is a
core of exactly **three**, which is precisely the size the blocking bounds
allow. So the remaining object has a name:

> **A uniquely 5-colourable unit-distance graph**, with a pivot whose three
> free-class representatives sit at radii and angles matching one of the 27 344
> blocking patterns, gives **χ(ℝ²) ≥ 6**.

The mechanism is visible one level down, and it recovers the oldest fact in the
subject. The triangular lattice **is** uniquely 3-colourable — its colouring is
the Eisenstein residue modulo (1 − ω) — and measuring a pivot there gives
pressure 2 and a core of size **one**, at squared distance 3. That is the
classical rhombus: √3 forces two points to agree at three colours *precisely
because the lattice's colour classes cannot move*.

It also closes the loop with criticality. A graph is uniquely k-colourable
exactly when the forced-same relation has k classes, and a k-vertex-critical
graph has no forced-same pair at all. **The two theorems are the two ends of one
axis**, and the object wanted sits at the far end from de Grey's G.

### Why the frontier is exactly at five

Two facts already proved combine into a ladder, with no search in it. A unit
circle is bipartite, so a neighbourhood containing an edge has pressure **2 for
free** and nothing more without ambient help. The pressure theorem then makes
the smallest possible core `k − 2`, and blocking reaches `r ≤ 3` and no further.

| k | free core | blocked by |
|---|---|---|
| 3 | **1** | the classical spindle. Rigidity is *free* here, since 2 = k−1 means the circle alone determines the centre's colour — which is exactly why the triangular lattice is uniquely 3-colourable and √3 forces agreement |
| 4 | **2** | counting. de Grey's construction |
| 5 | **3** | misalignment only, and only with reflections — the free pressure **exactly saturates** the blocking bound, with no slack anywhere |
| 6 | **4** | *nothing blocks four.* The method is not hard here, it is impossible |

So five is not the frontier because nobody searched hard enough. It is the last
value of k at which the free pressure of a unit circle and the largest blockable
core still meet — and they meet exactly.

### The distance to the target, in one number

> **Theorem.** A uniquely k-colourable graph has pressure **exactly k − 1 at every
> vertex.** Pressure is at most k − 1 always, since a vertex's own colour never
> appears in its neighbourhood; and if one had two free colours, switching between
> them would move it to a different class and give a genuinely different partition,
> not a permutation. ∎

That turns "how close is this graph to the target" into a per-vertex measurement
with a gradient — `unique_colouring_defect` reports `k − 1 − pressure(v)` — and
it is monotone in added points, so a search can climb it.

- **k = 3**: pressure 2 is exactly the free pressure of a unit circle, so unique
  3-colourability costs nothing. Measured on a triangular-lattice patch: defect
  **zero at all 37 vertices**.
- **k = 4**: pressure 3 is one above free. Sa reaches it at some vertices — that
  is where its forcing comes from.
- **k = 5**: pressure 4 is two above free, and nothing measured here reaches even
  3. Every graph tried comes back at a flat 2.

**The gap is 2 against 4**, in the same units as everything else in this file.

### One invariant, and the whole problem in one table

Strip the pivot out of the core condition and a single graph invariant is left:

> **ρ(W, k)** = the least size of a set that uses all k colours in **every**
> k-colouring — a rainbow-forcing set.

Every bound above is a statement about it. In a k-vertex-critical graph **ρ = n**
(colour W − u with k−1 and give u the kth; a set omitting u misses that colour,
so every vertex is needed — that *is* the inertness of de Grey's G). In a
uniquely k-colourable graph **ρ = k**, one representative per class.

And it converts straight into cores, with the pivot doing the work:

> **Theorem.** If S is rainbow-forcing and p is any point, **S ∖ N(p) is a core
> of p** — because N(p) together with it contains S. So the core has size
> `|S| − |S ∩ N(p)|`, and all that is needed is a pivot adjacent to `|S| − 3`
> elements of S.

|S| = 4 needs one such element, placed trivially; |S| = 5 needs two, which any
pair less than 2 apart supplies through its circle intersections; |S| = 6 needs
three concyclic at radius **exactly one**, and then p is their circumcentre. (One
trap, met on the first attempt: p must not lie in S. A minimal forcing set loses
the property when any element is dropped, so a circumcentre that happens to be a
forcing vertex removes its own target — measured on Sa as a perfectly good
circumcentre of degree 30 whose core was not one.)

Measured:

| graph | n | k | ρ |
|---|---:|---:|---:|
| triangular-lattice patch | 37 | 3 | **3** — equals k, as unique colourability requires |
| de Grey's Sa | 397 | 4 | **6** |
| de Grey's Y | 791 | 4 | **9** |
| de Grey's G | 1581 | 5 | **1581** — critical, so every vertex is needed |
| G ∪ (G+t), best overlap | 2938 | 5 | **> 400** |

**That is the whole gap in one number.** At three and four colours the forcing
sets are tiny, so cores of the blockable size are within reach of a well-placed
pivot. At five they explode past 400, and a core of three is not a matter of
searching harder — it is 400 away.

### And stacking copies can never close it

A forcing set meets every colour class of every colouring — equivalently, every
independent set I with `W − I` still (k−1)-colourable. That reading kills the
whole family of constructions this package kept returning to.

> **Theorem.** Let W split as A, B and a shared part, with `W − a − b`
> (k−1)-colourable for every non-adjacent a ∈ A, b ∈ B. Then every forcing set
> contains **all of A or all of B**, so **ρ ≥ min(|A|, |B|)**.
>
> *Proof.* For such a pair, colour W − a − b with k−1 and give a and b the kth;
> they are non-adjacent, so it is proper, and {a, b} is then a colour class. A
> forcing set must contain a or b. Over all pairs it is a vertex cover of the
> complete bipartite graph between A and B — and the only vertex covers of that
> are A and B. ∎

The hypothesis is exactly what a union of two k-critical graphs supplies:
deleting one vertex from each copy leaves both (k−1)-colourable. Verified on two
Moser spindles glued at two vertices — **all 22** non-adjacent cross pairs are
killable, and the forcing set comes back as one whole side plus both shared
vertices, ρ = 7 against the bound's 5.

For `G ∪ (G+t)`, |A| = |B| = **1357**. So ρ ≥ 1357 there, and the measured
"> 400" was not the loop running out of patience. **No stack of copies of a
critical graph can ever have a small forcing set, hence never a small core,
hence never a blockable one.** Every union in this package was dead before it
was built.

And it gets *worse* with more copies. With three pairwise-overlapping copies a
forcing set must hit every cross **triple**, so it is a vertex cover of a
complete tripartite 3-uniform hypergraph and has to swallow **two** whole parts.

### When ρ is n, and therefore when it can be small

Two disjoint reasons force ρ = n, and between them they cover everything here at
five colours.

- **k > χ(W).** Every vertex is then removable: χ(W − u) ≤ χ(W) ≤ k−1, so u can
  be left alone in the kth colour and every forcing set needs it.
- **W is k-vertex-critical.** Same conclusion, same colouring.

So ρ can only be small where **χ(W) = k exactly and W is not vertex-critical**.
Measured:

| graph | k | ρ | critical (sampled) |
|---|---:|---:|---:|
| three-hexagon gadget | 4 | **9** | **0 / 12** |
| three-hexagon gadget | 5 | > 123 | 12 / 12 |
| 19-vertex joint core | 4 | 19 = n | 12 / 12 |

The gadget has both properties at four colours and is useless at five for a
reason that has nothing to do with its geometry — it is 4-chromatic. G is
critical. The unions fall to the cross-pair theorem. **Nothing has both at
five.**

Which returns the same object from the other side: few realisable colour classes
is what makes ρ small, and a graph with exactly one 5-colouring up to permutation
has exactly five of them, giving ρ = 5. The target is a **uniquely 5-colourable
unit-distance graph** — which by the pressure theorem needs pressure 4, every
neighbourhood using four colours, at *every* vertex.

## Why the rigidity was never there: coset colourings

Let M be the ℤ-module generated by a graph's **edge vectors**. A group
homomorphism

> φ : M → ℤ/n with φ(d) ≠ 0 for every edge vector d

*is* a proper n-colouring — c(p) = φ(p), and adjacent points differ by an edge
vector. It is constant on cosets of ker φ, so it is as regular as a colouring
gets: the exact opposite of rigidity, since it composes with the automorphisms
of ℤ/n and with any translation of the module, so colour classes move freely.

Searched by SAT over the coordinates, so dimension is no obstacle, and verified
directly rather than on the solver's word:

| graph | dim | edge vectors | ℤ/5 |
|---|---:|---:|---|
| three-hexagon gadget | 8 | 18 | **exists** |
| 19-vertex joint core | 16 | 27 | **exists** |
| de Grey's G | 32 | 133 | **exists** |

**Every graph in this package is 5-colourable by cosets.** That is why every
pressure measurement came back at a flat 2, why no forced-same pair ever
appeared in 8549 queries, and why the forcing sets are enormous — it was one
fact, not five coincidences.

It also gives a **screen** that costs seconds: a graph whose edge-vector module
admits a ℤ/5 homomorphism can never be 6-chromatic, whatever its size.

And it names the design target sharply. φ must avoid every edge vector, so
writing D for their images in M/5M, a φ exists unless the hyperplanes d^⊥ cover
the whole dual — which needs D to meet every hyperplane, a **blocking set**. The
smallest is a projective line, **q + 1 = six points**. So: *six unit vectors
whose reductions mod 5 represent the six points of a projective line* rule out
every coset colouring at once.

### Where a coset colouring can and cannot be blocked

> **Theorem (rank 2 never blocks).** If the edge vectors span only rank 2 mod 5
> then M is a plane lattice, every unit vector has the same norm N, and the
> value of the quadratic form on a projective point is well defined up to
> squares — so all of them land in **one square class, three of the six points**
> of the projective line. Three hyperplanes cover 13 of the 25 points of
> (ℤ/5)², so a φ always survives.

Measured on the Eisenstein lattice at norms 1, 3, 7, 13, 21, 49 and 91: every
one gives exactly **three** classes, and which three is decided by whether N is a
square mod 5 — {(0,1),(1,0),(1,1)} when it is, {(1,2),(1,3),(1,4)} when it is
not. Even **24** unit vectors, at norm 91, give three.

**Rank 3 fails too**: the unit vectors lie on a conic of PG(2,5), and a conic has
10 exterior lines. From **rank 4** upward Chevalley–Warning makes every hyperplane
carry vectors of the right norm, so blocking becomes possible in principle.

**And counting is not enough.** A random φ survives m hyperplanes with
probability (4/5)^m, so coverage needs roughly `m > 7.2 r`. That predicts the
gadget (r = 4, m = 18, threshold 29) correctly and de Grey's G (r = 16, m = 133,
threshold 115) **wrongly** — G clears the count and still admits a φ. Being many
points is not the same as being the right points on the quadric. Every
root-of-unity step set up to n = 105 admits one too: ℤ[ζ_n] has n directions
against rank φ(n), a ratio never above about 3.5.

### And the cheapest blocking set is unreachable

Blocking asks the directions to meet every hyperplane — a blocking set — and in
PG(r−1,5) the smallest is a projective **line**, six points. So the cheapest
route is six unit vectors whose directions are the six points of a line.

> **Theorem.** No projective line over F₅ has all six of its points carrying a
> quadratic form value in one square class. Hence the directions of unit
> vectors — all satisfying Q(d) = N for a single N — never contain a complete
> line.

Checked exhaustively over every binary form (a,b,c) over F₅ and every target N:
the maximum is **five**, reached only when the form degenerates to rank one. By
hand: anisotropic gives 3 of 6 (the square classes split the line evenly);
hyperbolic gives 2 (two points are isotropic, carrying Q = 0); rank one gives 5,
the sixth point being the radical.

Measured across pairs of imaginary quadratic fields — the only source here of
infinitely many unit vectors at *fixed* rank, since the modulus-one elements of
ℚ(√−d₁, √−d₂) are exactly the products of one from each factor, generated by
things like (3+√−7)/4 and the Moser rotation (5+√−11)/6 — the direction sets
reach **25 at rank 4** and contain **zero** complete lines, on every pair tried.

One trap, which produced a false positive before it was found: some generators
have denominators divisible by 5, such as (1+3√−11)/10. The module then stops
being 5-integral, clearing denominators makes *every* vector divisible by 5, and
the search reports "no homomorphism" for reasons of arithmetic bookkeeping and
nothing else. Generators must be filtered to denominators coprime to 5.

### Where a blocking set finally exists

The multiquadratic families fail quantitatively and finally. The projective
direction count of the modulus-one group mod 5 is a **product of small factors**
— 3 where 5 is inert in a quadratic factor, 2 where it splits — so
ℚ(√−d₁,…,√−d_t) of degree 2t gives at most **3^t** directions against
PG(2t−1,5)'s (5^2t−1)/4 points. Measured: 6 at (7,11), 9 at (7,23), **27** at
(7,23,43) against **97 656**. The ratio (3/25)^t collapses, and more generators
do not help — the group is finite and already exhausted.

The cause is the splitting: in a multiquadratic field the Galois group is (ℤ/2)^t
and the decomposition group at 5 is cyclic, so the residue degree is at most 2
and 𝒪/5 breaks into tiny fields. A *large* norm-one group needs the decomposition
group to be everything — a **cyclic** Galois group with 5 **inert**, i.e. ℚ(ζ_n)
with 5 a primitive root mod n.

**n = 7 qualifies**: 5 has order 6 in (ℤ/7)^×. Then 𝒪/5 = F₅⁶ and the norm-one
subgroup has (5⁶−1)/(5³−1) = **126** elements. Modulus-one elements come free
from Hilbert 90 — u = α/conj(α) has modulus one for every α, with denominator a
norm that only has to stay coprime to 5.

> **63 of those directions, at full rank 6, admit NO homomorphism to ℤ/5.**
> Verified by brute force over all 5⁶ = 15 625 maps rather than on the solver's
> word: **zero survive.**

It is the first step set here with no coset 5-colouring at all. What that is and
is not: it removes the structural colouring that every other graph in this
package had; it does **not** make any particular graph 6-chromatic, since an
unstructured 5-colouring may still exist and only a solver can say. It is a
necessary condition, met for the first time.

Two caveats, stated because they bound what the tool proves. The search embeds M
in ℤ^d by clearing a common denominator, so every φ it **finds** is genuine
(restriction is a homomorphism) while a "none found" is only a statement about
that embedding — it screens *in*, not *out*. And that clearing multiplies every
vector by one integer, which made the whole set look divisible by 4 and 6 until
the common content was divided out; the divisibility obstructions at n = 3, 4, 6
are real only after that correction.

### Where blocking begins: the prime 29

Grouping the modulus-one elements of ℚ(ζ₇) by denominator locates the
obstruction exactly. Taking all of them up to each bound:

| denominator ≤ | 1 | 2 | 4 | 8 | 11 | 16 | **29** |
|---|--:|--:|--:|--:|--:|--:|--:|
| directions | 7 | 21 | 31 | 33 | 35 | 37 | **87** |
| blocks ℤ/5 | no | no | no | no | no | no | **yes** |

Nothing under 29 ever blocks, however many elements are collected — and
**29 = 4·7 + 1 is the least rational prime splitting completely in ℚ(ζ₇)**.
Splitting is what supplies many independent modulus-one elements sharing one
small denominator, which is exactly what covering PG(5,5) needs. The
denominators appearing next — 43, 71, 113, 127 — are the following primes
≡ 1 mod 7, as that reasoning predicts.

Blocking is a covering problem, so it has an exact optimum, and
`minimum_blocking_set` computes it by MaxSAT rather than greedily — greedy
carries a ln(points) factor and reported 25 directions where the structure
allows far fewer. It recovers the projective line of six at rank 2, and
correctly reports that the triangular lattice's three directions cover nothing.

### The first blocked graph, and why it is only 3-chromatic

Every modulus-one element of denominator ≤ 29 as a step, two rounds:
**15 313 vertices, 30 276 edges**, realising all 87 directions.

```
87 edge vectors; coset 5-colouring: NONE -- blocked
3-colourable: True
```

The first half is a first for this package — the structural colouring that
every other construction here admitted is gone. The second half is a theorem,
not bad luck.

> **Theorem.** A unit-distance graph over a field K contains a triangle **iff**
> ζ₆ ∈ K.

Three mutually unit-apart points are equilateral, so one is another rotated 60°
about the third: multiplication by a primitive sixth root of unity, which must
therefore lie in K. Conversely ζ₆ builds one outright. The roots of unity of
ℚ(ζ₇) are exactly **μ₁₄**, and 6 ∤ 14, so **no graph over ℚ(ζ₇) has a single
triangle however large it grows.** Clique number 2, measured average degree
3.95 — and every high-chromatic unit-distance graph known leans on triangles
throughout.

### Blocking and folding do not trade off

A first reading of the accompanying table suggests they do:

| graph | n | m | \|S\| | rank | corank | χ | blocked |
|---|--:|--:|--:|--:|--:|--:|:--:|
| Moser spindle | 7 | 11 | 7 | 4 | 3 | 4 | no |
| triangular patch | 121 | 320 | 3 | 2 | 1 | 3 | no |
| de Grey S | 39 | 18 | 8 | 4 | 4 | 3 | no |
| de Grey Sa | 397 | 1974 | 15 | 4 | 11 | 4 | no |
| de Grey Y | 791 | 3938 | 33 | 8 | 25 | 4 | no |

Low rank everywhere, high corank, never blocked. But folding needs *relations*
among the steps — ℤ-independent steps build a tree, and a tree is bipartite —
while blocking needs hyperplanes covering the dual, and **both are monotone
increasing in the step set**: another step adds a hyperplane and can only add
relations. They never trade against each other. The known graphs miss blocking
for want of steps, not for want of rank, and the object to build is the
**union** of a folding set and a blocking one.

That names the field between them. Triangles need ζ₆, so 3 | n; blocking needs
a completely split prime, which is where ℚ(ζ₇) came from. The smallest
cyclotomic field with both is

> **ℚ(ζ₂₁) = ℚ(ζ₃, ζ₇), degree φ(21) = 12**

— carrying the whole Eisenstein lattice and its triangles while inheriting
ℚ(ζ₇)'s split primes. In it 5 has order 6 in (ℤ/21)^×, so 5 splits into two
primes of residue degree 6 rather than staying inert. And blocking is
**inherited upward for free**: a φ on the degree-12 module would restrict to a
homomorphism on the ℚ(ζ₇) submodule, nonzero on all 87 directions, and none
exists.

### The wall at five colours is exactly one colour wide

Squeezing the circle into two colours does not hand every auxiliary the same
list. At k = 5, with the circle in colours {0,1}:

- an auxiliary seeing **both** gets the list {2,3,4} — size 3, and *every*
  such auxiliary gets the **same** list, so colouring that group is ordinary
  3-colouring;
- an auxiliary seeing only one gets a list of size 4, which is slack.

So the condition is not on the whole auxiliary graph but on the **confined
set**, and only in the orientations that actually arise:

> The **local** mechanism behind pressure 3 — the one the k = 4 construction
> actually runs on — is the confined set's list problem, in every proper
> 2-colouring of the circle.

Stated that way deliberately. An earlier draft of this section claimed the
condition was *necessary*, and it is not: a 3-colourable confined set does not
by itself give a colouring of the whole graph, and pressure 3 could in
principle come from deeper structure. What follows explains the measurements
and bounds what the local mechanism can do — it is not a proof that the
pressure must be 2.

Two circles meet in at most two points, so no auxiliary ever sees three circle
points: the list sizes are exactly 3 and 4, with nothing in between.

A first guess that the auxiliary graph's chromatic number was the obstruction
is **wrong**: measured on de Grey's Sa at its best pivot, χ(A) = 4 already,
degeneracy 4, and the pressure at five colours is still 2. The lists are what
differ. Measuring the confined set instead, over all 32 orientations (the
circle is five hexagons, each a 6-cycle with two proper 2-colourings):

| χ(confined set) | 2 | 3 | 4 |
|---|--:|--:|--:|
| orientations | 4 | **28** | **0** |

**Maximum 3, needed 4, in all 32.** The wall is exactly one colour wide, and
uniformly so — not a few awkward orientations but every one of them.

### The confined set is indexed by a cut

Write `bits[a]` for the parity chosen on hexagon `a`; the circle point at
position `i` in it takes colour `(i + bits[a]) mod 2`. So the auxiliary
`u_i + v_j` drawn from hexagons `a` and `b` is confined exactly when

```
i + j + bits[a] + bits[b]   is odd.
```

Two consequences, both measured rather than assumed.

**Same hexagon:** the bits cancel, so those auxiliaries are confined in *every*
orientation. They are the points at **√3** from the pivot — and their graph has
maximum degree **one**: a matching, 9 disjoint edges on 18 points. They can
contribute 2 to a chromatic number and never more.

**Across hexagons:** only `eps_ab = bits[a] XOR bits[b]` matters, and `eps` is a
**cut** of `K_t`, so `eps_ab + eps_bc + eps_ac = 0` for every triple. The `2^t`
orientations give only `2^(t−1)` distinct confined sets, and no design can
choose the pair-parities independently — at t = 3 exactly 4 of the 8 patterns
arise, and every missing one breaks the triangle identity.

Measured maxima of χ(confined set) over all orientations:

| | auxiliaries | χ(confined) |
|---|--:|---|
| 2 hexagons | 48 | [1, 1] |
| 3 hexagons | 126 | [2, 3] |
| de Grey Sa, 5 hexagons | 150 | [2, 3] |

Four is what is needed, in every orientation, and nothing here reaches it.

### Why k = 4 works and k = 5 does not, in one number

Adding hexagons in de Grey's own angle family takes the auxiliaries from 126 to
798 and leaves χ(confined) at [2, 3] throughout. That asks for an explanation,
not more search, and the explanation is a single statistic.

**The confined set is 2-degenerate.** Worst orientation at each size:

| hexagons | auxiliaries | confined | max degree | degeneracy |
|--:|--:|--:|--:|--:|
| 3 | 126 | 72 | 4 | **2** |
| 4 | 240 | 132 | 4 | **2** |
| 5 | 390 | 210 | 4 | **2** |
| 6 | 576 | 306 | 4 | **2** |
| 7 | 798 | 420 | 4 | **2** |

The confined set grows nearly sixfold and the degeneracy does not move. And
degeneracy is exactly what the list-colouring argument turns on — peeling
low-degree vertices greedily gives

```
d-degenerate  =>  (d+1)-choosable
```

so a 2-degenerate confined set is **3-choosable**. That settles both colour
counts at once, in opposite directions:

- **k = 4.** Lists of size 2. A 2-degenerate graph need *not* be 2-choosable —
  an odd cycle is 2-degenerate with list chromatic number 3 — so the lists can
  fail, and de Grey's Sa reaches pressure 3.
- **k = 5.** Lists of size 3. Every 2-degenerate graph **is** 3-choosable, so
  the lists always complete and the pressure is 2. No number of hexagons
  changes that, because none of them changes the degeneracy.

> The frontier at five colours is not a search that has not yet succeeded. The
> mechanism the method runs on lives in the gap between *"2-degenerate is not
> always 2-choosable"* and *"2-degenerate is always 3-choosable"*, and at five
> colours it falls on the wrong side by exactly one.

Reaching pressure 3 **through this mechanism** needs a confined set of
**degeneracy ≥ 3** — a statement about the geometry of Minkowski sums of
hexagons, not about how hard the search is run. Reaching it some *other* way is
not excluded by any of this; nothing measured here does.

### How much of the design space that actually covers

Two auxiliaries from the **same** pair of hexagons differ by `D1 + D2`, and

```
|q - q'|^2 = |D1|^2 + |D2|^2 + 2 Re(t z) = 1,    t = r_a conj(r_b)
```

is a **line in t meeting the unit circle** — so for two hexagons the good
angles are a finite set. Enumerating all 31 × 31 pairs of hexagon differences:
48 solutions, **eight** distinct modulo 60°. Measured at each:

| degeneracy | 0 | 1 | 2 | ≥3 |
|---|--:|--:|--:|--:|
| angles | 1 | 2 | 3 | **0** |

One case is free and useless: `D2 = 0` — same `v`, adjacent `u` — holds at
every angle, but confinement of `u_i + v_j` turns on `i + j` and adjacent `i`
differ by one, so **every always-available edge joins a confined auxiliary to
a non-confined one**. That is why Pythagorean rotations (cos and sin both
rational, so the whole thing lives in ℚ(√3)) leave the confined set with no
edges at all — measured, degeneracy 0 across 40 configurations.

**What this does not show.** The enumeration is complete for *two* hexagons
only. With three, an edge can join auxiliaries from **different pairs** — `q`
from (a,b) and `q'` from (a,c) — which is two free angles against one
equation: a **curve**, not a discrete set. de Grey's θ/2 is exactly such a
case, absent from the eight, and χ(confined) is 1 at two hexagons and jumps to
3 at three.

> So the flat degeneracy measured as hexagons are added describes **one curve
> through a continuous space**, not the space. It explains why de Grey's family
> stops where it does. It does not close the route.

### Scanning the surface, not the curve

With three hexagons the confined set's edges come from equations relating
**two** free angles, so the space is a surface. It can still be scanned,
because a configuration is rich exactly where many equations hold at once.
Fixing `α₁`, an edge between `q = u_i + v_j` from pair (0,1) and
`q' = u_i' + w_l` from (0,2) needs

```
|(w^i - w^i') + e^(i a1) w^j - e^(i a2) w^l| = 1
```

which for `A` = the first two terms and `C = A conj(w^l)` reads
`Re(C e^(-i a2)) = |A|²/2`, so `a2 = arg C ∓ arccos(|A|/2)` whenever `|A| ≤ 2`.
Each `(i, i', j, l)` gives up to two values, 2160 in all, and the rich
configurations are the **histogram peaks**.

The method reproduces de Grey's own configuration as a check: at `α₁ = θ/2`
the value `α₂ = θ` carries multiplicity 144, and the configuration measures
126 auxiliaries, 306 edges, confined degeneracy 2 — agreeing exactly with the
field arithmetic.

**Scanned over 900 values of α₁ and the six richest α₂ at each — 5400
configurations — nothing exceeds degeneracy 2.**

Two measurement bugs are worth recording, because both made the scan *blind*
rather than wrong-looking: rounding auxiliary coordinates to 1e-6 while testing
distances against 1e-7 found 50 of de Grey's 306 edges (keys are now rounded
to 1e-9, where the arithmetic really carries 1e-15); and the largest histogram
peak is always `α₂ = α₁`, the same six points, so coincident offsets are
dropped before ranking rather than after.

### What the confined set looks like

On de Grey's own three-hexagon configuration every confined auxiliary has
degree **0, 2 or 4 — never odd**, with 42 of 72 isolated, 36 edges, and
components of cycle rank 1 or 4.

That looked like a law. It is not: over 1200 configurations sampled across the
three-hexagon surface, **1276 orientations carry an odd confined degree** and
degree 5 occurs. The even degrees come from the reflective symmetry of that
particular offset family, not from the construction.

What does hold, across 1200 configurations and 9600 orientations:

| largest confined degree | largest **degeneracy** | needed |
|--:|--:|--:|
| 5 | **2** | 3 |

Never 3 — which is what lists of size 3 would need to fail.

The peak method extends a hexagon at a time: given the offsets already placed,
an edge joining an auxiliary that uses the **new** hexagon to one that does not
reads `|r_new w^m − X| = 1` with `X` the old auxiliary minus the circle point
it pairs with — the same line-meets-circle problem. Over three free angles, 360
four-hexagon configurations: **still 2.**

So the cap holds over a *complete* enumeration at two hexagons, 9600
orientations at three, 360 configurations at four, and de Grey's own family out
to seven. Stated as a measurement, not a theorem: the surface is continuous and
only its histogram peaks were evaluated.

### And the obvious generalisation is vacuous

The confined set is one subgraph. The honest object is everything outside the
circle, where the squeeze gives `L(v) = k − |colours v sees|` — 3, 4 or 5 at
k = 5 — and greedy completes the list colouring unless some subgraph has
`deg(v) ≥ L(v)` throughout. That maximal surviving subgraph, the **L-core**,
would be the exact greedy obstruction. Smallest core over every orientation:

| | Sa (397) | Sb (397) | Y (791) | gadget (127) |
|---|--:|--:|--:|--:|
| k = 4 | 366 | 366 | 730 | 108 |
| k = 5 | 360 | 360 | 718 | 108 |

**Non-empty everywhere and at both colour counts** — these graphs have average
degree around ten and almost nothing peels. So the greedy criterion separates
nothing: it does not explain why the pressure is 3 at four colours and 2 at
five, and it yields no new sufficient condition.

That is why the degeneracy result above is stated about the confined set's
*mechanism* and not about the pressure. The mechanism is genuinely sufficient
when it fails — an odd cycle on one shared 2-list is infeasible outright, which
is what k = 4 runs on — and at five colours it never fails. Whether pressure 3
could come from somewhere else stays open.

### Degree efficiency: what blocking actually costs

Blocking and folding are both monotone in the step set, so they do not trade
off directly. **Density per vertex does**, and one statistic shows it:

```
efficiency = average degree / unit vectors available in the module
```

A rank-2 module is discrete, so a patch of it is **closed** — nearly every
point has all its neighbours present. A module of rank ≥ 3 is dense in the
plane, so no bounded region contains a closed piece, every ball is boundary,
and most available steps land outside it.

| graph | n | m | steps | avg degree | efficiency | blocked |
|---|--:|--:|--:|--:|--:|:--:|
| triangular patch | 121 | 320 | 6 | 5.29 | **88.2%** | no |
| de Grey S | 39 | 18 | 16 | 0.92 | 5.8% | no |
| de Grey Sa | 397 | 1974 | 30 | 9.94 | 33.1% | no |
| de Grey Y | 791 | 3938 | 66 | 9.96 | 15.1% | no |
| de Grey G | 1581 | 7877 | 134 | 9.96 | 7.4% | no |
| ℚ(ζ₇) ball | 15 313 | 30 276 | 174 | 3.95 | **2.3%** | **YES** |
| ℚ(ζ₂₁) ball | 16 015 | 32 722 | 178 | 4.09 | **2.3%** | **YES** |

Efficiency falls monotonically as the structures get richer, and blocking
appears only at the bottom. Note what does *not* fall: de Grey's three sizes
all sit at average degree **9.96**, unchanged. What collapses is the fraction
of the available density a finite construction manages to collect.

> The cost of blocking is not fewer neighbours. It is needing forty times as
> many steps to get the same ten.

### A blocked graph cannot be a lattice

Adding steps never hurts blocking or folding, but it does **dilute a ball**: at
radius 2 over 178 steps only about one step-application in ninety lands back
inside, which is why the blocked graphs here come out with average degree near
4. Density wants **few, highly related** steps — a lattice, where every
translate lands on a point. And there the obstruction is absolute.

> **Theorem.** A unit-distance graph whose edge vectors generate a *discrete*
> subgroup of ℝ² can never block.

A discrete subgroup of ℝ² is ℤ^r with r ≤ 2, so the dual is PG(1,5): six
points whose hyperplanes are single points, so covering needs **six**
directions. The unit circle meets a plane lattice in at most six points — the
hexagonal case — giving at most **three**. Three against six, never. The bound
is tight in the wrong direction: the triangular lattice attains three, and is
exactly the classical coset 3-colouring.

So a blocked graph must live in a module **dense** in the plane, where "every
point of a region" is not a finite construction and a ball's boundary
dominates its interior. Measured, in ℚ(ζ₂₁) with 6 Eisenstein steps and 174
blocking ones: 16 015 vertices, 32 722 edges, **1056 triangles**, blocked —
and 3-chromatic. Both halves achieved at once for the first time, and still
far from rigid.

### Blocking needs the whole module

The cheapest cover — a pencil of six hyperplanes — needs its six directions to
lie in one 2-dimensional subspace, where they would have to hit all six points
of a PG(1,5) while sharing one norm. No square class allows that, so a blocking
set of unit vectors is **never minimal**, and the question becomes how much rank
it takes.

Constructed rather than sampled (random subsets of ℚ(ζ₇)'s steps have full rank
almost always): span a submodule, then collect every unit step inside it — the
whole unit circle of that submodule, which is what a graph built there would
have.

| rank | submodules with ≥6 unit steps | richest | blocks |
|--:|--:|--:|:--:|
| 3 | 2 | 12 directions | no |
| 4 | 34 | 4 directions | no |
| 5 | 159 | 5 directions | no |
| 6 | the whole module | **87 directions** | **yes** |

Proper submodules simply do not hold enough unit vectors — six to fourteen,
where covering needs far more. **So blocking requires the full rank**, and the
full rank is exactly what makes the module densest in the plane and a ball's
interior vanish against its boundary.

> That closes the loop with the efficiency table: blocking costs 2.3% not by
> accident, but because nothing short of the whole module blocks at all.

### The caveat that bounds all of it

Blocking is necessary and very far from sufficient, and it is worth being
blunt about how far. Attach a **pendant** edge in a blocking direction — one
new vertex of degree one — and that direction joins the edge module and helps
cover the dual, while the chromatic number changes by nothing at all: a
degree-one vertex extends any colouring of the rest greedily. **Any graph
whatever can be made blocked without becoming one colour harder.**

What blocking buys is therefore not a bound but a *measurement that means
something*. While a coset colouring exists the colour classes slide freely and
every rigidity statistic reads flat — which is exactly what happened to every
forcing search here. Removing it is a precondition for the forcing machinery
to have anything to detect. The chromatic number still has to be established
by a solver, on a graph that is blocked **and** rigid, and only the second half
is hard.

### What pressure 3 actually needs

de Grey's Sa reaches **pressure 3 at four colours** — k−1, the classical spindle
regime — with cores down to 7. And it does so with **no forced pair anywhere**:
the forced-different graph on a pivot's circle is exactly its 30 unit-distance
edges, every degree 2, no odd cycle, no forced non-edge at all. The ambient graph
kills every 2-colouring of the circle *jointly*, without pinning any single pair.

That removes the object this package kept failing to find. Pressure 3 at five
colours needs no rainbow and no forced non-edge — it needs ambient rigidity, and
pressure measures exactly that, and never decreases when points are added. What
it has not yielded to is stacking: translated copies at the most frequent
difference vector (|d|² = 1/3, folding 224 vertices), a 60-degree rotation about
a centre folding 789 of 1581, and de Grey's own dihedral group applied to G — all
still flat 2, up to 5533 vertices.


## The missing object, stated precisely

Reading the theorems together narrows the gap more than any one of them does.

The pressure theorem says a core of size `r` at `p` needs `pressure(p) ≥ k − r`.
At five colours a core of **three** therefore needs pressure ≥ 2 — and the
measured pressure is exactly 2, everywhere, on every graph here.

> **Pressure does not exclude a core of three.** It never did. The flat 2 was
> read as a wall when it is only a floor, and the floor is exactly met.

But *exactly met* is the worst place to be, and it would be easy to read that
as more room than it gives. Every colouring must use all k colours on
`N(p) ∪ T`, and the circle supplies at least `pressure(p)` of them — so in the
colourings that squeeze the circle to its minimum, `T` must supply
`k − pressure(p)` colours, **all new**. At k = 5 with pressure 2 and r = 3 that
means three targets carrying three distinct colours, none of them the circle's
two, in every such colouring. Three vertices forced to be rainbow *and*
colour-disjoint from the circle is close to the original problem restated.

Room appears only where pressure **exceeds** k − r, and that is exactly what
four colours has: Sa's pressure 3 against a requirement of far less, which is
why its cores go down to 7 and ρ comes out at 7 on 397 vertices. The spindle
method has always run on that slack.

Blocking a core of three is solved: 27 344 misalignment configurations do it,
and the counting bound (`r ≤ 2`) together with the pressure floor (`r ≥ 3`)
pins the size at exactly three from both sides.

So neither pressure nor blocking is what is missing. What is missing is:

> A 5-chromatic graph `W` with a pivot `p` that is **not essential** —
> `χ(W − p)` still 5, so the criticality corollary cannot give `p` a colour of
> its own — carrying an actual **core of three**: `min |c(N(p) ∪ T)| = 5` for
> some `|T| = 3`.

The first half is easy and known: `W = G ∪ (G + t)` works, because removing any
vertex leaves a whole 5-chromatic copy. The second half has never been found on
anything tested here. **That single measurement is the gap.**

Note what it is *not*. Not a bigger search over de Grey's G, which the
criticality corollary settles in one line. Not more pressure, which is already
at the required floor. Not a better blocking configuration, of which there are
27 344. **One core, of three, at five colours.**

### And all of it is one number after all

A core of size `r` at `p` means `min |c(N(p) ∪ T)| = k`: every colouring uses
all `k` colours on that set. That is exactly what **ρ**, the rainbow-forcing
number, measures — the least size of a set using all `k` colours in every
colouring. So a core of size `r` at `p` forces

```
rho <= deg(p) + r
```

On de Grey's G, whose maximum degree is 60, a core of three needs **ρ ≤ 63**.
Its actual value there is now **open**: the reading ρ = n rested on G being
vertex-critical, which is retracted above. What is measured is ρ ≤ 1580, and
that 1201 vertices are not forcing. The cross-pair bound still holds in form
but needs a critical graph to apply to, which G is not.

So the folded union is no longer excluded in advance. It was simply checked,
and came back empty:
its exact 60° rotation centre folds **789 of 1581**, giving 2373 vertices and
11 832 edges, and the forward core construction finds **no core in 60 steps at
any of four pivots**, with pressure 2 throughout. Confirmation, not discovery.

> So the ladder, the pressure, the cores, the criticality corollary and the
> cross-pair bound are all the same statement about ρ, and the target is:
>
> **a 5-chromatic unit-distance graph with ρ ≤ 63.**

For comparison ρ = k exactly when the graph is uniquely k-colourable. So what
is wanted is a unit-distance graph that is *nearly uniquely 5-colourable* — and
every graph measured in this package has ρ = n or ρ ≥ n/2. That is the whole
distance to the goal, in one number.

### The frontier in one measurement

The ρ reading is only worth anything if it comes out right where the
construction actually succeeds. At four colours de Grey's Sa reaches pressure 3
with cores down to 7 at a pivot of degree 30, which **predicts ρ(Sa, 4) ≤ 37**.
Measured on the same 397 points, one colour apart:

| | ρ |
|---|---|
| Moser spindle, k = 4 | **7** = n (4-vertex-critical, so ρ = n) |
| de Grey Sa, k = 4 | **7** — on 397 vertices |
| de Grey Sa, k = 5 | the same construction runs past **358** without closing |

Seven vertices out of 397 use all four colours in every colouring. That is the
rigidity the whole spindle method runs on, and the prediction is met with room
to spare. **Add one colour to the same graph and the construction that closed
at seven does not close at 358.**

Stated carefully: the 358 bounds the *greedy path*, not ρ itself, since a
different set of that size might still be forcing. What it shows is that the
counterexample construction — the same one, on the same points — is in
completely different regimes at four and five colours.

> That is the frontier, in one measurement, on one graph.

### Small ρ does not need a small critical subgraph

One direction is immediate: a k-chromatic subgraph uses all k colours in every
k-colouring of the whole graph, so

```
rho(G, k) <= |H|   for every k-chromatic subgraph H of G.
```

At four colours that would explain everything — Sa contains Moser spindles,
seven vertices, 4-critical — and it is the obvious reading of ρ(Sa, 4) = 7.

**It is the wrong reading.** The seven vertices the construction returns are
`S = [0, 3, 4, 6, 7, 8, 9]`, which induce **five edges**, have chromatic number
**3**, and include **three isolated vertices**. Not a spindle, not critical, not
4-chromatic. Re-derived from a CNF built from scratch rather than through the
same code: for each of the four colours, no proper 4-colouring of Sa leaves
that colour off S. And minimal — all seven single deletions break it.

> So the forcing is **ambient**. It is carried by the other 390 vertices, not
> by anything inside the set, and the converse of the theorem is false.

That matters for what is left. A core of three at a degree-60 pivot needs
ρ ≤ 63. Had small ρ required a small k-chromatic subgraph, this would be asking
for a 5-chromatic unit-distance graph on 63 vertices, against a published record
of around five hundred — hopeless. **It does not.** It asks for ambient rigidity
of exactly the kind Sa already exhibits at four colours, with a set that is
nearly edgeless.

### One pivot cannot reach five, and nor can seventy-six per cent of the graph

ρ(Sa, 4) = 7 has a mechanism behind it: the pressure at the pivot is 3, so
N(p) always carries three colours, and `p` together with a small piece of its
circle already forces all four. **One pivot suffices at four colours.**

At five it cannot. Pressure is 2 everywhere, so `{p} ∪ N(p)` forces three
colours at most — the pivot's own plus the circle's two — and three is not
five. Measured on de Grey's G, where each test is a single SAT call returning
in under two seconds:

| set | vertices | forcing |
|---|--:|:--:|
| one hub's closed neighbourhood | 61 | no |
| adjacent hub pairs | 38 | no |
| the 133 highest-degree closed neighbourhoods, united | **1201** | **no** |

Seventy-six per cent of the graph, and a 5-colouring still exists leaving one
colour off all of it — against **seven vertices out of 397** at four colours.

**What that proves, and what it does not.** A superset of a forcing set is
forcing, so a non-forcing set contains *no* 5-chromatic subgraph: any
5-chromatic subgraph of G must use a vertex outside those 1201. It does *not*
bound ρ itself, since some other set of that size might force.

It also bears on the criticality question from an unexpected side. If G had a
small 5-chromatic subgraph it would sit in the dense part, and this says it
does not — consistent with G being vertex-critical after all, and with the
published smaller 5-chromatic graphs being separate constructions rather than
subgraphs of this one.

### ρ only goes down when structure is added

If G sits inside W, every proper k-colouring of W restricts to one of G, so a
set forcing in G forces in W:

```
rho(W) <= rho(G)   whenever G is contained in W.
```

**That reverses the strategy.** The way to a small ρ is a *bigger* graph, not a
smaller one.

It was hidden while the cross-pair theorem was thought to apply to unions of
copies of G. The theorem is correct, but its hypothesis asks that `W − a − b`
be (k−1)-colourable for every non-adjacent cross pair — which for copies of G
means `G − a` must be 4-colourable, exactly the vertex-criticality now
retracted. So the bound `ρ ≥ 1357` for `G ∪ (G + t)` is **withdrawn**, and with
it the conclusion that every union here was dead before it was built.

Verified where ρ is computable by brute force: the Moser spindle has ρ = 7 at
four colours, and embedding it in a 13-vertex rotated union keeps ρ at 7 with
the same witness — never rising.

### Unioning copies drives ρ down to k + 1

Monotonicity gives `ρ(W) ≤ ρ(G)`, but ≤ is not <, and the strategy rests on the
difference. Measured at four colours, where ρ is cheap — Sa unioned with
rotated copies of itself through the half-Moser angle, each graph strictly
containing the last:

| copies | n | m | minimal ρ |
|--:|--:|--:|--:|
| 1 | 397 | 1974 | **7** |
| 2 | 619 | 3324 | **5** |
| 3 | 829 | 4638 | 5 |
| 7 | 1645 | 9558 | 5 |

**It drops**, from 7 to 5, and then holds. Five is one above the floor, since
ρ ≥ k always. So enlarging works — and works immediately: one extra copy
captures whatever rigidity there is, and further copies add nothing.

The five vertices are the interesting part. On `Sa ∪ rot(Sa)` they are
`[107, 208, 502, 568, 618]` — five points spanning seven edges in three
overlapping triangles, and **3-chromatic**: 208 can take 107's colour and 618
can take 502's. A 3-chromatic set of five points forces all **four** colours,
because of the 614 vertices around it. The same ambient mechanism as the seven
on Sa, now sharper: **ρ = k + 1 with a locally unremarkable set.**

> If five colours behave the same way, ρ falls to **6** — an order of magnitude
> below the 63 a core of three needs. That is the experiment the strategy now
> turns on, and at five colours it costs one hard UNSAT: assuming every
> selector is exactly the 4-colourability instance.

## Corrections to my own claims, kept rather than edited away

- **"Pressure > 2 at k = 5 *requires* the confined set to be 4-chromatic."**
  Not proven in either direction. A 3-colourable confined set does not by
  itself colour the whole graph, and pressure 3 could in principle come from
  deeper structure. It is the *local* mechanism — what the k = 4 construction
  runs on — and the degeneracy result bounds that mechanism, not the pressure.
- **"The auxiliary graph being 4-chromatic is the condition."** Measured wrong:
  χ(A) = 4 already on de Grey's Sa and the pressure at five colours is still 2.
  Only auxiliaries seeing *both* circle colours share a list.
- **Auxiliaries taken as cross-hexagon sums only.** A point one away from `u`
  and `v` is `u + v` for *any* pair on the circle; an intermediate fix then
  wrongly required `|u − v| = 1`. Both dropped most of the set, including the
  always-confined √3 points.
- **A hexagon family with `theta` alongside `theta/2`.** `CH² = CT` exactly, so
  it produced a coincident hexagon — the same mistake that caused a false
  positive earlier — and the guard skipped it silently. It now says so.
- **"Every union in this package was dead before it was built."** The
  cross-pair bound needs `G − a` to be 4-colourable, which is the retracted
  criticality. Withdrawn — and monotonicity says unions can only lower ρ.
- **"de Grey's G is 5-vertex-critical."** **False**, and it was load-bearing:
  G − 1420 is still 5-chromatic (UNSAT, 1581 s). Separability was measured and
  read backwards — it is a consequence of criticality, not a proof. Every
  conclusion drawn from it is withdrawn; the measurements themselves stand and
  are now unexplained.
- **"ρ(Sa,4) = 7 because Sa contains a Moser spindle."** The seven vertices
  induce five edges and are 3-chromatic with three isolated points. The
  forcing is ambient, and the converse of the subgraph theorem is false —
  which is the one piece of good news in this section.
- **"Every confined degree is even."** True on de Grey's configuration and
  false in general — 1276 of 9600 sampled orientations have an odd degree, and
  degree 5 occurs. It is a symmetry of that offset family.
- **A spurious degeneracy 4 at "60°."** Two hexagon offsets agreeing modulo
  60° are the *same* six points; a 1e-9 tolerance on the points missed a
  near-coincidence at 1e-7. The coincident-hexagon trap for the third time,
  now in floating point. The guard compares offsets, generously.
- **"The degeneracy result closes the route."** It covers de Grey's angle
  family, which is one curve; for three or more hexagons the space is
  continuous.
- **"Blocking is what a 6-chromatic candidate needs."** A pendant edge in a
  blocking direction blocks the graph and changes χ by nothing. Blocking is a
  precondition for rigidity to be *measurable*, never a bound.

## What has been ruled out so far

Recording failures is the point of a search log; these are real constraints on where
a proof can come from, not just wasted cycles.

**Plain balls are 4-colourable.** Every set of points reachable from the origin by at
most *s* unit steps drawn from {ω^a σ^m}, ω = e^{iπ/3} and σ = e^{i·arccos(5/6)}, was
found 4-colourable, up to **30 715 vertices and 315 222 edges** (radius 4). Since a
graph containing a non-4-colourable subgraph is itself non-4-colourable, this is
conclusive for those sets: they contain no 5-chromatic subgraph at all.

**No forced monochromatic pair at the origin.** The spindling argument needs a pair
(p, q) that *every* 4-colouring paints alike. Searching from the origin over all
spindle-able distances — d² ∈ {1/3, 5/9, 1, 7/3, 3, 13/3, 7, 71/9, 31/3} — found
none, in balls up to **82 357 vertices and 931 158 edges**. Forcing is monotone
under taking supergraphs, so each of these rules out every subgraph too.

**Why all of that failed — the field was wrong.** A ball of 82 357 vertices came
back 4-colourable. If de Grey's 1581-vertex graph were inside it, the ball could not
be. So it is not inside, and no amount of extra radius or compute was ever going to
put it there: the *generating set* was wrong, not too small.

Reading the construction in [de Grey 2018](https://arxiv.org/abs/1804.02385) says
why. His rotations are 2·arcsin(1/4) and 2·arcsin(1/8) — the spindles at distance 2
and 4 — plus π/2 ± arcsin(1/8). Their sines are √15/8, 3√7/32 and 3√7/8. So the
construction does **not** live in ℚ(√3, √11); it needs **√7 and √15 = √3·√5** too.

The filter that kept rotations inside ℚ(√3, √11) — the ones indexed by Eisenstein
norms — excludes precisely the two angles de Grey used. `rotation_joining(4)` and
`rotation_joining(16)` raise `ValueError` against the small field, which was correct
behaviour serving a false assumption. `hn/geometry.py` now carries `DEGREY_FIELD =
Field((3, 5, 7, 11))`, dimension 16, in which all of his rotations are exact and a
point at distance 2 rotated by 2·arcsin(1/4) moves by exactly 1. The general
multiquadratic `Field` meant this was a parameter change rather than a rewrite.

**No forced disjunction either.** The weaker three-copy hypothesis was then tested
across the six most central pivots of each ball and every spindle-able distance,
using the single-query separation test. In every case a 4-colouring existed that
separates the pivot from *all* targets at that distance at once — so no disjunction
over any subset is forced, and neither spindle applies. Balls up to 25 675 vertices
and 253 842 edges, searched this way.

Limitations worth naming rather than burying: the σ-exponent stays at |m| ≤ 1 because
|m| ≤ 2 overruns the vertex budget at these radii, pivots are drawn from the centre
outwards rather than exhaustively, and the target distances are capped at d² ≤ 40.

## Method

**The field.** Every coordinate lives in ℚ(√3, √11), represented as four rationals
over the basis {1, √3, √11, √33}. √3 builds the triangular lattice (ω = (1+i√3)/2);
√11 builds the Moser hinge (σ = (5+i√11)/6, the rotation by arccos(5/6) that carries
two points at distance √3 to distance exactly 1). Both have modulus 1, so any
product ω^a σ^m is a unit vector and sums of them are the reachable points.

Which rotations are available is a question the field answers: a spindle at squared
distance d needs sin = √(4d−1)/(2d) to be representable. Over the integers that
allows d ∈ {1, 3, 7, 19, 25, 37, 61, 91, …} — the centred hexagonal numbers
3k²+3k+1 always work, with sine (2k+1)√3/(2d) — and over the rationals more besides.

**Scale.** `Fraction` arithmetic is exact but caps out near 10⁴ points. Everything
in one search shares a denominator, so multiplying through turns generation into
int64 numpy: ~80 000 vertices and ~930 000 exactly-confirmed edges in well under a
minute. An overflow guard refuses the integer path rather than returning a wrong
answer, and `tests/test_fast_agrees.py` pins the fast path to the slow one.

**Spindling, automated.** If every k-colouring of G paints p and q alike, and ρ is
the rotation about p taking q to distance 1 from itself, then G ∪ ρ(G) has *no*
k-colouring: both copies force colour(p) onto q and ρ(q), which are adjacent. This
is why searching balls for a 5-chromatic subgraph finds nothing while spindling one
of those same balls can succeed.

**A third copy weakens what has to be forced.** Demanding an outright forced pair is
a lot to ask. Ask instead that every k-colouring tie p to *q₁ or q₂* — a disjunction,
much weaker, and therefore much likelier to be found. Three rotated copies then
suffice by pigeonhole: each copy forces one of the two, so two copies force the same
one, and their images of it share p's colour while being adjacent.

That last step needs the rotated images to be pairwise at distance 1, and the
geometry permits it at exactly one radius. Two such points fit on any circle of
radius ≥ 1/2; **three fit only on the circumcircle of a unit equilateral triangle**,
radius 1/√3, with the rotations at 120° and 240°. So d² = 1/3 is special, and cores
of size 3 or more are unusable — a fourth point pairwise at distance 1 does not exist
on a circle.

**One query per family.** Give each target q a selector asserting p and q take
disjoint colours and assume them all at once. SAT means some colouring separates p
from every target, killing every disjunction over that set in a single query rather
than O(|Q|²). UNSAT means the disjunction is forced, and the solver's UNSAT core,
shrunk by re-solving and single deletions, is the minimal forcing subset — whose size
says which spindle applies.

The same test run on a unit triangle shows exactly why k=4 is hard: with 3 colours
the triangle exhausts the palette and the centre is forced to repeat one of the
corners; with 4 colours a spare colour remains and nothing is forced at all.

## Running it

```bash
pip install python-sat numpy pytest
python -m hn.cli demo                                  # rebuild and certify chi >= 4
python -m hn.cli verify certificates/moser_spindle_no3coloring.json \
    --drat-trim /path/to/drat-trim
pytest tests/ -q

python scripts/search_forced.py        # forced pairs, two-copy spindle
python scripts/search_disjunction.py   # forced disjunctions, three-copy spindle
#   HN_K=4 searches for chi >= 5; HN_K=5 is the open problem
```

For real verification, build the checker from
[marijnheule/drat-trim](https://github.com/marijnheule/drat-trim) (`gcc -O2 -o
drat-trim drat-trim.c`). Without it, certificates still report the SAT result, marked
explicitly as unverified.

## Layout

| file | what |
|---|---|
| `hn/field.py` | exact arithmetic in ℚ(√d₁,…,√d_k) |
| `hn/geometry.py` | points, rotations, the exact distance-1 predicate |
| `hn/graph.py` | unit-distance graphs, k-core and component reductions |
| `hn/generate.py` | vertex sets: lattice balls, reachable sets, rotation families |
| `hn/fast.py` | the int64 path — generation, complete edge finding, overflow guard |
| `hn/coloring.py` | SAT encoding, k-colourability, UNSAT-core minimisation |
| `hn/spindle.py` | forced pairs and disjunctions; two- and three-copy spindles |
| `hn/certify.py` | certificate creation and independent verification |

## A gradient, and what it measured

The search was blind for most of a day. Asking a SAT solver "is this pair forced?"
returns yes or no, so a sweep is a run of independent coin flips with no sense of
getting warmer — which is why 29 930 separation queries taught us nothing beyond
"not that one either".

The solver was measuring the missing quantity all along. A pair separated with zero
conflicts is wide open; one costing thousands is nearly forced, since almost every
colouring ties those vertices together and the solver had to work to find an
exception. `SeparationDifficulty` reads it off.

**Controlled against size.** A larger formula can cost more conflicts for no reason
but its size. Adding a translate 100 units away doubles the vertex count and adds no
constraint between the copies: score 88 → 129, a factor of 1.5. A genuine tightening
of the same size reaches 4216, a factor of 48. The score tracks constraint.

**Calibrated at k=4**, where the answer is known. Nested balls of de Grey's graph:

| radius | vertices | score |
|---|---|---|
| 1.0 | 48 | 0 |
| 1.5 | 210 | 4 |
| 2.0 | 534 | 19 |
| 2.5 | 915 | **2932** |
| 3.0 | 1302 | **forced** |

It does not saturate. It sits flat, explodes by two orders of magnitude, and crosses
on the next step — at d² = 1/3, the radius the geometry and the gradient had each
singled out independently.

**Localisation was wrong.** Tightening the whole graph doubles the vertex count every
round, so the plan was to tighten only a ball around the hardest pivot — forcing being
a local question. Run side by side, the tight radius saturates while the global one
keeps climbing:

| round | radius 1.6 | radius 3.0 |
|---|---|---|
| 1 | 88 | 67 |
| 2 | 262 | 2110 |
| 3 | 2604 | 22 697 |
| 4 | 2882 | — |

The pair under test is local; what hardens it is not. So the compute wall that
localisation was meant to dodge is still there, and the dodge does not work. A blind
sweep could not have told these two runs apart at all.

## The upper bound, and why searching it is closed

Everyone repeats that lowering the upper bound needs an infinite object and so
cannot be searched. That is false for *periodic* colourings — Isbell's 7-colouring
is periodic, one hexagon repeated — and a periodic colouring is determined by a
finite fundamental domain.

So: tile by a lattice, cut the fundamental domain into cells, and forbid two cells
sharing a colour whenever a unit distance between them is *possible*, over every
lattice translation. Closed cells make that conservative, so SAT yields a genuine
colouring while UNSAT rules out only that lattice at that resolution.

It calibrates. A first attempt with cells of diameter 0.35 came back UNSAT even at
k=7, which the encoding must satisfy — too coarse, since cells of diameter d forbid
the whole band [1−d, 1+d] rather than the circle. Measuring the slack in the known
colouring fixed the resolution: with hexagons of circumradius 0.45 the same-colour
distances avoid [0.8953, 1.2256], so cells under 0.105 fit. On the index-7
sublattice with 28 cells a side, 784 cells and 78 792 constrained pairs, a
7-colouring is found in 0.7 s.

**And it cannot reach six, by a theorem rather than by a compute limit.** Cell-based
colourings are map-type colourings with polygonal regions, and those need at least
seven colours: Woodall (1973) and Townsend (1981) give six,
[arXiv:2502.01958](https://arxiv.org/abs/2502.01958) raises it to seven for maps
whose boundaries are not arcs of unit circles, with arbitrary polygons as a
corollary. Any discretisation into cells lands squarely inside that.

Which is worth knowing precisely because of the exception the theorem is careful to
state. Boundaries that *are* arcs of unit circles are not excluded. So a 6-colouring
of the plane, if one exists, cannot have polygonal colour classes — it needs curved
boundaries of radius exactly 1, or non-measurable ones. No grid will ever find it.

Recorded here so the next person does not spend the afternoon on it, as this did.

## What the literature already had

Read after the fact, which was the wrong order:

- **Forced monochromatic pairs are standard.** Parts calls them *mono-pairs* and
  builds 5-chromatic graphs from "cycles connecting two or more mono-pairs and one
  unit edge", closed by a pigeonhole argument. The core of the spindle lemma here
  is established work. What is not written down in that form is the conflict-graph
  and independent-transversal statement, which systematises an idea that already
  existed rather than introducing one.
- **The order bound is classical.** (ℤ/n)* has exponent 2 exactly when n | 24.
- **No consensus on the answer.** Polymath16 participants split between 5 and 6.
- **Scale.** Parts derived N = 6906 as a bound tied to 6-colourability, which says
  something about the size of the object being hunted.

## An unexplored field

One lead from that thread has not been followed. Philip Gibbs notes that 5/3 is
monochromatic with the origin in any homomorphic 5-colouring, since 5/3 = η² + (−η)²
in ℤ[ω₁, ω₃]. Spindling that point needs |5/3|² = 25/9, whose rotation
(34.9152°, cos 41/50) has sine 3√91/50 — so it lives in **ℚ(√7, √13)**, a field that
appears nowhere in the sources read here. Everything in the literature sits in
ℚ(√3, √11), with √5 and √7 added for de Grey's rotations.

The caveat that keeps this a lead rather than a result: Gibbs' observation concerns
*homomorphic* colourings, a restricted class, and forcing there does not imply
forcing in general.

## Slack: why the difficulty jumps where it does

A unit-distance graph in the plane has clique number **3**. Four points pairwise at
distance 1 would be a regular tetrahedron, which does not fit in two dimensions.

Define the **slack** of a colouring problem as s = k − 3, and let **f(k)** be the
number of vertices in the smallest unit-distance configuration containing a pair
forced monochromatic under k colours.

At **s = 0** a triangle exhausts the palette, so forcing is *local*: two triangles
sharing an edge already do it. f(3) = 4, which meets the trivial bound f(k) ≥ k + 1
(with n ≤ k, colour everything differently) **exactly**.

At **s ≥ 1** no local configuration exhausts anything — a spare colour always
remains — so forcing has to be assembled combinatorially across many vertices, and
f jumps. The same four vertices force nothing at k = 4, and neither does the
seven-vertex Moser spindle.

| k | slack | f(k) | settled |
|---|---|---|---|
| 3 | 0 | **4** (tight) | trivial |
| 4 | 1 | measured by `scripts/measure_fk.py` | 1961 |
| 5 | 2 | > 1581 here | 2018, 57 years later |
| 6 | 3 | — | open |

**Measured, and then corrected.** Minimising from a ball known to carry forcing
gives a 359-vertex configuration in which the pivot is forced monochromatic with
*something* at d² = 1/3. The first reading of that was that spindling it would give
a 5-chromatic unit-distance graph on ≤ 717 vertices, which is wrong: the forced
object is a **disjunction over 34 targets**, not a pair. Spindling a disjunction of
size r needs r+1 rotated images pairwise at distance 1, and a circle holds at most
three. So that configuration is not spindle-able at all, and no bound on the
5-chromatic record follows from it.

**What that correction reveals is more useful than the claim was.** Forcing at k = 4
is not absent — it is the wrong *shape*. It appears as a wide disjunction, and the
spindle needs a narrow one: core size 1, or size 2 at exactly d² = 1/3. That gap
between "forcing exists" and "forcing is usable" is the real obstruction, and it
explains why de Grey's and Parts' constructions are elaborate rather than a single
spindle of some forced pair.

**Why the core stays wide: symmetry.** On the measured configuration the forced
core is 34 of 36 targets, and all 36 form a closed orbit under the 60° rotation
about the pivot — every one maps to another. A colouring argument cannot single
out targets that a symmetry permutes, so symmetry sets a floor on how narrow a
core can get.

That rotation is not a full automorphism: 280 of 359 vertices land back in the
graph, not all of them. And that partial asymmetry is exactly what bought the two
exclusions taking 36 down to 34. **Asymmetry buys exclusions; symmetry blocks
them.** Getting to 2 needs 32 more.

Which indicts the tightening used here for most of a day. Unioning a graph with
copies of itself rotated about the *same* pivot **raises** symmetry, pushing the
core the wrong way — the effort score climbed while the quantity that decides
anything could not move. de Grey's construction is asymmetric on purpose, two
different rotations about an off-centre pivot, and this is what that is for.

**Asymmetric tightening narrows it.** Acting on that, tightening about pivots
*other* than the forcing one — breaking its symmetry instead of reinforcing it —
moves the number that symmetric tightening never touched:

| round | vertices | forced core |
|---|---|---|
| start | 359 | **34** |
| 1 | 634 | **22** |
| 2 | 987 | **13** |

Twelve exclusions, then nine, where a day of symmetric tightening held it at 34
while the effort score climbed through 88, 2110, 22 697.

**And it saturates at 11.** Carried on, the descent goes 34, 29, 24, 21, 17, 14,
12, 11 and then stops dead. The first suspicion was exhausted candidates — each
round took the eight highest-degree vertices as pivots, which after seven rounds
are the same eight — so the search was rerun with **32 candidates sampled across
the whole degree range**. Two further rounds returned 11 from all 32.

That refutes the exhaustion explanation and settles the number: **11 is the floor
of this construction, not an artefact of where it was looking.** The mechanism is
real — two thirds of the way down from 34, where symmetric tightening moved nothing
at all — and its reach stops well above the 2 a spindle needs.

**Why 11, structurally.** Two points on the circle of radius 1/√3 are at distance
1 exactly when they are 120° apart. The eleven surviving targets sit at

```
 38.55  98.55  158.55  218.55  278.55  338.55     six, step 60°
    69.6  129.6  189.6  249.6  309.6              five, step 60°
```

— two interleaved six-fold orbits, 31.05° apart. Within an orbit there are
120° pairs; **between the orbits there are none**, since 31.05 + 60k is never 120.
The two families place no constraint on each other, so an independent choice always
exists and no number of rotated copies makes the spindle lemma fire. The floor is
the orbit structure, not the copy count.

Finding that also turned up a restriction the lemma never needed. It says every
chosen image carries the pivot's colour, so *any* adjacent pair among them is a
contradiction — including one arising from two **different** targets. The
implementation compared only images of the same target. Lifting that adds 30
cross-target conflict edges against 33 same-target ones here, so it is a real
strengthening; it still does not block, for the orbit reason above.

**The mechanism is real; the metric for it is not.** The obvious measure — what
fraction of vertices the rotation about the forcing pivot still preserves — does
not predict which tightening works, and is anti-correlated in round 2: the winning
candidate (core 13) was the *most* symmetric of the five by that measure, and the
least symmetric stayed at 22. Particular positions win repeatedly; an aggregate
quantity does not explain why. Recorded as a mechanism without a metric rather
than dressed in an explanation that fits.

Scope: this is at k = 4, which is settled territory. It is a proof of mechanism,
not the prize. Applying it at k = 5 needs forcing to exist there first, and on de
Grey's graph 29 930 queries found none — there is no core to narrow until one
appears.

It also replaces the search objective. Chasing effort asks "how hard is this pair to
separate"; the quantity that actually has to move is **the size of the minimal
forced core**, which must reach 2. Here it sits at 34.

The gap f(k) − (k+1) is what the plane's missing K₄ costs, and the history of the
problem is that gap widening. The anchors are pinned in `tests/test_slack.py`; the
clique bound is elementary and f(3) = 4 is immediate, so what is worth anything
here is f(4) as an actual number, which decides whether a search reaching ten
thousand vertices is close to s = 2 or missing orders of magnitude.

Forcing and chromatic number are not the same thing — the Moser spindle is
4-chromatic while forcing nothing at k = 4 — and the framing is about f, not χ.

## Every generator in the package symmetrises

A minimal forced core has to reach 2. Asymmetric tightening floors it at 11, and
the reason turned out to be upstream of the tightening: **the search space was
symmetric by construction.** `hex_ball` is the six-fold triangular lattice by
definition. `walk_ball` steps along the sixth roots of unity. de Grey's own `Sa`
is, in his words, "all points obtained by rotating S about the origin by multiples
of 60 degrees and/or negating their y-coordinates" — a twelve-element group
applied deliberately.

That is not a flaw in those constructions; symmetry is what makes them tractable
to state and to verify. But it means every configuration examined here was a union
of orbits, and the 11 surviving core targets confirm it: they sit at 38.55 + 60k
and 69.6 + 60k, two interleaved six-fold orbits 31.05 degrees apart, with no
cross-orbit pair at 120 degrees. A core of 1 or 2 cannot be a union of six-fold
orbits. The object was excluded by the generator, not by the mathematics.

`scripts/asym_grow.py` drops the group: it accretes points at distance 1 from a
Moser-spindle seed and never symmetrises. The first version attached at random and
traded one structural defect for another — 699 vertices with 1070 edges, average
degree 3.06 against de Grey's 9.97, and a sparse graph forces nothing because it
colours with room to spare. Sampling a pool of candidate attachments and keeping
whichever lands adjacent to the most existing points fixes that: 499 vertices,
2453 edges, average degree 9.83, matched to de Grey's graph while staying
asymmetric.

Still no forcing at k = 5. But the comparison is finally a fair one — the earlier
negatives were measuring graphs too loose to force anything, not graphs too
symmetric.

## The ceiling was the field, and the number is 24

The narrowing programme needed a forced core of 2 and stalled at 11. It was
not the search that failed.

Targets at one distance from the pivot put every rotated image on a single
circle, and a point of a circle lies at unit distance from at most two of its
points, so the unit-distance graph on the distinct images has maximum degree 2
-- a disjoint union of paths and cycles. Bipartite pieces never trap, since one
side is an independent set meeting every copy, so blocking lives entirely in
the odd cycles, and an odd cycle `C_q` traps at most `(q+1)/2` targets.

Which odd cycles exist is a question about the field. A rotation of order `n`
needs `zeta_n` inside `K(i)`; for multiquadratic `K` that Galois group is an
elementary abelian 2-group, which forces `(Z/n)*` to be one too, and that
happens exactly when **n divides 24**. The odd divisors of `24 = 2^3 * 3` are 1
and 3. Over any multiquadratic field, on any circle, with any number of copies,
**the only odd cycle available is the triangle** -- and a triangle traps 2.

de Grey's field is `Q(sqrt3, sqrt5, sqrt7, sqrt11)`. Every configuration
examined here lives in a multiquadratic field. The ceiling arrived with the
arithmetic, before any searching started.

One measurement corrected the first version of this. Conflict-graph degrees
grow 4, 10, 22, 46, 94 as copies are added and look like an escape; that is
multiplicity. At squared distance 1/3, 384 images collapse onto 12 distinct
points, and among distinct points the degree is 2 at every copy count tried.
`tests/test_transversal.py` pins it.

### Niven closes it, elementarily

There is a second bound, and it is the one that actually bound every search
here. Trapping needs the images to close into a cycle, so the angle must be
commensurable with `2 pi`; and `cos t = 1 - 1/(2 d^2)` is rational exactly when
`d^2` is. **Niven's theorem** allows only `0, ±1/2, ±1` as rational cosines of
rational multiples of `pi`, so a rational squared distance admits a
finite-order rotation at just four values:

| `d^2` | angle | cycle | traps |
|------:|------:|------:|------:|
| 1     | 60°   | 6  | 0 |
| 1/2   | 90°   | 4  | 0 |
| 1/3   | 120°  | 3  | **2** |
| 1/4   | 180°  | 2  | 0 |

At every other rational `d^2` the images form a path -- bipartite -- and
nothing is trapped at all, over any field, with any number of copies.

Twelve call sites in this package skip a target whose squared distance is
irrational: four in `hn/spindle.py` and eight across `scripts/`. **Every search
run here was inside a space capped at 2 before it started**, and no choice of
field would have rescued it. The escape needs targets at an *irrational*
squared distance, `1 / (4 sin^2(pi t/n))` for an `n` with an odd divisor of at
least 5 -- which is also, independently, where the Galois bound above stops
applying.

### A hierarchy of magic circles

Closing the images into a cycle at all pins the radius. Adjacency joins images
`beta = 2 arcsin(1/2d)` apart, and that is `t` rotation steps only when
`cos beta = 1 - 1/(2d^2) = cos(2 pi t/n)`, so

    d = 1 / (2 sin(pi t / n)).

Each rotation order owns a radius, and targets anywhere else see a path and
cannot be trapped at all.

| order | radius | traps |
|------:|-------:|------:|
| 3  | 0.5774 = 1/sqrt3 | 2 |
| 5  | 0.8507 | 3 |
| 7  | 1.1524 | **4** |
| 9  | 1.4619 | 4 |
| 11 | 1.7747 | **5** |
| 13 | 2.0900 | 5 |
| 15 | 2.4049 | 4 |

`n = 3` is the circle where two points are adjacent exactly when they are 120
degrees apart -- the one every spindle argument in the literature is built on,
and the only one `n | 24` permits.

The capacities are computed, not read off the independence number. An earlier
version of this table claimed `(n+1)/2` -- 8 at order 15, 11 at order 21 -- and
that is wrong from length 9 upwards. The copies are the `n` rotations of the
cycle, so their image sets are the `n` shifts of the target set `T`, and a
colouring escapes when some independent set `S` of `C_n` meets every shift,
which happens exactly when `S - T = Z_n`. Blocking means no independent `S`
covers, and that is stricter than `S` merely being large. **Capacity does not
grow with the cycle: it peaks at 5 and comes back down**, because a longer
cycle also has larger independent sets and covering gets easier faster than
trapping does. Past length 15 the enumeration stops being affordable, and
`trapping_bound` raises rather than guessing.

`hn/cyclotomic.py` builds points where the arithmetic is different: `Z[zeta_n]`,
where `|z| = 1` forces `z conj(z) = 1` exactly, so Kronecker makes `z` a root of
unity and the unit steps are precisely the `2n` elements `±zeta_n^k`. `n = 3` is
the Eisenstein lattice everything here was built from: 6 steps, capacity 2.

It is not yet a construction. `Z[zeta_15]` gives 30 unit steps and keeps unit
triangles, so the plane's clique number survives, and a three-step walk is
3901 vertices at average degree 11.4 -- denser than de Grey's graph. It is also
**3-chromatic**, because `sqrt11` is not in `Q(zeta_15)` and so no Moser spindle
lives there. Richer arithmetic on its own buys nothing; what the theory asks for
is the compositum, a field holding both the spindles and an odd-order rotation.

## Ramification decides which magic circles exist

A unit step satisfies `u conj(u) = 1`, so at any prime `P` that complex
conjugation fixes, `2 v_P(u) = 0` and `v_P(u) = 0`. Every point reachable from
the origin by unit steps then has `v_P >= 0`, and so does every squared
distance inside one connected component.

The order-`n` magic radius squared is `1/((1 - zeta_n)(1 - zeta_n^-1))`. For
`n` a prime power that is the ramified prime above `p`, so the radius has
valuation `-2` and lies outside reach. Measured, and it is what every
disconnection all day was saying: the pivot's component held 1306 of 5017
vertices and **none** of the 11-gon, and 36481 points reached nothing one step
from a pentagon vertex.

The way out is for the prime to split with conjugation swapping the factors --
and that is exactly what the classical construction is. 3 splits in
`Q(sqrt(-11))` since `-11 = 1 mod 3` is a residue, `(1 + sqrt(-11))/2` has norm
3, and its quotient by its conjugate is `(-5 + sqrt(-11))/6`: **the Moser
rotation, up to sign**. `|1 - rho|^2 = 1/3` exactly, so `1 - rho` sits on the
classical magic circle two steps from the origin. The spindle's angle was never
chosen -- it is the split-prime quotient at 3.

`15` is not a prime power, so `1 - zeta_15` is a unit, the radius is an
algebraic integer, and nothing forbids it. Counting points on each circle in
one ball, at depth 3:

| order | prime power? | points found | capacity |
|------:|:------------:|-------------:|---------:|
| 3  | yes, ramified | 4  | 2 |
| 5  | yes, ramified | **0** | 3 |
| 15 | no, unit      | **30** | **4** |

### How close it gets, and where it stops

Over `Q(zeta_15)(sqrt(-11))` -- degree 16, the smallest field with `omega` for
the triangles, `zeta_5` for the rotations, `sqrt5` for the radius and
`sqrt(-11)` for the spindle -- the order-15 circle carries two full orbits of
the order-15 rotation, and at `k = 3` the pivot's colour **is** forced onto an
orbit. Off-centre balls take the minimal core from the full orbit of 15 down to
**5**. The capacity is 4.

Five is not four, and the shortfall is structural rather than a matter of
searching harder. The core comes out as five *consecutive* positions of the
15-cycle, and an arc of five never blocks: `S = {0, 5, 10}` is independent and
meets every 5-arc. Of the 1365 four-subsets of `C_15` exactly **45 block**, in
three classes up to rotation -- gaps `(1,2,10,2)`, `(2,3,4,6)`, `(2,6,4,3)` --
and none of the 45 is forced on either orbit.

### One order, a family of circles -- and the same answer on all of them

Adjacency on the circle need not be one step of the rotation. At `t` steps the
radius is `1/(2 sin(pi t/n))` and the cycle is `C_{n/gcd(n,t)}` taken in the
order `0, t, 2t, ...`, so the *same angular arc* lands on a different subset of
the cycle. At `n = 15` the family is:

| t | radius | cycle | capacity | points found |
|--:|-------:|:-----:|---------:|-------------:|
| 1 | 2.4049 | C₁₅ | 4 | 30 |
| 2 | 1.2293 | C₁₅ | 4 | 30 |
| 3 | 0.8507 | C₅  | 3 | **0** — the order-5 radius, ramified |
| 4 | 0.6728 | C₁₅ | 4 | 30 |
| 5 | 0.5774 | C₃  | 2 | 4 — this is `1/sqrt3`, the classical circle |
| 7 | 0.5028 | C₁₅ | 4 | 30 |

The theory falls out of the table twice over: `t = 5` reproduces the classical
`1/sqrt3` circle with its capacity of 2, and `t = 3` is the order-5 radius with
**no points at all**, exactly as ramification predicts.

The four capacity-4 circles have genuinely different radii and different point
sets, and **all four give the same minimal core: the angular arc `[0,1,2,3,4]`**.
Relabelling by the cycle turns that arc into `[0,1,2,8,9]`, `[0,1,4,8,12]` and
`[0,7,9,11,13]` — four different shapes, none of them forced down to four.

So the forcing is *angularly local*: the constraint reaches the pivot from its
own neighbourhood, which sees a bounded range of directions, and a blocking
subset needs its targets spread around the circle. That is a shape mismatch, not
a size one, and it does not look like something a bigger ball fixes. Bigger, in
fact, is worse: raising the cap to 120000 vertices left 20 points on the circle
and **no complete orbit**, because a truncated breadth-first frontier cuts
orbits in half.

At `k = 4` there is no forcing at all in this field: 55363 vertices at average
degree 12.6, separable. Forcing at four colours needs a construction of de
Grey's kind, not a ball.

## A whole avenue, closed with a reason

Every search here chases one binary fact, and until it turns up there is
nothing to show. There is a continuous quantity giving the same conclusion,
and every finite graph reports a value for it.

Let `S` be measurable, avoiding distance 1, of upper density `d`, and `G` a
finite unit-distance graph on `n` vertices. Average over rigid motions: the
expected size of `sigma(V) ∩ S` is `n d`, and that intersection is independent
in a copy of `G`, so it never exceeds `alpha(G)`. Hence

    m_1(R^2) <= alpha(G) / n

for every finite unit-distance graph, and five measurable classes covering the
plane force one of density at least `1/5`. So a single `G` with
`alpha/n < 1/5` would give `chi_m(R^2) >= 6`, from one independent-set
computation.

**It cannot happen.** Croft's 1967 construction is a measurable 1-avoiding set
of density about `0.2293`, and averaging runs the other way too: that forces
`alpha(G)/n >= 0.2293` for *every* finite unit-distance graph. Since
`0.2293 > 1/5`, no graph ever brings the ratio under `0.2`, and the most this
argument can yield is `1/0.2293 = 4.36` — `chi_m >= 5`, which is already known.

### Every fractional method is blind to this problem

The same set closes more than one avenue. Translate Croft's independent set of
density `m_1` over the isometry group: the copies cover every point equally, so
weighting them is a *fractional colouring* of total weight `1/m_1`. Hence

    chi_f(R^2) <= 1 / m_1 <= 1 / 0.2293 = 4.36 < 5.

The fractional chromatic number of the plane is **below five**. Every LP and
SDP relaxation is bounded by `chi_f`, so no relaxation can see `chi >= 5`, let
alone `chi >= 6`. That is why de Grey's 2018 result had to be an integral
combinatorial argument, and why the LP line had stalled for decades without it
being a shortage of computation. **The only route to six is integral.**

### And what that says about the object being searched for

Turned around, the same bound is a fact about the target rather than about the
method. `alpha(G)/n < 1/5` would give `chi(G) >= 6` outright, with no measure
theory at all: a 5-colouring splits `n` vertices into five independent sets and
one of them has at least `n/5`. Croft's bound says no unit-distance graph ever
qualifies.

So **a 6-chromatic unit-distance graph has independent sets of at least
`0.2293 n`**, and five of them cover at least `1.1465 n` -- more than the whole
graph, with room to spare. Whatever stops such a graph from 5-colouring, it is
never counting. It has to be structure, and the classes have to fail to fit
together for reasons no cardinality argument can see.

Recorded rather than dropped, because knowing why an avenue closes is worth as
much as a search that fails quietly inside it. The ceiling here is a
construction, not a shortage of computation, and `tests/test_density.py` pins
it. The Moser spindle's `2/7 = 0.2857` and the published `m_1 <= 0.2470`
(Ambrus, Csiszárik, Matolcsi, Varga, Zsámboki 2023, by Fourier methods rather
than from a graph) sit between the two.

## The gate: every 6-chromatic unit-distance graph is blocked

The blocking machinery was, until now, a thing running alongside the chromatic
number rather than under it. One line puts it underneath.

A **coset colouring** is a homomorphism `phi` from the edge module `M` to
`Z/n` that is nonzero on every edge vector. The vertices of a connected
unit-distance graph all lie in `v0 + M`, so colouring `v` by `phi(v - v0)` is
well defined, and adjacent `u, v` get `phi(u - v) != 0` — different colours. So
a coset colouring with `n` colours *is* an `n`-colouring, and `chi <= n`.
Contrapositive, at `n = 5`:

> **Every 6-chromatic unit-distance graph is blocked.**

Blocking is therefore the **first necessary condition** for `chi >= 6` — a gate
every candidate has to pass before any colouring argument begins. `phi` is
propagated along a spanning tree of `Sa` in `tests/test_homcol.py` and checked
on all 3938 edges, so the gate is walked, not asserted.

And the only 5-chromatic graph anyone has **fails it**. de Grey's `G` lives
over `Q(sqrt3, sqrt5, sqrt7, sqrt11)`, and a multiquadratic field provably
never blocks. Not "has not been shown to block": cannot.

### The first 5-chromatic blocked graph

Let `w` run over the unit steps of `Q(zeta_7)` of denominator exactly 29 — 29
being the least rational prime that splits completely there, which is what
makes those steps the cheapest blocking ones — and set

    U = G u {w.G},    each w.G the rotation of G about the origin by w.

`chi(U) >= 5` because `U` contains `G`. For the blocking, two levers:

* it is **monotone** in the direction set, so a blocking *subset* suffices;
* it is **invariant under a global rotation**, because `phi -> phi . u` is a
  bijection between the homomorphisms out of `uM` and those out of `M`
  preserving which vectors go to zero.

`G` contains a rotated copy `r.Y` of `Y`, and `Y` keeps all six hexagonal unit
edges at the origin — de Grey deletes only `(1/3, 0)` and `(-1/3, 0)`, and the
origin is measured here to have 60 neighbours, `(+-1, 0)` and
`(+-1/2, +-sqrt3/2)` among them. So `U`'s directions contain
`r.{w . zeta_6^k}`, which blocks exactly when `{w . zeta_6^k}` does — and that
set, 50 generators times the six sixth roots, **blocks**:

    300 directions, all projectively distinct, no coset 5-colouring

So `U` is **5-chromatic and blocked**: the first graph that is both, and the
first to pass the gate. Stated for what it is — `U` is not claimed to be
6-chromatic, and nothing above says it is.

It is also, as it turns out, the only object here that meets **every**
necessary condition known for a sixth colour. Its directions contain
`r·{w·zeta_6^k}`, which is the denominator-29 set up to a global rotation, and
that set blocks at `n = 2, 3, 4` **and** `5` while surviving the Cayley screen
at every modulus from 5 to 20. Rotation preserves blocking at every `n`, so
`U` inherits all of it.

### Every triangle of unit steps is a 60-degree pair

`|u| = |v| = |u - v| = 1` expands to `u.vbar + ubar.v = 1`; with `t = u.vbar`
that is `|t| = 1` and `t + tbar = 1`, so `t` is a primitive sixth root of
unity. Hence `u = zeta_6^{+-1} v` for **every** triangle of unit steps.

Two consequences. A `zeta_6`-closed step set carries triangles whatever its
denominators — and `zeta_6` is a unit at 29, so the denominator-29 steps are
closed under it and support 300 such pairs. The blocking directions are not a
rigidity-free fringe. And the rhombus is pinned: its tip is `u + zeta_6 u`,
with `|u + zeta_6 u|^2 = 2 + 1 = 3` — always `sqrt3`, whatever `u` is.

### Why the blocking keeps dying: the spindle is unique

A Moser spindle is two rhombi sharing an apex, tips at `sqrt3` by the identity
above, joined by a unit edge. If the arms are the rhombi on `u` and on `v`,
the closing edge has length `|u - v| sqrt3`, so

    |u - v|^2 = 1/3,   u.vbar + ubar.v = 5/3,   v = rho^{+-1} u,
    rho = (5 +- sqrt-11)/6.

**No freedom at all.** The arm rotation is a single algebraic number and the
spindle is unique up to isometry. With rotation invariance, one measurement
settles the whole class — and the measurement says the spindle admits a coset
5-colouring. Its 14 projective directions are not enough. So

> **No 7-vertex 4-critical unit-distance graph blocks.**

Which is exactly why the blocking keeps evaporating. Hang 144 rotated copies
of the spindle together and the union is 367 points, 671 edges, 427
directions, `chi = 4`, **blocked** — every edge a spindle edge, not a pendant
in sight. Then strip it to a 4-critical subgraph and what is left is

    7 points, 11 edges, 7 directions — one spindle, and not blocked.

The copies are each 4-chromatic alone, so the core keeps one and discards the
rest. That is the obstruction, and it sharpens the whole programme to a single
question:

> **Can a critical unit-distance graph block?**

A lower bound on what that would take: blocking at `n = 5` is covering
`PG(r-1,5)` by the hyperplanes `d-perp`, and the cheapest cover of a
projective space over `F_q` is a pencil of `q+1` hyperplanes through a fixed
codimension-2 subspace. So a blocked graph has **at least six edge
directions**, and six suffice only when their reductions mod 5 span a rank-2
space and occupy all six points of that projective line. Weak, but the right
shape: blocking is about how the directions sit mod 5, not how many there are.
The spindle has 14 and fails; the denominator-29 set has 300 and succeeds.

## Blocking is decided by the field — a necessary condition on `chi >= 6`

Every 6-chromatic unit-distance graph is blocked. So it is worth knowing which
fields *can* block, and the answer turns out to be arithmetic, decidable in
advance, and to rule out nearly everything anyone has built in.

Reduce mod 5. Let `K` be the field generated by the edge vectors, `sigma` its
complex conjugation, `F = K n R` the real subfield. When the edge module `M`
is 5-maximal, a coset 5-colouring is an `F_5`-linear functional on

    A = O/5 = prod_i F_{5^{f_i}},      sum f_i = [K:Q]

nonzero on every edge vector. Edge vectors are unit steps, `u.sigma(u) = 1`, so
they reduce into

    N = { x in A : x.sigma(x) = 1 },

and Hilbert 90 says every element of `N` lifts — `N` is exactly the set of
reductions available. Scaling the graph by `lambda` is a linear automorphism of
`A` permuting hyperplanes, so the question is scale-free: **a graph blocks only
if `N` meets every hyperplane of `A`.**

**The orbit-wise reduction.** `sigma` permutes the factors of `A`, and
`x.sigma(x) = 1` is componentwise, so `N` factors along the `sigma`-*orbits*:
`N = prod N_i`. A functional supported on one orbit and zero elsewhere takes
exactly that orbit's value. Hence

> if **one** orbit admits a functional nonzero on all of its `N_i`, the whole
> field admits one, and **no graph over `K` blocks** — whatever its size.

So blocking needs *every* orbit to block alone, and an orbit is one number: the
residue degree `f`, with `sigma` either fixing the prime (`f` even, `sigma` the
involution `Frob^{f/2}`) or swapping a conjugate pair. Every type of degree 2,
4, 6 and 8 was decided by exhaustion over `A`:

| orbit | `\|N\|` | hyperplanes missing | blocks |
|---|---|---|---|
| fixed `f=2` | 6 | 3 | no |
| fixed `f=4` | 26 | 13 | no |
| fixed `f=6` | 126 | 0 | **yes** |
| fixed `f=8` | 626 | 0 | **yes** |
| swapped `f=1` | 4 | 4 | no |
| swapped `f=2` | 24 | 12 | no |
| swapped `f=3` | 124 | 0 | **yes** |
| swapped `f=4` | 624 | 0 | **yes** |

A `sigma`-fixed prime of degree `f` lies over a prime of `F` of degree `f/2`; a
swapped pair lies over one of degree `f`. The four that block are exactly the
four of **residue degree at least 3 over `F`**. So:

> **Theorem.** Let `Gamma` be a unit-distance graph whose edge module is
> 5-maximal and whose edge vectors are integral at 5. If **either** `K/F`
> ramifies at some prime above 5, **or** some prime of `F = K n R` above 5 has
> residue degree at most 2, then `Gamma` has a coset 5-colouring, and
> `chi(Gamma) <= 5`.
>
> **Corollary.** Every 6-chromatic unit-distance graph has `K/F` unramified
> above 5 **and** residue degree at least 3 at every prime of its real
> subfield above 5 — every residue field of `F` at 5 has at least **125
> elements**.

The integrality hypothesis is automatic when every prime above 5 is
`sigma`-fixed: then `2 v_p(u) = v_p(u.sigma(u)) = 0` for every unit step.

**The ramified half has its own one-line proof.** Suppose `K/F` ramifies at a
prime `p` above 5. Then the residue extension is trivial, so `sigma` — which
generates `Gal(K/F)` — acts trivially on `O/p`, and a unit step's
`x.sigma(x) = 1` reduces to `xbar^2 = 1`. So `N mod p` is contained in
`{+1, -1}`: **two elements**. Compose `A -> O/p` with any `F_5`-functional
`psi` having `psi(1) != 0` and every element of `N` lands on `+-psi(1)`, never
on `0`. A coset colouring exists, so no graph over such a `K` blocks.

Ramification of 5 in `F/Q` by itself changes nothing — only the residue degree
does. Checked by exhaustion in three cases where `K/F` is ramified above 5:
`Q(zeta_5)` (11 of 156 hyperplanes miss `N`), `Q(zeta_15)` (3 of 97656) and
`Q(zeta_20)` (4 of 97656). Each was predicted unable to block before it was
computed.

### What it rules out

It **subsumes the multiquadratic result**. A multiquadratic field has
elementary abelian Galois group, so every decomposition group is cyclic of
order at most 2, so `f <= 2` and the residue degree over `F` is 1 — below the
bound, every time. de Grey's `G` was never going to block, and now that is a
corollary rather than a separate computation.

It closes the Moser spindle's own field for good. `Q(zeta_3, sqrt-11)` has
degree 4, and **every** degree-2 and degree-4 type fails, so no unit-distance
graph over it can ever block — at any size, by any construction. The same for
`Q(zeta_3, sqrt-m)` at `m = 1, 2, 5, 6, 7, 10, 11, 13, 15`, each measured
separately and each failing.

And it agrees with every measurement here. `Q(zeta_7)`, `Q(zeta_9)` and
`Q(zeta_21)` all have `ord(5) = 6` with `-1` in the group it generates, so one
`sigma`-fixed prime of degree 6, residue degree 3 over `F` — predicted able to
block, and the 300 denominator-29 directions do. `Q(zeta_7)` is the **smallest
cyclotomic field that can block**: 3 and 4 give degree 2, 5 ramifies, 7 is
next. That is the same 7 that put 29 — the least prime splitting completely in
`Q(zeta_7)` — at the bottom of the blocking denominators, reached here from the
opposite direction entirely.

### A verdict of my own that was wrong

`cyclotomic_can_block(n)` used to return "cannot block" for **every** `n`
divisible by 5, on the grounds that 5 ramifies there. That confuses
ramification of 5 in `K/Q` with ramification of `K/F`, and they are different.
Writing `n = 5^a m` with `5` not dividing `m`, the inertia group at 5 is
`(Z/5^a)* x {1}` and conjugation is `(-1, -1)`, so

    K/F ramifies above 5  <=>  conjugation lies in inertia  <=>  m | 2,

that is, only for `n = 5^a` and `n = 2.5^a`. `Q(zeta_35)` has
`f = ord(5 mod 7) = 6` with the primes `sigma`-fixed, hence residue degree
**3** over `F` — and the blanket verdict called it unable to block. Corrected,
and the two ramified cases settled by exhaustion agree with the new logic:
`Q(zeta_15)` and `Q(zeta_20)` both have residue degree 1 over `F`.

One caveat is kept explicit. When `5 | n` the local ring is `O_p/p^e` rather
than a field, so the orbit classification — which is about residue fields —
only rules out the functionals factoring through `O/p`. A `True` there means
"not obstructed by the residue field", and the result carries a
`residue_field_only` flag saying so.

`cyclotomic_can_block(n)` returns the verdict from the arithmetic of 5, and
`orbit_can_block(kind, f)` decides an orbit by exhaustion. Two that fall out:
`Q(zeta_11)` blocks (`f = 5`, swapped), and therefore so does `Q(zeta_33)` —
which contains `zeta_3` for the triangles **and** `sqrt-11` for the Moser
rotation, in degree 20.

## The condition at 3, and the smallest field where both can hold

The condition at 5 says which fields can block. There is a second one, at 3,
that says which fields can carry the construction at all — and they pull in
opposite directions.

A rhombus `0, w, zeta_6 w, w(1 + zeta_6)` forces its tip to the apex's colour
at three colours, and `|1 + zeta_6| = sqrt3` whatever `w` is. Chain `k` of them
in directions `w_1, .., w_k`: every partial sum is forced to the origin's
colour, and the chain **closes** into a 4-chromatic graph when the last point
is also adjacent to the origin —

    | w_1 + .. + w_k |^2 = 1/3.

Now let `q` be a prime above 3 fixed by conjugation. Every unit step has
`v_q(u) + v_q(ubar) = v_q(1) = 0`, and `v_q(ubar) = v_q(u)` because `q` is
fixed, so `v_q(u) = 0`. A sum of such steps has `v_q >= 0`, while `1/3` has
`v_q = -e < 0`. Hence

> **Theorem.** If every prime above 3 is fixed by complex conjugation, no
> rhombus chain ever closes — at any length `k`, over any such field.

This generalises the earlier no-Eisenstein-spindle result from `k = 2` to every
`k`, and it *explains* the Moser spindle rather than merely permitting it:
`rho = (5 + sqrt-11)/6` has valuation `+1` at one prime above 3 of
`Q(sqrt-11)` and `-1` at its conjugate, which is possible only because 3
**splits** there. A chain needs some prime of `F` above 3 to split in `K`.

Stated for what it is: this closes the rhombus-chain mechanism — how the Moser
spindle and every spindle-type construction here works. It does not claim no
4-chromatic unit-distance graph exists over such a field by some other route.

### The two conditions together

| field | degree | residue deg at 5 | blocks | 3 splits |
|---|---|---|---|---|
| `Q(zeta_9)`, `Q(zeta_18)`, `Q(zeta_21)`, `Q(zeta_27)`, `Q(zeta_42)` | 6–18 | `>= 3` | yes | **no** |
| `Q(zeta_24)`, `Q(zeta_48)` | 8, 16 | 2, 4 | **no** | yes |
| **`Q(zeta_33)`** | **20** | **10** | **yes** | **yes** |

`Q(zeta_21)` blocks — one `sigma`-fixed prime above 5 of residue degree 6 —
but 3 is totally ramified there with a single prime, so nothing closes. The
first field where a chain construction and a blocked one can be **the same
object** is

> `Q(zeta_33) = Q(zeta_3, zeta_11)`, degree 20,

which carries `zeta_3` for the triangles and `sqrt-11` for the Moser rotation,
has residue degree 10 over `F` at 5, and has its two primes above 3 swapped by
conjugation. That the arithmetic picks out exactly the field containing both
the triangle root and the spindle rotation was not put in by hand.

## Building in `Q(zeta_33)`: what worked, and what did not

Two things were built in the field the arithmetic picked out, and both results
are worth having.

**Chains of three exist.** Over `Q(zeta_33)` the closing condition
`|1 + a + b|^2 = 1/3` — equivalently `Re(a) + Re(b) + Re(a.bbar) = -2/3` — has
solutions, and two of them give

    10 points, 19 edges, chi = 4

genuinely new 4-chromatic unit-distance graphs: three chained rhombi rather
than two paired ones, not a Moser spindle. The `k = 2` control returns exactly
the two rotations `(5 +- sqrt-11)/6`, as it must.

**And they still do not block.** Their 14 directions have module rank 4 — the
steps landed inside `Q(zeta_3, sqrt-11)`, the spindle's own field, which the
residue-degree theorem closes permanently. The graph is new; the arithmetic it
sits on is not. 964 of the 986 steps enumerated do lie outside that subfield,
so the narrowed search is not vacuous; it is running.

### Rotating about a shared point gives a bouquet

The 367-point union of the spindle with its 144 blocking rotations is
4-chromatic and blocked, and its 4-critical core is one spindle. The obvious
repair: delete a hitting set for the spindles, so no copy is 4-chromatic alone
and any surviving core must span copies.

**The hitting set has size one.** All 145 copies are rotations about the
origin, so they share exactly that vertex; delete it and the remaining 366
points are 3-colourable. The union was never an interlocking structure — it is
a bouquet of spindles tied at a knot, and the whole chromatic number lives at
the knot.

That says what any repair has to do: rotate about a point the copies do *not*
share, so they overlap in many vertices — which is exactly what de Grey does,
rotating `Y` about `(-2, 0)` rather than about one of its own vertices.

### Errors caught here, kept rather than edited away

Two enumeration faults, both of the kind that returns a clean zero and reads
like a theorem:

* the first `Q(zeta_33)` pass enumerated `a/abar` over a box of low powers of
  `zeta_33` and reported **no Moser spindle at all** — impossible, since
  `sqrt-11` lies in `Q(zeta_11)`. The Gauss sum
  `g = sum (k|11) zeta_11^k` needs `zeta_33^{3k}` up to `k = 10`, far outside
  that box. Fixed by constructing `g` directly and generating the unit steps
  *multiplicatively* — they form a group, so products of a few generators
  reach further than any box.
* the subfield filter returned "0 of 986 steps lie outside `Q(zeta_3,
  sqrt-11)`", because the membership test was written as a pivot check on a
  stacked matrix rather than the rank comparison it should have been. The
  correct test says 964 of 986.

## The three-rhombus chain, parametrised — and the field becomes a choice

The chain closes when `|1 + a + b|^2 = 1/3`. With `c = 1 + a`, `m0 = |c|^2 =
2 + r` and `r = a + abar`, the closing quadratic `cbar.b^2 - T.b + c = 0` has
`T = -2/3 - m0` and discriminant

    Delta = T^2 - 4 m0 = r^2 + (4/3) r - 8/9,

which depends **only on the real part of `a`**. Since `K = F(sqrt-3)`, that is
a square in `K` exactly when `3(8 - 12r - 9r^2)` is one in `F`; and `a` itself
exists as a unit step with `a + abar = r` exactly when `3(4 - r^2)` is. So the
whole question is a curve:

    y1^2 = 3(4 - r^2),      y2^2 = 3(8 - 12 r - 9 r^2).

The first conic carries the rational point `(r, y1) = (1, 3)`, so it is
rationally parametrised, and substituting collapses the pair to **one quartic**:

    r  = (m^2 - 6m - 3)/(m^2 + 3),        y1 = 3 + m(r - 1)
    V(m) = -39 m^4 + 540 m^3 - 666 m^2 - 324 m + 297,    y2 = sqrt(V)/(m^2+3)
    a  = (r + y1 sqrt-3 / 3)/2
    b  = (-8/3 - r + y2 sqrt-3 / 9) / (2(1 + abar))

Two conditions then come for free. `r + 2 = 3(m-1)^2/(m^2+3)` and
`2 - r = (m+3)^2/(m^2+3)` are both non-negative, so `|r| <= 2` at *every* real
embedding and the unit step `a` always exists. And a closing chain forces some
prime above 3 to be moved by conjugation — proved above — so **building** the
chain settles the condition at 3 instead of having to impose it. Hence

> **Every real `m` with `V(m) > 0` gives a chain of three rhombi in the plane,
> over the field `K = Q(m, sqrt(V(m)), sqrt-3)`.**

Checked numerically to machine precision and then exactly, in a degree-12
tower.

### The field stops being given and starts being chosen

Take `m` cubic. Then `F = Q(m, sqrt V)` is a totally real sextic, and the
residue degrees at 5 are read straight off the factorisation of
`prod_i (X^2 - V(m_i))` mod 5. Among the totally real cubics with `V` totally
positive, **40 give a field where every prime above 5 has residue degree at
least 3** — blocking-capable by the theorem above. Built over
`T^3 - 9T^2 + 14T + 8`, the chain is

    10 points, 16 edges, chi = 4, 20 directions, MODULE RANK 6

the first 4-chromatic graph here whose directions escape the spindle's rank-4
field. **The arithmetic no longer obstructs.** What does is the direction
count: 20 hyperplanes do not cover `PG(5,5)`.

### Why a longer chain does not help

The rhombus on `w` has edges `w`, `zeta_6 w`, and diagonal
`w(1 - zeta_6) = w.zeta_6bar`. So a step contributes exactly its
**`zeta_6`-orbit** — three projective directions — and splitting `w` into
`w.zeta_6 + w.zeta_6bar`, always possible since `zeta_6 + zeta_6bar = 1`,
stays inside that orbit. Measured: 20 directions at `k = 3`, and still 20 at
`k = 32` across 77 points and 161 edges. Length is free and buys nothing.

New directions need new **orbits** — which the chain does allow, because only
the closing sum is constrained and `w_1, .., w_{k-2}` are free.

### Errors caught in this stretch

* the pair sweep compared `(a+abar) + (b+bbar) + (a.bbar+abar.b)` against
  `-4/3`, but expanding `|1+a+b|^2 = 3 + 2[...]` makes the target `-8/3`. It
  was searching for `|1+a+b|^2 = 5/3`. The `k = 2` control it carried was
  right either way — `2 + (a+abar) = 1/3` really is `-5/3` — which is exactly
  why the fault survived. The parametrisation above has no target constant to
  get wrong.
* the cubic sweep used `3 * V(m)` where `V(m)` already carries the factor 3,
  so it searched `Q(sqrt(3V))` instead of `Q(sqrt V)` — different fields
  entirely. Caught by writing the derivation out numerically first.

## Blocking that carries weight: one chain, every edge load bearing

The field is chosen first. Over the degree-12 tower
`K = Q(m, sqrt V(m), sqrt-3)` with `m` a root of `T^3 - 9T^2 + 14T + 8`, every
prime of the real subfield above 5 has residue degree at least 3 — so the
residue-degree theorem does not close it — and measured directly, **the 450
directions of its 75 `zeta_6`-orbits block, at rank 12**. The material is there
before any graph is drawn. That ordering matters: blocking is monotone in the
direction set, so establishing that the field *can* supply a blocking set
comes before hunting for a graph that uses one.

Then the chain is grown into it. The move is to **replace a step `w` by three
unit steps summing to `w`**: any `u` leaves `z = w - u` to be split as
`v1 + v2`, which is a hash lookup once the pairwise sums of the enumerated
steps are tabulated. The total sum is untouched, so `|sum w_i|^2 = 1/3`
survives — asserted exactly at every round — so the closing edge survives and
`chi` stays 4, while each replacement brings new orbits into the direction set.

    round 22:  47 steps, 25 orbits, 150 directions, rank 12,  BLOCKS
    graph:     141 points, 236 edges, chi = 4, 152 directions, BLOCKED

### Why this one is different

| construction | blocked | what carried it |
|---|---|---|
| 107 points, `chi = 4` | yes | **pendants** — which change no chromatic number |
| 367 points, `chi = 4` | yes | every edge a spindle edge, but a hitting set for all 145 spindles has size **one**: a bouquet tied at a knot |
| **141 points, `chi = 4`** | **yes** | **one chain** — no copies to discard, no pendants to strip |

Every edge here is an edge of a rhombus. But the next sentence, that every
rhombus carries one link of the forcing, was a prediction — and it is wrong.

### Retraction: the core collapses here too

The 4-critical core of the 141-point graph is

    10 points, 16 edges, 20 directions, rank 6 — and it does NOT block,

which is a three-rhombus chain: the profile the `k = 3` chain had before any
growth. The long chain is not minimal, and the blocking is still not load
bearing in the strict sense. The claim is withdrawn.

The reason is visible once looked for. Every partial sum
`B_j = (1 + zeta_6)(w_1 + .. + w_j)` is forced to the origin's colour, so if
**any** two of them land at distance 1, that pair already contradicts a
3-colouring and a shorter chain closes inside the long one. With 47 steps
there are over a thousand such pairs, and some hit. Greedy deletion keeps the
short one and discards the rest — the same collapse the bouquet showed,
reached by a different route.

The repair is a condition imposed during growth rather than hoped for: accept
a replacement only when no pair `(B_i, B_j)` other than the closing
`(B_0, B_k)` sits at distance 1. With that guard the chain grows through 5, 7,
9 steps at ranks 8, 10, 12 with no premature closing.

Recorded because a prediction about a critical core has now failed twice here,
in two different constructions. They are not worth making.

## The gate, strengthened: every finite abelian quotient

A coset colouring is the bottom rung of something much larger. For **any**
finite abelian `G` and homomorphism `phi` from the edge module `M` to `G` with
`0` not in `phi(D)`, colour `v` by the colour of `phi(v - v0)` in a proper
colouring of the Cayley graph `Cay(G, phi(D))`. Adjacent vertices differ by an
element of `D`, `phi` sends it to a nonzero connection element, so their
colours differ:

    chi(Gamma) <= chi(Cay(G, phi(D)))     for every such phi

> **Every 6-chromatic unit-distance graph has Cayley chromatic number at least
> 6 in every finite abelian quotient.**

Blocking is the case `G = Z/5`: there `phi(D)` missing `0` makes the Cayley
graph the complete graph `K_5`, whose chromatic number is exactly 5, so *any*
valid `phi` settles it. That is why "no `phi` to `Z/5`" was the right
condition — and why it is only the first of a family. Above 5 the Cayley graph
is no longer complete and its chromatic number has to be computed, so the
screen keeps biting after blocking stops.

`periodic_screen(vecs, n)` runs it by CEGAR: solve for a `phi`, colour its
Cayley graph, and if that needs more than five colours, exclude this `phi` and
solve again. `exhausted` with no `phi` at all is exactly blocking, for that
`n`. A hit at any modulus proves the graph 5-colourable outright, and no
amount of blocking at 5 can save it.

Run on the 300 denominator-29 directions, which block at 5:

| modulus | verdict |
|---|---|
| 5 | no homomorphism at all — blocking |
| 6 | **exhausted**: all 48 homomorphisms need more than 5 colours |
| 7 … 20 | 60 homomorphisms tried at each, every one needing more than 5 |

So the set survives the stronger gate at every modulus from 5 to 20. That is
worth more than blocking alone: it is the condition an actual 6-chromatic
candidate would have to meet, and it is checkable one modulus at a time.

## Why every known construction stops at five

A rotation `rho` has `|rho| = 1`, so `rho.rhobar = 1`, and its chord is
`|1 - rho|^2 = 2 - (rho + rhobar)`. Hence

    the chord is rational  <=>  rho + rhobar is rational
                           <=>  rho is a root of X^2 - (2-c)X + 1
                           <=>  [Q(rho) : Q] <= 2.

**A rotation with rational chord lives in a quadratic field.** Compose several
and the field generated is a compositum of quadratic fields — multiquadratic —
whose Galois group is elementary abelian. Every decomposition group at 5 is
then elementary abelian with cyclic quotient by inertia, so the residue degree
is at most 2, below the bound; and if `K/F` happens to ramify above 5, the
collapse of `N` to `{+1,-1}` settles it instead. Either way:

> **A multiquadratic unit-distance graph has a coset 5-colouring, so its
> chromatic number is at most 5.**

And **every rotation in de Grey's construction has a rational chord** —
computed from the angles, not read off the paper:

| rotation | `\|1 - rho\|^2` | field |
|---|---|---|
| hexagonal, 60° | `1` | `Q(sqrt-3)` |
| Moser, `2 arcsin(1/(2 sqrt3))` | `1/3` | `Q(sqrt-11)` |
| `Sb`, `2 arcsin(1/4)` | `1/4` | `Q(sqrt-15)` |
| `Ya`, `pi/2 + arcsin(1/8)` | `9/4` | `Q(sqrt-7)` |
| `Yb`, `pi/2 - arcsin(1/8)` | `7/4` | `Q(sqrt-7)` |

whose compositum is `Q(sqrt-3, sqrt-11, sqrt-15, sqrt-7)` — exactly the
`Q(sqrt3, sqrt5, sqrt7, sqrt11)` the paper names, since
`sqrt-15 = sqrt-3 . sqrt5`. The field was not assumed; it falls out of the
five chords.

> **Corollary. A 6-chromatic unit-distance graph must use a rotation whose
> chord is irrational.**

> **Corrected later**, in "Exoo–Ismailescu rebuilt, and the denominator that
> the project filtered out": the residue-degree theorem behind both statements
> needs edge vectors integral at 5, which this summary dropped. An edge vector
> with 5 in its denominator can block every coset 5-colouring, so neither the
> claim nor the corollary holds as written.

That is the whole explanation of the barrier at five. Rational chords are what
one naturally reaches for — they are the rotations carrying a lattice point to
another at a rational distance — and each buys exactly one square root. No
amount of ingenuity inside that habit can pass five, because the obstruction
is arithmetic and is fixed the moment the rotations are chosen.

It also says what the chain construction above is for, and that it delivers.
The chain's rotation has `a + abar = r = (m^2 - 6m - 3)/(m^2 + 3)` with `m` a
root of `T^3 - 9T^2 + 14T + 8`, so its chord is

    2 - r = (919 + 222 m - 36 m^2)/397,    of degree 3 over Q.

**Irrational** — the first construction here meeting the condition the
corollary demands of any 6-chromatic candidate.

## A 4-critical blocked unit-distance graph

The chain closed into a **necklace** is, as a graph, a cycle of `k` diamonds —
each rhombus `B_{j-1}, M_j, M'_j, B_j` is `K_4` minus the edge
`B_{j-1}B_j`, five edges — plus the closing edge `B_k B_0`. So it has

    3k + 1 points and 5k + 1 edges, and nothing else,

provided no two of its points accidentally land at distance 1. That proviso is
checkable by counting, and the 47-step necklace meets it exactly: **142 points,
236 edges**, and `5·47 + 1 = 236`. Every rhombus edge is forced to be present
by construction, so a total of `5k+1` leaves room for no others.

Then 4-criticality is a **proof**, not a prediction — which matters, since a
prediction about a critical core has failed twice here.

* `chi >= 4`: in each diamond `B_{j-1}, M_j, M'_j` is a triangle and `B_j` is
  adjacent to both middles, so `c(B_j) = c(B_{j-1})`. All tips share a colour,
  and `B_k` is adjacent to `B_0`.
* Remove a middle `M_j`: its diamond becomes the path `B_{j-1} M'_j B_j` and
  stops forcing. Colour the tips before `j` with 0, those from `j` on with 1,
  `M'_j` the third colour, and each intact diamond's middles the two colours
  its tips leave. `B_k = 1` differs from `B_0 = 0`. Proper.
* Remove a tip `B_j`: the two diamonds meeting there become triangles, the two
  runs of tips colour independently, and the same assignment works. Removing
  `B_0` deletes the closing edge outright.

**Every vertex is critical.** Verified combinatorially for `k = 3, 4, 5, 7`,
vertex by vertex, with the geometry stripped out.

> **A closed rhombus necklace whose step orbits block is a 4-critical blocked
> unit-distance graph.**

And the 47-step necklace over the degree-12 tower — grown with the
premature-closing guard, so no shorter necklace sits inside it — reaches 25
orbits and 150 directions, and **those block**.

**Confirmed by machine, independently of the proof.** Greedy deletion over all
142 vertices could remove none of them:

    4-CRITICAL CORE: 142 points, 236 edges, 152 directions, rank 12, BLOCKED

The core is the whole graph. Unlike the two earlier predictions about critical
cores — both of which failed — this one was a proof before it was a
measurement.

That is blocking which cannot be separated from the chromatic number. Not
pendants, which change nothing. Not a bouquet, where one vertex carried
everything. A graph in which deleting **any** vertex drops it to 3-colourable,
and whose direction set admits no coset 5-colouring at all.

A shorter one exists too: 43 steps, 23 orbits, 138 directions, also blocking —
130 points.

### And the stronger gate bites, on this very object

The necklace's directions block at 5, and survive `n = 6` and `n = 7`. At
`n = 8` they do not:

| modulus | verdict |
|---|---|
| 5 | no homomorphism at all — blocking |
| 6 | 60 homomorphisms tried, all needing more than 5 |
| 7 | no homomorphism at all — blocking |
| **8** | **5-colourable**: the seventh homomorphism has Cayley chromatic number **4** |

So colouring `v` by that colouring of `phi(v - v0)` is a proper 5-colouring of
*any* graph over those directions. **The necklace can never be the substrate
of a 6-chromatic graph**, however blocked and however critical it is.

### What that failure actually was

Diagnosed rather than left as a curiosity. The failing homomorphism sends the
directions onto `S = {1,2,3,5,6,7}` — everything but **4** — and
`Cay(Z/8, S)` has chromatic number 4, because `{0,4},{1,5},{2,6},{3,7}` are
independent pairs. Missing any element other than 4 leaves `chi = 8`.

And `{0, 4}` is a **subgroup** of `Z/8`. Avoiding it is the same as a
homomorphism to `Z/8 / {0,4} = Z/4` nonzero on every direction — a coset
colouring with **four** colours. So the `n = 8` result was no subtlety at all:
it said the necklace is 4-colourable, which it is.

That sharpens the screen and makes it *cheaper*. A coset colouring mod `n`
gives `chi <= n`, so

> **A 6-chromatic unit-distance graph must block at `n = 2, 3, 4` and `5`.**

Four SAT calls, no Cayley chromatic numbers. Measured:

| direction set | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| necklace, 138 | blocks | blocks | **coset colouring** | blocks |
| denominator-29, 300 | blocks | blocks | blocks | blocks |
| the field's 450 | blocks | blocks | blocks | blocks |

So the necklace's gap is **not forced** — the field has the material, and a
longer necklace can reach it. That is what the growth now targets.

And the target is cheap. Running the CEGAR cover over all four moduli at once
— hold a set of orbits, ask `has_homomorphism` at each `n` for a `phi` it
misses, add the orbit killing the most escapes collected so far — the loop
closes at

    26 orbits, 156 directions, blocking at 2, 3, 4 and 5

against the 25 that `n = 5` alone needed. The four conditions together cost
almost nothing over the one.

> **Corrected below.** Every "blocking at 2, 3, 4 and 5" in this section and
> the next was read off a search over `Z^d` rather than over the module the
> directions generate — see *A flaw in the blocking test*. On the right module
> `n = 4` falls away; `n = 5`, the gate, survives, and `n = 2, 3` were free
> all along for a 4-chromatic graph.

## A field with de Grey's spindle *and* the arithmetic a sixth colour needs

A distance `d` spindles iff `K` holds a rotation with `|1 - rho|^2 = 1/d^2`.
With `t = rho + rhobar = 2 - 1/d^2` one gets `4 - t^2 = (4d^2 - 1)/d^4`, so
since `K = F(sqrt-3)` the condition is

> a distance `d` spindles **iff `3(4d^2 - 1)` is a square in `F`**.

At `d^2 = 3` — de Grey's distance, the rhombus tip — that is **33**. And
`F = Q(m)(sqrt V)` contains `Q(sqrt33)` exactly when `V = 33 ·` (a square in
`Q(m)`): if `sqrt V` lay in `Q(m, sqrt D)` then `V = a^2 + b^2 D + 2ab sqrt D`
forces `ab = 0`, and `b = 0` would make `V` a square outright.

**33 is positive**, so `V = 33 s^2` is compatible with `V` totally positive,
which the chain needs. Searching the totally real cubics with `V` totally
positive, filtered by the necessary condition that `N(V)/33^3` be a rational
square, turns up two — and one of them is the object:

    x^3 - 10x^2 + 26x - 11        totally real, irreducible, 5 INERT
    V(m) = 1947 - 4653 m + 1848 m^2          totally positive
    V = 33 s^2   with   s = 13 - 26 m + 5 m^2        checked exactly

So `F = Q(m, sqrt33)` and `K = F(sqrt-3) = Q(m, sqrt-3, sqrt-11)`, degree 12,
which holds **the Moser rotation `(5 + sqrt-11)/6`**. And 5 is inert in `Q(m)`
— the cubic has no root mod 5 — and inert in `Q(sqrt33)`, since `33 = 3` mod 5
is a non-residue; so its residue degree in `F` is `lcm(3,2) = 6`, and `K/F` is
unramified there because `-3` is a unit at 5. **Both conditions of the
residue-degree theorem hold.**

> **The first field carrying both the rotation that takes four colours to five
> and the arithmetic that a sixth requires.**

Stated for what it is: this is where a 6-chromatic graph *could* live, not one
that does. And de Grey's **other** rotations are not available here — chord
`1/4` needs `sqrt-15`, hence `sqrt5`, and `F = Q(m)·Q(sqrt33)` has `Q(sqrt33)`
as its only quadratic subfield because `Q(m)` is cubic. His construction
cannot be transplanted whole. Only the final spindling step is available; the
4-chromatic forcing structure beneath it would have to be rebuilt from the
rotations this field *does* have, namely those with `3(4d^2 - 1)` in `(Q*)^2`
or in `33(Q*)^2`.

## Why four-colour forcing needs distance, not a gadget

The rhombus forces at three colours because the triangle `B, M, M'` uses all
three, so `B'` — adjacent to `M` and `M'` — has only `B`'s colour left. The
analogue at four would need a vertex whose neighbourhood uses three colours in
every 4-colouring. **It does not exist locally**, and the reason is short.

A triangle inside `N(v)` would make `{v} u T` a `K_4`, which no unit-distance
graph in the plane contains. More: `N(v)` lies on the unit circle about `v`,
where two points are adjacent exactly when they are 60° apart — so each has at
most two neighbours there, and **every cycle is a full hexagon**, of even
length. Hence

> **The neighbourhood of any vertex of a unit-distance graph is bipartite.**

So it is 2-colourable, and at four colours `v` always keeps at least two
choices. No vertex's colour is ever forced by its neighbourhood alone.

Measured on the necklace, which makes the point concretely: **0 of its 8385
non-adjacent pairs are monochromatic in every 4-colouring**. Each diamond
gives its tip two choices at four colours and nothing downstream closes them.
The necklace cannot be spindled to five.

This is why de Grey's graph is 1581 vertices and not 10. Four-colour forcing
has to be assembled out of long-range structure; it cannot be bought with a
local gadget the way the rhombus buys it at three.

## A 4-critical graph blocking at every modulus up to five

The first necklace blocked at 5 and failed at 4 — which is what its apparent
`n = 8` Cayley failure had been saying all along. Aimed at the 26-orbit cover
that blocks at all four, the growth closes one orbit further out:

    51 steps, 27 orbits, 162 directions, blocking at 2, 3, 4 and 5
    graph:  154 points, 256 edges  =  3·51+1 and 5·51+1, exactly
    chi = 4, 164 directions, BLOCKED
    4-critical core: the whole graph, 154 points, still BLOCKED

The counts matching `3k+1` and `5k+1` **exactly** is the hypothesis of the
necklace theorem — no two points accidentally at distance 1 — so
4-criticality is proved, and greedy deletion confirms it by removing nothing.

> **A 4-critical unit-distance graph whose directions admit no coset colouring
> with 2, 3, 4 or 5 colours.**

> **Corrected below.** Same flaw, same outcome: on the module the verdict is
> `2, 3, 5`. And of those four, only `n = 4` and `n = 5` were ever content —
> a 4-chromatic graph blocks at 2 and 3 for nothing, since a coset colouring
> mod `n` is a proper `n`-colouring. So what this graph really carries is the
> gate, and the `n = 4` claim is withdrawn.

And the Cayley screen, which caught the previous necklace at `n = 8`, now
passes there and at **every modulus from 5 to 20** — no homomorphism at all at
5, and at each of the other fifteen, sixty homomorphisms tried and every one
needing more than five colours.

Every necessary condition this work has produced, met at once, by a graph that
is critical rather than padded.

What it is **not** is 6-chromatic. It is 4-chromatic, and it cannot be
spindled upward: by the bipartite-neighbourhood theorem no vertex's colour is
ever forced by its neighbourhood at four colours, so the rhombus trick that
carried 3 to 4 has no analogue here. That is the honest shape of the gap.

## Four-colour forcing, extracted — and the same question at five

A single pair forced monochromatic in every 4-colouring is too much to hope
for: the necklace has none of its 8385, and forcing is never local. What works
is weaker and **disjunctive**. For a vertex `u`, say `(*)` holds at `u` when

> every 4-colouring puts `u`'s colour somewhere on `u`'s `sqrt3`-sphere.

One SAT call: fix `c(u) = 0`, forbid colour 0 on the whole sphere, and UNSAT
is the property. Measured:

| graph | points | vertices with `(*)` |
|---|---|---|
| `Sa` | 397 | **none** |
| `Y` | 791 | **the origin** only, on a sphere of 12 |

A full survey of all 791 vertices of `Y` finds **exactly one**: the shared
centre of the two copies. de Grey's forcing lives at a single point.

So the forcing appears exactly where de Grey puts it — in `Y = Sa u Sb`, not
in `Sa` alone — and the jump from 397 to 791 points is what buys it.

And `sqrt3` is forced to be the distance: `|1 - rho|^2 = 1/3` means the
rotation angle has `cos = 5/6`, so a point at `sqrt3` moves to distance
`sqrt3 · (1/sqrt3) = 1` — adjacent. The Moser rotation spindles this sphere
and no other.

### The same question one colour up

That is the shot at six, so it was asked directly of the only 5-chromatic
graph there is:

> is there a vertex `u` of `G` such that **every 5-colouring** puts `u`'s
> colour somewhere on `u`'s `sqrt3`-sphere?

**No** — not at one of `G`'s 1581 vertices. And since `sqrt3` was only forced
at four colours by the rotation that spindles it, the sweep was widened to the
**60 commonest squared distances** in `G`, each tested at every centre with a
sphere of four or more points — spheres up to 12, over a thousand centres per
distance:

    0 hits.

> **de Grey's `G` does not force at five colours on any common sphere.**

Consistent rather than surprising: `G` is 5-chromatic, so its 5-colourings are
plentiful where `Y`'s 4-colourings are tight. The step that bought forcing at
four was `Sa` to `Y`, 397 points to 791. The analogue at five has to be taken
above 1581 — and that is the honest frontier.

## de Grey's last step is one forced pair, and his field is his chain

Read `G = Ya u Yb` backwards. It is two copies of `Y` turned about
`p = (-2, 0)` by `pi/2 +- arcsin(1/8)`; the angle between the copies is
`2 arcsin(1/8)`, so a point at distance 4 from `p` has its two images exactly
one apart. For that to be a contradiction rather than a case split, `Y` has to
supply not a sphere condition but a **pair**:

> `p` and `q` at distance `d`, monochromatic in every proper 4-colouring.

Given one, the rest is a line. Let `rho` be a rotation about `p` with
`|rho| = 1` and `|1 - rho| = 1/d`, so `|q - rho(q)| = d.(1/d) = 1`. In
`H u rho_p(H)` both copies read `c(q) = c(p)`, because `rho` fixes `p` and
carries `H` onto a copy of itself — but `q` and `rho(q)` are adjacent.

> **Theorem.** If `H` carries a pair at distance `d` monochromatic in every
> proper `k`-colouring, and the field admits a rotation closing `d`, then
> `H u rho_p(H)` is `(k+1)`-chromatic.

That lemma is not new here — it is `hn/spindle.py`, stated in its docstring
and automated by `ForcedPairFinder`, from earlier in this work. What is new is
reading de Grey's own last step through it: *which* pair, *how many* edges the
union adds, and the arithmetic of which distances a field can close.

The edge counts say this is exactly what `G` is. `Y` has 791 vertices and 3938
edges; `G` has `1581 = 2.791 - 1` and `7877 = 2.3938 + 1`. **The union adds one
edge.** That edge is the whole step from four colours to five.

And the pair is there. Scanning `Y`'s six distance-4 pairs, exactly one is
monochromatic in every 4-colouring — `(2,0)` and `(-2,0)` — and the solver has
to work for it, 282 seconds of unsatisfiability proof against milliseconds for
the five that are not forced. `(-2,0)` is de Grey's pivot.

### Which distances a field can close

Pure arithmetic. With `D = d^2` the rotation satisfies
`rho + rhobar = 2 - 1/D` and `rho.rhobar = 1`, so it is a root of
`t^2 - (2 - 1/D)t + 1`, of discriminant `-(4D - 1)/D^2`. In real coordinates
that is `cos = 1 - 1/(2D)`, rational always, and

    sin = sqrt(4D - 1) / (2D),

so the only question is whether `sqrt(4D - 1)` lies in the field, with
`D >= 1/4` forced since `2d sin(theta/2) = 1` is unachievable below it. Now run
de Grey's own chain of distances through it:

| step | what it is | `D` | `4D - 1` | radicand |
|---|---|---|---|---|
| 1 | the triangular lattice | 1 | 3 | **3** |
| 2 | the Moser spindle | 3 | 11 | **11** |
| 3 | `Sb = rho(Sa)`, `sin = sqrt15/8` | 4 | 15 | **15** |
| 4 | `Ya u Yb`, `sin = 3 sqrt7 / 8` | 16 | 63 | **7** |

`Q(sqrt3, sqrt5, sqrt7, sqrt11)` — his field, in the order he introduces it.
**His field is not a choice. It is what his chain demands, term by term.**

### The chain doubles, and it stops where the field stops

`Y`'s forced pair is not anywhere: `(2,0)` and `(-2,0)` both sit on the ring of
radius `sqrt(D) = 2` about the union's pivot, antipodally, so at squared
distance `4D`. That is where the rigidity is — the ring is exactly the set the
rotation moves by one. So the chain doubles, `D -> 4D`, and the step from `D`
needs `sqrt(16D - 1)`:

| `D` | 1 | 4 | 16 | 64 | 256 |
|---|---|---|---|---|---|
| `16D - 1` | 15 | 63 | 255 | 1023 | 4095 |
| radicand | 15 | 7 | **255 = 3.5.17** | 1023 | 455 |

de Grey has `3, 5, 7, 11`, so he closes 1, 4 and 16 and stops: 255 asks for
`sqrt17`, which he has not got. **His construction is exactly as long as his
field allows.** Adjoining `sqrt17` reopens it — `closable_distance(64, (3,5,7,11,17))`
is true.

### What this buys, and what it does not

Two things, immediately. First a pruning: the forced pair lives on the ring,
so a union is scanned in sixty pairs instead of eighty thousand. Second a
statement of what is missing at five colours, which can now be measured rather
than guessed.

Measured, with the pair test (fixing `c(a) = 0` is free, colours being
interchangeable, so each query is one pair of assumptions against a single
incremental solver, and unsatisfiable means forced):

| graph | | pairs at a closable distance | forced |
|---|---|---|---|
| `Sa`, over `K` | 397 pts, 1974 edges | 4200 | none |
| `X`, the Moser closure over `K` | 597 pts, 1476 edges | 1762 | none |
| `Y`, de Grey's | 791 pts, 3938 edges | 6 at distance 4 | **one** |
| `G`, de Grey's, at **five** colours | 1581 pts, 7877 edges | 21358 | none |

At four colours the pattern is de Grey's own: neither seed forces, the union
does. At five colours `G` carries nothing, and not narrowly — every one of the
21358 queries came back satisfiable inside a 40000-conflict budget, and not one
needed the budget. **A sixth colour is not one rotation away from `G`.**

### What `K` can and cannot close

`K = Q(m, sqrt-3, sqrt-11)` is not of the form `F(i)`, and the difference
bites. A rational `r` has `sqrt(r)` in `K` only through one of `K`'s three
quadratic subfields `Q(sqrt-3)`, `Q(sqrt-11)`, `Q(sqrt33)` — the cubic `Q(m)`
admits none, 2 not dividing 3 — so `D` is closable over `K` exactly when the
squarefree part of `1 - 4D` is `-3` or `-11`. Integer `D` up to 139:

    1, 3, 7, 19, 25, 37, 61, 69, 91, 127, 135

`K` keeps de Grey's first two steps and **loses both of the last two**: 4 needs
`sqrt-15`, 16 needs `sqrt-7`, and `K` has neither. What it offers instead is
`1 -> 3 -> 7`, and the Moser closure over `K` realises `D = 1/3, 1, 3` and `7`
— 588 pairs at `sqrt7`. So the port of de Grey's architecture to the one field
that blocks at every modulus up to five is not a translation. It is a different
chain, and finding its forced pair is the open end of this work.

### The doubling step is almost nowhere over `K`

Where `Y`'s forced pair sits is not incidental: `(2,0)` and `(-2,0)` are
antipodal on the ring the rotation moves by one. That is an observed
mechanism, not a theorem — nothing forbids a forced pair off the ring — but it
is *reproducible*. Rerunning de Grey's own step under the ring pruning, with
no knowledge of his answer, turns up exactly one forced pair in
`Sa u rho_4(Sa)`, at `d^2 = 16 = 4D`, after 239 seconds of unsatisfiability
proof. The pruning recovers his construction from scratch.

If that is the mechanism, then asking which `D` admit the step — `D` and `4D`
both closable — separates the two fields sharply. Over `Q(sqrt3, sqrt5, sqrt7,
sqrt11)` they are common:

    1/2, 1, 4, 17/2, 37/4, 61/4, 86, 397/4, 721/4, 271, ...

and de Grey walks two in a row, `1 -> 4 -> 16`, stopping because 16 is not on
the list. Over `K`, searched to `D = 40000`, there are **two**: `7/12` and
`397/4` — and neither chains again (`7/3 -> 28/3` needs `sqrt327`;
`397 -> 1588` needs `sqrt6351 = sqrt(3.29.73)`).

The reason is a congruence. If `4D - 1 = 3u^2` then `16D - 1 = 3(4u^2 + 1)`,
and `3(4u^2+1) = 11v^2` forces `u^2 = 8 mod 11`, a non-residue, while
`4u^2 + 1 = w^2` has only `u = 0`. If `4D - 1 = 11u^2` then
`16D - 1 = 44u^2 + 3`, never `11v^2`, and `3v^2` only along the Pell equation
`v^2 - 132 u^2 = 1`, whose fundamental solution `23^2 - 132.2^2 = 1` gives
`D = 397/4` and whose next solution is far out of reach of any unit-distance
graph.

So the field that blocks at every modulus up to five is precisely the field
where de Grey's mechanism has almost nowhere to stand. That is not a defeat —
it is the sharpest statement yet of why the two halves of this problem resist
being solved at once, and it says exactly what a construction over `K` must
do instead: find its forced pair off the ring.

## A flaw in the blocking test, an audit, and one withdrawn claim

A coset colouring is a homomorphism `phi` from the group the points generate
to `Z/n`; it is proper exactly when `phi` is nonzero on every edge vector, and
blocking is the assertion that no such `phi` exists. The group that decides it
is `M`, the subgroup the edge vectors generate. `has_homomorphism` was
searching over `Z^d` — the ambient lattice of whatever coordinates the vectors
happened to be written in.

Those are different questions, and the difference points the wrong way. `M`
sits inside `Z^d`, possibly properly, and `Z/n` is not injective, so a
homomorphism `M -> Z/n` need not extend to `Z^d`. The search therefore finds
*too few* functionals and can report blocking that is not there. One dimension
is enough to see it:

    d = 1, D = {2}, n = 2.  Over M = 2Z the map phi(2) = 1 escapes.
    Over Z every psi has psi(2) = 0, and the test says "blocked".

Dividing out the global content kills that example, which is why the pipeline
has done it from the start — but content 1 does not give `M = Z^d`.

### When the two tests can disagree, decided without SAT

`M` sits in its saturation with finite quotient `T`, and restriction
`Hom(Z^d, Z/n) -> Hom(M, Z/n)` is onto exactly when `Ext^1(T, Z/n) = T/nT`
vanishes — when no invariant factor of `M` shares a prime with `n`. An
invariant factor is divisible by `p` exactly when the rank of the direction
matrix drops mod `p`. So

> the two tests agree at `n` **iff** `rank_p = rank_Q` for every prime `p | n`,

two Gaussian eliminations and no solver. And only one direction of a verdict
was ever at risk: *"does not block"* exhibits a `phi`, and a `phi` on `Z^d`
restricts to `M`, so **every escape ever found is genuine**. Only blocking
claims needed rechecking.

### The audit

| object | | over `Z^d` | over `M` | |
|---|---|---|---|---|
| denominator-29 set | 300 vectors, index 1 in `Z^12` | 2,3,4,5 | 2,3,4,5 | unchanged |
| de Grey's `G` | 133 vectors, rank 16 of 32 | 2,3,4 | 2,3,4 | unchanged |
| `U = G u w.G` | via a rotated denominator-29 subset | 2,3,4,5 | 2,3,4,5 | unchanged |
| the 133-point 4-critical necklace | 140 directions, rank 12 | 2,3,4,5 | **2,3,5** | **corrected** |
| `G u rho(G)` over `Q(m)(sqrt3,5,7,11)` | 268 directions, rank 32 | 2,3,4,5 | **2,3,4** | **withdrawn** |

`U`'s blocking survives because monotonicity survives the module version: if
`S` is inside `D` then every `phi` on `<D>` restricts to `<S>`, so a blocked
subset blocks the whole set. Its blocking comes from a rotated copy of the
denominator-29 set, which is saturated, so nothing moved.

The last row is the one that exposed the flaw. It was about to be reported as
a second 5-chromatic graph blocking at every modulus up to five, got cheaply
by bolting a rotated copy of de Grey's graph onto itself over a field where 5
has residue degree 3 — `m` a root of `x^3 - 10x^2 + 26x - 11`, irreducible mod
5. The chromatic number is untouched (the two copies share one vertex, and
`1581 = 2.791 - 1` all over again), and over `Z^32` the directions blocked at
five. Over their own lattice they do not. **The claim is withdrawn.**

### The gate verdict was never at risk

The rank test says so without any Hermite reduction. The necklace's 140
directions have rank 12 over `Q` and

    rank mod 5 = 12,   rank mod 3 = 4,   rank mod 2 = 2.

Full rank mod 5 means `M` is 5-saturated, so restriction
`Hom(Z^12, Z/5) -> Hom(M, Z/5)` is onto and the two tests agree at `n = 5`
**by theorem**. The collapse to rank 2 mod 2 and rank 4 mod 3 is exactly where
the ambient test was saying nothing — and `n = 4` is where the verdict moved.

### What survives

The gate is `n = 5`, and everything that mattered there stands. `U` still
meets the sharpened gate in full. The necklace still blocks at the gate, and
at 2 and 3 — what it loses is `n = 4`, so it no longer meets the sharpened
gate in full. For a 4-chromatic graph, admitting a coset colouring mod 4 is
what one should expect; the interesting thing was always `n = 5`.

## What a forced pair actually is, and how far five is from six

`H` carries a pair `a, b` monochromatic in every proper `k`-colouring exactly
when `H` plus the single edge `ab` has no proper `k`-colouring. Same sentence
twice, and the second reading is the useful one:

> **A forced pair is one extra edge that pushes `chi` from `k` to `k+1`.**

The spindle turns that virtual edge into geometry. If the field closes
`|a - b|` then `rho` about `a` sends `b` to distance exactly 1 from itself,
both copies of `H` in `H u rho_a(H)` read `c(b) = c(a)`, and the virtual edge
is a real one. So the whole of *get from `k` colours to `k+1`* is:

> find a `k`-chromatic unit-distance graph that becomes `(k+1)`-chromatic when
> **one** edge is added, at a distance the field can close.

At `k = 3` that graph is the unit rhombus and the edge joins its two tips at
`sqrt3` — `4D - 1 = 11`, and the closing rotation `(5 + sqrt-11)/6` is Moser's.
At `k = 4` it is `Y` and the edge joins `(2,0)` to `(-2,0)` at 4 — `4D-1 = 63`.
At `k = 5` nothing of the sort is known, and this says how far off it is:
**`G` is 5-chromatic, and adding any one of its 21358 pairs at a closable
distance leaves it 5-colourable**, with no query even reaching a
40000-conflict budget.

### Density is not obtained by closing under rotations

`S` works as a seed because it is 39 points in a small region whose twelve
images overlap heavily — 397, not 468. Doing the same to `G`:

    W = G closed under the 12-element dihedral group
    18966 points out of 12 x 1581 = 18972, so six shared and nothing else
    94548 edges, mean degree 10.0 -- G's own 10.0
    0 forced pairs at five colours among the 400 busiest vertices

`G` spans too much for its own rotations to interlock. **The move that builds
the fourth floor does not build the fifth.**

## Most of the sharpened gate is free — and one of my own claims is a tautology

A coset colouring mod `n` **is** a proper `n`-colouring. So if one exists then
`chi <= n`, and contrapositively

> `chi(Gamma) > n`  implies  `Gamma` blocks at `n`,

with no arithmetic involved at all. Two consequences, and the first is against
this work's own framing.

**The sharpened gate is a tautology.** "Every 6-chromatic unit-distance graph
blocks at `n = 2, 3, 4, 5`" is true because such a graph has no proper
`n`-colouring for any `n <= 5`, so it certainly has no coset one. As a
necessary condition on the target it says nothing that `chi >= 6` does not
already say. Measured on `G`, which is 5-chromatic:

    blocks at 2, 3, 4  -- every modulus below its chromatic number, free
    fails at 5         -- the barrier theorem, the one place content existed
    blocks at 6        -- content of a different kind

**What survives is the reading as a design criterion**, and there it is sharp.
For a graph of chromatic number `c`, blocking at `n < c` is free and blocking
at `n >= c` is real. So of the necklace's four moduli:

    chi = 4:   n = 2, 3 free;   n = 4 and n = 5 are the content.

The corrected verdict `[2, 3, 5]` therefore says the necklace **keeps the
gate** — which is the thing a sixth colour needs — and loses the one other
real fact it had. The headline "blocks at every modulus up to five" was, for a
4-chromatic graph, two free facts and two real ones, and one of the real ones
went.

## The anatomy of `Y`: one shared point and six edges

Measured, not inferred. `Sa` and `Sb = rho_4(Sa)` have 397 points each and
their intersection is **a single point, the origin**; the union is 793, and
`Y` is that less the two vertices de Grey drops. The edges decompose exactly:

    2 x 1974 = 3948   the two copies
          + 6         edges joining the exclusive half of Sa to that of Sb
         - 16         taken away with the two dropped vertices
        = 3938

So `Y` is two copies of `Sa` glued at one point by six edges — and that is the
entire difference between a graph with no forced pair and one that has one.
`Sa` alone: 4200 pairs at a closable distance, **none forced**. `Y`: the pair
`(2,0), (-2,0)`, which lies **wholly inside the `Sa` copy**. Six edges landing
elsewhere are what make two points of `Sa` unable to differ.

That is a sharp and cheap filter for the same search over another field. A
rotation about the origin producing no cross edge glues two lumps at a point
and can force nothing new, so only the ones that *bite* deserve a pair scan.
Sweeping the 3030 units of `K = Q(m, sqrt-3, sqrt-11)` as rotations of `Sa`
finds them in quantity, and harder than his:

> **one shared point and thirty cross edges**, against his six.

## Solving for the rotations that bite, instead of sampling for them

Sampling finds none. All 3030 known units of `K`, used as rotations of `Sa`,
leave `Sa` and its image sharing only the origin with **zero** cross edges
once the float filter is checked exactly — the counts a `1e-9` tolerance
reported were artefacts of evaluating long elements of `K` in double
precision. Nor should sampling work: a cross edge is an exact algebraic
condition, and de Grey did not stumble on `2 arcsin(1/4)` either.

Set it up properly. A cross edge is `p, q` in `Sa` with `|p - u.q| = 1` and
`|u| = 1`. Put `w = u.q`, `A = |q|^2`, `P = |p|^2`; then `|w|^2 = A` and
`|w - p|^2 = 1` expand to

    w.pbar + wbar.p = A + P - 1 =: 2R,      |w.pbar|^2 = A.P,

so `w.pbar` is a root of `t^2 - 2Rt + A.P`, and the rotation exists over a
field exactly when `sqrt(R^2 - A.P)` does. Both `R` and `A.P` are real and the
radicand is negative, so over `K = F(sqrt-3)` with `F = Q(m, sqrt33)`:

> the rotation lies in `K` **iff** `(R^2 - A.P) / (-3)` is a square in `F`.

And `F = Q(m)(sqrt33)` turns that into two square-root tests in the cubic
`Q(m)`, each screened for free by the norm being a rational square, then read
off the cubic's three real embeddings to fifty digits and confirmed exactly.
`Sa` is `zeta_6`-invariant and `u.Sa = (zeta_6.u).Sa`, so `q` runs over 66
orbit representatives rather than 397 points — **26136 pairs, 81 rotations**.

And they bite far harder than his six. Exact cross-edge counts:

    156, 156, 126, 126, 108, 108, 108, 108, 60, 60, 48, 48, 48, 48, 36, 36, ...

## `rho_4` is not one rotation among many

The unions built over `K` are built on *his* `Sa`, so `(2,0)` and `(-2,0)` sit
in every one of them — and the filtered pair scan never asked about them,
because `D = 16` is not closable over `K`. It should have. A hit would have
been spindleable one quadratic step away, over

    K(sqrt-15) = Q(m, sqrt-3, sqrt-11, sqrt-15),

and adjoining a quadratic to `Q(m)` cannot drop the residue degree at 5 below
3 — the prime with `f = 3` either stays inert, splits into primes still of
`f = 3`, or ramifies with `f = 3`. So that field still blocks, while carrying
`m` through `u`, which is exactly what de Grey's own field can never supply.

Asked of all 81: **the pair is free in every one.** A couple of seconds of
satisfiability each, against the 282 seconds of unsatisfiability proof `rho_4`
costs. Among every rotation `K` supplies that bites `Sa` at all, not one
forces his pair.

Which says what those six cross edges are worth. The `K` unions reach **156**
of them against his six and force nothing — at any `K`-closable distance, or
at his. **Six edges in the right place beat a hundred and fifty-six anywhere.**

## How far five is from six, measured rather than guessed

A pair is forced when the solver cannot separate it. Short of that, the *work*
it takes to separate it is a distance to forcing: a pair that comes apart in
no conflicts is wide open, one that costs tens of thousands is nearly pinned.
The scale is calibrated by the case where the answer is known.

| | conflicts |
|---|---|
| `Y` at four colours, the five distance-4 pairs that separate | 157, 4756, 12227, 16483, 35365 |
| `Y` at four colours, the pair that is **forced** | > 2 000 000 (budget exhausted) |
| `G` at five colours, **all 21358 pairs together** | **6410 in total**; 2044 cost anything; dearest 29 |

So `G`'s hardest pair at five colours costs **29** conflicts, against the
**157** of `Y`'s *easiest non-forced* pair at four, and the two million of the
one that is forced. The entire five-colour scan costs a fifth of what one easy
four-colour pair does. Five is not near six, and this is how near it is not.

And the forcing does not shrink. The ball of radius 1.5 about the segment
joining `(2,0)` and `(-2,0)` holds **773 of `Y`'s 791** vertices and the pair
comes apart — so the eighteen vertices *furthest* from the pair are load
bearing. Forcing is not a local phenomenon that a smaller gadget could carry,
and whatever forces at five colours will not be small either.

## `Y` has no slack: every vertex is essential to its forcing

A stratified sample of 61 of `Y`'s 789 non-pair vertices, taken by degree, and
**every one of them is essential** — delete it and `(2,0)`, `(-2,0)` take
different colours. Thirty-nine answered inside a 60000-conflict budget; the
other twenty-two answered the *same way*, in seconds, once the budget was
removed.

That last part is a lesson in itself. `solve_limited` gave up long before its
bound: vertices 184 and 359 came back undecided at 60000 conflicts and then
separated in three and seven seconds with no budget at all. **Every budgeted
"undecided" has to be rerun to the end before it means anything.** The error
points the safe way, though — giving up early produces *more* undecideds,
never fewer, so a scan reporting **zero** hard queries decided every one of
them, and the negatives elsewhere in this work stand untouched.

> **Corrected.** The sample was not the population. A full sweep of all 789
> non-pair vertices, no budget anywhere, finds **vertex 630 — degree 4 —
> slack**: `Y` minus it still forces the pair. So `Y` is **not**
> vertex-critical, and "every vertex is essential" holds only of the 61
> sampled. Sampling by degree took every thirteenth vertex and this one fell
> between. What survives is that `Y` is *nearly* minimal — one exception in
> the first 630 checked — not that it is minimal.

## The counting threshold, and exactly where it applies

A direction set of rank `r` over `F_n` blocks once it carries enough
independent lines: a random functional survives `L` of them with probability
`((n-1)/n)^L`, and there are `n^r` functionals, so

    L  >  r . log n / log(n/(n-1)),    which at n = 5 is 7.21 r.

It is a heuristic — the directions are structured, not random — so it is worth
checking against every verdict already established. **At the gate it calls
both known cases right:**

| | rank | lines mod 5 | threshold | predicts | truth |
|---|---|---|---|---|---|
| denominator-29 set | 12 | 120 | 87 | blocks | **blocks** |
| de Grey's `G` | 16 | 54 | 116 | escapes | **escapes** |

Below the gate it is inapplicable, and the reason is sharp. `G`'s directions
contain a vector **congruent to zero** mod 2, mod 3 and mod 4 — the set holds
`n` times one of its own members — so every `phi` kills it and the blocking is
trivial rather than statistical. At `n = 5` there is no zero residue among its
109, and the count governs.

So the prediction the threshold makes for the sixteen-rotation set — 880 lines
against a threshold of 693, and no zero residue — is worth the wait.

## Translates, never tried here — and the measurement finally moves

Every union in this attack was `G` with a **rotated** copy, because the
spindle needs a rotation. But the union that carries the *forcing* does not.
`G u (G + t)` is just as good a unit-distance graph, just as 5-chromatic, and
translates cost nothing arithmetically:

    a cross edge is  t = p - q - v    for p, q in G and v a unit direction,

so every difference of `G` shifted by any unit vector is a candidate. No
square has to lie in the field, no rotation has to exist. The family is vastly
larger than the rotations and it had gone completely unexamined.

The histogram over `G` — 21 million distinct translates, keyed on floats and
confirmed exactly — puts the best non-trivial one at 451 from a sample of 150
of the 1581 points. Measured exactly it gives **3026 points, 136 shared,
16820 edges, 1442 cross**. And then the solver cost, which is the quantity
that matters:

| | pairs | conflicts | dearest |
|---|---|---|---|
| `G` alone at five colours | 21358 | 6 410 | **29** |
| **two disjoint copies of `G`** (control) | 42716 | 13 301 | **28** |
| the translate union, 1442 cross | 44011 | **723 211** | **3 689** |
| a second translate union, 1440 cross | 43871 | 499 237 | 3 278 |

**The control is what makes this evidence.** Two copies of `G` a thousand
apart, where no cross edge is geometrically possible, cost exactly twice
`G`'s own total and their dearest pair does not move at all. So the 54-fold
jump in the total and the 132-fold jump in the dearest pair come **entirely
from the 1442 cross edges**, not from the instance being bigger.

For scale: `Y`'s *unforced* distance-4 pairs at four colours cost 157 to 35365
conflicts, and its forced one over two million. A single translate lifts `G`
out of "nowhere near" — 29 — and into the bottom of the band where forcing
actually lives.

Still no forced pair. But this is the first thing in the whole search that has
moved the measurement rather than the vertex count.

### The law I fitted to it, and how it broke

Three points — the control at 0 cross and 28, a glide reflection at 393 and
151, a translate at 1442 and 3689 — fit `dearest = 28 . exp(cross/L)` with
`L = 233` and `295`. At `L = 300` that put the forced regime, where `Y`'s own
forced pair sits above two million, at about **3350 cross edges**, between
depth 2 of the translate stack (2864) and depth 3 (4306).

**It broke at the first test.** Depth 3 carries 4306 cross, where the law
predicts 48 million; the dearest pair in 50 sampled came to **131 260** —
short by a factor of 365. And the variance at fixed crossing is large:

    1440, 1442, 1442, 1444 cross  ->  3278, 3689, 3184, 2279
    1528 cross                    ->  17 555

A fivefold jump for six percent more crossing. So the cross count does **not**
determine the cost — structure does — and no single exponent fits.

What survives is worth keeping: the cost grows steeply with crossing, the
control at zero cross is flat, and the best union so far — 1528 cross, dearest
17 555 — sits near the **top** of `Y`'s unforced band of 157 to 35 365 rather
than at its bottom. That is the closest anything in this search has come.

## Sampling, and a metric that finally has a calibration

Deciding whether a pair is forced by SAT costs the solver real time, and a
union of two copies of `G` offers about 44 000 candidate pairs. Separating
them all runs to roughly 700 000 conflicts and twenty minutes, so the census of
399 translates would have taken a week.

But a pair is forced only if it agrees in *every* proper 5-colouring, and one
colouring in which it differs settles it for good. Sampling colourings
collapses the candidate set geometrically — 44 414 pairs to zero in a dozen
samples — and only the survivors need a SAT call. The filter is sound in the
direction that matters: a separation comes with an explicit colouring as its
witness, so no forced pair can be filtered away. Twelve seconds a union
instead of twenty minutes, and a *stronger* verdict, because "zero survive"
proves there is no forced pair rather than reporting that none was found
inside a budget.

Diversity has to come from randomising the solver's saved phases. Cadical
ignores `set_phases` and hands back the same model every time — one distinct
colouring in eight tries — so it samples nothing; glucose honours them and
gives eight out of eight.

The point of a curve is that it can be calibrated, and `Y` is the calibration:
it forces a monochromatic pair in every proper 4-colouring, so its curve
cannot reach zero. The filter found that pair unaided — one of 10 647
candidates, by sampling alone, with nothing told to it about where to look.

| graph | colours | head rate ÷ (1/k) | tail rate | floor | forced |
|---|---|---|---|---|---|
| `Y` | 4 | 1.36 | **0.882** | 2 | **1** |
| `Y` | 5 | 1.07 | — | 0 | 0 |
| `G` | 5 | 1.16 | — | 0 | 0 |
| stack, depth 1 | 5 | 1.19 | — | 0 | 0 |
| stack, depth 2 | 5 | 1.21 | — | 0 | 0 |
| stack, depth 3 | 5 | 1.18 | 0.525 | 0 | 0 |

The head of the curve says almost nothing: it is dominated by generic pairs,
and generic pairs behave generically — everything at five colours sits between
1.07 and 1.21. The tail is the signal. `Y` at four colours decays at 0.882 a
sample and settles on a floor; everything at five crashes to zero by sample 7
to 9. Depth 3 of the stack is the first object in this search to grow a tail
at all, and it still reaches zero.

> **Corrected.** An earlier reading averaged ratios across marks that are not
> evenly spaced — 1, 2, 3, 5, 10, 20, 40 — mixing a two-sample step in with
> one-sample steps, which understated the rates (0.228 and 0.229 for depths 1
> and 2, against 0.238 and 0.242 computed from consecutive marks alone). The
> conclusion drawn from them, that stacking does nothing, was too strong. The
> head does barely move; the tail at depth 3 is new.

## Symmetrise first, then rotate

`Sa` is the 12-element dihedral orbit of `S` about the origin, and it is
exactly invariant: 397 points in, 397 out. That invariance is what makes
rotating it useful. Every point stays on its own ring, so a rotation by any
angle whose chord at that radius is 1 produces cross edges in bulk. Rotating a
lopsided graph produces nothing, which is exactly what the earlier scan found:
987 rotations of de Grey's field bite `G`, and the best of them contributes
four cross edges.

So close under the symmetry *before* rotating — and about the centre the graph
was actually built around. About the origin, `G`'s dihedral orbit is 6 × 1581
points and 6 × 7877 edges, to the last unit: six copies that never touch. `G`
does not live at the origin. It lives at `(-2,0)`, the pivot its two copies of
`Y` were turned about, and about that point its orbit is 13 873 points and
73 782 edges — five thousand points short of twelve disjoint copies.

`Y`'s own closure is smaller than it looks: 791 points become 1189, not 9492.
Rotating `Sb` by 60° gives `Sb` back, because `Sa` is 60°-invariant and `Sb`
is a rotation of `Sa`; so the whole dihedral orbit of `Y` is
`Sa ∪ Sb ∪ Sb'`, with `Sb'` the mirror rotation, and 397 + 396 + 396 = 1189.

## The ring has to pay twice

Turning a symmetric core about its own centre is the step from `Sa` to `Y`,
and the ring it turns on is spent twice over.

Once on the rotation itself. To make a ring of squared radius `D` produce
cross edges, the turn has to carry each of its points to distance exactly 1
from where it started, which means `cos t = 1 - 1/(2D)` and
`sin t = √(4D-1)/(2D)`. The field has to contain `√(4D-1)`.

Once on the spindle afterwards. The pair the construction is aiming at is the
*antipodal* pair on that ring, which sits at squared distance `4D`, and
forcing it is worth nothing unless it can then be spindled — which needs
`√(16D-1)`.

Asking for both at once is a savage filter, and the reason to trust it is that
it retrodicts de Grey's own choice. `Sa` has five rational closable rings —
1, 1/3, 5/9, 4 and 3. `D = 1` is the 60° turn and fixes `Sa`, so it is
degenerate. Of the rest **exactly one** passes both tests, and it is `D = 4`,
the ring de Grey used, whose antipodal pair `(2,0)`, `(-2,0)` is precisely the
pair `Y` forces.

It is not the populous ring. `D = 1` carries 30 of `Sa`'s 397 points and
`D = 4` carries six. Ranking rings by population — the obvious heuristic, and
the one tried here first — picks the wrong one. Population is not the
property; being spendable twice is.

Applied to `G*`, which is 5-chromatic and exactly dihedral and so plays `Sa`'s
role one level up: of 572 rings, twelve are rational and closable, and two pay
both costs — `D = 4` again, and `D = 17/2`, which `Sa` has no analogue of
(`ρ` by `√33`, spindle by `√135 = 3√15`). `D = 16`, where de Grey's own
doubling chain stopped for want of `√17`, still does not qualify.

So the five-colour step has two candidates rather than a search space, and the
pairs worth interrogating are the ones on the ring — a few dozen, not the 360
million of a 27 000-point union. Along the way, `rotation_joining(4)` turns
out to be exactly de Grey's rotation, `cos 7/8` and `sin √15/8`, and it
carries `Sa` onto `Sb`, all 397 points of it.

## Forcing is monotone, and that cuts both ways

A pair forced in `H` is forced in every graph containing `H`: each proper
colouring of the larger restricts to one of the smaller, so a pair that agrees
in all colourings of `H` agrees in all colourings of anything built on it.
Growth can create forcing; it can never destroy it.

That answers a question for free. The dihedral closure `Y*` contains `Y`, so
it forces `Y`'s pair, and since `Y*` is invariant under the 12-element group
which carries `(2,0)`, `(-2,0)` to the antipodal pairs of the ring of squared
radius 4, all six of those are forced. A sampling-and-SAT sweep of `Y*` had
spent twenty-five minutes on the first of them before it was stopped.

The other direction stings more. A union that forces nothing has *no subgraph
that forces anything* — so every negative recorded here is a negative about
everything inside it, which is a large claim cheaply bought, and it means a
search that only ever grows its graphs is spending its budget in the wrong
place. `G*` has 13 873 points and forces nothing; `Y` has 791 and forces a
pair. Size is not the variable.

## The recipe, driven only by the criterion, rebuilds de Grey's graph

Start from `Sa` and the criterion alone, knowing nothing about the answer. Of
`Sa`'s five rational closable rings exactly one pays both costs, `D = 4`, and
that single choice fixes everything downstream with no further decisions:

| step | what the criterion forces |
|---|---|
| `ρ` | `rotation_joining(4)` about `Sa`'s centre — `cos 7/8`, `sin √15/8` |
| `Z` | `Sa ∪ ρ(Sa)`, which forces the ring's antipodal pair |
| pivot | one end of that pair, `(-2,0)` |
| `σ` | `rotation_joining(16)` about the pivot — the pair is at squared distance `4D = 16` |
| `W` | `Z ∪ σ(Z)` |

Then, measured rather than assumed:

- `ρ(Sa) = Sb`, all 397 points;
- `A = (2,0)` and `σ(A)` come out at squared distance **exactly 1**;
- `Y ∪ σ(Y)`, turned by `π/2 - arcsin(1/8)` about the pivot, **is `G`** —
  1581 of 1581 points, and `|G| = 1581`.

That last line is the validation. de Grey builds `G` as `Ya ∪ Yb`, two turns of
`Y` through `π/2 ± arcsin(1/8)`; those two angles differ by `2 arcsin(1/8)`,
which is precisely `σ`. His pair of turns and this recipe's single turn
describe the same configuration, so the recipe lands on his graph exactly
rather than on something merely like it.

`W` therefore contains a rotated copy of `G` and cannot be 4-coloured — a
proof by containment, with no solver anywhere in the argument.

What this licenses is the same template one level up, with `G*` in `Sa`'s place
and its two doubly-usable rings, `D = 4` and `D = 17/2`, in place of the one.
What it does not license is optimism about the outcome: the criterion says
where to look, not that anything is there.

## The same template, one level up, does not climb

`G*` is the dihedral closure of `G` about its own pivot: 13 873 points,
5-chromatic because it contains `G`, and exactly invariant under the
12-element group. That is `Sa`'s role one level up, and the recipe that
rebuilds de Grey's graph from `Sa` can be run on it unchanged.

Geometrically it runs through without a hitch. `ρ = rotation_joining(4)` about
the pivot gives `Z = G* ∪ ρ(G*)` at 27 673 points; the ring's antipodal pair
sits at squared distance 16; `σ = rotation_joining(16)` about one of its ends
puts `A` and `σ(A)` at squared distance **exactly 1**, as it must; and
`W = Z ∪ σ(Z)` is 55 345 points and 314 283 edges.

And `W` is 5-colourable. 842 373 conflicts to find the colouring — sixteen
times what `G*` alone costs and an order of magnitude past anything else here
— but found. So `Z` does not force the antipodal pair on de Grey's own ring at
five colours.

The contrast with the level below is the measurement, not the verdict:

| | `Z₄ = Sa ∪ ρ(Sa)` at 4 | `Z₅ = G* ∪ ρ(G*)` at 5 |
|---|---|---|
| points | 793 | 27 673 |
| antipodal pairs surviving 14 samples | **6 of 6** | — |
| pairs agreeing in every sample | 75 (35 closable) | — |
| spindle lands at distance 1 | yes | yes |
| spindled union colourable | **no** (that is `χ ≥ 5`) | **yes**, 842 373 conflicts |

The recipe is not what fails. It reproduces de Grey's graph exactly one level
down and its geometry is exact one level up. What fails is that `G*`'s
5-colourings are not rigid the way `Sa`'s 4-colourings are — and closing under
the symmetry, which was the one structural lever available, did not make them
so.

## The census closes at zero, and it was the wrong thing to count

All 399 translates are done — each unioned with `G`, each interrogated for a
forced pair at five colours — and the answer is zero. Not a budget running
out: every candidate pair in every union was separated by an explicit sampled
colouring, so each negative is witnessed, and monotonicity extends it to every
subgraph of every union. Cross counts ran from 1552 at the top of the census
down to 126 at the bottom; the candidate sets, 43 000 to 46 000 pairs a union,
collapsed to nothing in eight to fifty-five samples. 4773 seconds for the lot.

And the quantity being counted was never the right one. `Z₄ = Sa ∪ ρ(Sa)` has
3954 edges where two loose copies of `Sa` have 3948 — so de Grey's forcing at
four colours is produced by **six** cross edges. The translate unions here
carried 1552 and the stacks 4306.

The obvious repair is that six edges can only close a structure that is nearly
closed already, so the real quantity would be how constrained a single copy
is. That is wrong too, and the measurement is flat about it:

| graph | colours | pairs agreeing in all 14 samples |
|---|---|---|
| `Sa` | 4 | **0** |
| `Sa ∪ ρ(Sa)` | 4 | **75** (35 closable, incl. all 6 antipodal) |
| `Y` | 4 | 20 (1 closable — the forced one) |
| `Sa` | 5 | 0 |
| `G` | 5 | 0 |
| `G*` | 5 | 0 |

`Sa` alone is completely free. Nothing is nearly closed. Six edges take a
graph with no agreeing pairs at all to one with seventy-five, and the pair de
Grey needs is among them.

So forcing does not accumulate; it arrives. No quantity measured anywhere in
this search — cross edges, agreement, size, solver conflicts — rises towards
it beforehand, which means there is no gradient for a search to climb. That is
a real constraint on method, not a complaint: candidates have to be *proposed*
by structure and then tested, and the ring criterion is exactly such a
proposer.

That `G*` is free at five colours is therefore not evidence against it. `Sa`
is equally free at four, and forcing still arrived when it was rotated.

## Past where the chain stopped: adjoining √17

de Grey's doubling chain is `D = 1, 3, 4, 16`, and it halts at 16 because the
next step wants `√(16·16 - 1) = √255 = √(3·5·17)` and his field has no `√17`.
That field was forced on him by the chain, not chosen, and nothing obliges a
search to stay inside it.

`G*`'s ring `D = 16` pays its two radicals separately. `ρ` needs
`√(4·16-1) = √63 = 3√7`, already present — `ρ` is `31/32 + i(3/32)√7`, the
same turn that served as the *spindle* at `D = 4`. Only the spindle of the
antipodal pair, at squared distance 64, needs `√255`. So the extension buys
one step at the very end, and the copy it produces is still an exact
unit-distance graph: its coordinates live in a 32-dimensional field instead of
a 16-dimensional one.

| ring | field | shared | `Z` | `W` | edges | conflicts | 5-colourable |
|---|---|---|---|---|---|---|---|
| `D = 4` | `K` | 73 | 27 673 | 55 345 | 314 283 | 842 373 | yes |
| `D = 17/2` | `K` | — | 27 745 | 55 489 | — | 169 887 | yes |
| `D = 16` | `K(√17)` | 6 937 | 20 809 | 41 617 | 221 353 | 136 144 | yes |

All three land exactly — `A` and `σ(A)` at squared distance 1 every time — and
all three colour. The structural programme on `G*` is exhausted, cleanly.

One more intuition dies on the way out. `D = 16` is the ring whose rotation
overlaps `G*` with itself most — 6937 shared points, a third of a copy,
against 73 for `D = 4` — and its spindled union is the *easiest* of the three
to colour. Overlap does not predict hardness either.

## The control, and what it rules out

A negative is worth what the method behind it is worth, so the pivot scan was
run again on `Sa` at four colours, unchanged: every pivot orbit, every
rational closable ring about it, the rotation carrying that ring to distance
1, the same twelve sampled colourings and the same bucketed agreement filter.

The first union it looked at:

```
pivot 0, D = 4: 793 points, 6 cross edges, 211 agree, 28 closable
```

Pivot 0 is the origin, `D = 4` is de Grey's ring, and 793 points with six
cross edges is `Sa ∪ Sb` exactly. **The scan rediscovers his construction
unprompted, as the first thing it tries** — and finds others besides: 256
agreeing pairs at pivot 7 with `D = 3`, 209 at `D = 1/3`, 148 at `D = 5/9`,
and 22 from a union carrying a *single* cross edge.

Against which, at five colours: 12 199 (pivot, ring) candidates over `G`, 44
distinct rings that actually bite, exact cross counts to 812, and **zero**
agreeing pairs in all but one — which had one, at a distance that cannot be
closed.

So the search is not what fails. The negative at five colours is a statement
about five colours.

## Why, as far as the measurements reach

The first explanation was criticality: `Y` is nearly 4-critical, so its
colourings are scarce, while `G` is 5-chromatic on 1581 vertices where about
500 suffice. It does not survive contact with the numbers.

> **Withdrawn.** `χ(Sa) = 4`, measured — not 3-colourable in 21 conflicts,
> 4-colourable in 909 — so the four-colour core does sit at its threshold, as
> `G` does at five. But `Sa` reaches `χ = 4` on **397** vertices where the
> Moser spindle does it on **seven**, fifty-seven times over, while `G`
> reaches `χ = 5` on 1581 where ~500 suffice, three times over. The graph
> that pins 211 pairs is the more redundant *by a factor of twenty*.
> Redundancy is not the variable.
>
> **And `G` is not redundant at all.** Nineteen of thirty probed vertices are
> essential — `G - v` is 4-colourable — with an explicit verified colouring
> for one of them (0 monochromatic edges of 7869, every one of 1580 vertices
> with exactly one colour). `G` is *nearly vertex-critical*, and its
> five-colour relation is still exactly its edge set. That is a stronger
> statement than the one withdrawn: the negative does not rest on `G` being
> loose. (No contradiction with the literature — 1581 is already de Grey's own
> reduction of his 20 425, and the smaller graphs come from fresh search, not
> from deleting vertices here.)

What is left is the plane's, not the graph's. A unit-distance graph in the
plane has no clique larger than a triangle. Against four colours a triangle
spends three and leaves one spare; against five it spends three and leaves
two. Every local structure that propagates constraint — triangle, rhombus,
spindle — works against twice the slack, and the measurements agree
uniformly: at four colours six cross edges pin 211 pairs and one cross edge
pins 22, while at five colours 812 pin nothing.

That is an explanation fitted to the negatives, not a theorem. What it
predicts is that no amount of ingenuity with rotations of *these* graphs
closes the gap: the search has to find local structure that constrains five
colours, and the plane does not obviously provide any.

## The one measurement that explains all the negatives

A graph constrains its colourings through two relations: pairs forced to
*agree*, and pairs forced to *differ*. Edges supply the second for free, so
the question is what survives beyond them. Both halves can be filtered the
same way, run in opposite directions — a pair that differs in some sampled
colouring cannot be forced to agree, and a pair that agrees in one cannot be
forced to differ — and every survivor put to the solver.

| graph | colours | pairs | chance | candidates | ratio | verified same / different |
|---|---|---|---|---|---|---|
| `G` | 5 | 1 248 990 | 5 898 | 6 502 | **1.10** | **0 / 0** |
| `Y` | 5 | 312 445 | 1 475 | 1 659 | **1.12** | **0 / 0** |
| `Sa` | 5 | 78 606 | 371 | 496 | 1.34 | **0 / 0** |
| `Sa` | 4 | 78 606 | 79 | 1 548 | **19.63** | counts only |
| `Y` | 4 | 312 445 | 314 | 6 871 | **21.92** | counts only |

"Chance" is `((k-1)/k)²⁴`, the rate at which an unconstrained pair survives
twenty-four samples. At four colours these graphs sit **twenty times above
it**; at five they sit *at* it, and the verification there is absolute rather
than statistical: all 6502 of `G`'s candidates were checked individually, and
every one is free.

> **Careful.** The four-colour excess is *correlation*, not forcing. Forty of
> `Sa`'s 1548 candidates were drawn at random and put to the solver with no
> budget: all forty are free. Two controls fix the distinction — a rhombus at
> three colours has one non-edge and it *is* forced-same, so the test does
> find constrained pairs; the Moser spindle at four is 4-critical with ten
> non-edges and not one is constrained, so criticality alone guarantees
> nothing. What the table shows is that at four colours these colourings are
> strongly correlated and at five they are indistinguishable from
> independent. Forcing is absent from both bare cores — de Grey's pair lives
> in `Sa ∪ ρ(Sa)`, not in `Sa`. Correlation is the raw material; the rotation
> converts it. At five colours there is nothing to convert.

> **`G`'s colour relation at five colours is exactly its edge set.** No pair
> of non-adjacent points is constrained in either direction, in any proper
> 5-colouring.

That single fact accounts for every negative in this work at once — the
translates, the stacks, the pivot unions, the symmetric closures, the
spindles. None of them failed for want of size, crossing, or the right field.
They failed because there was nothing to work with: **a rotation of a graph
with no colour relation has no relation to combine.**

The rainbow search agrees. `forced.py` identifies a 5-rainbow — five points
pairwise forced onto different colours — as what `χ ≥ 6` would need, and since
the plane's clique number is 3 it must use pairs forced apart at distances
other than 1. There are none: of `G`'s 2840 triangles, four have a vertex
avoiding all three colours across 24 samples, and not one survives
verification. Not even a *4*-rainbow exists.

## The same graph, the same samples, only the colour count changes

Every comparison above moves two things at once: a four-colour graph against a
five-colour one, of different sizes and different constructions. One pair of
numbers does not.

| | candidates | chance | ratio |
|---|---|---|---|
| `Sa` at **four** colours | 1548 | 79 | **19.6×** |
| `Sa` at **five** colours | 496 | 371 | **1.34×** |

Three hundred and ninety-seven points, 1974 edges, twenty-four samples, one
graph. The only thing that differs is whether the solver is handed four
colours or five, and the correlation in its colourings collapses by a factor
of fifteen.

That is the finding in its cleanest form. It is not about de Grey's
construction, not about size, symmetry, the field or the crossing: **a
unit-distance graph whose colourings are strongly correlated at four colours
has colourings indistinguishable from independent at five.** The fifth colour
is enough slack to decouple them, and a rotation of a decoupled graph has
nothing to combine.

## The ladder, and why six is hard

One measurement, twenty-four samples, three graphs, each at its own chromatic
number and above it:

| | `k=3` | `k=4` | `k=5` | `k=6` |
|---|---|---|---|---|
| triangular lattice (`χ=3`) | **11 014×** | 2.2 | 0.9 | — |
| `Sa` (`χ=4`) | — | **16.8×** | 1.4 | 1.1 |
| `G` (`χ=5`) | — | — | **1.1×** | 1.0 |

At `k = χ` the correlation collapses as `χ` rises: **eleven thousand at three
colours, seventeen at four, one at five.** One colour above `χ` it is always
about one — 2.2, 1.4, 1.0 — so a single spare colour already decouples
everything, at every level.

The lattice validates the instrument on the way past. Its proper 3-colouring
is unique up to permuting the colours, so two points differ in every colouring
exactly when they lie in different classes, and for three classes of ~133 that
is `1 - 3·C(133,2)/C(400,2) = 66%`. Measured: 52 212 of 78 679 non-edge pairs.
Also 66%. The filter recognises the uniqueness with nothing told to it.

> So the difficulty of `χ(ℝ²) ≥ 6` is not a matter of finding the right graph.
> The plane's unit-distance graphs **stop constraining their own colourings
> somewhere between four colours and five**, and by five there is nothing left
> to constrain with. Every negative recorded in this work is a corollary of
> that table.

(`Sa` at four reads 16.8 here and 19.6 under a different seed; the order of
magnitude is the reading, not the digit.)

### The instrument check the ladder rests on

A correlation ratio inflates whenever the sampler returns near-identical
colourings — twenty-four copies of one colouring make every pair look
constrained — and that failure is real here: about one random subgraph of `Sa`
in forty returns twenty-four colourings differing on two per cent of vertices
and scores **700**. Renaming does not catch it, since colourings differing in
three vertices are still distinct tuples. Distance does.

| graph | `k` | ratio | samples disagree on | independent |
|---|---|---|---|---|
| triangular lattice | 3 | 11 014 | **0.000** | 0.667 |
| `Sa` | 4 | 24.0 | 0.652 | 0.750 |
| `Sa` | 5 | 1.2 | 0.725 | 0.800 |
| `Y` | 4 | 25.4 | 0.598 | 0.750 |
| `G` | 5 | 1.1 | **0.762** | 0.800 |

The lattice's zero is not a failure but the answer: its 3-colouring is unique
up to permutation, so twenty-four samples *must* coincide, and that is exactly
where its ratio of eleven thousand comes from. The rest are well mixed — `Sa`
at 87% of independent, `G` at 95% — so their ratios measure the graphs, not
the solver. (The spread sitting below the independent value throughout is the
thing being measured, not a defect: proper colourings of a constrained graph
are not independent random assignments.)

## The gadget is the unit, and at five colours it is seventy times bigger

The triangular lattice is rigid at three colours because every triangle spends
all three. At four a triangle leaves one spare, so something bigger must
propagate — and the smallest 4-chromatic graph is the Moser spindle: seven
vertices, two unit rhombi hinged at a point with their far diagonals at
distance 1. Counted directly:

| graph | `χ` | points | spindles | per point |
|---|---|---|---|---|
| triangular lattice | 3 | 400 | **0** | 0.00 |
| `Sa` | 4 | 397 | 576 | 1.45 |
| `Y` | 4 | 791 | 1152 | 1.46 |
| `G` | 5 | 1581 | 2304 | **1.46** |

The lattice has none, exactly as it must — it is 3-colourable and a spindle is
not — which is what makes the counter worth trusting elsewhere. And de Grey's
whole family carries one density, 1.45 per point, inherited rather than
achieved: `G` is built from `Y` and `Y` from `Sa`.

That is where `Sa`'s correlation at four colours comes from, and why `G`'s does
not survive at five. **A spindle is a 4-*critical* gadget**: it pins four
colours and leaves slack against five. The unit that would pin five is a
5-critical subgraph, and the smallest known is about five hundred vertices
against the spindle's seven.

Which says what scale would be needed rather than that none would do. To carry
5-critical gadgets as densely as `Sa` carries spindles takes roughly `500/7`
times as many points — `397 × 71`, about **28 000**.

`Z = G* ∪ ρ₄(G*)` has **27 673**, and is being measured.

## How the forcer grows, and the scale six would need

The spindling lemma turns a forced pair at `k` colours into a
`(k+1)`-chromatic graph, so the object whose size matters is the **forcer**:
the smallest graph carrying a pair monochromatic in every proper
`k`-colouring, at a distance the field can close.

| `k` | forcer | vertices | minimal? |
|---|---|---|---|
| 3 | the rhombus | **4** | yes — a triangle has no non-adjacent pair to force |
| 4 | `Y` | **791** | nearly — 787 of 789 non-pair vertices are individually indispensable |

A factor of **198**. If the next step costs the same, a five-colour forcer
runs to about **157 000 vertices** — and everything built in this session tops
out at 55 345.

Two points are two points, and this is an extrapolation rather than a theorem.
But both ends are measured rather than assumed, and it agrees with everything
else here: with the gadget count, where the 4-critical spindle is seven
vertices and the 5-critical unit is about five hundred; with the ladder, where
correlation at `k = χ` falls 11 014 → 24 → 1.1; and with the relation, which
at five colours is empty. Four independent measurements pointing at the same
wall.

## The mechanism, confirmed — and three escapes that fail

Three families, ten objects, spindles per point against correlation at four
colours:

| object | spindles/point | ratio at 4 |
|---|---|---|
| triangular lattice | 0.000 | 2.2 |
| lattice + 1 hinge | 0.000 | 2.0 |
| lattice + 4 hinges | 0.004 | 1.2 |
| `Sa` thinned to 0.4 | 0.289 | 1.3 |
| `Sa` thinned to 0.6 | 0.534 | 1.7 |
| `Sa` thinned to 0.7 | 0.705 | 2.4 |
| `Sa` thinned to 0.8 | 0.984 | 5.0 |
| **`Sa`** | **1.451** | **24.0** |
| `Y` | 1.456 | 25.4 |

Monotone end to end, across a lattice, a spindled lattice and thinned copies
of de Grey's core.

> **Corrected — the mechanism is confounded.** Thinning removes *edges* along
> with spindles, and average degree tracks the ratio better than spindle
> density does: Spearman +0.782 against +0.685, with the two predictors
> themselves correlated at +0.867. Worse, the data contain a control the first
> reading walked past. At average degree ≈ 5.6 there are five objects whose
> spindle density runs 0.000 → 0.534, and their ratios are 2.2, 2.0, 1.5, 1.2,
> 1.7 — **no trend at all.** At fixed degree, spindles do nothing.
>
> What survives is weaker: within a fixed colour count, correlation rises with
> the density of *constraints*, and edges are the constraints. Across colour
> counts neither predictor works — `Sa` and `G` have average degrees 9.94 and
> 9.96 and ratios 24.0 and 1.1. Only the ladder isolates the variable, by
> holding the graph fixed.

(An aside worth keeping: **`lattice + 1 hinge` is 4-chromatic and contains no
Moser spindle at all.** The count and the chromatic number are different
questions.)

Three ways out of the ladder were tried, and all three fail:

- **Size.** `G*` at 28× the minimum, sixty samples over 96 million pairs:
  ratio **1.14**, against `G`'s 1.10 at 3×. The five-colour size curve is flat
  throughout — 1.10, 1.18, 1.22, 1.14.
- **Base.** The lattice sits a thousand times above `Sa` in correlation, and
  its unique 3-colouring makes `(0,0)`,`(3,0)` a forced pair at distance 3
  with `√35 = √5·√7` in the field. Spindling it works perfectly — 721 points,
  not 3-colourable — and inherits **nothing**: ratio 1.6 at four colours,
  *below* `Sa`'s 24.
- **Hinges.** Stacking spindles on the lattice lowers it: 2.0 → 1.5 → 1.4 →
  1.2, because the added copies carry no gadgets of their own (0.004 per
  point against `Sa`'s 1.451).

## And the last piece: `G` contains one 5-critical subgraph, not thirty

The mechanism needs gadgets *at density*. `Sa` carries 576 Moser spindles on
397 points — ten memberships per point — and for 500-vertex gadgets to match
that, a graph would need about `n/50` distinct ones: thirty-two inside `G`.

It has one. If `G - v` is **4-colourable**, then `v` lies in *every* 5-critical
subgraph, and that is the cheap side of the query, since the solver need only
exhibit a colouring. Thirty random vertices, 400 000 conflicts each:

| | count |
|---|---|
| essential (`G - v` is 4-colourable) | **19** |
| dispensable (`G - v` still 5-chromatic) | **0** |
| undecided within the budget | 11 |

Nineteen proved, none refuted. So every 5-critical subgraph of `G` contains at
least those nineteen and, by proportion, at least 63% of its 1581 vertices —
about a thousand. The density of five-colour gadgets in `G` is therefore at
most **0.001 per point**, against the spindle's **1.45**. A factor of fourteen
hundred, and the account closes.

(The budget is safe-sided in the direction used: it can fail to find a
colouring but never invent one, so every "essential" verdict is sound. The
eleven undecided are genuinely open and counted neither way.)

## Criticality and correlation run opposite

The criticality diagnosis was withdrawn once already, on the weak grounds that
the graph which pins pairs is the more redundant. Measuring criticality
directly finishes it off, and the other way round from the intuition this
search began with:

| graph | `k` | one deletion drops `χ`? | correlation at `k` |
|---|---|---|---|
| `Sa` | 4 | **0 of 40** — never | **24.0** |
| `Y` | 4 | **0 of 40** — never | 25.4 |
| `G` | 5 | **19 of 30** — usually | **1.1** |

`Sa` and `Y` are as far from critical as a graph can be, and carry the
strongest colour relations measured here. `G` is nearly vertex-critical and
carries none.

In hindsight that is the right way round. **Critical means minimal** — exactly
enough constraint to force the chromatic number and nothing to spare — so
everything not load-bearing is free. It is the *redundant* structure that
correlates.

Which says the five-colour search cannot be fixed by finding a leaner
5-chromatic graph. Leaner is the wrong direction. What forced de Grey's pair
was `Sa`'s **surplus** — 397 vertices where seven suffice, 576 spindles packed
into them — and a five-colour object with that surplus has to be built, not
trimmed.

## The fifth colour does not shift the curve — it flattens it

With the gadget story confounded, what survived was that correlation rises
with the density of *constraints*. The controlled version: thin two graphs by
the same procedure and seeds, match their average degrees, and change only the
number of colours.

| average degree (`G` / `Sa`) | `G` at **five** | `Sa` at **four** |
|---|---|---|
| 9.96 / 9.94 | **1.1** | **24.0** |
| 8.06 / 8.11 | 1.0 | 5.0 |
| 6.88 / 7.14 | 1.0 | 2.4 |
| 6.07 / 5.66 | 0.9 | 1.7 |
| 4.95 / 4.87 | 0.8 | 2.1 |
| 3.83 / 3.84 | 0.6 | 1.3 |

At four colours the ratio climbs from 1.3 to 24 across that range. At five it
is **flat at one throughout**, drifting *below* one at low degree — which is
what noise around independence looks like.

> Constraint density is not weakened at five colours. It **stops operating**.
> Adding edges correlates a four-colour graph's colourings and does nothing at
> all to a five-colour graph's.

That is the cleanest form of everything above: same thinning, same seeds, same
samples, matched degrees, only the colour count differing. And it explains,
without any appeal to gadgets, why every construction here failed — translates,
stacks, pivot unions, symmetric closures and spindles all add points and
edges, and at five colours that is the one lever measured to do nothing.

(The four-colour series is monotone down to a ratio of about 2 and inverts
once below it — 1.7 at degree 5.66 against 2.1 at 4.87 — where the values sit
near one and the noise is the size of the signal. The 251.4 at degree 8.97 is
the sampler artefact recorded above, not a data point.)

## Surplus does not predict either

The last hypothesis standing was **surplus** — that correlation needs a graph
many times larger than the minimum for its own chromatic number. `Sa` is 397
points where seven suffice for `χ = 4` (57×) and scores 24; `G` is 1581 where
~500 suffice (3×) and scores 1.1. `Z` at 55× was being measured to test it, at
something like twelve hours.

The data already in hand refute it. Ordering the four-colour objects by
surplus:

| surplus | ratio | degree |
|---|---|---|
| 34× | 1.7 | 5.66 |
| 57× | **24.0** | 9.94 |
| 103× | 2.0 | 5.60 |
| 113× | **25.4** | 9.96 |
| 257× | 1.2 | 5.61 |

Surplus runs 34 → 257 and the ratio goes 1.7, 24.0, 2.0, 25.4, 1.2 — no order
at all. The degree column separates them perfectly: the two ratios above 20
are exactly the two objects at degree 9.95, the three below 2.1 exactly the
three at 5.6.

So within four colours degree predicts and surplus is irrelevant; at five
colours degree is flat. Neither axis leads anywhere, and the `Z` measurement
was stopped rather than spend half a day placing a point on one that carries
no signal.

## Two extrapolations to the scale six would need — and they disagree by eighty

**Through degree.** The four-colour ratio rises roughly exponentially above
degree eight — 5.0 at 8.11 and 24.0 at 9.94, a slope of 0.857 per degree. At
five colours the curve is flat at one up to the highest degree reached here,
10.64. If the five-colour curve shares that slope and simply starts later,
reaching a ratio of 24 takes **degree ≈ 14.3**. Average degree grows very
slowly with size in this family — 9.96 at `G`'s 1581 points, 10.64 at `G*`'s
13 873, 11.36 at `Z`'s 27 673, fitting `9.96 + 0.49 ln(n/1581)` — which puts
degree 14.3 at about **1.3 × 10⁷ points**.

**Through the forcer.** The rhombus forces at three colours on 4 vertices and
`Y` at four on 791 — a factor of 198 — and one more such factor gives about
**157 000**.

They disagree by a factor of **eighty**, and that is the honest headline. Each
is an extrapolation from two or three points; one assumes the five-colour
curve shares the four-colour slope, the other that the forcer keeps growing by
the same factor, and nothing supports either. Quoting one figure alone would
be quoting the assumption rather than the data.

What they agree on is the only part worth having: the scale is somewhere
between a hundred thousand and ten million points. Everything built in this
session tops out at 55 345.

## The correlation ratio is not a proximity-to-forcing metric

Breaking the degree ceiling looked like the way out. `Sa`'s family tops out
near degree 11, and Erdős's lattice — the integer grid scaled by `1/√r`, where
two points are at distance 1 exactly when `dx² + dy² = r`, so the degree is
`r₂(r)` and is chosen by arithmetic — reaches 17.6 easily. At five colours it
scores **8.7, 12.7, 26.4** where `G` at comparable degree scores 1.1.

They are all **2-chromatic**, and not by accident.

> `dx² + dy² = r` forces `dx + dy ≡ r (mod 2)`. For **odd `r`** exactly one of
> `dx`, `dy` is odd, so every edge flips the parity of `i + j`. For
> **`r ≡ 2 (mod 4)`** both are odd, so every edge flips the parity of `i`. And
> `r ≡ 0 (mod 4)` is a scaled copy. A parity is a proper 2-colouring every
> time — **the integer lattice at a single squared distance is always
> bipartite**, however large its degree.

Confirmed on eight values of `r` up to degree 19.5 and 25 281 points: `χ = 2`
throughout, including the four with even-parity representations that looked
like they should break it.

So a graph can score **26.4** at five colours and be incapable of forcing
anything whatsoever. The ratio measures *local correlation*, which rises with
degree in any graph at all; it does **not** measure proximity to forcing.

> **This qualifies much of what is above.** The ladder, the degree curve and
> the gadget table all describe that weaker quantity, and the two were run
> together without being separated. What is unaffected is everything verified
> rather than sampled: `G`'s colour relation at five colours is exactly its
> edge set — all 6502 candidates put to the solver, every one free — and `Y`
> remains the only graph here with a forced-*same* pair, which is what the
> spindle actually needs.

## The exact test, and what it retires

Forced-same is not a quantity to estimate. The pair `(i, j)` takes one colour in
every proper `k`-colouring exactly when the graph with the edge `(i, j)` added is
not `k`-colourable, and colour symmetry collapses the `k(k-1)` ordered colour
pairs to a single normal form: assume `c(i) = 0` and `c(j) != 0`, and read the
answer off one SAT call. `hn.homcol.forced_same` is that call.

> **Corrected.** The preceding section ranks graphs by their best non-edge pair
> over thirty-two sampled colourings, and a hill-climb was built on that ranking:
> scan unions `G u rho(G)`, take the maximum over every non-edge pair, watch it
> climb. It did climb — 21, 22, 23, 24, 25 over the first fifty unions — and the
> climb is free. The null control is two copies of `G` translated a thousand
> apart, so **not one edge crosses between them**: same vertex count, same clause
> count, same search, zero coupling. It reaches 23 and 24, against 21 and 22 for
> `G` alone. Twice the vertices is four times the pairs, and the maximum of four
> times as many `Binomial(32, p)` draws is higher for nothing; a solver under 5%
> phase randomisation also returns less diverse colourings on a bigger instance,
> lifting every pair at once. The winning union had **one** cross edge and its
> winning pair lay entirely inside the original copy, which was the tell. The
> graded metric comparing equal-sized graphs at equal sample counts stands; the
> ranking across sizes does not.

The exact test is also faster than what it replaces, by three orders of
magnitude. `G`'s entire closable-pair census — every non-edge pair whose squared
distance is rational *and* whose spindle rotation the field admits — is 21344
pairs, and it runs in twenty-one seconds:

| | |
|---|---|
| pairs at a rational squared distance | 44341 |
| distinct rational squared distances | 106 |
| closable ones | 36 |
| pairs they carry | 21344 |
| **forced at five colours** | **0** |

Complete, not sampled. And forcing survives adding vertices, so no subgraph of
`G` carries a forced closable pair either.

## de Grey's forced pair, identified

Everything above had guessed at what his forcer forces. The exact test made
asking cheap, and the answer is what the ring criterion predicted, to the letter.

`Sa` is the twelve-element dihedral closure of `S` about the origin, and it
carries a ring of squared radius 4. The rotation that makes that ring bite itself
is `cos 7/8, sin sqrt(15)/8` — it wants `sqrt(4D - 1) = sqrt(15)`. `Sb` is its
image, and `Y = Sa u Sb` less two points forces the ring's **antipodal** pair
`(-2, 0), (2, 0)`, at squared distance `4D = 16`, whose spindle rotation wants
`sqrt(16D - 1) = sqrt(63) = 3 sqrt(7)`. Both radicals are in the field. That is
the ring paying twice, and it is why `D = 4` and nothing else.

| | |
|---|---|
| `Y` at four colours | **forced same** (UNSAT, 210 s) |
| `Y` at five colours | free — the control |
| `Sa` at four, before the bite | free — the rotation makes it |
| `Y`'s vertices at squared distance 16 from `(-2,0)` | 1 |

The spindle then turns about **one end** of that pair, and `build_G`'s pivot is
`(-2, 0)`: an end, not a midpoint. Its two rotations are `pi/2 -+ arcsin(1/8)`,
and their composition has `cos 31/32, sin sqrt(63)/32` — the `D = 16` spindle
rotation exactly. The two symmetric turns are one spindle written evenly. That
identity is checked in the tests, not asserted here.

## The bottleneck is arithmetic, not scale

The template wants a ring that pays twice. De Grey's field admits twenty-two such
`D` among the rationals with numerator and denominator at most forty — `2/7`,
`2/5`, `4/9`, `1/2`, `8/7`, `4`, `17/2` and more.

`G` populates one of them. Scanning every vertex of `G` and four thousand
midpoints as the centre, the richest doubly-usable ring anywhere in `G` holds
**twelve points** — which is exactly the ring `Sa` already had. Four times the
size bought no extra coupling, so the bite at five colours runs with the same
twelve cross edges that sufficed at four, against a harder problem. Scale was
never the binding constraint. Ring population is.

Two levers raise it. The first is the bite itself: `|1 - rho| = 1/2` for `D = 4`,
so a radius-2 point and its image are one apart, and so are that image and *its*
image. Iterating `rho` threads a path along the circle, one per point of the
orbit already there, and `arccos(7/8)` is not a rational part of a turn, so the
paths never close. Centred at `G[0]`, where `G`'s ring actually is, the ring goes
12, 24, 36, 48 and the antipodal pairs 12, 18, 24 — at 2372 extra vertices a
step. Centred at the origin it goes 1, 2, 3, 4, 5, 6, because `G` has a single
vertex at radius 2 from the origin. The centre is not a detail.

The second lever is free. **Adjoining a radical costs no vertices**: the points
of `G` do not move, and only the rotations need the new square root. Censusing
every rational ring at every vertex centre and asking which single adjunction
turns the biggest one usable:

| adjoin | unlocks | points on the ring |
|---|---|---|
| `sqrt(17)` | `D = 5/3` | **48** |
| `sqrt(13)` | `D = 1/3` | 36 |
| `sqrt(71)` | `D = 5/9` | 24 |
| — (de Grey's own field) | `D = 4` | 12 |

`D = 5/3` is the one to want, and its arithmetic is clean: `4D - 1 = 17/3`, whose
squarefree part is `51 = 3 * 17`, so it is the **bite** that wants the new
generator; `16D - 1 = 77/3`, squarefree part `231 = 3 * 7 * 11`, already there,
so the **spindle is free**. Four times de Grey's coupling for one generator and
not one extra point. Built over `Q(sqrt3, sqrt5, sqrt7, sqrt11, sqrt17)` about
`G[0]`, the bite is `cos 7/10, sin sqrt(51)/10`, ring points and their images
come out at distance exactly one as they must, and one step gives 4741 points
with **72 antipodal pairs** against de Grey's 6.

None of them is forced yet, and every union so far is still 5-colourable. What
has changed is that the mechanism is now measured rather than guessed at, and the
quantity it turns on has a name and a lever.

## The gateway: a weaker property, and one call to test it

Everything above tests the **strong** hypothesis — a *named* pair, the same
colour in every proper colouring — because that is what the spindling lemma
consumes. De Grey's construction does not start there.

| stage | statement | where it holds |
|---|---|---|
| weak | in every 4-colouring of `Sa`, at least **one of three** antipodal pairs on the `D = 4` ring is monochromatic | `Sa` **alone** — no bite, no union |
| strong | the bite `rho_4` sharpens "one of three" to the named pair `(-2,0),(2,0)` | `Y = Sa u rho_4(Sa)` |
| — | spindle at squared distance 16 | `G` |

The weak stage costs **one SAT call**, because forbidding a pair from being
monochromatic is exactly adding the edge. So "some antipodal pair of this ring
is always monochromatic" is "the graph with all those antipodal edges added is
uncolourable" — one call for the whole ring, not one per pair. `Sa` carries it
on two rings: `D = 1` (fifteen pairs, spindle distance 4) and `D = 4` (three
pairs, spindle distance 16, de Grey's own).

Strong implies weak, so a graph failing the weak property cannot carry a forced
pair on that ring at all. That makes it the cheap gateway every scan here should
have started from — and it settles `G` in under a minute:

| | |
|---|---|
| centres scanned (every vertex of `G`) | 1581 |
| `(centre, ring)` candidates with ≥2 antipodal pairs | 1664 |
| of those, with a closable spindle to follow | 688 |
| **carrying the weak property at five colours** | **0** |
| cost | 55 seconds, complete |

This subsumes the negatives above rather than adding to them. The exhaustive
closable-pair census, the bitten unions, the thickened rings — all were hunting
the strong property, which cannot hold where the weak one does not. `G` is not
a graph whose bite can be sharpened; it is a graph with nothing to sharpen.

What it does not touch is other graphs. The weak property belongs to a graph,
and `G` is one graph.

## The by-distance gateway, and where G's edge is

Antipodality was de Grey's *symmetric* choice, not the lemma's requirement. The
weak property needs a set of non-adjacent pairs all at the **same distance**,
such that forbidding every one kills colourability — they share the distance, so
they share the spindle rotation's angle; only its centre moves from pair to
pair. Dropping antipodality widens the family from 119 antipodal pairs over 8
rings to 21344 pairs over 36 closable distance classes.

`Sa` carries it on three of its fourteen closable classes — `D = 4/9` (393
pairs), `D = 4` (273), `D = 16/9` (156) — and size is not what decides: those sit
5th, 6th and 9th of fourteen, while the four biggest classes all colour. De Grey
took `D = 4` and then only the three **antipodal** pairs out of its 273, the
tightest of the three statements. His lemma has no slack: drop any one of the
three and `Sa` colours again.

The same test on `G` at five colours is enormously more expensive, and the
expense is the measurement. `G` clears its four biggest classes in seconds each
and then stops dead on the fifth — which is also the *smallest* of the five:

| class | pairs forbidden | cost |
|---|---|---|
| `D = 1/3` | 6510 | instant |
| `D = 7/3` | 3648 | instant |
| `D = 3` | 3216 | instant |
| `D = 5/9` | 2448 | 20 s |
| **`D = 4/9`** | **1558** | **over an hour, on four solvers, unresolved** |

Balls about the hub locate it. Every one of them **colours**, so the full graph
almost certainly does too — what the curve says is where `G` stands:

| ball | cost |
|---|---|
| 1200 pts | 21 s |
| 1400 pts | 72 s |
| 1500 pts | 1279 s |
| 1540 pts | 4969 s |
| 1560 pts | 6893 s — with **1548 of the 1558** pairs forbidden, and it still colours |

> **Settled, and negative.** A calibrated TabuCol reached zero conflicts in
> **151 s** where the four CDCL solvers had run for over three hours between
> them, and the colouring was verified from scratch against the exact geometry
> rather than the search's own bookkeeping: 1581 of 1581 vertices, five colours,
> **0** of 7877 unit edges monochromatic, **0** of 1558 pairs at `4/9`
> monochromatic, and every edge and pair re-derived exactly. So `G` does **not**
> carry the weak property there. The same tool settled `Y`'s `4/9` — 577 s of
> cadical — in seconds. CDCL explodes on both sides of a phase transition;
> local search is asymmetric, landing on satisfiable instances and plateauing on
> unsatisfiable ones.

So the verdict is negative and the reading is not the verdict. `G` is nowhere
near the property on nine of its closable classes and right at the edge on one. That one is `D = 4/9`, which is **doubly usable with no adjunction at
all** — `4D - 1 = 7/9` wants `sqrt(7)` for the bite, `16D - 1 = 55/9` wants
`sqrt(55) = sqrt(5) sqrt(11)` for the spindle, both already in de Grey's field —
and it is one of the three classes `Sa` carries at four.

Solver difficulty is not a principled distance to the property, and is not
offered as one. It is a real asymmetry across classes of the same graph at the
same size, and it points at one class rather than the others.

## A correction worth its own section: the gap is 10.8, not 700

An earlier section compares `Sa`'s 1.45 spindles per point against `G`'s 0.002
5-critical subgraphs per point and calls the gap a factor of seven hundred.
Those are not comparable quantities. A spindle is 7 vertices and a 5-critical
subgraph is about 500, so "gadgets per point" measures different things on the
two sides. What a saturation argument needs is how much of the graph the gadgets
**cover** — gadget-vertex incidences per point:

| | incidences per point |
|---|---|
| `Sa`, 576 spindles × 7 / 397 | 10.16 |
| `Y`, 1152 × 7 / 791 | 10.19 |
| `G`, 2304 × 7 / 1581 | 10.20 |
| `G`, 3 five-critical × 500 / 1581 | **0.95** |

The gap is **10.8**, and `764 = 71.1 × 10.8` where `71.1 = 500/7` is the gadget
size ratio — which the section above already accounts for separately. Counting
it again inside the density comparison overstates the gap sixty-five fold.

What does **not** change is the conclusion drawn from it. Density is inherited,
never raised: 1.45, 1.46, 1.46 spindles per point across `Sa`, `Y`, `G`, and
4.972, 4.979, 4.982 edges per vertex, with both bites of `G` landing at 4.987 —
union and rotation preserve density to two parts in a thousand. Neither closes a
gap of 10.8 any more than one of 764. Only the size of the thing to be closed
changes, and it changes by a lot.

## What each operation can and cannot do

Three moves are used throughout: **bite** (union with one rotated copy at the
angle that makes a ring touch its image), **thicken** (iterate the bite), and
**close** (the twelve-element dihedral closure about a centre). Their effect on
edge density separates cleanly, and sharply:

| | points | edges | per vertex |
|---|---|---|---|
| `Sa` | 397 | 1974 | 4.972 |
| `Y = Sa u rho_4(Sa)` | 791 | 3938 | 4.979 |
| `G` | 1581 | 7877 | 4.982 |
| `G` bitten at `D = 4`, both ways | 3953 | 19715 | 4.987 |
| `G` bitten at `D = 4/9`, both ways | 4741 | 23643 | 4.987 |
| `G` bitten at `D = 5/3` (needs `sqrt17`) | 4741 | 23727 | 5.005 |
| `G` bitten at `5/3`, `1/3` and `4` at once | 8693 | 43443 | 4.997 |
| `G` with the ring thickened, `m = 9` | 22929 | 114419 | 4.990 |
| `G` bitten at `D = 5/3`, `m = 4` | 14221 | 71277 | 5.012 |
| `G` closed about the origin, then bitten | 37932 | 189114 | 4.986 |
| `G` closed about its hub | 11047 | 59919 | **5.424** |
| `G*` about `(-2,0)` ∪ `rho_16(G*)` | 27673 | 157140 | **5.678** |

Bites and thickenings land between 4.98 and 5.01 without exception — twelve
thousand points added at `m = 9` and the density moves by eight parts in ten
thousand. Only the dihedral closure raises it at all, and by fourteen per cent.

The sharper point is that density is not the discriminator anyway. `Sa` sits at
4.972 and **carries** the weak property at four colours; `G` sits at 4.982 and
does not at five. What separates them is 5-critical coverage — 10.16
gadget-vertex incidences per point against 0.95 — and no operation available
here touches that.

Which gives the structural conclusion these measurements converge on: **the
bite sharpens a weak property into a named pair, it does not create one.** `Sa`
carries it alone and `Y` names it. `G` lacks it, and every union, thickening and
closure of `G` lacks it too. So the target is not a bigger union of `G`. It is a
five-chromatic graph carrying the weak property *on its own* — and nothing in
this repository produces one.

## `D = 4/9` is this family's distance, at four colours and at five

Sorting the closable classes by how the weak-property call behaves separates one
of them from all the rest — and the same one every time:

| graph | colours | class | verdict |
|---|---|---|---|
| `Sa` | 4 | `4/9`, 393 pairs | **does not colour** — carries it, 8 s |
| `Sa u Sb` | 4 | `4/9`, 786 pairs | **does not colour** — carries it, 8 s |
| `Y` | 4 | `4/9`, 779 pairs | **does not colour** — carries it, 8 s |
| `Sa` | 5 | `4/9`, 393 pairs | colours, in 10 s — the other thirteen classes take none |
| `Sb` | 5 | `4/9` | the same |
| `G` | 5 | `4/9`, 1558 pairs | over an hour against four solvers — the only one of thirty-six to resist |

`Y` at five is the same story one size up: 0 of its 23 closable classes carry
the property, `4/9` takes **577 s**, the next slowest (`5/9`) takes 10, and the
rest take none. The dihedral closure about the hub — 11047 points, the densest
object here — clears all 34 of its small and middle classes in five seconds or
less and then takes **372 s** on `16/9`.

And `build_Y` deletes exactly two points from `Sa u Sb`: `(1/3, 0)` and
`(-1/3, 0)`, whose squared distance is `4/9`. They are the antipodal pair of the
`D = 1/9` ring about the origin. Whether that is *why* he deleted them is his
business; what is measurable is that it costs nothing — `Y` still carries the
property on that class afterwards, with 779 pairs instead of 786.

> **Corrected.** An earlier reading of this listed three carrying classes and
> called their distances an arithmetic progression of step `2/3`. There are
> **four**, and the fourth breaks it. The complete list for `Sa` at four
> colours:
>
> | class | pairs | distance |
> |---|---|---|
> | `D = 4/9` | 393 | `2/3` |
> | `D = 4` | 273 | `2` |
> | `D = 16/9` | 156 | `4/3` |
> | `D = 16` | **3** | `4` |
>
> The fourth is de Grey's own — three pairs, which are exactly the antipodal
> pairs of the `D = 4` ring — and it is the *smallest* carrying class of the
> four, which is why it is the one worth building on. The distances are `2/3`
> times `1, 2, 3, 6`, not an arithmetic progression.

Solver time is not a metric and is not offered as one. What is stated is that
one class behaves unlike the other thirteen, or thirty-five, on five different
graphs — and that it is the class that carries the property at four colours.

## The family is empty at five colours, and the last heuristic goes with it

`Sa` carries the weak property at four colours on four closable classes —
`4/9`, `16/9`, `4`, `16`. At five it carries **none of them, singly or
together**: glucose found a 5-colouring of `Sa` avoiding all four at once, in
850 s, with the colour symmetry broken 60-fold by a triangle. The two-class
statement needs no run of its own — `{4/9, 16/9}` is a subset of the four, so a
colouring avoiding the larger set avoids the smaller.

That settles the whole family at five colours:

| | |
|---|---|
| `Sa` | 0 of 14 classes, and 0 for all four together |
| `Sb` | 0 of 14 |
| `Y` | 0 of 23 |
| `G` | 0 of 1664 `(centre, ring)` candidates; class `4/9` satisfiable |
| `G` closed about its hub | 0 of the 35 classes reached |

> **Retired.** Sixty-seven runs of the calibrated TabuCol on `Sa` with
> `{4/9, 16/9}` — sixty at 600k moves and seven at three million — never reached
> zero, and touched **one** conflict twenty-four times. Against a
> known-unsatisfiable calibration that plateaus at 37, that reads as barely
> unsatisfiable. It was a hard **satisfiable** instance.
>
> The same thing had already happened on `G`'s class `4/9`: three hours of CDCL,
> a TabuCol plateau at seven, answer SAT. Twice is a pattern. **A local search
> that does not land is not evidence of anything** — landing proves
> satisfiability, and nothing else proves anything. Every reading here was
> hedged as evidence rather than result, which was the right hedge; the hedge is
> now the finding.

What survives is the structural account, which never rested on a heuristic: the
bite sharpens a weak property into a named pair rather than creating one, `Sa`
carries one at four and `G` carries none at five, and no union, thickening or
closure of `G` changes that.

## Leaving the family: two-distance lattices, and why they fail

Everything above descends from one seed. `S` gives `Sa` gives `Y` gives `G`, and
the closures, bites and thickenings of `G` are still `G`. The family carries the
weak property at four colours and nothing at five, and the bite cannot create
one — so the search has to leave the family, not grow it.

The mechanism does not require the family. The weak property is exactly "the
graph on this point set with **both** distance 1 and distance `d` joined is not
`k`-colourable", and the point set is free. Single-distance lattices were ruled
out here earlier and correctly — they are always bipartite — but that argument
says nothing about two distances, where the parity classes interleave and the
bipartition dies.

Forty-five instances: triangular patches of 127, 301 and 517 points with second
distance squared in `{3, 4, 7, 9, 16, 19, 25, 37}`, square patches of 113, 253
and 441 with `{2, 4, 9, 16, 25, 34, 37}`, every second distance closable over de
Grey's own field. **Every one 5-colourable, every one in under a second.**

The reason is structural: a two-distance lattice graph is a **Cayley graph of
`Z^2`**. Both edge sets are translation-invariant, so a periodic colouring by
cosets of a sublattice avoiding both generator sets works, and five colours
leave far too much room. Homogeneity is exactly what makes them easy.

Which says where the structure has to come from. `Sa` is inhomogeneous: it is a
dihedral closure under a rotation of **infinite order**, so its points sit at
many scales — radii `1`, `1/3`, `5/9`, `4/3`, `5/3` — and no translation
preserves it. A lattice has one scale and order-six rotations.

## Sumsets of unit vectors: the first lever that moves density

Density is what the structural account says matters, and nothing in the family
moves it — `Sa`, `Y` and `G` all sit at 4.98 edges per vertex, and bites,
thickenings and closures leave it there. The bitten lattices came out at 2.5 to
3.1, which is why they coloured instantly. So build for density directly.

Take unit vectors `v_1 … v_m` of the field and form every integer combination
`Σ a_i v_i` with `|a_i| ≤ c` inside a radius. Each point has an edge to every
point differing by one generator, so the degree is `2m` before any coincidence,
and every algebraic relation between the generators folds the set onto itself
and adds more. Generic generators give a hypercube, which is bipartite; the
generators used here are `rho_D` for `D = 2, 3, 4, 7, 9, 14, 16`, all of
**infinite order**, so the point set keeps the many scales the lattices lacked.

| generators | points | edges per vertex |
|---|---|---|
| `1, w` | 49 | 2.45 |
| `1, w, r4` | 323 | 3.84 |
| `1, w, r7` | 323 | 4.43 |
| `1, w, r4, r7` | 593 | 4.72 |
| `1, w, r3, r4, r7` | 2819 | 5.47 |
| `1, w, r2, r3, r4, r7` | 13483 | 6.20 |
| `1, w, r3, c3, r4, c4` | 13639 | 6.25 |
| `1, w, r3, r4, r7, r9` | 13505 | 6.37 |
| **`1, w, r3, r4, c4, r7`** | **13579** | **6.59** |

against 4.97 for `Sa`, 4.98 for `Y` and `G`, 5.42 for the densest dihedral
closure and 5.68 for `G*` union its spindle. Five generators already beat
everything in the family, and the trend has not turned over.

Every one of them is 5-colourable, so this is a **lever, not a result**. What it
is, is the first lever in this work that moves the quantity the structural
account points at, rather than the size.

## Corrected: density was not the missing ingredient

The sumset programme above rests on an inference, and the inference is wrong.
The observation was real — `Sa` and `G` both sit at 4.98 edges per vertex and no
operation in the family moves it. The step from there to *build for density* is
what fails, and there was already a result in this repository cutting against
it: **single-distance lattices are bipartite**. The densest unit-distance graphs
known are lattice-like, and lattice-like is 2-chromatic, so pushing density
pushes toward the achromatic end.

Measured, with chromatic numbers computed by SAT:

| edges per vertex | χ | graph |
|---|---|---|
| **5.49** | **3** | walk `rho_4 e≤2`, 7000 points |
| 4.98 | 4 | `Y`, 791 |
| 4.97 | 4 | `Sa`, 397 |
| 4.91 | 4 | walk `rho_3, rho_4`, 7000 |
| 3.57 | 4 | walk `rho_3, rho_4, rho_7`, 7000 |
| 2.84 | 3 | triangular lattice, 517 |

The **densest** graph built here, at seven thousand points, is **3-chromatic** —
below `Sa`, which reaches four with 397 points and less density. Density and
chromatic number are decoupled, and the densest entry is nearly the least
chromatic.

The synthesis meant to rescue it fails too. Closing a dense sumset patch under
the dihedral group — density from one direction, inhomogeneity from the other,
which is exactly what `Sa` is with a 39-point seed — gives 499, 1339 and 2587
points at 4.32 to 4.56 per vertex, and **all three are 4-colourable**. A denser
seed than de Grey's produces a *less* chromatic closure.

So what `Sa` is doing is not density. It is the **seed**: de Grey's 39 points are
chosen, and choosing a seed for density makes the closure worse. That is the
honest residue of this stretch — the measurements stand, the framing does not.

## The seed can be improved, and it buys nothing at five

Greedy growth from `Sa`, keeping the D6 symmetry, runs 397 points at 4.97 edges
per vertex out to **853 at 6.77** — the densest symmetric graph here, past the
whole family at 4.98 and past its densest dihedral closure at 5.42 — with χ
staying at 4 throughout. Density alone has already been shown not to matter, so
the question is what growth does to the **weak property**, and either answer was
worth having: destroying it would have said the seed's power is fragile and
hand-chosen.

It does not destroy it. The grown seed at 541 points carries the property at
four colours on **five** closable classes where `Sa` carries four:

| class | pairs | |
|---|---|---|
| `D = 16` | **3** | de Grey's own, the tightest statement — **unchanged** |
| `D = 20/3` | 102 | **new** — `Sa` does not carry this class |
| `D = 16/9` | 228 | |
| `D = 4/9` | 513 | |
| `D = 4` | 699 | |

So a constructed seed beats a hand-chosen one, keeping everything `Sa` has
including the three-pair statement de Grey built on. That is the first
construction here to improve the **seed** rather than rearrange what the seed
produces.

**And it buys nothing at five.** Measured two independent ways: CDCL finds 0 of
the 10 smallest closable classes carrying the property, and the calibrated
TabuCol **lands** on seven of the fourteen classes in under a second each —
landing proves colourable, which is the one direction local search can prove.
`Sa` carries nothing at five and neither does the grown seed. The bitten grown
seed (`rho_4`, 1081 points, 6390 edges) is 4-colourable too.

Which decouples the two levels. The seed can be improved at four colours,
measurably, while keeping de Grey's own statement intact — and whatever the
barrier at five is, it is **not** a deficiency of the seed that better
seed-building fixes.

## Every neighbourhood is bipartite, which bounds the single-point attack

The sharpest question about a new point is exact: `G + p` is not 5-colourable
precisely when, in **every** proper 5-colouring of `G`, the neighbours of `p`
already use all five colours. Colour symmetry collapses it to one SAT call with
assumptions — if any colour is free on `N(p)`, colour 0 is free in some
colouring.

Asked exhaustively rather than by sampling — a point with `k` neighbours is the
intersection of `C(k,2)` pairs of unit circles, so counting multiplicities finds
the high-concurrence points that random pair-sampling throws away — `G` has
34424 distinct new intersection points, the largest neighbourhood is **13**, and
every one is placeable.

> **Caught.** The first run reported a 60-neighbour point. Vertices *already in
> the graph* appear as intersections of their own neighbours, and `G`'s hub has
> degree 60, so it is the intersection of `C(60,2) = 1770` pairs. Such a point is
> trivially placeable — its own colour is free — and is not a candidate at all.

Graded rather than yes/no, it is far worse than "no". The minimum, over proper
5-colourings, of the number of colours on `N(p)`:

| | |
|---|---|
| `G`, 400 top candidates | **all minimum 2** |
| `Y`, 400 | all minimum 2 |
| `Sa`, 400 | six at 1, the rest at 2 |

including all twenty-four candidates with thirteen neighbours. A blocked point is
not one colour away. It is **three**.

**And the reason is a theorem.** `N(p)` lies on the unit circle about `p`, and
two of its points are adjacent exactly when they subtend 60°. A cycle there is a
sequence of ±60° steps returning to its start, so with `a` steps of `+60` and `b`
of `−60` it needs `60(a − b) ≡ 0 (mod 360)`, that is `6 | (a − b)`. Since `a + b`
shares parity with `a − b`, the cycle length is **even**. No triangles, no odd
cycles:

> **The neighbourhood of any point in a planar unit-distance graph is bipartite**,
> so `χ(N(p)) = 2` for every `p`, always.

That is what the measurement sees. To block `p`, the ambient graph must force
five colours onto a set that is 2-chromatic on its own — and the ceiling of
thirteen was never the binding constraint, since a five-point neighbourhood
forced rainbow would win just as well, and nothing here is forced past two.

> **Refined.** "The forcing is global with no local witness" is true of `G` at
> five and **wrong in general**. At four colours, `Sa` has 187 of its 397
> vertices forced to **3 of 4** — the hub included — against `G`'s 1581 forced
> to 2 of 5:
>
> | | forced | available | ratio |
> |---|---|---|---|
> | `Sa` at four | **3** (187 vertices) | 4 | 0.75 |
> | `G` at five | 2 (all 1581) | 5 | 0.40 |
>
> So the lens was never useless: it discriminates sharply between the level
> where the mechanism works and the level where it does not. `Sa` is one colour
> from rainbow at nearly half its vertices; `G` is three colours away at all of
> them. The bipartite theorem is untouched — `χ(N(p)) = 2` always, and what
> varies is how many colours the *ambient* graph forces onto it. And it agrees
> with the class count from the other side: 1 class and 3 pairs at four, 5
> classes and 229 pairs at five.

## G does carry a weak property at five colours — over five distances

Single classes were tested and all of them colour. `Sa`'s four carrying classes
were tested together and they colour too. The **graded** question was never
asked of `G`: forbid its closable distance classes *cumulatively* and count how
many it takes. Twenty-seven, added smallest first. Minimised by dropping the
largest first and keeping every drop that still fails to colour, twenty-seven
comes down to **five**, and they are irreducible:

| class | pairs |
|---|---|
| `D = 15/16` | **1** |
| `D = 16` | 12 |
| `D = 17/2` | 24 |
| `D = 9` | 48 |
| `D = 7` | 144 |
| | **229** |

> **In every proper 5-colouring of `G`, at least one of those 229 pairs is
> monochromatic.**

That is the first positive statement about `G` at five colours here, and it is
not a density artefact: 229 pairs is **2.9%** of `G`'s 7877 edges, and the mean
degree moves from 10.0 to 10.3.

Verified rather than reported. Three solvers — cadical, glucose, minisat — all
return UNSAT on the 229-pair formula. Each of the five classes, dropped in turn,
lets the graph colour again, so every one is needed; `D = 15/16` contributes a
**single pair** and is still necessary. All 229 distances were re-derived from
the exact coordinates and all five are closable, so the statement is about the
plane and not about an indexing slip.

What it is **not** is de Grey's hypothesis. His is **one** distance and **three**
pairs, which two rotated copies consume. Five distances and 229 pairs would need
a multispindle over five rotation angles and 229 centres, and nothing here builds
one. But the gap between the levels finally has a figure attached:

| | classes | pairs |
|---|---|---|
| `Sa` at four | 1 | 3 |
| `G` at five | 5 | 229 |

## A spindle that survives a disjunction — and its ceiling

The ordinary spindle needs a *named* forced pair. Every disjunction found here
refuses it, because each pair carries its own centre of rotation and one
rotation cannot serve them all. Making the pairs **share an endpoint** removes
that objection: a hub `u`, a ring about it, and a set `W` on that ring with `u`
monochromatic with *some* point of `W` in every proper `k`-colouring. The
rotation about `u` fixes `u`, so every rotated copy repeats the disjunction
about the same hub over a shifted `W`, and the colour class of `u` must hold an
**independent transversal** inside one circle.

Write ring points as exponents of the rotation; two are at unit distance
exactly when their exponents differ by one. For `W = {0, −2}` and three copies:

    a0 = 0   ->  a1 in {1, -1}, both one away.            dead
    a0 = -2  ->  a1 = 1  ->  a2 in {2, 0}, both one away.  dead

Exhaustive search settles the general case: the stack closes **iff the
separation `m` is even**, and then at exactly `m + 1` copies. The hand argument
had claimed only `m = 2`; it checked consecutive copies and missed that the
non-consecutive constraints bite as well.

The separation is a distance, and rational for every `m` — carrying the chord
recurrence as a rational pair keeps it exact, since one of the two parts always
vanishes. So the hypothesis is searchable. Searched, and absent: Sa at four and
`G` at five offer 1962 and 7800 candidates and not one carries it, and widening
`m` from 2 to {2,4,6,8} added **no candidate at all** — the only chords either
graph realises belong to the ring `D = 1/3`, whose rotation has order six so
that its two-, four- and eight-step chords collapse onto the single value one.

And the ceiling is hard. Over every shape inside a window of thirteen exponents
and up to twenty-six copies — 66 shapes of width three, 220 of width four —
**not one closes**. Only widths one and two do, and width two is barely weaker
than the forced pair it was meant to replace.

## de Grey's ring lemma is a palette cap, not an antipodal pair

The six points of Sa's `D = 4` ring carry no edges among themselves — adjacent
ones are two apart, antipodal ones four — so all 187 partitions of six things
into at most four blocks are a priori available. Forcing the ring to each in
turn says exactly what Sa claims, with no paraphrase in the way:

| | patterns available | surviving | colours used |
|---|---|---|---|
| `k = 4` | 187 | **10** | one 1-block, nine 2-block |
| `k = 5` | 202 | **202** | up to and including all five |

So at four colours the ring takes **at most two colours**, and in every
two-block pattern the minority class is a pair at chord 2 or chord 4 — never
one of the six at chord 2√3. "Some antipodal pair is monochromatic" is a
corollary, and a much weaker one.

That explains how the bite works on so little contact. `Y` is Sa with one
rotated copy, and the two share **one vertex and exactly six edges** — the
matching that carries each ring point to its image, which the rotation places
at distance one. Thirteen stacked copies share 72. Only a statement as strong
as the cap survives that little contact.

## The cap is the design target, and nothing has it at five

The cap screens in one SAT call: can the ring show all `k` colours at once?
Swept over every centre and ring of `G` at five — 3943 rings of six points or
more — 1557 came back capped, and **every one of them is a unit ring capped at
4**. That is not rigidity. A neighbourhood showing all five colours would leave
its centre uncolourable, so "the unit ring is capped at `k − 1`" is nothing but
the statement that the graph colours. The degree-60 hub is capped for that
reason and no other; no ring at any other distance is capped at all.

Beside it, two more censuses, both empty:

| | result |
|---|---|
| forced-same pairs | 0 of 21344 closable non-edge pairs |
| forced-**different** pairs | 0 of 39923 non-adjacent pairs at `d² < 4` |
| apex centres | 0 — no neighbourhood forces five colours |

`G` has no rigidity of any kind at five colours.

Thickening does not help, and the reason is structural rather than empirical.
Rotations about the origin commute with the sixty-degree rotation, so every
image `ρᵗ(Sa)` is dihedrally symmetric and the stack keeps the whole structure
while the ring grows six points a level. At four colours `Sa` alone already
carries the gateway — the calibration reproduces de Grey exactly. At five it
fails at every certified level, and it must: ring points are adjacent only
between consecutive levels at the same hexagon position, so the ring plus its
antipodal pairs is **three disjoint ladders**, bipartite at any thickness.

> **A bug worth stating.** Each bite multiplies the shared denominator by
> eight, so thirteen of them push it past `8¹³` and the int64 squared form
> **wraps silently**. `fast_edges_complete` returned a graph with *no edges at
> all* and the sweep went on reporting "colourable, gateway fails" for three
> more levels — true of the empty graph and vacuous about the real one.
> `IntBasis.overflow_headroom` existed for exactly this and no script here was
> calling it. It is conservative, so it certifies rather than detects: levels 0
> through 8 are certified and the rest are not reported.


## The bite, derived by hand

The census taken on the ring alone gives "at most two colours". Taken on the
**centre and the ring together** it gives much more: at four colours **ten of
715** patterns survive, the seven points take at most two colours between
them, and the centre is **never alone**. Writing `A` for the ring points that
share the centre's colour, the ten sort into three shapes:

| shape of `A` | how many |
|---|---|
| all six ring points | 1 |
| the complement of an **adjacent** ring pair | 6 |
| an **antipodal** pair | 3 |

The bite fixes the centre and adds six edges, `h_j — ρ(h_j)`. So a 4-colouring
of `Sa ∪ ρ(Sa)` picks one pattern per copy, the two agree on the centre's
colour, and corresponding ring points must differ — which says exactly that
`A` and `B` are **disjoint**. That single condition does all the work:

    A = all six                 ->  B empty                       impossible
    A = complement of {j, j+1}  ->  B inside an adjacent pair      impossible
    A = antipodal {j, j+3}      ->  B = {j+1,j+4} or {j+2,j+5}     the only case

**Six of a hundred** shape pairs survive, every one of the last kind. And in
that case the remaining four ring points all carry the second colour, so the
other two antipodal pairs are monochromatic as well — **all three are**.

That is de Grey's forced pair, and naming one of them is a convenience rather
than a fact about the union: rotating by sixty degrees carries `Sa ∪ ρ(Sa)` to
itself and permutes the three pairs cyclically, so no single pair can be
singled out, and none needs to be.

At five colours the same census returns **855 of 855**. There is nothing to
eliminate, so the derivation has no input — which restates the target
exactly: *a graph whose joint centre-and-ring census at five colours is small,
with shapes closed under complement the way these are.*


## The pruning costs two of the three pairs, and the bite is a one-off

`Y` is `Sa ∪ ρ(Sa)` **less** the two vertices `(±1/3, 0)`. Removing vertices
can only *add* surviving patterns, so `Y` is strictly weaker than the union it
comes from — and the census says by exactly how much. Numbering the centre 0
and the ring by angle, so the antipodal pairs are `(1,4)`, `(2,5)`, `(3,6)`:

| | patterns | pairs forced | centre alone |
|---|---|---|---|
| `Sa ∪ ρ(Sa)`, 793 points | **3** | **all three** | 0 |
| `Y`, 791 points | 7 | one, namely `(3,6)` | 4 |

The three the union keeps are precisely those where the centre sits beside an
antipodal pair — exactly what the hand derivation predicts, since `A` and `B`
must be disjoint *and non-empty*. The four `Y` gains all have the centre
alone, and each loses one of the other two pairs. Confirmed both ways: seven
SAT calls on the union, and an `UNSAT` proof per pair.

So the pruning is a saving of two vertices that costs two thirds of the
conclusion. For a write-up that needs one forced pair, de Grey's choice is
right; for a search, the unpruned union is the stronger object.

And biting again buys nothing. The census tightens *across* the first bite —
`Sa` 10 of 715, `Y` 7 — which invites the obvious question, and seven SAT
calls answer it, because survivors only ever shrink. Biting `Y` takes it to
1187 points and the ring from twelve to eighteen, and **all seven patterns
survive**. The chain has exactly the levels de Grey used.


## There is no gradient

Every proxy built here reads nothing at five colours, in every graph, at every
size. The natural reading is that the graphs are far away and a better search
would close the distance. Peeling `Sa` and censusing as it shrinks says
otherwise.

| points | edges | `k = 4` | `k = 5` |
|---:|---:|---:|---:|
| 7 | 0 | 715 of 715 | 855 of 855 |
| 107 | 143 | 715 | 855 |
| 207 | 537 | 715 | 855 |
| 307 | 1145 | 715 | 855 |
| 347 | 1445 | **577** | 855 |
| 397 | 1974 | **10** | 855 |

At four colours the census sits at the **ceiling** until 307 points — 77 % of
`Sa` — and then falls off a cliff in the last ninety vertices. Nothing in the
first three hundred points hints at what the last ninety do. The same vertices
censused at five never leave the ceiling at all.

So a reading of 855 carries no information about distance. It is what an
incomplete object reads, and an object *one orbit short of complete* reads it
too — which is what orbit-irreducibility says from the other side: all 39 of
`Sa`'s dihedral orbits are needed, and none can go.

The consequence for method is sharp. Every search strategy that optimises a
proxy — climb the agreement, thicken the ring, buy contact, pack gadgets —
assumes the proxy improves as the object improves. **It does not.** There is
nothing to climb until the object is essentially complete, and by then there is
no climb left to do.

And it is not simply that four colours are scarce only near the end. Walking
the same peel and asking for the chromatic number refutes that: `χ` reaches
**4 at 207 points** and the census does not move for another 140 vertices.

| points | `χ` | census at four |
|---:|---:|---:|
| 157 | 3 | 715 of 715 |
| **207** | **4** | 715 |
| 307 | 4 | 715 |
| 347 | 4 | **577** |
| 397 | 4 | **10** |

Being `k`-chromatic is necessary for rigidity at `k` and nowhere near
sufficient — a factor of about **1.9** in size separates the two, and
everything in that gap is structure the chromatic number does not see.

Naming that structure is worth more than another search, and it is measurable
on the same peel. A Moser spindle is assembled from rhombi — two points at
squared distance 3 with two common unit neighbours — so rhombus memberships per
point count the raw material:

| points | edges/v | rhombi | memberships | `χ` | census |
|---:|---:|---:|---:|---:|---:|
| 157 | 1.83 | 4 | 0.10 | 3 | 715 |
| **207** | 2.59 | 29 | **0.56** | **4** | 715 |
| 307 | 3.73 | 136 | 1.77 | 4 | 715 |
| 347 | 4.16 | 214 | **2.47** | 4 | **577** |
| 397 | 4.97 | 444 | **4.47** | 4 | **10** |

Across the gap density rises by **×1.9** and saturation by **×8**. Saturation
is what moves, and it moves with a threshold: rigidity begins near 2.5
memberships per point and collapses near 4.5.

`G`'s five-colour saturation is at most **0.001** gadgets per point, measured
independently. Against a threshold of 2.5–4.5 that is a shortfall of **two to
four thousand** — not the factor of 1400 estimated earlier from gadget sizes
alone, and now measured on both sides rather than argued from one.

And it predicts, which is what makes it a criterion. At 207 points a *random*
subset never leaves the ceiling; the same number of points chosen greedily for
rhombus membership reads **341 of 715**:

| points | random: sat / census | chosen: sat / census |
|---:|---:|---:|
| 207 | 0.56 / **715** | 3.88 / **341** |
| 257 | 1.01 / 715 | 4.81 / 259 |
| 307 | 1.77 / **715** | 4.38 / **245** |
| 347 | 2.47 / 577 | 4.55 / 174 |
| 397 | 4.47 / 10 | 4.47 / 10 |

This is the first design variable in the whole of this work that moves the
quantity of interest at all — density does not, size does not, symmetry does
not, contact does not. It is also **not sufficient**, and the table says so: at
257 points the chosen subset has *higher* saturation than complete `Sa` and a
census of 259 rather than 10. Saturation carries the census most of the way
down; the collapse still needs the complete object.

> **Corrected, by the next measurement.** Read as a statement about the
> *chromatic number*, that is false, and two things say so. First, de Grey's
> own chain has **identical local statistics** throughout — `Sa`, `Y` and `G`
> all sit at 4.98 edges/point, 4.47 rhombus memberships and 0.57 spindles per
> point, to two decimals — while `χ` goes 4, 4, **5**. Y and G are unions of
> rotated copies joined by six edges and a spindle, so nothing local changes
> and the whole gain is global. Second, growing *for* saturation drives it far
> past the family's and buys nothing: starting from `Sa` and adding the points
> that complete the most rhombi reaches **1397 points at 13.65 memberships**,
> three times `G`'s, and it is 4-colourable at every one of twenty-five steps. So saturation
> tracks the census *within* a design at fixed `k`, and predicts which subset
> of a given graph is rigid. It does not predict `χ`.

> **Corrected: the spindle count was counting hinges.** `count_spindles` checks
> three distances — `|a−d|² = |a−g|² = 3`, `|d−g|² = 1` — which is the *hinge*,
> three of a Moser spindle's seven vertices, and never that either rhombus is
> present. Requiring both rhombi and seven distinct vertices:
>
> | | hinge triples | true spindles |
> |---|---|---|
> | `Sa` | 576 (1.45/pt) | **228 (0.57/pt)** |
> | `Y` | 1152 (1.46/pt) | **452 (0.57/pt)** |
> | `G` | 2304 (1.46/pt) | **904 (0.57/pt)** |
>
> Sixty per cent of what was counted is not a spindle. The qualitative account
> survives — the triangular lattice has none by either count — but every figure
> scaled from 1.45 per point is 2.5 times too generous. That is
what `G`'s readings mean: 5-chromatic puts it at the analogue of 207 points,
not 397, and the 13 356-point unions of translates are *wide* rather than
*tight*, which is the same place.


## What the object would have to be, and how big

With no gradient to climb, searching *towards* the object is not a strategy:
it has to be built complete and tested once. That is what de Grey did, and the
shape of the design is now clear enough to scale.

The pattern is one of ratios. His construction rests on a graph whose census
falls — `Sa`, 397 points — and what makes it fall is saturation with the gadget
one level down, the **Moser spindle**, which has seven points: 576 of them
inside `Sa`, 1.45 per point. Carrier over gadget is `397/7 ≈ 57`.

One level up the gadget is a 5-chromatic unit-distance graph. The smallest
known has **509** vertices; `G` has 1581. At the same ratio:

| gadget | carrier | SAT variables |
|---|---|---|
| 509 (smallest known) | ≈ 29 000 points | ≈ 145 000 |
| 1581 (`G`) | ≈ 90 000 points | ≈ 450 000 |

The lower figure sits at the edge of what a modern solver settles; the upper is
past it. But neither can be attempted, because a graph of that size chosen
*without* a design simply colours — and nothing here gives the seed.

So the bottleneck is not solver time and not the size of the search space. It
is that the five-colour gadget is seventy times larger than the three-colour
one while the construction's ratio stays fixed.


## The global quantity is redundancy

The chain's local statistics do not move, so the gain is global. Asking each
graph at its **own** number of colours says which global quantity it is —
criticality is only meaningful at `k = χ − 1`, the colour relation only at
`k = χ`, and asking either at the wrong `k` measures nothing. (The first
attempt did exactly that: it asked criticality at `k = 4` on `Sa`, which *is*
4-colourable, and reported 60 of 60 essential for the empty reason.)

| | vertices essential at `k = χ − 1` |
|---|---|
| `Sa` (χ=4) | **0** of 60 — and *provably* zero |
| `Y` (χ=4) | **0** of 60 — and provably zero |
| `G` (χ=5) | **42** of 60, measured exactly |

`Sa`'s zero needs no budget and no sampling. It carries 228 Moser spindles and
its busiest vertex lies in only **72** of them, so deleting any single vertex
leaves at least 156 intact — and a graph containing a spindle is 4-chromatic.
`Y` likewise: 452 spindles, busiest vertex in 144, at least 308 survive.

So the two ends of the chain are opposites. `Sa` and `Y` are `k`-chromatic with
**maximal redundancy** — no vertex matters, because the property is carried
hundreds of times over. `G` is `k`-chromatic with **none** — the property is
carried once, by a 5-critical subgraph spanning most of the graph.

That is where rigidity comes from, and it is neither density nor saturation:
many *independent* `k`-chromatic subgraphs each constrain a colouring and the
constraints accumulate into a census of ten. `G`'s five-chromaticity is used up
exactly once, in being five-chromatic at all, and nothing is left over.

> **Refined, again by the next thing checked.** Redundancy is not sufficient
> either. The packed union of ten overlapping copies of `G` — 13 356 points —
> is redundantly 5-chromatic by exactly this definition: removing any vertex
> leaves at least seven copies untouched, so **no vertex is essential**. Its
> census is **855 of 855**. Redundant and entirely loose.
>
> And there is a measure that makes the gap look small and is wrong:
>
> | | gadgets per point | incidences per point |
> |---|---|---|
> | `Sa` at four | **0.57** | 4.02 |
> | `G` at five | **0.00063** | 0.63 |
> | ten copies of `G` | **0.00075** | 0.75 |
>
> By incidences the shortfall is a factor of six; by gadget **count** it is a
> factor of nine hundred. The count is the one that matters — each gadget
> imposes one constraint whatever its size, so a single critical subgraph
> spanning every vertex scores high on incidences and constrains nothing beyond
> making the graph `k`-chromatic. Counting incidences would have turned 900 into
> 6 and made the problem look nearly solved.

> **Withdrawn, immediately.** Packing density looked like an independent
> confirmation of the size estimate: `Sa` holds 228 spindles of seven vertices
> in 397 points — four incidences per point — so the same packing with
> 509-vertex gadgets needs `228 × 509 / 4 ≈ 29 000`, agreeing with
> `397/7 × 509`. It is not independent. `228 × 509 / (228 × 7 / 397)`
> simplifies to `509 × 397 / 7`; the spindle count cancels. One identity
> written twice, and the agreement was guaranteed.


## The obstruction, stated exactly

Rigidity is not density, not saturation, not redundancy. What is left is the
**number** of independent critical subgraphs per point — each imposing one
constraint on a colouring whatever its size. Along the census curve:

| points | spindles | per point | census at four |
|---:|---:|---:|---:|
| 107 | 0 | 0.000 | 715 |
| 207 | 2 | 0.010 | 715 |
| 307 | 31 | 0.101 | 715 |
| 347 | 56 | **0.161** | **577** ← leaves the ceiling |
| 397 | 228 | **0.574** | **10** ← collapses |

Now the arithmetic that closes the account. A union of copies of one gadget has
gadget density **exactly `1 / (new points per copy)`** — `c` copies costing `p`
new points each occupy `c·p` points and contain `c` gadgets. That ratio is an
**invariant of the overlap** and does not depend on the number of copies at
all, which is why every "add more copies" experiment here moved nothing.

| | new points per copy | density |
|---|---:|---:|
| `Sa`'s spindles | **1.74** | 0.574 |
| `G` translated, measured | **1336** | 0.00075 |
| a 509-gadget at `Sa`'s overlap *fraction* | 126.6 | 0.0079 |

The third line is the one that hurts. Even reproducing `Sa`'s overlap
proportion exactly — three quarters of every gadget already present — a
509-vertex gadget costs 127 new points and lands **seventy times below the
threshold**, because density counts *gadgets* and a bigger gadget costs more
points for the same single constraint.

So the demand is this: to reach the threshold with a 509-vertex gadget, each
copy may contribute **at most six new points** — 98.8 % of every gadget must
already be in the graph. `G`'s translates contribute **1336**.

That demand rested on a choice no argument settles: does the threshold travel
in **gadgets** per point (0.161) or in vertex **incidences** per point (1.13)?
Carried to a 509-vertex gadget those differ by 143 — the first demanding 98.8 %
overlap, the second only `Sa`'s own 75 %. A third level decides it. At three
colours the critical gadget is the **triangle**, and a peeled lattice patch
gives:

| | gadgets/pt | incidences/pt | census |
|---|---:|---:|---|
| `k=3`, leaves the ceiling | **0.150** | 0.450 | 284 of 365 |
| `k=4`, leaves the ceiling | **0.161** | 1.13 | 577 of 715 |
| `k=3`, collapses | 1.092 | **3.277** | 1 of 365 |
| `k=4`, collapses | 0.574 | **4.02** | 10 of 715 |

The **onset** travels in gadget count — 0.150 against 0.161, seven per cent
apart, where the incidence reading differs by a factor of 2.5. The **collapse**
goes the other way. Both are reported; the onset is what decides whether
rigidity appears at all, and it is the harder verdict. So ≈ **0.155 critical
subgraphs per point**, which for a 509-vertex gadget means at most 6.5 new
points per copy and **98.7 % overlap** — now supported by two levels rather
than assumed from one.


## Rigidity is not local, and that retires the rest

The peel showed census and gadget density moving together, which invited the
reading that gadget density is what rigidity tracks — the last local quantity
standing after density, saturation and redundancy had each been ruled out. The
constructive test refutes it outright.

Growing `Sa` for rhombi raises the **spindle** count per point fourfold, from
0.574 to **2.309** at 997 points — four times the density at which `Sa`'s own
census collapsed. The census over the same growth:

| | points | spindles/pt | census |
|---|---:|---:|---:|
| `Sa` | 397 | 0.574 | **10 of 10** |
| step 4 | 557 | 1.413 | **10 of 10** |
| step 8 | 717 | 1.863 | **10 of 10** |
| step 11 | 837 | ≈2.0 | **10 of 10** |
| step 15 | 997 | **2.309** | — |

Not one pattern eliminated, at any step. Adding vertices can only *remove*
survivors, so ten of ten means the tightening is exactly zero.

The peel's correlation was not causal. Removing vertices destroys a graph and
everything falls together; it does not follow that raising one of the things
that fell will raise the others — and it does not.

What is left is the operation that *does* tighten: **the bite**, which takes
the census from ten to seven and, unpruned, to three — and which adds six edges
and one shared vertex, changing no local statistic at all. `Sa`, `Y` and `G`
agree to two decimals on every local quantity measured here while `χ` goes
4, 4, **5**.

So: **rigidity is a global property.** No local statistic predicts it, produces
it, or improves it, and the only levers that move it are the bite and the
spindle. Every constructive strategy attempted in this work optimised something
local — which is exactly why all of them read flat.


## The search over operations is closed, by a theorem

Rigidity is global, so the question becomes which global operations exist. That
turns out to be answerable, not merely searchable.

**The principle.** An operation can tighten the census of a set `S` only if it
maps `S` to itself — a matching that reaches outside `S` says nothing about
`S`'s own patterns. So the family to try is the **stabiliser** of whatever is
being constrained, and both stabilisers here are small enough to exhaust.

| constrained | stabiliser | tried | result |
|---|---|---|---|
| a ring about `c` | rotations about `c`, mirrors through `c` | 104 bites about the origin — **the complete family** | best census **3**, at `D = 4` alone |
| | (confirmation) bites about every other centre | 2440 | census 10, every one |
| | (confirmation) mirrors, up to 13 matched points vs the bite's 6 | 60 | census 10, every one |
| a pair `{u,v}` | half-turn about the midpoint, both axis mirrors | 81 pairs, whole stabiliser | **0 forced** |

**And the enumeration is complete by necessity, not by choice.** Beckman and
Quarles (1953): *every unit-distance preserving map of `ℝⁿ` into itself, `n ≥ 2`,
is an isometry.* So gluing a congruent copy of a unit-distance graph into the
plane has no choice about what it is. There is no exotic transformation waiting
to be found; "copy the graph and glue it on" is a family with exactly the
members above.

**Nor is the floor a property of the arithmetic.** The field decides which
rings are closable and hence which bites exist at all — and `Sa`'s points lie
in the base field, so the closure is the *same graph* in every field containing
it, with only the operations changing:

| field | closable rings (of 279) | bites | best census |
|---|---:|---:|---:|
| `ℚ(√3,√5,√7,√11)` | 52 | 104 | **3** |
| `ℚ(√2,√3,√5,√7,√11)` | 80 | 160 | **3** |
| `ℚ(√2,√3,√5,√7,√11,√13)` | 97 | 194 | **3** |

Ninety operations that did not exist in de Grey's own field, and not one beats
his `D = 4`.

**So the chain reaches the floor of its own family and the floor is three of
ten.** All 344 second bites of `Sa ∪ Sb` about the origin — where the knife
edge means the answer can only be 3 or 0 — come back 4-colourable. De Grey's
pruning to `Y` trades two of three forced pairs for two vertices; his spindle
converts the one that remains. Nothing goes further, and *nothing can*.

What this does **not** close is the other kind of construction: adding points
that are the image of nothing — the seed search. That is open, and its size is
measured: **one in 2³⁹** for a single level.

## Testing the universe instead of the seed

The seed search is the right target, but it was being asked the wrong way. The
cap is **monotone in both directions**: it survives adding vertices, since more
vertices means fewer colourings; and an *un*capped ring stays uncapped in every
subgraph that still contains it, since a colouring restricts. So for a ring of
at least `k` points, testing the **maximal** point set answers the question for
every subgraph containing that ring — and a universe that fails rules out its
entire subgraph lattice at once.

Built as richly as the arithmetic allows: every point a short unit walk reaches
using every unit direction the field offers, inside a disc.

| radius | directions | points | edges | rings | capped at 5 |
|---:|---:|---:|---:|---:|---:|
| 2.6 | 30 | 1 717 | 8 292 | 7 | **0** |
| 3.5 | 66 | 24 003 | 122 848 | 21 | **0** |
| 5.0 | 90 | **70 021** | **374 333** | 27 | **0** |

Seventy thousand points and three hundred and seventy thousand edges, swept in
**105 seconds**, because "not capped" is the satisfiable answer. That is by far
the largest object tested in this work, and its negative is not about it alone.

> **Corrected by its own calibration.** The *method* is sound; the *universe*
> is the wrong shape. Run at **four** colours — where a capped ring is known to
> exist, `Sa`'s `D = 4` being capped at two — the same 24 003-point universe
> returns **0 of 22**. An instrument that misses the known case says little
> about the unknown one.
>
> The reason is exactly what invalidated the first seed search, and should have
> been checked the same way: the universe is built by unit **walks**, reaching
> radii 0.027 and upward, while de Grey's points sit at 0.168, 0.264, 0.292,
> 0.333 — arising as **intersections of unit circles**, which no walk produces.
> Only **19 of `Sa`'s 397 points** lie in it.
>
> So the five-colour negative is true of every subgraph of a 70 021-point
> universe that keeps one of its rings — a wide claim about a large space — but
> it is **not** a claim about the space that matters, and calling it "a whole
> class of constructions" was wrong. The fix is to build the universe from
> unit-circle intersections instead.

**Rebuilt the right way, the instrument calibrates.** Points that carry caps
arise as *intersections* of unit circles about points already present — which
is exactly where de Grey's radii of 0.168 and 0.333 come from. Seeding with
`Sa` and closing once under that generator takes 397 points to 3 001, and at
**four** colours it finds **three** capped rings where the walk universe found
none:

| ring | points | |
|---|---:|---|
| `D = 7/3` | 19 | capped |
| `D = 4` | 10 | capped — de Grey's own |
| `D = 16/9` | 4 | capped |

It detects the known case and two besides. Pointed at **five** colours on the
same 3 001 points: **0 of 7**.

That negative means what the earlier one did not. The universe demonstrably
*contains* capped configurations at four, so its silence at five is a statement
about five rather than about the universe — and by monotonicity it covers every
subgraph that keeps one of those rings.

Four times larger, the same shape of answer. Closing `Sa` once at radius 4.5
gives **8 953 points and 47 724 edges**:

| | rings of ≥ k points | capped |
|---|---:|---|
| `k = 4` | 11 | **4** — `D = 5/9` (36 pts), `D = 7/3` (36), `D = 4` (18, de Grey's), `D = 16/9` (6) |
| `k = 5` | 11 | **0**, in 915 s |

The tell is in the timings: deciding five-colourability took **164 seconds**
against **one** at four colours, so the graph is far tighter at five — and caps
nothing.

> The `int64` guard earned its keep here too. Asking for *two* rounds of
> intersections reached 14 001 points with headroom **1.30**, and the run
> aborted rather than handing back a graph with no edges and calling it
> colourable — which is exactly the silent failure that cost three levels of
> meaningless measurements earlier in this work.


## Two mistakes in the instrument, and what they were hiding

Both were mine, both were structural, and both were found by making the
instrument reproduce de Grey's own construction instead of only searching
past it.

**The carrier was never measured.** Every cap test at five colours had to run
on some graph, and which graph was never examined. de Grey's `H` is
4-chromatic and his lemma caps a centre together with its ring at two colours
out of four: carrier and question carry the same number. That is the
condition, not a coincidence. A cap says the colouring has no room, and
having no room is what chromaticity means — with five colours on a graph that
needs four, one colour is free at every vertex and any ring can be shown all
five by permuting inside the slack.

Sa's universe, 8953 points and 47724 edges, is 4-colourable. So every
five-colour cap test run on an Sa- or Y-seeded universe was asking a
4-chromatic graph to have no room for a fifth colour, and would have returned
zero whatever the geometry inside it. G-seeded universes are the only ones in
this work that have ever been in the right regime; the 26002-point one is
5-chromatic exactly, which is de Grey's own situation one level up.

**The ring scans discarded de Grey's own bite.** Every ring scan here filtered
to distance squared at most 4, on the reflex that a wider ring cannot have two
adjacent points on it. For a palette bound that is harmless. For a bite it
throws away the mechanism, because a bite does not stitch ring points to each
other — it stitches each point to its own image under the turn, and
`2·rho·sin(theta/2) = 1` has a solution at every radius at least 1/2.

de Grey's `Sb` turns `Sa` by cosine 7/8, the bite on the ring of radius 2,
which the filter keeps. His `G` turns `Y` by a relative `2·arcsin(1/8)`,
cosine 31/32 — and that is the bite on the ring of radius **four**, distance
squared 16, thrown away by every scan in this work.

## de Grey's construction is two bites whose radius doubles

Written out, `S -> Sa -> Y -> G` is: close under the order-twelve group once,
then bite twice. The chromatic number does not move on the first bite — `Sa`,
`Y` are both 4-chromatic — it moves on the second.

The two bites are not independent. The first turns `Sa` about the origin on
the ring of radius 2. The second turns `Y` about `(-2, 0)`, which is not an
arbitrary new centre: it is a point **on that first ring**, since the origin's
radius-2 orbit in `Sa` contains it. And the second radius is twice the first.

So the step is: bite a ring of radius `rho` about `c`, move the centre to a
point of that ring, bite a ring of radius `2·rho` about it.

Whether the radius can double again is a question about the field. Since
`sin(theta) = sqrt(4·rho² − 1) / (2·rho²)`, the quantity `4·rho² − 1` must be
a square times one of K's radicands:

| radius | `4·rho² − 1` | in `Q(sqrt3, sqrt5, sqrt7, sqrt11)`? |
|---|---|---|
| 2 | 15 | yes — de Grey's `Sb` |
| 4 | 63 = 9·7 | yes — de Grey's `G` |
| 8 | 255 = 3·5·17 | **no** — needs `sqrt17` |
| 3 | 35 | yes, unused |
| 5 | 99 = 9·11 | yes, unused |

The obvious third level does not exist over his field. Two others do, and
neither has been tried.

## Symmetrising only works at the right centre

Closing `G` under the order-twelve group about the origin gives 18966 points
and 94548 edges. That is twelve times G's 7877 exactly, at 9.97 edges per
vertex, which is G's own density exactly: twelve disjoint copies wearing one
name, 5-colourable for free. `S` straddles the origin so closing it there
fuses; `G` is `Y` turned about `(-2, 0)` and lives in the upper half plane
around `(-2, 2)`, so closing it there does not. Y's centroid, by contrast, is
exactly the origin.

The centre `G` actually has is that pivot. It is a vertex of `G`, 129 of its
1581 points sit at rational distance from it across 31 rings, and the twelve
images about it are centred two apart on a circle of radius two while each has
radius about 2.8 — so they overlap. The closure there, `Gp`, has 13873 points
and 73782 edges, 10.64 per vertex, contains `G` so it needs five colours, and
colours with five so it needs no more.

`Gp` is the first object in this work that is tight and symmetric at the same
time. It has 30 rings of at least six points about its pivot, twelve of them
biteable in K — eleven being turns nobody has applied.

## The cap is abundant at four and absent at five

The contrast is not one of degree.

At four colours, on Sa's universe, at de Grey's exact strength — a centre
together with its ring held to two colours out of four — the scan returns one
every couple of minutes: twenty of them, across five distinct centres, on the
rings `D = 4/9` and `D = 16/9` and `D = 4`. His configuration is not a needle
in a haystack.

At five colours, in every carrier whose chromatic number matches the question:

| carrier | points | test | capped |
|---|---|---|---|
| G's universe | 35132 | 579 balls, 40 centres | 0 |
| `Ga` (origin closure) | 18966 | whole graph, every ring | 0 |
| `Gp` (pivot closure) | 13873 | whole graph, every biteable ring | 0 |

The one exception is free and carries no information: a ring at distance
exactly 1 from its centre is capped because the centre touches all of it, so
the ring loses the centre's colour. Centre-plus-ring is then back to `k`.

Biting `Gp` on each of its wide rings genuinely stitches — `D = 4` adds 30744
cross edges, `D = 16` adds 17712, `D = 9` adds 23496 — and every union is
5-colourable.

## Forcing is absent even at four, which leaves only the cap

A forced pair needs no bite: the spindle converts it straight into another
colour. It had never been censused broadly here because the test was run on
whole graphs, which is expensive. But forcing travels upward exactly as
capping does — "`u` and `v` agree in every `k`-colouring" is "the graph plus
the edge `uv` does not colour", and adding vertices only removes colourings —
so a ball suffices, and the census becomes cheap.

Calibrated on the smallest thing that forces: a rhombus of two unit triangles
forces its two tips to agree at three colours and not at four, and the Moser
spindle built from two of them is not 3-colourable. On Y's universe at four
colours, 740 ball tests around the twenty-five busiest vertices return **zero**
forced pairs.

That is the expected shape rather than a surprise. The local mechanism for
forcing needs the neighbourhood of one endpoint to contain a `K(k−1)`, so that
together with the other endpoint it closes a `K(k)` and leaves exactly one
colour. Unit-distance graphs in the plane have clique number 3, so the
mechanism is available at three colours and nowhere above it. de Grey did not
use a forced pair, and at four colours there is none to use.

## Honest odds

Polymath16 worked on this for years. The chance that this finds a 6-chromatic
unit-distance graph is small, and nothing here did.

What it provides instead is threefold. A correct, fast, fully certifying
search, against which any future claim can be checked in one command. A
**positive structural account**: de Grey's bite reduced to a computed census
and a finite elimination that fits on a page, with the consequence — all three
antipodal pairs forced, not one — that his own presentation understates; plus
the fact that his pruning of two vertices costs two thirds of that conclusion.
And a **precise statement of what is missing**, arrived at by eliminating every
alternative one at a time:

> Rigidity is **global**. Density, rhombus saturation, redundancy and gadget
> density were each measured, each correlates along a peel, and each was then
> tested constructively and moved nothing — spindle density four times `Sa`'s
> eliminates *not one* of its ten patterns. The only operations that have ever
> tightened a census here are the bite and the spindle, and the bite is a
> one-off.

The negatives are recorded to the same standard as anything else: what was
measured, at what cost, and what it does and does not license. Several of them
are corrections to claims made earlier in this same work — the spindle count
that was counting hinges, "saturation is the discriminator" retracted by the
next measurement, a packing estimate whose independent-looking agreement was
algebraically guaranteed, and a silent `int64` overflow that returned a graph
with no edges while the sweep called it colourable. Those are in the record
beside the originals, not in place of them.


### What the second pass added

The account above was written before the instrument was turned on de Grey's
own construction rather than only pointed past it. Doing that produced three
things.

**His mechanism, derived rather than described.** An exhaustive scan of `Sa`
— every vertex as a centre, every radius, palette bounds on the whole graph —
returns his lemma with his numbers: the origin together with the regular
hexagon of radius 2, `(±2, 0)` and `(±1, ±√3)`, held to **two colours out of
four**. The bite that ring admits is `cos 7/8, sin √15/8`, which is `Sb`. The
cap and the bite are the same ring seen twice. A decision procedure — sample
colourings, demand the chosen set miss enough colours in each, verify exactly,
feed the refuting colouring back — reconstructs that set from nothing but the
graph and the question, in 175 seconds. Read as a rule, the whole construction
is two bites whose radius doubles, `2` then `4`, with the second centre a
point of the first bitten ring.

**Two mistakes in the instrument, and one claim withdrawn.** Every five-colour
cap test had been run on carriers that only need four colours, where a spare
colour at every vertex makes any ring showable — those zeros were not evidence.
Every ring scan had filtered to `D ≤ 4`, which discards de Grey's own second
bite at `D = 16`. And a hill-climbing search over free-form sets was withdrawn
after it failed its own calibration: a tight set is an isolated minimum, so a
local search cannot reach one.

**The gap, narrowed to one sentence.** At four colours de Grey's exact
configuration is *abundant* — 37 found in one sweep of `Sa`'s universe, across
five centres. At five colours nothing is there, in any carrier whose chromatic
number matches the question, and not only the cap: the **weak disjunction**,
which is all the bite actually consumes — *some antipodal pair of the ring is
monochromatic in every colouring* — is absent too, across 150 biteable rings
in `G`, 195 in the best carrier built here, and 12 in the pivot closure. So
the shortfall at five is not that a cap is too much to ask for. The weakest
statement the machine can run on is not there either.

And the recursion stops for a reason that is not about colourings at all. The
third level wants radius 8, whose turn needs `√17` — adjoining it costs
nothing and makes the angle exact — but no point of `G` lies 8 from the
radius-4 ring, because `G` only reaches 4.8 from its pivot. Every operation in
the catalogue fills a carrier in; none reaches out. The missing operation is
not another turn.

The tests recompute the numbers rather than quoting them, so a reader who
doubts any figure above can run it.


### What the third pass added: a scale instead of a verdict

Everything above asks one question — is this graph 5-colourable — and gets
one answer, about a hundred times. Two failures were indistinguishable: no
number said which had come closer. The third pass replaced the question.

The **circular chromatic number** χ_c is the least p/q for which a graph maps
to the circular clique K(p/q), whose vertices are p points of a circle with
adjacency at circular distance ≥ q. It satisfies χ_f ≤ χ_c ≤ χ and
⌈χ_c⌉ = χ, so it is a real number living inside the integer, and for every
5-chromatic graph here it lies in (4, 5]. The target is unchanged in
substance and sharper in form: **χ_c > 5 is exactly χ ≥ 6**.

The encoder is checked against known values before any of its answers are
believed: C₅ = 5/2, C₇ = 7/3, C₉ = 9/4, K₃ = 3, K₄ = 4, K₅ = 5, all exact.

**What the scale reads on the graphs here.**

| graph | χ | χ_c |
|---|---|---|
| Moser spindle | 4 | 7/2 = 3.5 |
| `Sa` (397 pts) | 4 | **4 exactly** |
| `Y` (791 pts) | 4 | **4 exactly** |
| `G` (1581 pts) | 5 | ≤ 9/2 = 4.5 |

`Sa` and `Y` refuse all 28 ratios below 4, down to 3 — they sit at the
ceiling of their interval with no circular slack. The Moser spindle does
not, so tightness is a property of the particular graph, not of
unit-distance graphs in general. And `G`, which contains `Y` and is
5-chromatic, maps to K(9/2): it is **not** tight at five. Its 4-chromatic
components reach four exactly; the construction on top of them lands about
half a colour short.

Read as a statement about de Grey's step, it buys roughly half a colour of
χ_c. That is enough to move ⌈·⌉ from 4 to 5, which is why one application
settled χ ≥ 5 — and it is also why one more would not settle χ ≥ 6, since
4.5 plus a half-step is 5.0, whose ceiling is still 5. Two steps reach 5.5,
and ⌈5.5⌉ = 6. The reading is a reading, not a theorem; what is measured is
the table.

**Why every growth run in this project decayed.** They scored each candidate
point by how many sampled 5-colourings it kills, on graphs with more than
half a colour of slack at five. Almost no single point can kill a
5-colouring of a graph that loose, so the score was nearly always zero — a
blind walk in the costume of a hill climb. Scoring at the ratio where the
graph is actually tight, on the same graph and the same pool: the best
candidate kills **120 of 120** homomorphisms where the old score managed
**7 of 250**.

That is a better question, not a working search. Twelve rounds later the
sampler still returned a full 120 fresh homomorphisms each time, because the
number of homomorphisms of a 1600-point graph to K(9/2) is astronomical and
killing the sampled ones leaves the space untouched. **The kill rate is no
more evidence than the solve time was.** Only an UNSAT proves anything.

**Where the difficulty actually sits.** At p/q a neighbour forbids 2q−1 of
the p positions, so blocking a point takes ⌈p/(2q−1)⌉ neighbours. Across
every ratio in (4, 5] with q ≥ 2 — 14/3, 19/4, 24/5, 29/6, 34/7, 39/8, 44/9
— that number is **3**. Only at 5/1 itself, where 2q−1 = 1, does it jump to
**5**. The whole local difficulty is concentrated in the last rung, which is
a precise way of saying why approaching five from below is cheap and
crossing it is not.

**The cone, and exactly which ratios it refuses.** Put a graph H on the unit
circle about a point v. Then v forbids a window of 2q−1 positions and H is
confined to the complementary arc of p−2q+1. Inside an arc, circular
distance equals linear distance, and the graph on M positions with i ~ j iff
|i−j| ≥ q is the complement of a unit interval graph — hence perfect, so its
chromatic number equals its clique number ⌊p/q⌋−1, and it contains a clique
that size. So H maps into the arc exactly when χ(H) ≤ ⌊p/q⌋−1, and v ∪ H
refuses every ratio in (4, 5) precisely when χ(H) = 4.

What the plane forbids is supplying that H: two points of a unit circle at
distance 1 subtend 60°, so a unit circle always carries a subgraph of
disjoint hexagons — bipartite, χ ≤ 2. The criterion needs 4 and the circle
gives 2. The obstruction is old; naming the exact ratios turns it into a
design criterion, and it says where to look next: one centre can never
confine a graph to an arc too short for its chromatic number, but three can,
since three windows cover the circle at every ratio in (4, 5].

**The milestone this opens.** Since χ_c is monotone under subgraphs,
χ_c(ℝ²) ≥ χ_c(H) for any unit-distance H, so every ratio a concrete graph
refuses is a quantitative statement about the plane proved by one UNSAT.
χ ≥ 5 gives only χ_c(ℝ²) > 4, and the fractional number — the quantity
Polymath16 pushed, to 3.8992 — sits below χ_c, so the interval between them
is a gap the usual instruments do not read. **Is there a unit-distance graph
with χ_c > 9/2?** It is strictly stronger than χ ≥ 5, strictly weaker than
χ ≥ 6, and it is a finite question.

#### The instruments the rung is out of reach of

Four routes can be closed off with numbers rather than opinion.

**Counting.** A map to K(p/q) needs p independent sets covering every vertex
exactly q times, so it needs an independence ratio of at least q/p. Turn it
round: a graph refuses K(p/q) by counting alone when its ratio is below q/p.
Screening every graph in this project by a greedy lower bound takes 25
seconds and puts all of them between 0.25 and 0.30, against the 2/9 = 0.2222
that K(9/2) would need. More decisively, for ratio ρ counting proves a
refusal only below 1/ρ, and the best unit-distance graphs known reach
ρ ≈ 0.2565, capping counting at 3.898 — below 4, which χ ≥ 5 already gives.
**For every unit-distance graph anyone has built, counting is strictly
weaker than the chromatic number.** That is χ_f(ℝ²) ≥ 3.8992 sitting below 4,
seen from a new angle.

**Fractional and measurable arguments.** Same counting, so the same cap: a
measurable homomorphism needs m₁ ≥ q/p, and m₁ ≤ 0.2544 caps the route near
3.93.

**Rigidity.** A uniquely k-colourable graph has χ_c = k, which would have
explained `Sa`'s and `Y`'s tightness — and would have made the next rung
demand a uniquely 5-colourable unit-distance graph. Measured instead:
around one 4-colouring of `Sa`, **zero of 100** same-class pairs are forced
to agree and 4 of 100 cross-class pairs forced to differ. Rigidity ≈ 2%
(`Y`: 0.8%), and both are tight anyway. The classical route is not the
mechanism, which also means the next rung does not need that object.

**Greedy growth**, at the slack ratio and at the tight one alike, as above.

#### The Minty side, which is what is left

Goddyn, Tarsi and Zhang: **χ_c(G) = min over orientations of max over cycles
of |C| / min(|C⁺|, |C⁻|)**. Brute force over every orientation and every
cycle reproduces the solver's independent answers five for five — C₃ = 3,
C₅ = 5/2, C₇ = 7/3, K₄ = 4, and the Moser spindle at 7/2, the one value no
textbook supplies. Two computations sharing no code.

A (p,q)-colouring then hands over a good orientation free: order by position
and orient upward. Seeded that way, `Sa` scores 4.0001 and `G` scores
**4.5003** on 1581 points — exactly the promised bound. Two hundred local
moves, each flipping an edge of a cycle witnessing the current score,
**could not get below 4.5**. Evidence, not proof, and of the weak kind; but
it agrees with the UNSAT run, and together they say **`G` is tight at 4.5**.
So the growth aimed at crossing K(9/2) was aimed at the right rung all
along; what failed was the method of crossing.

The witnesses are small. In `G`'s best orientation the cycles realising 4.5
are **80 of length 9**, each with exactly two edges on the minority side,
plus two of length 18 that are doublings — diameters 1.79 to 2.72. Compact,
local, enumerable.

Working out which cycles bind: triangles and 4-cycles are free, lengths 5–9
need ≥ 2 on the minority side, 10–13 need ≥ 3. A cycle with **one** edge the
minority way is a directed path closed by an edge back, so

> **χ_c(G) ≥ 5 iff every acyclic orientation contains a directed path on
> five vertices whose two ends are adjacent.**

And Gallai–Roy says every orientation of a 5-chromatic graph already has a
directed path on five vertices. The entire gap is whether some such path has
its ends a unit apart.

#### Why the ends never are

Measured on `G`: 42 distinct end pairs of directed 4-edge paths, distances
0.29 to 2.39, and **none in [0.99, 1.01]**; seven within 10% of a unit, the
nearest 2.8% off. The hole is forced, by one line of arithmetic. The
orientation comes from the position, adjacency makes each step at least
q = 2, and with p = 9 the highest a 4-edge path reaches is 0 → 2 → 4 → 6 → 8
— four steps of ≥ 2 already total 8, so every directed 4-edge path has
positions **exactly** 0, 2, 4, 6, 8. Its ends sit at 0 and 8, circular
distance 1, below the threshold. They cannot be adjacent, ever. The same
arithmetic explains why 600 starting vertices in a graph with 7877 arcs give
42 end pairs and not thousands.

So a homomorphism's orientation is built, by the wrap-around, to avoid
precisely the configuration that would raise the ratio to 5. `G` comes
within 2.8% of a closure it cannot make while a map survives. Killing every
map is still the only way through — this says exactly what a map is doing to
stay alive.

#### What the map looks like inside

`G`'s nine position classes size 69 to 269; their consecutive unions are all
independent and total 2n exactly, densities averaging 2/9 to the last digit.
All nine share a centroid near (−2.0, 2.0) and a spread near 0.7: **not
stripes, not translates, not a lattice — interleaved throughout.** The map
carries no geometric shape, which is a plain reason why adding rotated
copies of things never aimed at anything.

#### Searching the space of lemmas

De Grey did not grow a graph from a seed. He found a seven-point gadget with
a provable cap — the centre and its radius-2 ring take at most **two of four**
colours — and the construction is that lemma plus a way to make copies
conflict. He found it by insight in 2018. Nobody has enumerated small
configurations and asked each one for its cap, which is a search over
**lemmas** rather than over graphs, and it turns out to cost one solver call
each: *does any proper k-colouring spread all k colours across S?* UNSAT is a
cap.

Two things had to be fixed before it meant anything. A set with a common
neighbour outside it is capped **for free** — if w is adjacent to all of S
then w needs a colour none of them has — so the first version reported every
neighbourhood in `G` as a discovery. And "capped below k" is nearly
worthless when de Grey's cap is 2 of 4; the palette is now measured by
bisection so each lemma's *strength* is reported.

**The control decides whether the silence elsewhere counts.** On `Sa` at four
colours, over 1200 of 39566 candidates:

```
CAP 2 of 4  |S|=6  ring r²=4.000 about vertex 0      ← de Grey's hexagon
CAP 2 of 4  |S|=7  that ring with its centre         ← his gadget, exactly
CAP 3 of 4  |S|=6  ring r²=1.708 about vertex 1      ← one he doesn't mention
```

Found unaided, with his configuration and his exact value.

**The question.** `G` at five colours: 35500 configurations, **0 caps**.
`wide0` at five colours: 40000 configurations, **0 caps**. The carriers were
chosen to make it easy to pass — the cap is monotone *downward* in the
carrier, so `wide0`'s 11047 points give every configuration a better chance
than `G`'s 1581 — and still nothing. 75500 configurations, two carriers,
zero caps at five.

#### A 24-point graph with χ_c = χ = 4

`Sa` refuses every circular clique below 4. Activation literals turn that
into a shrinking problem — one literal per vertex implying it takes a
position, solve under assumptions, and the UNSAT core is a smaller vertex set
— which converged in two rounds from 397 points to 91. Greedy deletion took
it to vertex-minimality:

> **24 points, 54 edges, degrees 3–8, 17 triangles, diameter 3.221.**
> Not 3-colourable, 4-colourable, and refusing **all 92** circular cliques
> below 4 with denominator ≤ 12. So χ_c = χ = 4 exactly.

Exact field coordinates in `data/tight_four.json`; `tests/test_tight_four.py`
recomputes the edges, the chromatic number and the refusals from that file
alone — 11 tests, 3.4 seconds. The Moser spindle is 4-chromatic on **seven**
points but has χ_c = 7/2, half a colour loose. **Seventeen extra points is
what tightness costs at four colours in this family.**

Its shape says where tightness lives: 17 triangles in 24 points, mean degree
4.5, diameter 3.2 — a compact heavily linked cluster, not a long chain of
gadgets. That is the opposite of how the chromatic number is built up.

#### Iterating the step on it moves nothing

| construction | result |
|---|---|
| all 11 rotations by 30° (exact in the field) | 157 pts, 456 edges — **χ_c = 4** |
| one spindle, angle set by cos θ = 1 − 1/(2d²) | 88 of 552 pairs admit one — **χ_c = 4** |
| accumulated spindles, 24→38→75→145→269→519→1037 | 3683 edges — **χ_c = 4 throughout** |

`G` reaches 4.5 with 1581 points; this reaches 1037 and stays at 4. At
comparable size de Grey's construction gets half a colour further, so **his
choice of gadget and angles does work that repetition does not reproduce.**
It is the old invariant with a finer reading: gadget density is 1/(new points
per copy) and does not depend on the count, so doubling the graph doubles
copies and points together — and on the circular scale that shows as χ_c
staying at *exactly* 4 with no drift across a fortyfold size increase.

#### A synthesis, and its falsification

Caps are abundant at four and absent at five; `Sa` and `Y` are tight at four
and `G` is loose at five. Read as one fact, that says **caps live where the
carrier is tight** — a cap asserts the colouring has no room, tightness
asserts the same globally.

The tight cores falsify it. Both the 91-point and the 24-point cores refuse
every ratio below 4, so they are tight in exactly that sense:

| | tight at 4 | caps at 4 |
|---|---|---|
| `Sa`, 397 pts | yes | **yes** (2 of 4) |
| core, 91 pts | yes | **no** (0 of 2214) |
| minimal core, 24 pts | yes | **no** (0 of 124) |

**The two properties live in different parts of the graph.** The 24 points
carry all of `Sa`'s tightness and produce no cap; `Sa`'s caps come from the
other 306. A small compact core carries the tightness, the bulk around it
carries the cap — which is consistent with the monotonicity already recorded,
but the synthesis read a correlation off two points instead of drawing the
consequence.

What survives is weaker and honest: `G` has no caps at five *and* is not
tight at five, so tightness may still be **necessary** — nothing here tests
that. What is measured is that it is nowhere near sufficient. So "find a
carrier tight at five" is not the programme either; it would need bulk as
well, and bulk is exactly what the density measurements say cannot be added
without diluting.

### What the fourth pass added: a grammar, and a graph that is not his

The third pass ended by saying that "try another seed" is not a five-minute
experiment but a project: reproducing a construction of de Grey's scale in a
different field. Starting that project meant reading his design instead of
guessing at it, and the design turned out to have a grammar nobody seems to
have written down.

#### Every construction is a word in spindle letters

A rotation by `2·arcsin(1/(2r))` about a centre `c` sends **every** point at
distance exactly `r` from `c` to a point at distance exactly `1` from itself —
the chord is `2r·sin(t/2) = 1` by construction. Its cosine is `1 − 1/(2r²)` and
its sine is `√(4r²−1)/(2r²)`, so the rotation is rational in `r²` together with
one radical, `√(4r²−1)`.

Call `(c, r²)` a **spindle letter**. A construction is a word in those letters,
and the field it needs is the compositum of their radicals. `G` is a three-letter
word:

| letter | where | `4r²−1` | radical |
|---|---|---|---|
| `r² = 3` | inside `S` itself | 11 | `√11` |
| `r² = 4` | about the origin | 15 | `√3·√5` |
| `r² = 16` | about `(−2,0)` | 63 = 9·7 | `√7` |

and `ℚ(√3,√5,√7,√11)` is exactly what those three radii force. It is not a
choice; it is a consequence. His remark that the construction could not have
been found inside `ℚ(√3,√11)` is the same statement read backwards: the letters
at `r² = 4` and `16` are the two that leave it.

The alphabet is much larger than the three letters he used. On `Sa` the circles
whose spindle costs **no new radical at all** include `r² = 5/9` (twelve points,
`√11`), `7/3` (nine, `√3`), `7` and `13/3` (four each, `√3`), beside `r² = 1`
and `1/3`, which are the 60° and 120° rotations and therefore already symmetries
— a letter that is already a symmetry of the carrier buys nothing. His own
`r² = 4` has the **smallest** hinge in the whole list, six points, and is one of
the two he paid a new radical for. So the hinge count is not what he was
maximising.

#### A pure hinge cannot raise the chromatic number

It cannot be, and the reason is three lines. Let `H'` be a disjoint copy of `H`
and let the only new edges be `v–v'` for `v` in some set `R`. Take any proper
`k`-colouring `c` of `H` with `k = χ(H) ≥ 2`, and let `σ` be a fixed-point-free
permutation of the `k` colours — a `k`-cycle will do. Then `c' = σ∘c` is proper
on `H'`, and `c'(v) = σ(c(v)) ≠ c(v)` for **every** `v`, not merely for `v ∈ R`.
So the glued graph is still `k`-colourable however large `R` is.

Measured on de Grey's own two glues, the interfaces are:

| | shared points | cross edges | of which hinge |
|---|---|---|---|
| `Sa ∪ ρ(Sa) = Y`, 397 + 397 | 1 | 6 | 6 |
| `Ya ∪ Yb = G`, 791 + 791 | 1 | **1** | 0 |

`G` is 1582 points joined by **a single edge and a single common vertex**. The
six-edge hinge below it cannot have raised anything — and indeed `Y` is still
4-colourable. The fifth colour comes from the one incidental edge.

#### So `G` is a Moser spindle with `Y` as the rhombus

`Ya` and `Yb` are two rotations of `Y` about `(−2,0)`, so they share that pivot;
the one cross edge joins the two images of a single point `q` of `Y` at distance
4 from it. For that to be a contradiction, `Y` has to satisfy

> in **every** 4-colouring of `Y`, the pivot and `q` take the same colour,

and the solver says it does, on both sides. Then the two copies agree at the
pivot, hence agree at the two images of `q`, and those are one apart. That is
the whole of it. The rhombus forces its tips equal at three colours; `Y` forces
`(pivot, q)` equal at four.

#### The glue manufactures the forcing, and `G` has none at five

`Sa` has **no** forced-equal pair at four colours. Not a sample — a proof: forty
explicit 4-colourings separate all 78 606 of its pairs, and a pair that differs
in some colouring is not forced. `Y` has six. The glue made them out of a
carrier that had none.

The same instrument, run on `G` at five colours, separates all **1 248 990**
pairs with forty explicit 5-colourings. So `G` has no forced-equal pair at five,
and **cannot be spindled to six by de Grey's own step**. That is a complete
negative result about the most natural next move, and it cost ten seconds
because the filter is `O(colourings × solve) + O(n)` rather than `O(n²)`.

#### The cap is not what does the work

de Grey glues along a *capped* circle: the six points at distance 2 from the
origin show at most 2 of 4 colours. Enumerating the colour patterns that circle
can carry gives exactly ten classes — one monochromatic, nine with two colours —
which is his lemma, recovered without being told it. And it is not a small
gadget: shrinking `Sa` while keeping the lemma still needs **358 of the 397**
vertices. It is a global property, and nothing to transplant.

But gluing along circles that are **not** capped manufactures forcing just as
well, and usually more of it. That is what reopens the route.

#### A 5-chromatic graph that does not live in his field

Glue `Sa` to its image under the 60° rotation about a **vertex** — the glue
circle is that vertex's own unit circle, twenty points, uncapped, and the
rotation costs nothing. The 570-point union has **eight** forced-equal pairs at
four colours, all at squared distance `64/9`. Spindle one of them:

```
Z  =  1139 vertices,  6475 edges,  χ(Z) = 5
```

The spindle at `64/9` has `cos = 119/128` and `sin = 384√247 / 16384`, so `Z`
lives in `ℚ(√3, √11, √247)` with `247 = 13·19`. And it is not a redrawing of
anything of his: squared distances are invariant under every isometry of the
plane, so the field they generate is an invariant of the graph rather than of
this picture of it — and `Z` has a squared distance with a `√741 = √3·√247`
component, which `ℚ(√3,√5,√7,√11)` cannot express. **No congruent copy of `Z`
lies in the field `G` needs.**

The 39-point seed is de Grey's. The assembly — glue at an uncapped unit circle
about a vertex, then spindle the pair that appears — is not, and the result is
28% smaller than `G`. It is recorded exactly in `data/five_247.json` and
re-derived from those coordinates alone, edges included, by
`tests/test_five_247.py`.

That does not move `χ(ℝ²) ≥ 6`. What it does is show the ladder is a mechanism
and not a coincidence: **cap or no cap, a glue manufactures forced pairs, and a
forced pair at a spindleable distance is one rotation away from another colour.**
The rung at six now has a single explicit requirement — a 5-chromatic carrier
with a forced-equal pair at five colours — and an instrument that can certify
its absence in seconds rather than assert it from a sample.

#### The four-colour side improves without limit; the five-colour side does not move

The glue is a densifier that keeps the carrier rigid, so it can be iterated.
Each level takes the highest-overlap glue available and re-measures:

| level | n | mean degree | free@4 | forced pairs at 4 | its spindle | free@5 | forced at 5 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 (`Sa`) | 397 | 9.94 | 0.00% | **0** | — | — | — |
| 1 | 570 | 11.36 | 0.00% | **8** | 1139 | 10.18% | **0** |
| 2 | 679 | 12.05 | 0.00% | **27** | 1357 | 10.10% | **0** |
| 3 | 782 | 12.54 | 0.00% | **52** | 1563 | 10.88% | **0** |
| 4 | 880 | 12.86 | 0.00% | … | | | |

The size growth decays — ×1.44, ×1.19, ×1.15, ×1.13 — because the available
overlap climbs with each level (224/397, 461/570, 576/679). Mean degree rises
monotonically, rigidity at four holds exactly, and the forcing grows roughly as
the square of the level. **None of it reaches five.** The spindled graphs hold
at 10–11 % slack with no trend and no forced pair, at every level.

Spending all the forcing at once does not help either. `H` has eight forced
pairs; spindling all eight simultaneously — the copies sharing the whole of
`H`, not translated copies of a finished graph — gives 5122 vertices at mean
degree 11.71, still 5-colourable, slack still 12–14 %. The gadget-density
arithmetic predicts exactly that: each extra spindle costs 569 new points per
5-chromatic copy where the measured threshold is 6.5.

#### Rigidity is not the gate, and there is no gradient after all

Two data points made rigidity look like the requirement: `Sa` is at 0.00 % slack
at four and its glue forces eight pairs; `G` is at 17 % at five and its glue
forces none. So thin `Sa` and watch. Deleting the same glue's worth of vertices
and re-measuring:

| deleted | n | free@4 | forced pairs |
|---:|---:|---:|---:|
| 0 | 397 | 0.00% | 8 |
| 4 | 393 | 0.00% | 8 |
| 5 | 392 | 0.00% | 7 |
| 7 | 390 | 0.26% | 6 |
| 10 | 387 | **0.00%** | **0** |

**The forcing dies while the slack is still exactly zero.** Ten deletions out of
397 — two and a half per cent — and the carrier is every bit as rigid by the
measure that was supposed to explain it. Deleting the ten *lowest-degree*
vertices instead leaves all eight intact and even raises the mean degree. So
the slack is not what the glue is consuming, the 10 % at five colours is not the
obstruction, and the gradient this pass was steering by does not exist. The
earlier section's verdict stands after all, for a different reason than it gave.

*(The first version of this experiment re-chose the best-overlap glue for each
thinned carrier and picked near-symmetries — overlap 393 of 395, union 397 — so
its zeros measured that choice and nothing else. The numbers above hold the
glue fixed.)*

#### Three copies, which need only a palette of two

Everyone glues two copies, because the spindle is a two-copy device. If
`|s − q|² = 1/3` then `q` and its images under the 120° and 240° rotations about
`s` form an equilateral triangle of side exactly 1, so

    U = H ∪ rot120_s(H) ∪ rot240_s(H)

is `k`-colourable only if three `k`-colourings of `H`, agreeing on every shared
point, give `q` three different colours. That demands a palette of **two** where
the spindle demands a palette of one — strictly weaker — and `rot120` costs `√3`,
which every carrier here already has. Four copies are not available: four points
pairwise one apart do not exist in the plane, so three is the whole of the extra
room.

Scanned by overlap on `Sa` at four and on `G` and `Z` at five, every union is
colourable, up to overlap 1031 of 3417 on `Z`. The route is open; the carriers
are not good enough for it.

#### What the wall actually is

Putting the measurements together, the obstruction is not slack, not density,
not saturation and not redundancy:

> `Sa` is rigid at four because it packs **Moser spindles at 1.74 new points per
> copy**. The ladder needs that kind of packing at the colour it is attacking.
> At five the critical gadget is a 5-chromatic unit-distance graph, and the
> smallest known has **509 vertices**, where the measured onset needs at most
> **6.5 new points per copy**.

So `χ(ℝ²) ≥ 6` runs through a *small* 5-chromatic unit-distance graph — which is
an open problem in its own right, and the one place where a new construction in
a new field is worth something. That is what `Z` is for.

#### Why four works and five does not, as a number

The glue forces because agreeing on the shared set pins the rest of the
colouring. So what decides it is neither slack nor density but **how many
distinct patterns a small set of vertices can carry** over all proper
`k`-colourings. A carrier whose interface admits a handful of patterns leaves
the second copy no freedom; one that admits hundreds leaves the σ-argument
intact.

Eight high-degree vertices, 1500 colourings each, counted up to permutation of
the colours:

| carrier | `k` | distinct patterns | saturated? |
|---|---:|---:|---|
| `Sa` | 4 | **72** | yes — 63, 70, 72, 72 |
| `H` = `Sa` glued | 4 | **38** | yes, flat from 200 on |
| `G` | 5 | **653** | no — 44 % of colourings still new |
| the 807-vertex graph | 5 | **921** | no — 61 % still new |

At four colours the interface space closes. At five it does not close at all,
and is an order of magnitude larger.

**A correction, made the same day.** Reading those two rows together — 72 for
`Sa`, 38 for `Sa` glued — as *"the glue halves the interface"* was wrong: each
probe was the eight highest-degree vertices **of its own graph**, so they are
different vertex sets. Pinning the same eight points by coordinate and gluing
repeatedly gives

| | depth 0 | 1 | 2 | 3 |
|---|---:|---:|---:|---:|
| `Sa`, `k=4` | 72 | **59** | 59 | 59 |
| the 807, `k=5` | 689 | **693** | | |

so at four the glue contracts the interface **once** and then plateaus, and at
five it does not contract at all. There is no rate to extrapolate and no depth
at which this would close. What survives — and it is the part that matters — is
the absolute gap: **59 against 689**, with no saturation at five however many
colourings are drawn. And since the forcing does keep growing with depth (8
pairs, then 27, then 52), a fixed eight-vertex window sees the contrast between
four and five without seeing the mechanism that produces it.

*(Two corrections this needed, both of which produced confident nonsense first.
The encoding deliberately omits at-most-one clauses — sound for deciding
colourability, since the edge clauses already make adjacent colour **sets**
disjoint — but wrong for **reading** a colour off a model, where a vertex with
three true colour variables reads as whichever is smallest. Without
at-most-one, an eight-vertex interface reported **one** pattern over four
hundred colourings. And blocking alone does not diversify: forbidding the last
solution on forty random vertices lets the solver change one of them and leave
everything else, so a fixed interface still reports one pattern — out of
laziness, not constraint. Randomised polarity moves the whole assignment;
blocking guarantees it moves at all.)*

#### Two solver facts, so they are not paid for twice

Randomised decision polarity is what makes the forced-pair filter work. On a
tight graph, random **assumptions** are refused almost always — forty tries on
`Y` at four colours produced nothing at all — while random **phases** steer
minisat into a genuinely different corner every time, and that is how 78 606
pairs of `Sa` and 1 248 990 of `G` were separated.

But **cadical ignores `set_phases`.** Twenty "different" colourings of a
6607-vertex orbit came back identical and left 4 366 463 surviving pairs, which
is no filter at all. And **full phase randomisation actively hurts minisat on a
large loose instance**: fifteen minutes on a 6607-vertex 5-colouring that
cadical settled in **four seconds**, because the random hints fight the solver's
own heuristic. Blocking a random sample of the last solution works whatever the
solver does with hints, and took the same orbit from 4 366 463 surviving pairs
to zero.

#### Symmetry buys overlap and spends slack, and the slack is worth more

Overlap is the resource, and the deficit is plain: `Sa`'s best glue reuses 56 %
of it because `Sa` is a `D6` orbit, while the spindled 5-chromatic graphs top
out at 39 % because the spindle rotation is about a vertex and breaks the
symmetry. Restoring it costs 6607 points rather than twelve times 951, since
`Sa` inside is already fixed by the group — and it **loosens**:

| | mean degree | free@5 | forced pairs at 5 |
|---|---:|---:|---:|
| the 951-vertex graph | 10.87 | 11.46 % | 0 |
| its `D6` orbit, 6607 points | 11.19 | **13.14 %** | **0**, certified |

Twenty-four explicit 5-colourings separate all 21.8 million pairs of the orbit.

### What the fifth pass added: a group of our own

The fourth pass ended with the interface measurement and a wall. The fifth
begins with a question that turned out to be exactly right: *if the `D6` orbit
does not help, is that because `D6` is not our group?*

It is not. Measured, by searching for every isometry that maps a point set onto
itself — any such map fixes the centroid, so each candidate is determined by
where it sends one vertex, and the whole group is computable in `O(n²)`:

| | isometry group |
|---|---:|
| `Sa` | **order 12** (`D6`) |
| every graph built from it here | **order 1** |

The glue rotation is about **one** vertex and the spindle about **one** pivot,
and each destroys every symmetry `Sa` had. Taking the `D6` orbit of the
finished graph afterwards imposes `Sa`'s group on an object not shaped for it,
which is why it loosened rather than tightened.

#### Never break it instead of repairing it

Both operations respect conjugation — `g·rot_w·g⁻¹ = rot_{g(w)}` for `g` in the
group, and `g(Sa) = Sa` — so gluing at **every vertex of an orbit at once**
leaves the union invariant, and the forcing, being equivariant, arrives in
whole orbits. The same argument lets the spindle be applied at every pivot of a
forced orbit simultaneously.

| carrier | n | mean degree | free@4 | forced pairs |
|---|---:|---:|---:|---:|
| `Sa` | 397 | 9.94 | 0.00 % | **0** |
| sequential chain, 3 glues | 782 | 12.54 | 0.00 % | 52 |
| **symmetric glue, 6 centres at once** | **1021** | **13.34** | 0.00 % | **153** |

Three times the forcing at a comparable size — and a forced distance the
sequential chain never produced: nine pairs at `d² = 64/3`, whose spindle needs
`√759 = √3·√11·√23`, a third field.

#### The first 5-chromatic graph here that carries a symmetry

Spindling that forced orbit over all six of its pivots keeps the invariance:

```
n = 7141   m = 47682   mean degree 13.35   χ = 5   C6-invariant
```

The group is `C6` and not `D6`: the reflection is not in it, because the glue
centres form a rotation orbit only. Spending **both** forced orbits at once —
twelve spindles, still invariant, in `ℚ(√3,√11,√23,√247)` — gives
`n = 13261`, `m = 88548`, and `χ = 5` again.

What does change is the cost. Cadical, with a triangle pinned throughout:

| | | |
|---|---:|---|
| 6607 points, symmetrised *after* the fact, no spindle | **4 s** | |
| 7141 points, one symmetric spindle | **322 s** | 80× |
| 13261 points, two symmetric spindles | **1718 s** | 430× |

Four hundred times the work across the series for the same verdict. It is the
only monotone signal this project has produced — and it is worth exactly what
the earlier measurement says it is worth: `uniSa` maps in 13 seconds and
`deepgrow`, half the size and denser, in 2853, for the same answer. **Cost is
not evidence.** It is recorded as a direction, not a result.

#### Which orbit of centres, and a correction

The first symmetric carrier used the orbit of `Sa[25]` because that vertex had
come up earlier. Ranking all 38 distinct orbits by the density they produce
shows it is near the bottom:

| orbit | centres | n | mean degree |
|---|---:|---:|---:|
| `Sa[265]` | 12 | 2689 | **15.92** |
| `Sa[253]` | 12 | 2605 | 15.75 |
| … | | | |
| `Sa[25]` | 6 | 1021 | 13.34 |
| the origin | 1 | — | nothing: `rot60` about it is a symmetry of `Sa` |

and stacking the best orbits keeps climbing — 15.92, 16.81, 17.48, 17.75,
18.14, **18.32** at six orbits and 6043 points — with the full group of order
12 intact at every step.

**But density is not the lever, and reading it as one was wrong.** Every one of
those carriers, up to mean degree 18.32 — nearly double `Sa`'s 9.94 and half as
dense again as anything the sequential chain reached — is still 4-chromatic,
still exactly 0.00 % free, and cadical 4-colours each in **under a second**.
Raising the degree did not make the colouring problem harder at all. What makes
it harder is the spindle: the 7141-vertex graph, at mean degree 13.35, costs
322 seconds where a 6607-vertex graph of the same degree costs 4. The
constraint that matters is whether a forced pair has been spent, not how many
edges the carrier carries.


## Forced pairs compose, and the composite distance is a free parameter

Every construction so far took the field it was handed. The radii pick the
radicals: `r² = 3` needs `√11`, `r² = 4` needs `√15`, `r² = 16` needs `√7`, and
`ℚ(√3,√5,√7,√11)` is not a choice but a consequence. Composition inverts that.

Let `H` force `c(v) = c(q)` in every `k`-colouring, with `|v − q| = D`, and let
`τ` be the isometry carrying `v` to `q` whose rotational part is `R_φ`:

    τ(p) = q + R_φ (p − v)

`τ(H)` is a congruent copy, so it forces `c(τv) = c(τq)`, that is
`c(q) = c(τq)`. Forcing is transitive, so `H ∪ τ(H)` forces `c(v) = c(τq)` — a
**new** forced pair — and with `u = v − q`,

    v − τ(q) = u + R_φ u,        |u + R_φ u|² = 2 D² (1 + cos φ)

The composite distance sweeps **all of `[0, 2D]`** as `φ` turns, and the angle
landing on a chosen `d` is

> `cos φ = d²/(2D²) − 1`, **rational whenever `d²` and `D²` are.**

Only `sin φ` carries a radical. With `D² = 64/9` and `d² = p/q`:

    cos φ = (9p − 128q)/(128q)          sin φ = 3√(p(256q − 9p))/(128q)

so `√(p(256q − 9p))` is the entire arithmetic cost, and **choosing `p` and `q`
chooses the field**. Four targets, each built from the 1021-point symmetric
carrier, with the composite distance checked as an exact field element, the
composite forcing confirmed by the solver rather than assumed, and the result
tested for four colours:

| `d²` | radical | field | n | spindle |
|---|---|---|---:|---|
| 1 | `√247` | `(3,5,7,11,247)` | **2041** | **not needed** |
| 1/3 | `√759 = 3·11·23` | `(3,5,7,11,23)` | 4061 | yes |
| 2 | `√476 = 4·7·17` | `(3,5,7,11,17)` | 4081 | yes |
| 4 | `√880 = 16·5·11` | de Grey's own, by another route | 4081 | yes |

`d² = 1` is the case worth keeping: the composite pair lands at distance
exactly 1, so it is an **edge**, and the chained union refuses four on its own —
half the vertices and no second rotation. Its radical is the one the ordinary
spindle at `64/9` already needs, so that case is a sanity check rather than a
new field. The others are new.

### And it carries density up to five

This was the standing bottleneck. Stacking orbits of glue centres takes a
4-chromatic carrier to mean degree 18.32, and every one of those carriers is
4-coloured in under a second — density alone never bought a chromatic number.
Meanwhile every 5-chromatic graph here topped out at 13.82, because a spindle
adds a sparse second copy.

The chain tuned to 1 is *two copies of the carrier*, so it keeps the carrier's
degree:

| carrier | n | deg | chained n | chained deg |
|---|---:|---:|---:|---:|
| symmetric glue | 1021 | 13.34 | 2041 | 13.35 |
| two stacked orbits | 3463 | 16.81 | 6925 | 16.81 |
| **ten stacked orbits** | 6235 | 18.47 | **12469** | **18.48** |

A 5-chromatic unit-distance graph at mean degree 18.48, built in 27 seconds,
and with no new radical — `cos φ = −119/128`, `sin φ = 3√247/128` is already
the carrier's own field.

## The null model was wrong twice, and the raw numbers say it better

`free@k` counts the vertices with a spare colour, and a vertex of degree `d` is
free with probability about `(k−1)(1−1/(k−1))^d` under a random neighbourhood.
The honest null is the **average of that over the degree sequence**, not its
value at the mean degree. The function is convex, so by Jensen the mean-degree
version is strictly smaller whenever the degrees are spread — and every graph
here is a union of rotated or translated copies, dense in the middle and thin
at the rim, spread by about 4.6.

Two claims made against the wrong null are withdrawn: *"the Minkowski sum with
the hexagon is 0.8× the law, worse than random"* and *"the tuned chain is
0.62×"*. Against the per-vertex null nothing is worse than random:

| graph | n | deg | null | measured (12 colourings) | ratio |
|---|---:|---:|---:|---:|---:|
| `five_247_c` | 803 | 10.12 | 30.92 % | 15.73 ± 1.00 % | 0.51 |
| `five_247_b` | 951 | 10.87 | 27.18 % | 13.62 ± 1.33 % | 0.50 |
| `five_247` | 1139 | 11.37 | 24.38 % | 12.50 ± 0.86 % | 0.51 |
| `five_tuned_1_1` | 2041 | 13.35 | 16.00 % | 13.86 ± 1.75 % | 0.87 |

The 5-chromatic graphs beat the null by a steady factor of two; the tuned
chain — two copies joined at a point rather than interlocked — by only 1.15.

The per-vertex null also relocates the target: at mean degree 18.48 the
mean-degree law predicts `free@5 = 0.49 %` and the per-vertex null **7.09 %**,
and the factor of fourteen is the **boundary**. `free@5` is dominated by the
low-degree tail — a vertex of degree 6 is free with probability 0.71 — and a
finite unit-distance graph is a patch of the plane, so it always carries a rim.
Stacking orbits adds rim as fast as it adds middle.

### And then the per-vertex null fails too

The densest object measured, 6925 points at mean degree 16.81, comes back at
`free@5 = 16.40 %` against a per-vertex null of 9.67 % — **1.7 times more free
than random**, which no structural mechanism can produce and which therefore
condemns the null rather than the graph.

The reason is the bipartite-neighbourhood theorem, showing up quantitatively.
A neighbourhood is a union of paths and hexagons, and a proper colouring of a
path uses **two** colours. So the colours on `N(v)` are not a random draw from
`k−1` values at all: they are drawn from a structured distribution that
strongly favours few colours, and a missing colour is far more likely than any
independence assumption allows. Against a baseline the geometry forbids, a
ratio means nothing.

What the measurements say without a baseline is simpler and stands on its own:

| graph | n | mean degree | free@5 |
|---|---:|---:|---:|
| `five_247_c` | 803 | 10.12 | 15.73 % |
| `five_247_b` | 951 | 10.87 | 13.62 % |
| `five_247` | 1139 | 11.37 | 12.50 % |
| `five_tuned_1_1` | 2041 | 13.35 | 13.86 % |
| `five_dense_2` | 6925 | 16.81 | 16.40 % |

> **`free@5` is flat at 12–16 % across a 66 % rise in mean degree.** It does not
> trend down; density does nothing to it. Meanwhile `Sa` reaches exactly
> **0.00 %** at 397 points with **minimum degree 4**. Rigidity is not a
> quantity these constructions move at all — which is what the bipartite
> neighbourhood said in advance, since a vertex cannot be pinned locally at
> four colours or more.

The cost is worth recording too: finding a single 5-colouring of that
6925-point graph took cadical **8921 seconds**, against 50 s for the 2041-point
one. The instances are hard. They are also satisfiable.

## Both relations are empty at five

Every filter here had hunted pairs forced to **agree**, since that is what the
spindle consumes. The bipartite neighbourhood points at the other one: a
neighbourhood has maximum degree 2 *and no 4-cycle* — four steps of ±60° cannot
sum to 360° — so a degree-4 vertex has at most 3 edges among the 6 pairs of its
neighbourhood, and pinning it at five colours needs at least **3 of those pairs
forced apart while non-adjacent**.

Forced-apart also closes the problem by itself, with no rotation at all: five
points pairwise forced apart use all five colours in every 5-colouring, and a
sixth forced apart from all five has nothing left. Six pairwise forced apart
cannot sit inside a 5-colourable graph, so five is the ceiling and the sixth
point is what a construction step would have to supply.

The filter is the mirror image and certifies the same way — a pair that agrees
in any exhibited proper colouring is **proved** not forced apart:

| graph | n | non-adjacent pairs | rounds to empty | survivors |
|---|---:|---:|---:|---:|
| `five_247_c` | 803 | 253 642 | 140 | **0** |
| `five_247` | 1139 | 512 342 | 85 | **0** |

Two seconds and seven seconds. At five colours **both relations are empty**:
nothing is forced to agree and nothing is forced to differ.

## All 25,493,370 pairs of the symmetric graph, settled

The 7141-point `C6`-invariant graph has no pair forced to share a colour in
every 5-colouring, and the statement is now complete rather than partial:

| | pairs |
|---|---:|
| eliminated by exhibited 5-colourings, through the group-expanded filter | 25 491 171 |
| candidates surviving | 2 199 |
| orbit representatives, since forcing is equivariant | **376** |
| forced | **0** |

Three things had to be right. A pinned triangle is free for *deciding*
colourability but destroys colour symmetry, so the usual shortcut — "`x = 0` and
`y = 1` is impossible, hence they agree" — is unsound there; the test forbids
the pair from sharing **any** colour instead. The per-pair selectors gating
those clauses must be **assumed false** for every pair not under test, since
left free they are unconstrained variables the solver may set true, switching
on all 376 constraints at once — and cadical ignores `set_phases`, so
assumptions are the only switch. That took the base solve from over twenty
minutes to 428 s. And the 71 410 at-most-one clauses, which matter only when a
colour is read off a model, were dropped, because nothing is read here.

## A floor that cannot be lowered by deletion

`data/five_247_c.json` — 803 points, `ℚ(√3,√11,√247)` — is **vertex-critical**.
All 803 vertices checked one at a time, none removable.

The check is cheap precisely because the answer is yes. Asking *"does `H − v`
still refuse four?"* asks the solver to **refute** `H − v`, which is the slow
direction; if `v` is essential then `H − v` is 4-colourable and a colouring
comes back fast. Only a removable vertex costs a slow proof, and there are
none: 803 calls in 1389 s, against 47 s for a single refutation of the whole
graph — itself affordable only because a pinned triangle deletes the `4!`
symmetric copies of every refutation, taking it from over five minutes to 47 s.

## Two results re-derived that were already here

Twice this pass a theorem was worked out from scratch that the earlier passes
had already recorded — the bipartite neighbourhood (and with it *rigidity is
not local*), and consuming a disjunction with a stack of rotated copies. The
credit belongs to those sections; what the new pass adds is narrower and is
stated as such: the consequence that local pinning needs `χ(N(v)) = k−1` and so
exists only at `k = 3`, the absence of 4-cycles in a neighbourhood, and the
**two-radius** form of the stack, where the options sit at different radii, the
pigeonhole becomes a statement about two circulants at once, and one of them
must contain an odd cycle — which is what pins one radius to `1/√3` in any
multiquadratic field and frees the other entirely.

## What a forced pair costs, and the 60° grading that explains the wall

Two measurements this pass changed what the search should be doing, and one
theorem explains both.

**A forced pair costs the whole near-critical graph.** `five_247_c` is
5-chromatic and *vertex-critical*, so `H = G − p` is 4-colourable and
`μ₄(H,p) = 4` for **every** one of its 803 vertices — 803 ceiling points, free
of charge. `μ₄ = 4` alone does not make `N(p)` a rainbow: five points carrying
four colours leave one pair sharing. Pairwise forcing needs the cap attained,
`|N(p)| = k`, and exactly two vertices have degree 4. At `p = 315`,
`N = {130, 461, 561, 757}` holds one 60° edge, so five pairs there are
genuinely forced apart. Forcing is monotone in the vertex set, so binary
searching the smallest distance-ordered **prefix** that still forces each pair
prices it in about ten solves:

| pair | `d²` | smallest prefix | radius |
|---|---:|---:|---:|
| (130,461) | 2.124094 | **802 of 802** | 6 |
| (130,561) | 1 | 2 of 802 | 0 |
| (130,757) | 0.028382 | **802 of 802** | 7 |
| (461,561) | 1/3 | **802 of 802** | 6 |
| (461,757) | 2.457427 | **802 of 802** | 6 |
| (561,757) | 1.304951 | **802 of 802** | 7 |

Unanimous. Only the 60° edge is local; every pair forced *for a reason* needs
all of it, and dropping the single farthest vertex already breaks the forcing.
So a forced pair exists only one vertex short of `(k+1)`-chromatic, and hunting
one at five as a stepping stone to a 6-chromatic graph is circular. De Grey did
not find a 5-chromatic graph by finding forced pairs at four; he built one, and
the forced pairs came with it. **Construction is the route; forced pairs are the
receipt.**

**The 60° grading.** Inside `N(p)` every edge is a 60° step round the unit
circle, so along any path the angle moves by ±60° per step and the parity of the
path length is the parity of the angle index. A component is therefore a
self-avoiding walk on `ℤ/6` with ±1 steps: at most **six** points, all in one
coset of 60°, and the only possible cycle is the hexagon itself. And a
2-colouring of `N(p)` *is* the parity:

| separation | `d` | under `\|c(N(p))\| ≤ 2` |
|---|---|---|
| 60° (index 1) | 1 | apart — they are adjacent anyway |
| 120° (index 2) | √3 | **monochromatic, automatically** |
| 180° (index 3) | 2 | **apart, automatically** |

Checked rather than assumed: over the candidate neighbourhoods of the 803-,
2041- and 6925-point graphs, index parity and bipartition class agree on every
pair in every component. **Zero violations.**

This kills the one construction whose hypothesis survives its own rotation. If
a 5-colouring of `G' = ⋃ρⁱG` has `|c(N_{G'}(p))| ≤ 2` then — because `ρ` fixes
`p` and carries `N_G(p)` into `N_{G'}(p)` — each copy's restriction is an escape
colouring *of that copy*, so every copy contributes its own escape conclusion
and the conclusions chain. A conclusion at separation `θ` links a point to the
one `θ` further round; the colour classes are cosets of `⟨θ₁,θ₂,…⟩`; and the
contradiction wanted is two monochromatic points 60° apart. But
`⟨120°⟩ = {0,120,240}` misses 60° — the classes are the two inscribed triangles
— and by the grading, 120° is the *only* separation the hypothesis can ever
hand over. `⟨120°,180°⟩ = ⟨60°⟩` would close, and so would `⟨90°,120°⟩ = ⟨30°⟩`;
180° comes out forced *apart*, and 90° means `d² = 2`, which no triangular
lattice ever realises — **2 is not a Loeschian number**, being inert in the
Eisenstein integers with norm 4. Measured: zero 90° pairs in all four graphs
scanned. The stack reproduces the escape hypothesis instead of contradicting it,
for every graph, not just ours.

**The burden, and a heuristic that was backwards.** Blocking `p` means refuting
every proper colouring of `N(p)` that uses at most four colours. That count is
exactly computable, and on the 803-graph it is 2 380 for `|N| = 5`, 10 820 for
6, 28 240 for 7 — about 4.5× per extra neighbour, so a 30-point neighbourhood
(five hexagons, by the grading) carries a burden of order `10¹⁴`. The floor is
`|N| ≥ 5`. **The cheapest blocked point has exactly five neighbours**, and every
`μ` scan here had filtered to `|N| ≥ 8` or `≥ 10`. The ceiling at four says the
same from the other side: it is attained at `|N| = 4`, the minimum.

So the small neighbourhoods were scanned for the first time — all 802 candidate
points of the 803-graph with `5 ≤ |N| ≤ 7`, exactly, no budget. `μ₅ = 2`
throughout (`μ₅ = 1` at 11 of them). Still two; but the target is now a number
instead of a hope.

## Niven picks the radius, and it is 1/√3 — where μ finally reads three

Everything that made the single-point attack hopeless was about the **unit
circle**, not about neighbourhoods. On a circle of radius `r`, two points are a
unit apart at the angle `θ` with `cos θ = 1 − 1/(2r²)`. A cycle inside that ring
needs `θ` to be a rational multiple of `π`, and `cos θ` is rational there, so by
**Niven's theorem** `cos θ ∈ {0, ±½, ±1}`:

| `cos θ` | `θ` | `r` | cycles |
|---|---|---|---|
| ½ | 60° | 1 | even only — **the unit circle** |
| 0 | 90° | 1/√2 | even only |
| −½ | **120°** | **1/√3** | **triangles** |
| −1 | 180° | ½ | degenerate |

Exactly one radius in the plane puts an **odd cycle** on a ring, and it is
`1/√3`, the circumradius of a unit triangle. That is why the number turns up in
every disjunction this project has built, and it is why the unit circle is a
dead end: `χ(N(p)) = 2` there always, in every graph, so a 2-coloured
neighbourhood exists locally whatever the rest of the graph does.

On the `1/√3` ring the local floor is **three**. Measured for the first time:

| graph | hubs tested | `μ₅` on the `1/√3` ring |
|---|---:|---|
| `five_247_c` (803) | 120 richest | **3** at every one |
| `five_247` (1139) | 120 richest | **3** at every one |

Every `μ` this project has measured before read 2. This is the first number
above the floor — and it is the *floor* that moved: from a blocked hub the
distance is two rungs, not three.

A blocked hub is also the **consumable** shape. If the ring about `h` shows all
five colours in every 5-colouring then in particular `h`'s own colour is on it,
so every 5-colouring gives `c(h) = c(v)` for some `v` at distance `1/√3` — and a
rotation by 120° about `h` carries that ring *to itself*, so every rotated copy
names a partner on the same ring, and two partners 120° apart are a unit apart.
That is the disjunctive spindle with the odd cycle the odd-cycle theorem
demands. Not there yet: scanning every vertex of the 803- and 1139-point graphs,
none is a ring hub — each can avoid its own `1/√3` ring, rings of up to 18
points and 4 072 internal unit edges included.

## Forcing a distance, and the pigeonhole that nearly passed for one

Every forcing test here had asked the strongest question — *is this pair
monochromatic in every 5-colouring?* — and the certificate pricing says why the
answer is always no. The weak question is a disjunction and costs one solve,
since forbidding a distance is just adding its pairs as edges:

```
H forces d at five colours   ⟺   χ(H; {1,d}) > 5
```

No single distance is forced on the 803-graph: all sixty of its richest classes
leave it 5-colourable, one with 3 500 extra edges. **Two are:**

```
χ( five_247_c ; { 1,  3/2 − √33/6,  7/6 − √33/6 } )  >  5
```

so every 5-colouring has a monochromatic pair at one of those two squared
distances. The second is exactly a distance the escape analysis had named, and
the two differ by `1/3`.

A first pass asked instead for the fewest individual **pairs**, called it `δ`
and got ten. That is pigeonhole, not strain: the ten pairs complete `K₆` on six
vertices that already carry five unit edges, and some six vertices of this graph
carry nine, so six pairs suffice in *any* unit-distance graph with a rich enough
vertex. Withdrawn. The distance-class version survives because **a two-distance
set in the plane has at most five points**, so a `{1,d}`-graph can never contain
`K₆`; three distances can (the regular hexagon realises `1, √3, 2`), so the
two-class result was checked — 8 793 edges, no 6-clique.

Restricting every candidate pair to those two classes makes the trap impossible
by inheritance, and the UNSAT core plus a greedy pass gives the narrowest
disjunction this project has: **133 specific pairs** (76 + 57) over 65 vertices.
Genuine, and too wide to consume — no vertex meets more than 11 of them, so no
rotation about a point can chain the conclusions.

## Closing a disjunction: the complete list, and the unit square

Rotating `m` copies of a graph about a hub `h` consumes a statement
`c(h) = c(v)` for `v ∈ S` of width `w` when every choice function forces two
named points a unit apart. There are exactly two ways to arrange that.

**By pigeonhole**, with `m > w` copies. The adversary picks which pair repeats,
so `2r·sin(kψ/2) = 1` must hold for every `k = 1 … m−1` at once, which needs
`sin(kψ/2)` constant. `m = 2` solves it at any radius — the ordinary spindle,
width 1, a forced pair. `m = 3` needs `sin(ψ/2) = sin ψ`, so `ψ = 120°` and
`r = 1/√3` exactly. `m ≥ 4` has no solution at all.

**By all-pairs clash**, with `m = 2` and `w = 2`. Demand instead that all four
cross pairs clash: `|v_s − ρv_s| = 1` forces `r₁ = r₂ = r`, the two cross
conditions force `cos(a₂−a₁+φ) = cos(φ−a₂+a₁)` so the partners are **antipodal**
about `h`, and substituting gives `1 + cos φ = 1 − cos φ`, hence `φ = 90°` and
`r = 1/√2`. **That is the unit square** — `h` its centre, the partners one
diagonal, their images the other, and the four clashing pairs are the square's
four sides. Checked exactly on `(0,0)`, `(1,1)`, `h = (½,½)`.

> If a vertex `h` of a unit-distance graph `G` centres a unit square whose
> diagonal `{v₁,v₂}` lies in `G`, and every 5-colouring gives `c(h) = c(v₁)` or
> `c(h) = c(v₂)`, then `G ∪ ρ₉₀(G)` has no 5-colouring.

Two copies where `1/√3` needs three, and the 90° rotation is *rational*, so the
field never grows. What it needs instead is a pair at distance `√2` — and **2 is
not a Loeschian number**, so no triangular-lattice carrier has one. Every
carrier here starts in the Eisenstein lattice, the spindle rotations preserve
norms, and `σ` divides them by three (`2 = 6/3` would need 6 Loeschian, which it
is not). Measured after one `σ` round on the 803-graph, 8 323 points: pairs at
`d² = 2`, **zero**; at `d² = 1/2`, **zero**; square centres, **zero**. Three
translates by a vector of length `√2` do supply them — every point of `G` then
centres a square — but the copies barely touch: 2 409 vertices with **169** edges
between them, for every in-field half-offset tried, and one exhibited colouring
frees 206 of the 803 triples outright.

So of the two closures only the `1/√3` one is native here. And both come down to
the same rung: two partners cover at most two colours, so `c(h)` is forced onto
them only if `c(N(h))` covers the other three — **`μ₅(N(h)) ≥ 3`**, rung one of
the original ladder, now known to be the bottleneck of *every* route rather than
of one of them.

## σ, and what the centre-closure could not build

`σ(x,y) = ((x − √11·y)/6, (√11·x + y)/6)` is exact over `ℚ(√3,√11)`, since
`cos θ₀ = √3/6` and `sin θ₀ = √33/6`. Two lines of algebra give
`|σu|² = (x²+y²)/3` and `|u − σu|² = x²+y²`: it lands every point on the circle
of `1/√3` times its radius while moving it by *exactly its radius*. Centred
anywhere, `σ_c(u) = c + σ(u−c)` therefore sits at `1/√3` from `c` and a unit
from `u` — a cross-ring edge, built rather than hoped for.

That is what the **centre-closure** could not do. Adjoining the centre of every
unit triangle saturates the 803-graph at 1 851 vertices (803 → 1363 → 1781 →
1839 → 1851, then nothing), with 10 312 pairs at `1/√3` and rings of up to 30
points — and it never helped, because a triangular lattice has only angles that
are multiples of 30° while `θ₀ ≈ 73.22°` is not one, so those rings had no edges
to the neighbourhood they were meant to squeeze. All 1 851 vertices escape their
own ring.

One `σ` round on the 803-graph gives 8 323 points, 49 987 edges, rings of up to
50 and **86 896** cross-ring triangles — thirteen times what it started with.
Still 5-colourable, still no ring hub, and `χ(N(h) ∪ ring) = 3` — as the
two-ring theorem requires, since every such configuration is a disjoint union of
copies of one graph. A single ring triangle is capped at 3 for a sharper reason:
a point of `N(h)` reaches the ring only at `φ ± θ₀`, and those two are
`2θ₀ = 146.44°` apart, never a multiple of 120°, so it has **at most one**
neighbour in any given triangle, keeps a list of two of its three colours, and
`N(h)` is 2-choosable.

## Exoo–Ismailescu rebuilt, and the denominator that the project filtered out

**A reduction.** Exoo and Ismailescu (arXiv 1909.13177, *Geombinatorics* 2020)
proved that every 5-colouring of the plane has two points of one colour at
distance 1 **or 2**. Their 426-vertex `{1,2}`-graph lives in this package's
field — vertices `[a,b,c,d] = (a√3/12 + b√11/12, c/12 + d√33/12)` — and
rebuilt from the paper's 23 points every count matches: `G` 205/966/423, `H`
214/1004/446 with its pair `A, B` (distance 5) monochromatic in every
5-colouring (checked: `H` plus `c(A) ≠ c(B)` is UNSAT), `K = H ∪ λ_A(H)`
426/2009/892 with `|BB'| = 1`, where `λ = (49 + 3√−11)/50`. So:

> one unit-distance graph with two points at distance 2 that no 5-colouring
> colours alike proves `χ(ℝ²) ≥ 6`

— every congruent copy would split every 2-pair, and a 5-colouring of the
plane would become a proper colouring of the `{1,2}`-graph. In the 803-graph,
though, distance 2 is *attractive* (`P(same) = 0.29` against 0.2 at random:
the tangency point is the pair's only common neighbour), no pair closer than 3
is forced apart, and the most repulsive distance is `2/√3` (0.09).

**A correction.** The residue-degree theorem assumes edge vectors **integral
at 5**, and the summary above ("a multiquadratic unit-distance graph has a
coset 5-colouring") dropped that hypothesis; an earlier note even filtered
generators like `(1 + 3√−11)/10` out as bookkeeping. On the module the edge
vectors actually generate, a denominator of 5 genuinely blocks: in `ℚ(√−11)`
the prime 5 splits into two primes swapped by conjugation, `λ` has valuations
`±2` there, and `λ + λ̄ = 49/25` puts `e/25` in the module beside `e`.

| graph | coset colouring mod 5 |
|---|---|
| 803-graph | exists (1 728 admissible functionals) |
| E–I's `H`, unit edges | exists |
| E–I's `K = H ∪ λ_A(H)`, unit edges | **none** — an edge vector lies in `5M` |
| `803 ∪ λ(803)` | **none**; none mod 2, 3, 4; no periodic 5-colouring found through ℤ/6…12, 20, 25, 50, 125 |

So the Moser field is **not** closed to a sixth colour. The Moser rotation's
denominator 6 kills the coset colourings mod 2 and 3; de Grey's closing
rotation (`cos 31/32`) has denominator `2⁵`; the step to six needs a 5.

**Where a forced pair can live.** A coset colouring colours `u, v` alike iff
`ψ(u − v) = 0`. Enumerating all `5⁸` functionals on the rank-8 edge modules,
the admissible ones span the whole dual in both the 803-graph and the tight
4159-graph, so their common kernel is exactly `5M`: **every pair whose
difference is not five times a module element is split by an explicit coset
colouring**, and no growth along the graph's own directions can force it —
the one-line reason every forced-same hunt at five colours here came back
empty. What survives is `u − v = 5e`, distance 5 along a unit direction:
precisely Exoo and Ismailescu's pair, closed by precisely their rotation `λ`,
whose 5 in the denominator lets the union escape the coset colourings. In this
field their construction is the only shape a forced-pair spindle to six can
take. (The tight graph's 24 such pairs are all split by ordinary colourings.)

**And a gate for pairs.** Colourings through larger quotients *can* split a
`5M` pair. Exoo and Ismailescu's coordinate module holds exactly 30 unit
vectors (rank 4); enumerating every `ψ : M → ℤ/n` nonzero on all of them, the
Cayley graphs over ℤ/10, 15, 20, 25, 30, 40 and 50 all carry 5-colourings that
split a `5e` pair — and `ψ = (1,5,3,6)` into ℤ/10 splits *their own* `A, B`,
colouring every unit-distance graph on the module. Their forcing lives
entirely in the 446 two-edges. The 803-graph's module fails the same way
(ℤ/10, 20, 50). So:

> a pair `(A, A+5e)` can be forced by unit distances on a module `M` only if
> no finite quotient `ψ : M → G`, nonzero on every unit vector, has a
> 5-colouring of `Cay(G, ψ(units))` with `c(0) ≠ c(ψ(5e))` — the **pair
> gate**, to a forced pair what blocking is to a sixth colour.

**The ladder, and what it says the next rung is.** A pair `u, u + ke` (`e` a
unit step) lies in `kM`, so every coset `k`-colouring colours it alike; the
rotation closing distance `k` has `cos = 1 − 1/(2k²)` and denominator `2k²`,
which carries every prime of `k`:

| k | pair | `4k²−1` | closing rotation | used by |
|---|---|---|---|---|
| 2 | `2e` | 15 | `(7 + √−15)/8` | de Grey's `Sb` |
| 4 | `4e` | 63 = 9·7 | `(31 + 3√−7)/32` | de Grey's last step, `(±2, 0)` |
| 5 | `5e` | 99 = 9·11 | `(49 + 3√−11)/50` | Exoo–Ismailescu |

At `k = 5` the rotation lives in the Moser field with no new square root.

**The apart gate, and one surviving direction.** Mirror of the pair gate: a
`2e` pair can be forced *apart* only if no periodic 5-colouring keeps it alike.
On E–I's module ℤ/4 keeps every `2e` pair alike — so neither half of a
unit-distance version of their proof can live there. On the 803-graph's module
the pair gate falls in all 31 directions (ℤ/10, ℤ/15) and the apart gate in 30
at once; the `√247` edge `e = (−3/16, √247/16)`, which lies in `8M` so that
every quotient with a 2-part of order at most 16 is blind to `2e`, held out
through ℤ/4…25 — and then fell to ℤ/15 (and 30, 35, …, 60). *(Corrected: it
was first reported here as the survivor.)*

**Why one module is never enough.** Pairs no coset colouring keeps alike come in
whole distance classes on the 803-graph (`d² = 7/3, 1/9, 4/9, 19/9, 4, …`), but
the *same* coset colouring splits all of them at once: inside one integral
multiquadratic module, gadgets for any set of distances can never combine to
six. The step out needs a 5 in a denominator — and Exoo–Ismailescu's `K`
supplies it, so a distance-2 gadget in an *integral* module would still finish
the proof. That gadget is what the growth runs are now hunting, along `√247`
in the 803 module and on the blocked `803 ∪ λ(803)`, with
`scripts/verify_gadget.py` waiting for any UNSAT.

**Twisted coset colourings.** `c(p) = ψ(p) + t·⌊φ(p)/L⌋ (mod 5)`, with `ψ` an
admissible coset colouring and `φ` any integer functional on the edge module:
a unit step moves the level by at most one, so `c` is proper on every graph
with the module's directions exactly when the class `{u : ψ(u) = t}` lies in the
half-space `φ ≥ 0` — one LP, by Farkas. It is neither periodic nor a coset
colouring, and it is what the gates could not see: it **splits** `5e` pairs and,
with `t = −2ψ(e)`, **keeps `2e` pairs alike**. On the 803-module all 31
directions are refuted for both; on Exoo–Ismailescu's, all 15. Built
explicitly on the gadget search's own grown graph (2 185 points, still on the
same 31 directions): **0 improper edges of 12 883, and `c(a) = c(b)`**. No
distance-2 gadget and no forced `5e` pair can live on these modules at any size;
the integral searches were stopped. What survives is a module with **no**
admissible `ψ` — a 5 in a denominator, as in `803 ∪ λ(803)`, where the search
continues — or one whose directions are so rich that every class positively
spans.

**Credit.** Homomorphic ("coset") colourings as an obstruction are Polymath16's:
Philip Gibbs showed in 2018 that `ℤ[ω₁,ω₃,ω₄]` and `ℤ[ω₁,ω₃,ω₄,ω₇]` have
homomorphic 5-colourings, so no 6-chromatic unit-distance graph lies in either
ring (summarised in A. P. Goucher's *cp4space* post "Royal Wedding and
Polymath16"). What this package adds on top — as far as a search of the public
record shows, which is not a guarantee — is the integrality-at-5 criterion and
the residue-degree theorem, the escape through a 5 in a denominator and its
identification with Exoo–Ismailescu's rotation, the `5M` kernel and the spindle
ladder, the pair and apart gates, and the twisted coset colourings.

## κ, the rotation every coset colouring is blind to

The last paragraph left one way to reopen an integral module: directions so
rich that every class `{u : ψ(u) = t}` positively spans. One rotation turns
out to supply that.

`ρ₇ = (1 + 4√−3)/7 = (2 + √−3)/(2 − √−3)` is the rotation of `ℚ(√−3)` at the
prime 7. Polymath16's `ω₇ = (13 + 3√−3)/14`, the chord-1 turn at radius `√7`,
is `ω⁻¹ρ₇`. Mod 5 we have `1/7 = 3`, so `ρ₇ = 3(8ω − 3) = 1 − ω = ω⁻¹`. Hence

    κ = ω·ρ₇ = ω²·ω₇ = (−11 + 5√−3)/14,        κ − 1 = 5·(ω − 3)/7

is an irrational rotation (`cos = −11/14`, about 141.8°, not a root of unity)
that is **congruent to 1 mod 5**. Every homomorphism `ψ : M → ℤ/5` therefore
has `ψ(κx) = ψ(x)` whenever `(ω − 3)x/7` lies in `M`.

**The check.** `five_rho7` is the graph `five_247_c` together with its `ρ₇`
and `ρ₇⁻¹` images about its densest vertex: rank 8, 93 directions and 960
admissible `ψ`.
- In 120 of its 186 unit vectors `u`, `κu` is again a unit vector.
- `ψ(κu) = ψ(u)` in all **115 200 of 115 200** cases `(ψ, u)` (`scripts/kappa_check.py`).
  In fact `κu − u ∈ 5M` for each of the 120, whatever `ψ` is (`tests/test_kappa.py`).
- Neither `ω` nor `ρ₇` alone is even scalar on a single `ψ`. Only the
  product is invisible.

**Why that kills the twisted colourings.** A twisted colouring needs a real
functional `φ ≥ 0` on a class `D_t`. Suppose `D_t` were κ-invariant, as in any
`ℤ[ω, 1/7]`-module. Then the cone of such `φ` would be invariant under a
rotation that turns every complex embedding by an irrational angle. The Haar
average of an orbit is 0, so the cone contains the orbit's span, and `φ`
vanishes on `D_t`. Once `D_t` spans, `φ = 0`.

A finite module is only partly κ-stable, so on `five_rho7` this is checked
rather than assumed:
- **Apart gate:** there are exact rational Farkas certificates
  `−e ∈ cone(D_{−2ψ(e)})` in **all 178 560 cases** (93 directions × 960 `ψ` ×
  2 orientations, `scripts/farkas.py`).
- **Every `t`:** an LP sweep finds no twisted colouring in any direction, for
  either gate. This part was checked in floating point, not with exact
  certificates.

**Everything cheap we tried is a coset colouring there.**
- **Periodic apart gate** through `ℤ/n`, for `n = 10, 12, 15, 16, 20, 25, 30,
  40, 50`:
  - The search is exhaustive up to 300 000 `ψ` and sampled beyond that.
  - 5-colourable quotients appear only for the `n` divisible by 5, about 750
    of them each, and they are the coset colourings again.
  - None keeps a `2e` pair alike.
- **Full quotients `M/qM`** do not help either (`scripts/quotgate.py`,
  `scripts/unitsinqm.py`):
  - A unit vector lies in `qM` for `q = 2, 3, 4, 7, 8`. The unit contents are
    3, 7, 8, 21 and 56.
  - So no 2M-, 3M-, 4M-, 7M- or 8M-periodic colouring exists at all.

Coset colourings split every `2e` pair and join every `5M` pair, which is
exactly what a distance-2 gadget or a forced pair needs. So `five_rho7` is the
first module here where neither route is refuted by any colouring we tried.
That is not a proof. Exotic colourings are not excluded, and only a finite
UNSAT certificate would settle the question.

Two searches now run there:
- the gadget, with `a, a + 2e` forced apart;
- the forced pair, with `a, a + 5e` forced alike.

**Credit, again.** `ρ₇` is Polymath16's `ω₇` up to a sixth root of unity, and
`five_rho7` lies in their family `ℤ[ω₁, ω₇, …]`. That family has homomorphic
5-colourings, which is why only gadgets and forced pairs are left to hunt. The
congruence `κ ≡ 1 (mod 5)` and its use against the twisted colourings are, as
far as we know, new here. Polymath16's third thread
(D. Mixon's blog, 2018-05-01, "Is 6-chromatic within reach?") asked whether
`ℤ[ω₁,ω₃,ω₄]` admits a homomorphic 5-colouring. We found no discussion in the
public record of whether *every* 5-colouring of such a ring must be
homomorphic, and none of strip colourings or of rotations `≡ 1 (mod 5)`. The
rigidity theorems below are, to our knowledge, new. That is a search result,
not a guarantee.

**Any pair, not one pair.** The growth scripts ask about a single chosen pair.
But **any** non-adjacent pair forced alike proves `χ ≥ 6`: chain translated
copies until the gap is at least ½, then close with the spindle rotation whose
chord at that radius is 1. Likewise, any pair at distance 2 forced apart proves
`χ ≥ 6` through Exoo–Ismailescu.

`scripts/backbone.py` asks about all pairs of a grown graph at once:
- it filters with SAT colourings in random vertex order, plus Kempe chains
  (a pair survives only if no single Kempe swap can split it, or join it);
- it then settles the survivors with one incremental disjunctive SAT query:
  "some survivor breaks".

### Every other structured colouring we could think of, on both modules

**The twisted colourings were periodic after all, through a non-cyclic quotient.**
`ψ + t⌊φ/L⌋` repeats when `φ` moves by `5L`, so it factors through
`ℤ/5L × ℤ/5`. That quotient is not cyclic, which is why the `ℤ/n` gates never
saw these colourings.

The **product gate** (`scripts/prodgate.py`) works as follows:
- It samples a homomorphism `χ : M → ℤ/n₁ × … × ℤ/n_k`, often with one factor
  an admissible `ψ`.
- It looks for **any** proper 5-colouring of `Cay(A, χ(U))` that joins a `2e`
  pair or splits a `5e` pair, not only staircases.

On `five_rho7` it drew 4 000 samples of each of:
- `ℤ/5 × ℤ/5 × ℤ/5`;
- `ℤ/5 × ℤ/10`, `ℤ/5 × ℤ/15`, `ℤ/5 × ℤ/20`, `ℤ/5 × ℤ/25`;
- `ℤ/5 × ℤ/5 × ℤ/2` and `ℤ/5 × ℤ/5 × ℤ/3`.

Every sampled quotient is 5-colourable, and **none refutes a direction**. The
cyclic pair gate (`ℤ/10 … ℤ/50`) has finished too: all 93 directions are open.
*(That gate sampled homomorphisms; its verdicts for `n ≥ 10` are withdrawn as
evidence in "The circular gate: colourings through a real character", below.)*

**Kempe swaps of coset colourings are trivial.** In `ψ + k`, take the Kempe
component of `x` in colours `α` and `α + δ`. It is `x + L ∪ x + p₀ + L`, with
`L = ⟨p − p′ : ψ(p) = ψ(p′) = δ⟩`. A swap can join a `2e` pair or split a `5e`
pair only if `[ker ψ : L] > 1` (`scripts/kempegate.py`). That index is 1 for
every admissible `ψ` and every `δ`:
- on `five_rho7`, 3 840 of 3 840 cases;
- on `five_247_c`, 6 912 of 6 912 cases.

So a swap only relabels a whole coset.

**The blocked module has no known colouring at all.** In `803 ∪ λ(803)` some
unit vectors lie in `2M`, `3M`, `4M`, `5M`, `8M` and `15M`; the 5 is the
denominator of `λ`. The sampled searches found no 5-colourable periodic
colouring:
- through `ℤ/n` for `n = 6 … 25`;
- through `ℤ/2 × ℤ/2`, `ℤ/3 × ℤ/3`, `ℤ/4 × ℤ/4`, `ℤ/6 × ℤ/6`, `ℤ/7 × ℤ/7`,
  `ℤ/3 × ℤ/9`, `(ℤ/2)³`, `(ℤ/3)³`, `(ℤ/2)⁴` and three mixed products.

Most of these quotients are not even loopless. Whether the unit-distance graph
on that module is 5-colourable at all is **open**; if it is not, it already
contains a finite 6-chromatic graph. A plain growth runs there, beside the
gadget search.

**Local search first.** `scripts/grow_ls2.py` repairs the previous colouring
with a small C tabu search (`scripts/tabucol.c`) and calls the CDCL solver only
when that fails. The pool stores only indices, which uses a third of the
memory.
- **Speed:** growth runs about fifty times faster than with CDCL alone. On
  `five_rho7` the gadget search reached 6 600 points in 36 seconds, where the
  CDCL-only run needed 21 minutes to reach 5 900.
- **Where the tabu search gives out:** about 6 600 points for the gadget and
  10 200 for the forced pair. From there on every step is a CDCL call, and
  that is where an UNSAT would show up.

## Exoo–Ismailescu's graph lives inside our modules

Their 30 unit vectors, turned by 150°, are all unit vectors of `five_247_c`,
and therefore also of `five_rho7` and of the blocked `803 ∪ λ(803)`.

- **`H` inside `five_rho7`.** Their `H` (214 points, 1 004 unit edges, 446
  two-edges) is turned by 150° and placed with `A` at `five_rho7`'s densest
  vertex. It then lies in `five_rho7`'s module (`scripts/mk_kseed.py`).
- **`K` inside the λ-closure.** `K = H ∪ λ_A(H)` lies in
  `five_rho7 ∪ λ(five_rho7)`: 426 points, `|BB′|² = 1`, 892 two-edges.
- **Re-verified in place.** `H` with its two-edges and `c(A) ≠ c(B)` is UNSAT
  for cadical, glucose and minisat.

On `five_rho7` no colouring we tried splits a `5e` pair, so E–I's pair
`A, A + 5e` is an open target there. A unit-distance graph in `five_rho7`
forcing `c(A) = c(B)` would close with `λ` into a **6-chromatic unit-distance
graph**.

**The skeleton search** (`scripts/grow_kw.py`) starts from `five_rho7 ∪ H`. A
two-phase tabu search (`scripts/tabu2.c`) looks for a colouring:
1. First it finds a proper colouring with `c(A) ≠ c(B)`.
2. Then, without breaking a unit edge, it reduces the number of **alike**
   two-edges as far as it can.

Rainbow points are inserted next to the alike two-edges that remain. The count
of alike two-edges is a progress measure. Its true minimum can only rise as the
graph grows; the tabu search gives an upper bound, not the exact minimum.

| points | alike two-edges |
|---|---|
| 2 606 | 26 |
| 2 806 | 39 |
| 3 006 | 47 |
| 3 206 | 50 |

**Two sobering measurements.**
- **Finite-scale rigidity is invisible.** Colour a graph by tabu search and
  compare it with the best coset colouring (`scripts/rigid.py`,
  `scripts/cosetfit.py`). The share of matching vertices is near chance, which
  is about 20–25%:
  - `five_rho7`: 24%;
  - the vertex-critical 803-graph: 27–30%;
  - a grown gadget graph of 10 141 points: 22%.

  So the coset colourings may be the only structured colourings of the infinite
  module graph, but finite graphs of these sizes do not feel them at all.
- **Greed is weak.** Pointed at a problem whose answer is known (E–I's
  two-distance forcing, 214 points), the greedy growth had not found it by 344
  points. At that size a single SAT call already needed 2.4 million conflicts.

  Structure beats greed, and the skeleton search is the attempt to give the
  growth some.

### Coarse rigidity, exactly

A twisted colouring needs a functional `φ ≠ 0` with `φ ≥ 0` on a class
`D_t = {u : ψ(u) = t}`. **Stiemke's lemma** says no such `φ` exists iff `D_t`
spans and has a strictly positive linear dependency.

`scripts/stiemke.py` finds such a dependency and proves it exactly:
1. An LP finds one.
2. All coefficients but a basis are rounded to rationals.
3. The basis coefficients are then solved **exactly**.
4. The identity `Σ λ_u u = 0`, with every `λ_u > 0`, is checked in the field
   coordinates.

| module | admissible `ψ` | `(ψ, t)` with an exact certificate |
|---|---|---|
| `five_rho7` (93 directions) | 960 | **3 840 of 3 840** |
| `five_247_c` (31 directions) | 1 728 | 1 664 of 6 912; the other 5 248 lie in a half-space |

So on `five_rho7` there is no twisted colouring at all, now exactly and not only
by the LP sweep. `ρ₇` is precisely what makes the difference.

**Theorem (coarse rigidity).** Call a 5-colouring of the unit-distance graph on
`five_rho7`'s module *coarse* if two things hold:
- `M ⊗ ℝ = ℝ⁸` is cut into polyhedral cells whose faces are all much larger
  than a unit step;
- on each cell `P` the colouring is `ψ_P + k_P`, with `ψ_P` admissible.

Then every coarse colouring is a coset colouring.

*Proof.* Look at a facet `F` between two cells.
1. Suppose `ψ_P` and `ψ_Q` differ. Their difference can be constant along the
   lattice points next to `F` only when:
   - `F` is rational, `F = {φ = 0}`, and
   - `ψ_Q − ψ_P = λφ`.

   In that case take a unit `u` with `φ(u) = s`. Crossing from level `−j` to
   level `s − j` changes the colour by `ψ_P(u) + δ + λ(s − j)`. Once `s ≥ 5`
   this runs through every residue, so some `j` gives a conflict. And every
   nonzero integer functional takes a value of at least 7 on some unit: none has
   `|φ(u)| ≤ 6` on all 186 units (Fincke–Pohst on the LLL-reduced value
   lattice).
2. So `ψ_P = ψ_Q`, and the colour jumps by some `δ ≠ 0`. Then any unit of
   `D_{−δ}` that crosses `F` is a conflict. Such a unit exists because
   `D_{−δ}` positively spans (the certificates).

Hence no facet carries a jump, and a single `(ψ, k)` covers the whole module. ∎

**What it does not say.** A colouring can still fail to be a coset colouring in
two ways:
- its interfaces are *thick*, a slab where the colouring looks nothing like
  `ψ + k`;
- it has no large-scale structure at all.

The 1-D wall gate looked for the first kind (`scripts/wallgate.py`):
- It searched along the 24 shortest integer functionals, with values up to 14
  on the units, for every admissible `ψ`, in windows 161 levels wide.
- No window colouring joins a `2e` pair; 12 functionals were checked before
  the time limit.
- None splits a `5e` pair; 8 functionals were checked.

### Rigidity propagates

**Theorem.** Let `ψ` be admissible on a module `M`, and suppose every class
`D_t = {u : ψ(u) = t}`, `t = 1…4`, positively spans `M ⊗ ℝ`. (For `five_rho7`,
that is exactly the 3 840 certificates above.) Suppose a proper 5-colouring
`c` agrees with `ψ + k` on some half-space `{φ < a}`. Then `c = ψ + k` on **all**
of `M`.

*Proof.* Let `δ₀ = min_t max_{u ∈ D_t} (−φ(u))`. It is positive, because a class
that positively spans contains some `u` with `φ(u) < 0`. Suppose `c` is already
known on `{φ < a + mδ₀}`, and take `x` with `φ(x) < a + (m + 1)δ₀`. For each `t`,
pick `u_t ∈ D_t` with `φ(u_t) ≤ −δ₀`. Then `x + u_t` lies in the known region and
has colour `ψ(x) + t + k`. Since `x` is adjacent to all four, the only colour
left for `x` is `ψ(x) + k`. Induct on `m`. ∎

**The ball version.** `M` is a lattice in `M ⊗ ℝ = ℝ⁸`, so a ball contains only
finitely many of its points. The inradius of `conv(D_t)` about 0 has a positive
minimum over all `(ψ, t)`, and the same induction runs outward from any
sufficiently large ball: if `c = ψ + k` on `B(z, R) ∩ M` with `R ≥ R₀`, then
`c = ψ + k` everywhere.

**Why it matters.** Coset colourings are **infectious**. They cannot sit next to
anything else, across a sharp wall or a thick one: coarse rigidity is just the
half-space case. So a 5-colouring of the `ρ₇` module that is not a coset
colouring must be *nowhere locally coset*: it disagrees with every coset
colouring on every ball of radius `R₀`.

Whether such colourings exist is exactly the question that remains. If none
does, every `5M` pair is forced alike and every `2e` pair is forced apart, and
`χ(ℝ²) ≥ 6` follows by the λ-closure.

On the 803 module the hypothesis fails for 5 248 of 6 912 classes. There the
twisted colourings are precisely the colourings that are one coset colouring on
a half-space and another beyond it.

**The scale, measured.** The balls here are balls of ℝ⁸, so they hold about
`r⁸` points. We coloured the complete combinatorial 2-ball of `five_rho7`
(every point two unit steps from a vertex: 17 305 points, 68 796 edges) by tabu
search (`scripts/ball2.py`). Its best coset fit is 22.8%, which is chance.

So nothing near that size is rigid. The propagation needs balls much larger
than anything a solver can hold, which is why every finite graph so far has
nowhere-coset colourings and forces nothing.

**How far away is that?** The propagation step needs a unit of every class
pointing back by at least the inradius of `conv(D_t)`.

In the trace form every unit has the same length, so `U` sits on a sphere in
ℝ⁸. On 160 sampled classes (qhull), that inradius is:
- 0.057 of the unit radius at worst;
- 0.092 at the median.

So the induction needs a starting ball of about **nine unit steps**. A ball of
ℝ⁸ that large holds astronomically many lattice points; the module even
contains vectors as short as `u/56`.

Hoffman's bound gives no shortcut. The mod-5 quotient `Cay(F₅⁸, U)` has least
eigenvalue −64.58, not `−d/4 = −46.5`, so it only bounds the independence ratio
by 0.258, not by 1/5.

The rigidity is real, but it acts far beyond what a solver can hold. A finite
gadget in this module, if one exists, has to come from a cheaper mechanism.

## Periodic through an ideal, and the weakest hypothesis that suffices

**The module is a ℤ[ω]-module, with more unit vectors than the graph uses.** The
186 edge directions of `five_rho7` are not closed under the 60° turn. Turning
them adds 12 more unit vectors, and all 12 lie in the module (exact lattice
test), so the module holds at least **198** unit vectors. All 960 admissible
`ψ` stay nonzero on the new ones. Enumerating every unit vector is out of reach:
the trace form has minimum about 0.0003, so Fincke–Pohst would visit about 10¹²
lattice points.

**Periodic through a small ideal, and why it proves little.** Because the
module is a ℤ[ω]-module, `αM` is a period lattice for every `α = a + bω`, of
index `N(α)⁴`, and `scripts/idealquot.py` asks kissat for a 5-colouring of
`Cay(M/αM, U)`. For `α = 3 + ω, 4 − ω` (norm 13) and `3 + 2ω, 2 + 3ω` (norm 19)
there is none, on `five_rho7` and on its `λ`-closure. **But the reason is
trivial.**
- The UNSAT core of the norm-13 torus has **15 vertices and 75 edges**
  (`scripts/torecore.py`). It is a small, dense graph that exists only because
  the short period folds far-apart points of the module together. It does not
  lift back to the module: 53 of its 75 edges wrap around the torus.
- The pair queries on `M/(ker ψ ∩ αM)` (`scripts/idealquot2.py`, 30 queries
  UNSAT) are artifacts in the same way. The split of `(0, 5e)` is already
  impossible on the radius-1 torus ball around the pair, about 385 vertices.
- The twisted-ideal family (`scripts/twistquot.py`) and the earlier small
  cyclic and product gates are exposed to the same folding.

These tests say little about rigidity. A periodic test only means something
when the period is long compared with the unit steps.

**Unbiased colourings do not feel the rigidity.** `scripts/ballr.py` builds the
complete combinatorial 2-ball in exact lattice coordinates (17 305 points) and
colours it by tabu search. For each of the 181 centres with nearly all their
neighbours present, it compares pairs `x + u, x + v`. The pairs with
`u − v ∈ 5M`, which every coset colouring colours alike, are alike **25.6%** of
the time. All other pairs are alike **25.6%** of the time. The grown skeleton
graphs of 3 k to 12 k points are worse still: their colourings split the `5M`
pairs more often than chance (15% alike against 30%), because the growth's own
colouring search is built to split things.

**The hypothesis can be weakened to one about measures.** Here is the precise
statement that would finish the proof. Let `M₁` be the `ρ₇`-module, `e` one of
its unit directions, and `λ` Exoo–Ismailescu's rotation, with `|5e − 5λe| = 1`.

> **Theorem.** Suppose every `M₁`-translation-invariant probability measure on
> proper 5-colourings of `Γ(M₁)` satisfies `c(0) = c(5e)` almost surely. Then
> `Γ(M₁ + λM₁)` is not 5-colourable, and so `χ(ℝ²) ≥ 6`.

*Proof.*
1. Let `c` be a 5-colouring of `Γ(M′)`, where `M′ = M₁ + λM₁`. Averaging its
   translates over Følner boxes of `M′ ≅ ℤ⁸` gives an `M′`-invariant measure
   `ν` on proper colourings.
2. Restricted to `M₁`, `ν` is `M₁`-invariant, so `c(0) = c(5e)` holds
   `ν`-almost surely.
3. Pushed forward by `c ↦ (y ↦ c(λy))`, `ν` is again `M₁`-invariant. A
   translation by `m` becomes a translation by `λm ∈ M′`. So
   `c(0) = c(5λe)` also holds `ν`-almost surely.
4. Hence `c(5e) = c(5λe)` holds `ν`-almost surely. But `5e − 5λe` is a unit
   vector of `M′`, so no proper colouring does this. ∎

So even colourings that are not coset colourings are harmless, provided they
split `5e` pairs only on a set of density zero.

**Why this is hard, measured against a known open problem.** Suppose that
every 5-colouring of `Γ(M)` splits every pair `x, x + m` with `m ∈ u₀ + 5M`.
Coset colourings do exactly this. Then no *measurable* 5-colouring `f` of the
plane can be proper along the finitely many unit directions of `M`.

*Proof.*
1. The coset `u₀ + 5M` is dense in the plane, so it contains vectors `mⱼ → 0`.
2. For almost every offset `y`, `m ↦ f(y + m)` is a proper colouring of
   `Γ(M)`, so `f(y + mⱼ) ≠ f(y)`.
3. Continuity of translation in `L¹` says the opposite, on a set of full
   measure. ∎

So the rigidity of any finite-direction module would already imply that the
**measurable** chromatic number of the plane is at least 6, which is open.
Conversely, a measurable 5-colouring along a module's directions would refute
its rigidity. The torus gate searched for pixelated ones along the 803
directions and found none.

**The blocked growth's colouring, looked at directly.** The plain growth on
`803 ∪ λ(803)` ran out of rainbow points: the colouring extends to every
candidate point without a conflict. We checked whether that colouring hides a
period (`scripts/periodfind.py`, `scripts/quot2d.py`).
- Of 1 890 translations `u ± v`, `ku`, the best leaves colours unchanged
  **87%** of the time. The median is 46%, against about 20% by chance.
- The 42 best translations span a rank-6 lattice `N`. But `N` contains 60 unit
  vectors, so no colouring factors through `M/N`.
- The colour classes show no regions in the plane, nor in any of the three
  conjugate planes.

The colouring is strongly organised, but not along anything we can name.


## The circular gate: colourings through a real character

A **circular colouring** is `c(x) = ⌊5·frac(φ(x))⌋`, for a character
`φ ∈ Hom(M, ℝ/ℤ)`. It is proper exactly when `frac(φ(u)) ∈ [1/5, 4/5]` for
every unit vector `u`: a shift by at least a fifth, and at most four fifths,
always changes the fifth of the circle a point lies in, even at the boundary.

Coset colourings are the special case `φ = ψ/5`. So the **feasible region**

    F(M) = { φ ∈ Hom(M, ℝ/ℤ) = T^r : frac(φ(u)) ∈ [1/5, 4/5] for all u ∈ U }

contains every coset colouring as a 5-torsion point. By the Stiemke
certificates those points are isolated, but `F(M)` can have other points.
- A point of `F(M)` with rational coordinates gives a periodic colouring
  through a cyclic quotient `ℤ/q`.
- An interior point gives irrational `φ` nearby, and so quasi-periodic
  colourings.

`scripts/circgate.py` maximises the least slack
`s = min_u dist(frac(φ(u)), outside [1/5, 4/5])` as a MILP (HiGHS). Its
variables are `φ`, one integer `n_u` per direction, and `s`.
`scripts/circverify.py` then re-checks a solution in exact rational
arithmetic: first on all unit vectors, then on every edge of a grown graph.

| module | result |
|---|---|
| the 803-graph, **4** colours (sanity) | infeasible, as it must be for a 5-chromatic graph |
| E–I's `H` | `φ` of order 4, slack 0.05, and it splits their pair |
| blocked `803 ∪ λ(803)` (144 units) | **slack 387/31250 > 0**. The graph grown there, 46 496 points and 453 731 edges, is coloured with **0** monochromatic edges |
| `five_rho7` (198 units) | **a non-coset `φ` of order 280**, slack 0. It colours the 12 226-point skeleton graph with 0 monochromatic edges, and `φ(5e) = 96`, so it colours Exoo–Ismailescu's pair **alike** |

Explicitly, in the LLL basis of `circgate.py`,

    φ = (11/20, 9/20, 39/56, 9/140, 1/5, 3/4, 69/280, 0).

`5φ` is not integral, so this is not a coset colouring: it is periodic through
`ℤ/280`, a modulus the cyclic gates, which stopped at 50, never reached.

**What this changes.**
1. **The blocked module was never a candidate.** A circular colouring 5-colours
   every graph along its directions. This is why its growth ran out of pressure,
   and why its colouring agreed with itself under small translations.
2. **The rigidity conjecture, in its strong form, is false.** `five_rho7` has a
   5-colouring that is not a coset colouring. The propagation theorem stands:
   this colouring is nowhere locally coset, exactly as the theorem requires.
3. **And the Exoo–Ismailescu pair does not survive a small move.** The order-280
   colouring colours `(A, A + 5e)` alike. It is not isolated in `F(M)`, though:
   the units tight at it have rank 7, and the local cone contains
   `θ = (0, 0, −1, 1, 0, 0, 1, −1)`. The moved character

       φ′ = φ + θ/300 = (11/20, 9/20, 2911/4200, 71/1050, 1/5, 3/4, 1049/4200, −1/300)

   keeps every one of the 198 units in `[1/5, 4/5]`, with exact slack 0. It
   colours the 12 226-point skeleton graph (`five_rho7 ∪ H`, grown) with **0**
   monochromatic edges, and `φ′(5e) = 1927/20`, so `frac = 7/20`. It therefore
   gives **`c(A) ≠ c(B)` on every translate of the pair.**
   - **No unit-distance graph inside the `ρ₇` module can force Exoo–Ismailescu's
     pair.** The skeleton search there was futile, and it has been stopped.
   - The measure hypothesis of the closure theorem fails too.
   - Neither the order-280 character nor `φ′` extends to the `λ`-closure: none of
     the 390 625 extensions is proper (`scripts/circextend.py`). Whether the
     `λ`-closure has any circular colouring at all is the MILP still running.

**A necessary condition, stronger than blocking.** Every unit-distance graph
whose edge module `M` has `F(M) ≠ ∅` is 5-colourable. So a 6-chromatic graph
needs a module with **no homomorphic circular 5-colouring**. No such module has
been found yet.

By compactness, the condition also bears on the other side of the problem. The
characters of the plane, seen as a discrete group, form a compact group, and
characters of a subgroup extend to the whole group. So: *if every finite set
of unit vectors admits a character with all values in `[1/5, 4/5]`, then the
plane has a homomorphic 5-colouring and `χ(ℝ²) = 5`.* A single finite set of
unit vectors with `F = ∅` is therefore a necessary first step toward six, and a
test that the upper bound 5 must also pass.

### The `λ`-closure against the circular gate

**Its unit set had been under-counted.** For every unit `e` of `M₁`,
`w_e = 5(1 − λ)e = 5e − 5λe` is a unit vector of `M′ = M₁ + λM₁`. These are
exactly the edges Exoo–Ismailescu's argument uses. The grown graphs had
carried only 420 unit vectors. The module has at least **594**: 198 of `M₁`,
198 of `λM₁` and 198 `w_e` (`scripts/lamclosure_units.py`). Among the 15 430
grown points the missing directions add only 59 edges. Even so, the growth was
restarted on the full set of 594 directions (`data/rho7_lamfull_kw1.json`).

**What a circular colouring of `M′` would have to look like.** Let
`φ = φ′|M₁` and `φ̃ = φ′∘λ|M₁`. The `w_e` force

    χ = 5(φ − φ̃)  ∈  F(M₁)

on top of `φ, φ̃ ∈ F(M₁)`.
- **Neither restriction can be a coset colouring.** If `φ = ψ/5`, then
  `χ ≡ −5φ̃`, so `φ̃` and `5φ̃` would both be circular colourings of `M₁`. The
  MILP with the multiples `{u, 5u}` of every unit is **infeasible**
  (`MULT=1,5 scripts/circgate.py`, HiGHS, 25 min). The same holds with `φ` and
  `φ̃` exchanged.
- **The one non-coset cell we know does not extend.** The cell of `F(M₁)` around
  the order-280 point is a short segment in the direction `θ`: 9 vertices, and
  every coordinate varies by at most 0.009. We tested 200 samples along it, each
  with all 390 625 extensions to the 420-unit closure. None is proper; the best
  slack is −0.14 (`scripts/cellsample.py`).

The complete answer, whether `F(M′) = ∅` on all 594 units, is the MILP still
running.

**A second non-coset colouring, and a correction to the old cyclic gates.** At
`Q = 280`, CP-SAT (`scripts/circsat.py`, the circular gate as a pure integer
model, `φ = y/Q`) found in 47 seconds

    φ₂₀ = (3/10, 13/20, 3/5, 1/20, 1/5, 3/4, 0, 1/4),

a character of order 20. The exact check shows:
- every unit lies in `[1/5, 4/5]`;
- the 12 226-point skeleton graph is coloured with 0 monochromatic edges;
- `frac(φ₂₀(5e)) = 1/4`, so it too splits Exoo–Ismailescu's pair on every
  translate.

It is periodic through `ℤ/20`, a quotient the old cyclic pair gate reported as
refuting nothing. That gate **sampled** homomorphisms `M → ℤ/n`, and there are
`20⁸ ≈ 2.6·10¹⁰` of them for `n = 20`. Its verdicts for `n ≥ 10` are therefore
withdrawn as evidence. Solver-based gates (MILP, CP-SAT) replace sampling from
here on.

## The gate in relation space, and the local criterion

**An exact, fast form of the circular gate.** Take as unknowns the values
`t_u = frac(φ(u)) ∈ [1/k, 1 − 1/k]`, one per direction. A vector `t` comes from
a character exactly when `Σ a_u t_u ∈ ℤ` for every integer relation
`Σ a_u u = 0`, and it suffices to impose this on a basis of the relation
lattice. `scripts/circrel.py` computes that basis by LLL (python-flint, on
`[I | W·D]`). On the `λ`-closure of `five_rho7` it has 289 relations, all with
coefficients in `{−1, 0, 1}` and at most 20 terms. The MILP then has small
integer ranges only. With `SOLVER=SCIP` it runs through OR-Tools' SCIP, an
independent branch and bound.

It maximises the least slack, so it returns `κ(U) = max_φ min_u ‖φ(u)‖` exactly,
the lonely-runner constant of the unit set:

| unit set | rank | directions | `κ` (SCIP, optimal) |
|---|---|---|---|
| E–I's `H` (unit edges) | 4 | 9 | 1/4 |
| its `λ`-closure | 4 | 27 | 1/4 |
| rotation words `ω, σ, λ, ρ₇` (`rot_1110`) | 4 | 81 | 1/4 (34 s) |
| rotation words, `σ^{±2}` (`rot_2110`) | 4 | 135 | 1/4 (74 s) |

**Why always 1/4: the Moser field is 4-colourable.** Circular colourings are
*adelic*. The feasible region `F(M)` is a finite union of polytopes with
rational vertices, so if it is nonempty it contains a character of finite
order `N`. Such a character is a sum of `p`-adic pieces, one for each prime
`p | N`. One place can already colour a whole field:

> **Local criterion.** Let `K` be a CM field, so that `z z̄ = 1` for its unit
> vectors, and let `v` be a place of the real subfield `K⁺` that does not split
> in `K`. The norm-one group `T(K⁺_v)` is then compact, so every unit vector of
> `K` is a `v`-adic unit. If some additive `φ(z) = frac_p(L(z))`, with `L` a
> `ℚ_p`-linear form on `K_v`, keeps `‖φ‖ ≥ 1/k` on all of `T(K⁺_v)`, then
> `⌊k·φ(z)⌋` properly `k`-colours the unit-distance graph on **all of `K`**.

**Theorem.** The unit-distance graph on the Moser field `ℚ(√−3, √−11)` has
chromatic number exactly 4.

The proof takes three steps.
- **The field sits in a 2-adic field where 2 does not split.** Since
  `33 ≡ 1 mod 8`, `√33` is a 2-adic integer. The field therefore embeds in
  `ℚ₂(ω)`, the unramified quadratic extension, and 2 does not split there.
- **The unit vectors fall into six classes.** Its norm-one units are, mod 4,
  exactly the six sixth roots of unity.
- **One character separates them.** `φ(α + βω) = frac₂((α + 2β)/4)` takes only
  the values 1/4, 1/2 and 3/4 on them.

The Moser spindle gives the lower bound.

`hn/adelic.py` and `scripts/moser2adic.py` check this exactly, at both places
above 2:
- 1 236 unit vectors, words in `σ`, `λ`, `ρ₇` and `ω`;
- all 1 004 unit edges of Exoo–Ismailescu's `H`, with none monochromatic;
- `tests/test_moser_field.py`.

**No 5-chromatic unit-distance graph lives in the Moser field**, whatever its
size. The project's 5-chromatic graphs needed `√247`, and de Grey's needed
`√−15`. That is now forced.

**The denominator principle, explained.** A unit vector with `p` in its
denominator can exist only where a place above `p` splits, because a
non-split place has a compact torus. So:

| field | effect |
|---|---|
| `ℚ(√−3)` | 3 ramifies and does not split: `x mod √−3` 3-colours it |
| Moser field | 2 does not split: the 2-adic 4-colouring |
| `+ √−15` (de Grey's `2e` rotation, denominator 8) | 2 now splits, and the 4-colouring dies |
| `+ λ` (Exoo–Ismailescu, denominator 25) | 5 splits: no field-wide coset colouring mod 5 |

**Which places could colour.** At an odd place with residue field `𝔽_p`,
unramified and non-split, the character reduces to `Tr(ȳu)/p` on the `p + 1`
norm-one residues. Deeper levels always smear a full coset of `(1/p)ℤ` across
the forbidden arc: the tangent lines `ū√d·𝔽_p` take `(p + 1)/2` distinct
values, so no single kernel contains them all. The traces met by the norm
class `n` are the `c` with `c² − 4n` a non-residue or 0. That leaves a finite
list:

| `k` | residue primes `p` with a local `k`-colouring |
|---|---|
| 4 | 3, 7 |
| 5 | 3, 5, 7, 19 |
| 6 | 3, 5, 7, 11, 17, 19 |

The primes up to `10⁵` were checked by computer. Beyond that the Weil bound
excludes the rest.

**The field of `five_rho7` passes the local test for six.** In
`K = ℚ(√−3, √−11, √−247)`:
- 2 splits (because of `τ = (119 + 3√−247)/128`), and so do 3, 5, 7, 13 and 19.
- The non-split places lie above 11 (residue field `𝔽₁₁`) and above the primes
  with `(−3/p) = (−11/p) = (−247/p) = −1` (29, …).
- None of these is 3, 5, 7 or 19.

So no single place 5-colours `K`, and the universe of the search is not
excluded by this test. The criterion is only sufficient, though. A finite
module can still be coloured by characters that combine several places, split
places among them, as the order-280 and order-20 colourings of `five_rho7` do.
Deciding that is the relation-space MILP's job, now running on the 594 units of
the `λ`-closure (HiGHS at `k = 5, 6`, SCIP at `k = 5`).

Related prior work: G. E. Moorhouse, *On the chromatic numbers of planes*
(draft, 2010), studies `χ(K²)` over fields by reducing to finite fields.

## The field of `five_rho7` is 5-colourable: reduction at 11

**The result.** Every unit-distance graph with vertices in
`K = ℚ(√−3, √−11, √−247)` is 5-colourable, however large. That is the field of
the 803-graph, of `five_rho7`, of the blocked module and of the `λ`-closure. So
**no search in this field could ever have found six**, and every growth run
there was bound to colour.

**Why.** Characters were the wrong family to look at. The right family is
*reduction at a place*.

The place above 11 does not split:
- `K⁺_v = ℚ₁₁(√33)`, which is ramified, with uniformiser `π = √33`;
- `K_v = K⁺_v(√−3)`, which is unramified over it, because −3 is a non-residue
  mod 11.

So every unit vector `u` of `K` is a `v`-adic unit of relative norm 1. It
reduces to the norm-one group `N₁` of `𝔽₁₂₁` (12 elements).

Every point `z ∈ K` has a residue `g(z) ∈ 𝔽₁₂₁`, taken relative to its coset of
`O_v`, and `g(z + u) = g(z) + ū`. So any colouring of the **finite plane**
`Cay(𝔽₁₂₁, N₁)`, 121 points of degree 12, colours the whole field. Its
chromatic number is **5** (SAT; 4 is UNSAT).

In coordinates: write `z = α + β√−3` with `α = a + bπ` and `β = c + dπ`. The
residue is `(units digit₁₁(a), units digit₁₁(c))`, where
`a = a₀ + a₂√741` and `c = b₀ + b₂√741/3`, with `√741 ∈ ℤ₁₁`.

`hn/adelic.py` (`reduce11`, `finite_plane_11_colouring`) and
`scripts/reduce11.py` check it exactly, at both places above 11:

| graph | points | edges | directions | monochromatic |
|---|---|---|---|---|
| `five_rho7` | 2 403 | 12 223 | 186 | 0 |
| 803-graph | 803 | 4 065 | 62 | 0 |
| `λ`-closure growth (`rho7_lamfull_kw1`) | 15 630 | 119 788 | 594 | 0 |
| blocked growth | 32 312 | 282 909 | 134 | 0 |
| `λK` growth | 15 930 | 122 336 | 392 | 0 |
| `ρ₇`/`H` growth | 15 942 | 134 255 | 186 | 0 |

In each case every unit vector reduces into `N₁`.

**What the old gates missed.** On a finitely generated module this colouring
is periodic through a finite 11-group. The periodic gates tried Eisenstein
ideals (`ℤ[ω]`, norm 13) and sampled cyclic quotients. The circular gate tried
arcs of the circle. None of them tried the residue field of the one prime that
does not split.

**The general rule.** Take any place `v` of `K⁺` that does not split in `K`:
- **Unramified, residue field `𝔽_q`:** `χ(K) ≤ χ(Cay(𝔽_{q²}, N₁))`, the
  unit-distance graph of the finite plane.
- **Ramified:** the norm-one residues are `±1`, so `χ(K) ≤ 3`. A field that
  holds a Moser spindle has no such place.
- **In a field containing `√−3`:** a non-split place has `(−3/q) = −1`, so
  `q ≡ 5 mod 6` or `q = 2`.

| `q` | 2 | 5 | 11 | 17 | 23 | 29 |
|---|---|---|---|---|---|---|
| `χ(Cay(𝔽_{q²}, N₁))` | 4 | 4 | **5** | … | … | … |

The row for `q = 3`, 7 and 13, which cannot occur next to `√−3`, gives 3, 4 and
(not needed). The last three columns are being computed.

**What a field for six must look like.** It must contain `√−3`, and every
place above 2, 5 and 11 must split, as must any other `q` whose finite plane is
5-colourable. It can have no ramified non-split place. The unramified
non-split primes of some candidates:

| field | ramified | unramified non-split `p < 400` |
|---|---|---|
| Moser `ℚ(√−3, √−11)` | 3, 11 (**11 does not split**) | **2**, 17, 29, 41, … |
| `ℚ(√−3, √−11, √−247)` | 3, 11, 13, 19 (**11 does not split**) | 29, 83, 107, … |
| `+ √−7` (de Grey's `4e` rotation) | 3, 7, 11, 13, 19 (all split) | 83, 173 |
| `+ √−7, √−15` | 3, 5, 7, 11, 13, 19 | none |

Adjoining `√−7` (or `√−2`) puts `√21` in the real subfield. That is a
non-residue unit at 11, so the place above 11 splits.

## Finite planes, the field screen, and a new universe: `ℚ(√−3, √−7, √−11)`

**Finite planes.** Take `q ≡ 5 mod 6`, the case that occurs next to `√−3`.
The finite plane is `G_q = Cay(𝔽_{q²}, N₁)`, where `N₁` is the norm-one circle
(`q + 1` points).

The graph has a clean structure:
- **Two unit steps and their difference.** Two unit steps whose difference is
  also a unit step differ by a sixth root of unity. So `G_q` is the union of
  `(q+1)/6` rotated triangular tori.
- **It is Ramanujan.** Its eigenvalues are Salié-type sums. We checked
  `|λ| ≤ 2√q` numerically for every `q < 400`.

Hoffman's ratio bound then gives `χ_f(G_q) > 5`, so **`χ(G_q) ≥ 6` for every
`q ≥ 53`**. For `q > 62` this follows from the Ramanujan bound alone; 53 and 59
were computed.

The same bound holds at every level `O/π^r`, because level-`j` eigenvalues are
`q^{r−j}` times level-`j` sums. The Delsarte LP gives nothing more.

| `q` | 2 | 5 | 11 | 17 | 23 | 29 | 41 | 47 | ≥ 53 |
|---|---|---|---|---|---|---|---|---|---|
| `χ(G_q)` | 4 | 4 | 5 | 6? | ≥ 7? | ≥ 6? | ≥ 6? | ? | ≥ 6 (Hoffman) |

How the uncertain entries were reached:
- **`q = 17`:** tabu search finds no 5-colouring in 2·10⁸ moves, and finds a
  6-colouring at once. The best independent sets found have 57 points, against
  `289/5 = 57.8`. SAT proofs are running.
- **`q = 23`:** tabu fails at 5 and at 6 colours.
- **`q = 29`, `q = 41`:** tabu fails at 5.

A marked entry is not yet a proof.

**The field screen.** `scripts/fieldscreen.py` lists the non-split places of a
multiquadratic CM field, using Kummer theory for its local Galois groups.

| field | non-split places, residue `𝔽_q` | verdict for six |
|---|---|---|
| Moser `ℚ(√−3, √−11)` | 2, 11, 17, 29, … | dead (2-adic: `χ = 4`) |
| `ℚ(√−3, √−7, √−15)` | 5, 41, … | dead (`χ ≤ 4`) |
| `ℚ(√−3, √−11, √−15)`, `ℚ(√−3, √−11, √−23)` | 11, … | dead (`χ ≤ 5`) |
| `ℚ(√−3, √−11, √−247)` (`five_rho7`) | 11, 29, … | dead (`χ = 5`) |
| **`ℚ(√−3, √−7, √−11)`** | **17, 41, 83, 101, …** | **open**, if `G₁₇` and `G₄₁` have `χ ≥ 6` |
| de Grey's `ℚ(√−3, √−7, √−11, √−15)` | 41, 101, 131, … | open, if `G₄₁` has `χ ≥ 6` |
| `ℚ(√−3, √−7, √2717)` | 59, 83, 89, … | open (Hoffman at every place) |
| `ℚ(√−3, √−7, √−11, √−247)` | 83, 173, … | open (Hoffman at every place) |

**A 5-chromatic graph in the new universe.** Take the carrier's forced pair at
`d² = 64/9`, with the composition rule of *Forced pairs compose*. It moves to
`d² = 16` with the rotation `cos φ = 1/8`, `sin φ = 3√7/8`, whose radical is
`√(16(256 − 144)) = 16√7`. De Grey's `4e` rotation `(31 + 3√−7)/32` then
spindles it.
- It is built by `TARGETS=16 BASE=3,7,11 scripts/tune.py`.
- Everything stays in `ℚ(√3, √7, √11)`, with complex coordinates in
  `ℚ(√−3, √−7, √−11)`, and needs no `√5` and no `√247`.
- `data/five_tuned_16_1_3_7_11.json` has **4 081 points and 27 242 edges**. It
  is not 4-colourable: CaDiCaL (284 s) and kissat independently.

**The first search there.**
- **Its `λ`-closure kills the cheap colourings.** Exoo–Ismailescu's rotation lies
  in `ℚ(√−11)`, so the closure stays inside the field. It has 378 unit vectors
  and rank 8. Every homomorphism to `ℤ/5`, `ℤ/7`, `ℤ/10`, `ℤ/12`, `ℤ/14` or
  `ℤ/15` that is nonzero on all units is gone. Sampling `ℤ/9`, `ℤ/11` and
  `ℤ/13` found none that is 5-colourable.
- **The circular gate is undecided.** SCIP finds no arc character in 30 minutes.
- **Growth is running.** It starts from the graph together with its `λ`-image,
  8 161 points (`data/F8_lam_kw1.json`).

**Two searches now run there.**
- **A plain search for a non-5-colourable graph.** The first `λ`-closure put
  5-denominators at only one of the two places above 5. At the other place
  every unit vector was integral, which left room for a 5-adic periodic
  colouring. Tabu search did not find one at level `5³` in 50 minutes, but
  the module is now closed under `λ̄` as well: 666 unit vectors, with both
  places above 5 non-integral. The run starts again from 12 241 points. The
  first run reached 27 614 points in the old module. There each hard step took
  kissat 8–30 minutes, and the translation agreement of its colourings was
  0.45–0.53, against 0.87 in the dead field.
- **A search for a distance-2 gadget.** Exoo and Ismailescu's theorem turns
  any unit-distance graph in which two points at distance 2 always get
  different colours into `χ ≥ 6`. The 2-edges of their `{1,2}`-graph lie in
  the Moser field, and so does every rotation needed to put a copy of the
  gadget on each of them, so the whole construction stays in
  `ℚ(√−3, √−7, √−11)`.

  The search imposes `c(A) = c(B)` for a pair at distance 2
  (`MODE=same grow_kw.py`). An UNSAT answer would be that gadget, and it would
  go through `verify6`-style checks before any claim.

Single-place residue gates on the 378-unit module show no local colouring:

| `p` | residues of the unit vectors | Cayley graph |
|---|---|---|
| 7 | 48, all of `𝔽₄₉*` | not 5-colourable |
| 11 | 54 | not 5-colourable |
| 13 | 84 | not 5-colourable |
| 19 | 120 | not 5-colourable |
| 23 | 132 | not 5-colourable |

These come from `scripts/ramgate.py`, which also handles ramified primes.

## Blind edges, the level-2 plane at 17, and fields with no local obstruction

**The growth scripts were half-blind.** `grow_kw.py` and its predecessors
looked for edges only along their own unit set `U`. The exact graph on the same
points can have more: unit-distance pairs along directions outside `U`, mostly
between rotated copies of the seed.
- In the `F8` checkpoint of 19 851 points (666 units) the exact graph has
  143 778 edges. The growth saw 136 954.
- The 6 824 unseen edges lie along 26 directions, and the saved colouring made
  **6 788 of them monochromatic**. The colourings were exploiting edges the
  growth could not see. This is where the high "translation agreement" came
  from.
- The exact graph is still 5-colourable: tabu repairs the colouring in 20 s.
  Kissat's model of the 21 627-point hard instance leaves 2 exact edges
  monochromatic, so that verdict needs a repair as well.

`scripts/grow_lean.py` now finds every unit-distance pair:
1. A KD-tree finds all pairs at float distance 1.
2. One pair per unknown direction is checked in exact arithmetic.
3. The direction joins `U` with its negative and conjugates, together with
   every edge along it.

On the `F8` λ-seed this recovers exactly the 6 811 missing edges, along 64
directions. The same script keeps its candidates in sorted numpy arrays: 0.66 GB
at 20 000 points, against 10.8 GB for the dict-based pool of `grow_kw.py` at
32 000 points.

**The level-2 plane at 17.** Every level gives an upper bound:
`χ(Γ(F8)) ≤ χ(Cay(O/17^r, T_r))` for each `r`. We tested whether level 2 could
be easier than level 1.
- **The graph.** Level 2 has 83 521 vertices of degree 306.
- **Tabu from the lift.** It starts from the best 5-colouring found for `G₁₇`
  (34 conflicts), which becomes 167 042 conflicts at level 2. In 2.4·10⁷
  moves it never improves on the lift.
- **No embedding of level 1.** A lift `σ(a) = τ(a) + q·s(a)`, with `τ` the
  Teichmüller lift, would embed `G_q` in level 2 and make the two levels
  equivalent. The conditions on `s` form a linear system over `𝔽_q`, and it
  is inconsistent for `q = 5, 11, 17` (rank 47 of 50, 239 of 242, 575 of 578).

So level 2 is not trivially equivalent to level 1. Still, no sign of a
5-colouring there either (`scratchpad` scripts `cay_tabu.c`, `lift_section.py`).

**Fields with no local obstruction at all.** Extending `F8` by one more
`√−d` removes the non-split places 17 and 41 for many `d`. With `d` squarefree
and below 400, `F8(√−d)` has no non-split place of norm below 53 (and none
ramified) for:

> 1, 2, 42, 43, 59, 66, 83, 86, 87, 103, 115, 118, 127, 154, 155, 166, 174,
> 185, 195, 203, 206, 213, 223, 230, 237, **247**, 251, 254, …

In these fields Proposition B leaves no local 5-colouring at any place or
level.

`F8(i)` needs no new coordinate field, because `i` is just the vector `(0, 1)`.
But `U ∪ iU` alone gives the Cartesian product `Γ(F8) □ Γ(F8)`, of chromatic
number 5. Only unit vectors outside `F8 ∪ iF8`, such as `(3+4i)/5`, mix the
two layers, and even they never close a triangle across them.

**`L16 = ℚ(√−3, √−7, √−11, √−247)` holds both families.** `five_tuned_16`
and `five_rho7` share their 402-point Moser-field carrier.
- **The union.** It has 6 080 points and 37 474 edges. CaDiCaL needs 107 s to
  5-colour it, far more than either graph alone.
- **Its units.** There are 918: the `F8` λ-closure, the 134 of the 247-module,
  and every edge direction of the seed.
- **Growth there.** Tabu first failed at 12 284 points, where kissat needed
  about 20 minutes (SAT). In `F8`, tabu first failed at about 22 000 points.

## Two-step distances, and the repulsive distance `2/√3`

**Where Exoo–Ismailescu's route stalls.** A unit-distance graph `H` with a
pair at distance `d` split by every 5-colouring, together with a finite
`{1, d}`-graph `W` that is not 5-colourable, gives `χ(ℝ²) ≥ 6`: put a copy of
`H` on every `d`-edge of `W`. `W` is known for four values of `d`:

| `d` | who |
|---|---|
| `(1+√5)/2` | Huddleston; 31 vertices (Parts, [arXiv 2010.12656](https://arxiv.org/abs/2010.12656)) |
| `2` | Exoo–Ismailescu ([arXiv 1909.13177](https://arxiv.org/abs/1909.13177)) |
| `√3`, `(√3+1)/√2` | Ágoston–Pálvölgyi |

**All four are two-step distances** `|u + v|`, with `u, v` unit vectors, at
turning angles 72°, 0°, 60° and 30°. Two points at such a distance share a
unit neighbour, so 5-colourings colour them alike more often than chance. The
gadget `H` asks for the opposite, and that is why it has never been found. In
the L16 seed (6 080 points, 10 random 5-colourings):

| `d²` | `d` | pairs | `P(same)` |
|---|---|---|---|
| **4/3** | **1.1547** | 17 138 | **0.079** |
| **16/3** | 2.3094 | 2 667 | **0.085** |
| **7** | 2.6458 | 4 398 | **0.090** |
| **28/3** | 3.0551 | 1 345 | **0.098** |
| 52/3, 19 | 4.16, 4.36 | 172, 146 | 0.121 |
| 9 | 3 | 1 490 | 0.141 |
| 19/3 | 2.5166 | 3 579 | 0.187 |
| 7/3 | 1.5275 | 21 580 | 0.241 |
| 3 | `√3` | 20 837 | 0.292 |
| 4 | 2 | 8 850 | 0.339 |

`2/√3` is the most repulsive distance. At `2/√3`, 5 202 pairs are split, with
both ends in one Kempe component, in 12 of 12 tabu colourings; at distance 2,
only 347 are. In L16, `2/√3` is **not** a two-step distance: `|1 + e^{iθ}|² = 4/3`
needs `cos θ = −1/3`, hence `√−2`, which L16 lacks. So these pairs have no common
neighbour.

**The witness side for `2/√3`.** None is known.
- The lattice `ℤ[ω]/√−3` with `{1, 2/√3}` forbidden is 5-colourable,
  periodically mod 6.
- The L16 seed's `{1, 2/√3}`-graph (54 612 edges) is 5-colourable (kissat,
  23 s).
- Minkowski sums in `ℚ(√−2, √−3)`, where the triangle `(1, 1, 2/√3)` exists,
  are 5-colourable up to 931 points.

**Several forbidden distances.** With more than one extra distance, the
lattice becomes a small witness:
- with `{1, 2/√3, √(19/3)}` forbidden, a 127-point patch is not 5-colourable;
- with `{1, 2/√3, √(7/3)}` forbidden, a 61-point patch is not 5-colourable.

A `W` with several distances needs one gadget per distance. So the search now
targets gadgets at the repulsive distances. It runs in four cloud workers,
following `notes/worker_jobs.md` and checked by `scripts/verify_pair.py`:
- two workers grow gadgets for eight `2/√3` pairs;
- one searches for a `{1, 2/√3}` witness;
- one maps the repulsion spectrum across the project's graphs.

**The witness side with several distances (24 September, evening).** Put the
repulsive distances together.
- **What the lattice alone gives.** `ℤ[ω]/√−3` forbids the distances
  `√(n/3)`. No single extra distance makes it 6-chromatic. With two or three
  extra distances it is often not 5-colourable (`scripts/lattice_witness.py`;
  16 such sets use only distances with `P(same) ≤ 0.141` in the L16 seed).
- **A verified witness.** `data/W_lattice_16_21_28_61.json` is a vertex-critical
  72-point lattice graph, with edges at `1, 4/√3, √7, √(28/3), √(61/3)`, that is
  not 5-colourable. CaDiCaL, Glucose, MiniSat and kissat agree, and drat-trim
  verifies kissat's DRAT proof.
- **What it would finish.** One unit-distance gadget per distance would give
  `χ(ℝ²) ≥ 6`, and the four distances are repulsive (0.085, 0.090, 0.098, 0.133).

**A cascade instead of four gadgets.** Once a unit-distance gadget forces
`2/√3` apart, every 5-colouring of the plane avoids `{1, 2/√3}`. A
`{1, 2/√3}`-graph can then act as the gadget for the next distance. It is far
denser, so forcing is easier. For example, the lattice with
`{1, 2/√3, √(7/3)}` forbidden is not 5-colourable on 61 points. So one
unit-distance gadget plus one two-distance gadget at `√(7/3)` would do.

**The L16 module has no cheap colourings.**
- **No local ones.** There are no non-split places below 53.
- **No circular ones, numerically.** Its units give 471 directions in rank 12.
  - The largest `min_u ‖φ(u)‖` found over characters `φ` is 0.024, against the
    0.2 a circular 5-colouring needs (4 000 restarts of coordinate ascent).
  - The same search finds 0.204 at once for the Moser module, which is
    4-colourable.
  - The exact MILP ran out of time, so this is evidence, not proof.

## An open question closed: the plane over `ℚ(√3, √11)` is 4-chromatic (25 September)

This is a by-product of the arithmetic built for six, and not a step towards
six. It settles a question that was open in print.

**Theorem.** `χ(ℚ(√3, √11)²) = 4`. So no 5-chromatic unit-distance graph has
all its coordinates in `ℚ(√3, √11)`.

**The question, as it stood.**
- `ℚ(√3, √11)` is the smallest field whose plane contains a Moser spindle
  (Moorhouse 2010, Prop. 1.4). So `χ ≥ 4`, and Moorhouse wrote: "We have not
  determined the exact value of `χ(K²)` in this case."
- Madore ([arXiv 1509.07023](https://arxiv.org/abs/1509.07023), Prop. 4.6)
  proved `4 ≤ χ ≤ 5` by reducing at a place over 11.
- Exoo–Ismailescu ([arXiv 1805.00157](https://arxiv.org/abs/1805.00157), DCG
  2020): "One interesting question is whether there exists a 5-chromatic unit
  distance graph which can be embedded in `ℚ[√3, √11] × ℚ[√3, √11]`."
- Polymath16, [thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/):
  "it would be nice to decide if the chromatic number of this subfield is 5."
- Voronov, in Polymath16 [thread 17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/)
  (July 2021): "it seems likely that `χ(Q(i, √3, √11)) = 4` ... But as far as
  I know, nobody has proved this yet."

**The proof.** Identify the plane with `K = ℚ(i, √3, √11)`. A unit vector is
then a `u ∈ K` with `u ū = 1`.
1. **The completion over 2.** `√33` is 2-adic, since `33 ≡ 1 (mod 8)`. So
   `ℚ(√3, √11)` has two places over 2, each with completion `ℚ₂(√3)`.
2. **The place is inert.** Neither `−1` nor `−3` is a square in `ℚ₂`, so `i`
   is not in `ℚ₂(√3)`. Adjoining it gives `ℚ₂(√3)(√−3)`, which is unramified.
3. **Unit vectors reduce to nonzero residues.** The local ring is
   `ℤ₂[√3][ω]`, with residue field `𝔽₄`. Every unit vector is a unit of it,
   with a nonzero residue.
4. **The colouring.** Colour `z` by the residue of `z` minus its 2-adic
   fractional part. A unit step adds a nonzero element of `𝔽₄`, so the
   colouring is proper, with four colours.

This is `hn.adelic.q311_colour`.

**Checks.** `tests/test_q311.py` verifies:
- 810 unit vectors, including `(3 + 4i)/5` and `(√33 + 4i)/7`, all reduce to
  nonzero residues at both places;
- the colouring is proper on random unit steps, on the Moser spindle, and on
  Exoo–Ismailescu's 214-point graph;
- on 3 012 further edges (`data/ei_rho7.json`) it has no monochromatic edge.

**Why Moorhouse and Madore stopped short.** They reduce the coordinates
`(x, y)`, which needs `x² + y²` to be anisotropic modulo `𝔪` or `𝔪²`.
- Over 2 with `√3` present that fails, since `1 + 1 ≡ 0 (mod 𝔪²)`.
- A unit vector like `(−1/2, √3/2)` has non-integral coordinates, yet it is
  integral in `ℤ₂[√3][ω]`.

Reducing `z = x + iy` instead — the Hermitian form, `notes/local_colourings.md`
§2 — sees the `ω` that the coordinates hide.

**Who came closest.** The 2-adic idea is not ours.
- In thread 3 (May 2018) David Speyer 4-coloured the *Moser ring*, the
  elements of `ℚ(√−3, √−11)` integral over `ℤ[1/3]`. He used exactly this
  reduction, with colours in `ℤ[ω]/2 = 𝔽₄`.
- Gibbs and Hubai found that all such colourings have period 8.
- Dúcz ([arXiv 2606.12325](https://arxiv.org/abs/2606.12325)) 4-coloured the
  Moser lattice and ring in 2026.

What is ours is the step from the ring to the whole plane. The place over 2
is inert in `ℚ(i, √3, √11)`, so every unit vector of the plane is a 2-adic
unit, not only those of the ring. Denominators of 2 are handled by cosets.

**A paradox from Polymath16, explained.**
- In [thread 13](https://dustingmixon.wordpress.com/2019/07/08/polymath16-thirteenth-thread-bumping-the-deadline/)
  Parts chained Exoo–Ismailescu's forced alike pairs at distance 8/3 into
  alike pairs at every distance `8/9ⁿ`, whose sum is 1. He offered this as a
  "funny proof" that this field needs five colours, and asked why the colour
  is lost in the limit.
- In the colouring above, every pair at distance `8/9ⁿ` is indeed alike,
  since `(8/9ⁿ)u` is divisible by 8 in `ℤ₂[√3][ω]`. The colouring is still
  proper.
- It is continuous for the 2-adic topology, not for the real one, which is
  where the limit was taken.

**By-products** (`notes/local_colourings.md` §8):
- `χ(ℚ(√d)²) ≤ 4` for `d ≡ 3 (mod 8)`. This extends Moorhouse's Theorem 8.1,
  which excluded `d ≡ 47, 59, 83 (mod 84)`, to `d ≡ 59, 83, 131 (mod 168)`.
- A real field `F` with `χ(F²) ≥ 5` needs `−1` to be a square in its
  completions over 2 (for small residue fields). It also needs no place with
  residue field `𝔽₃` or `𝔽₇`.

We found no proof of this theorem in any paper, Polymath16 thread or web
search; the closest is Speyer's colouring of the Moser ring. It has not been
refereed.

## Where this work stands in the literature

`notes/literature.md` compares the project with the published record, with
links. In short:
- **The route to six is known.** It is Exoo–Ismailescu's: a two-distance
  witness plus a unit-distance gadget, which Polymath16 calls "virtual edges".
- **Known witnesses.** Witnesses exist for `φ`, `2`, `√3` and
  `(√6 + √2)/2`. Colouring-guided growth at four colours is Heule's and
  Parts'.
- **Known reduction.** Reduction of field planes modulo a prime is Moorhouse's
  and Madore's.
- **Known 2-adic colourings.** Speyer 4-coloured the Moser ring by the
  Hermitian reduction over 2 (Polymath16, thread 3).
- **Not found anywhere else:**
  - the theorem above, which carries that reduction to the whole plane;
  - the CM whole-field results, with their thresholds for six;
  - the repulsion spectra, and the two-step explanation of why every known
    witness distance is attractive;
  - `2/√3` as the gadget distance;
  - the 72-point multi-distance witness;
  - skeleton and `MODE` growth.

## Galois conjugation carries gadgets: a witness that needs one gadget (25 September)

**The observation.** Let `K` be a CM field, as all of ours are, and `σ` any
automorphism of `K`.
- **A graph automorphism.** Complex conjugation commutes with `σ`. So if
  `u ū = 1` then `σ(u) · conj σ(u) = σ(u ū) = 1`: `σ` maps unit vectors to unit
  vectors. It is therefore an automorphism of the unit-distance graph `Γ(K)`,
  although it moves points wildly in the real plane.
- **It carries gadgets.** Let `H ⊂ K` be a unit-distance graph in which a pair
  at squared distance `d²` is split in every 5-colouring. Then `σ(H)` does the
  same for `σ(d²)`, which is again positive.

**What it changes.** In the Exoo–Ismailescu route a witness `W` needs one gadget
per non-unit distance. With Galois it needs one per Galois orbit of distances.
We found nothing like this in the literature or in Polymath16.

**A verified witness of the new kind.** `data/W_moser_orbit_9_33.json`:
- **Size.** 187 points in `ℚ(√−3, √−11)`, vertex-critical.
- **Edges.** 508 at distance 1, and 495 at the orbit `d² = (9 ∓ √33)/6`
  (`d = 0.7366` and `1.5676`).
- **Verdict.** Not 5-colourable. CaDiCaL, Glucose, MiniSat and kissat agree,
  and drat-trim verifies the proof.
- **Cost.** It needs one gadget, at 0.7366 or equivalently at 1.5676.

**Repulsion is Galois-invariant, as it must be.** Across 16–24 tabu colourings
of the L16 seed, conjugate distances have almost the same `P(same)`:

| orbit | `P(same)` |
|---|---|
| `(9 ∓ √33)/6` | 0.205 / 0.212 |
| `(14 ∓ 2√33)/3` | 0.079 / 0.086 |
| `(7 ∓ √33)/2` | 0.133 / 0.132 |
| `(31/6 ∓ 5√33/18)` | 0.135 / 0.139 |

So an orbit is repulsive or attractive as a whole.

**The tension that remains.**
- **The witness side.** A witness needs distances that constrain colourings,
  that is, attractive ones. The only orbit that makes a 865-point Moser ball
  refuse five colours is the two-step orbit `(9 ∓ √33)/6`.
- **The gadget side.** A gadget is easiest at a repulsive distance.

**Four cloud searches follow from this.**
- **`orbitw`.** Does the 18 524-point L16 growth graph (`data/L16_kw2.json`)
  refuse five colours once the edges of the repulsive orbit `(14 ∓ 2√33)/3`
  (`d = 0.9149`, `2.9149`) are added?
  - tabu fails with 366 conflicts, against 108 without them;
  - kissat is running for 12 hours.
- **`w0915`.** Grow a witness directly as a `{1, 0.9149, 2.9149}`-graph.
  `scripts/grow_lean.py` now accepts exact extra distances, via the JSON keys
  `dist2_exact` and `units2`.
- **`g0915`.** Grow a unit-distance gadget at 0.9149 or 2.9149, on four
  target pairs.
- **`g0737`.** Grow a unit-distance gadget at 0.7366 or 1.5676, which would
  complete the verified 187-point witness.


## Voronov's second case: the plane over `ℚ(√2, √3)` is 4-chromatic (25 September)

In Polymath16 thread 17, Voronov named two cases: `χ(ℚ(i, √3, √11)) = 4` and
"case (2, 3)". The 2-adic argument settles the second as well.

**Upper bound.** `ℚ(√2, √3)` has a single place over 2.
- It is totally ramified, with residue field `𝔽₂`.
- `i` is not in the completion, so the place does not split in `ℚ(i, √2, √3)`.
- The extension contains `√−3`, so it is unramified, with residue field `𝔽₄`.

So every unit vector reduces to a nonzero element of `𝔽₄`.

**Lower bound.** A chain of three unit rhombi along `(1, 0)`, `(0, 1)` and
`(−2/3 − √2/6, −2/3 + √2/6)`. Their sum has squared length exactly 1/3, so the
chain closes at distance 1. The graph has 10 vertices and 16 edges.
- Three solvers find no 3-colouring.
- drat-trim verifies kissat's proof.

See `notes/local_colourings.md` §10, `tests/test_q23.py` and
`certificates/chain23_no3coloring.json`.


## Correction: `χ(ℚ(√3, √11)²) = 4` is Fischer's theorem (1994) (25 September)

A literature search for the quadratic field `ℚ(√47)` turned up two papers of
K. G. Fischer that none of Moorhouse, Madore, Exoo–Ismailescu or Polymath16
cites:
- *A planar geometric graph of chromatic number four*, Congr. Numer. 104
  (1994) 73–79 (Zbl 0836.05030). It proves that `ℚ(√p, √q)²` has an additive
  4-colouring, into `ℤ/4`, for squarefree coprime `p ≡ 3`, `q ≡ 11 (mod 16)` with
  `pq ≡ 1 (mod 32)`, and concludes `χ(ℚ(√3, √11)²) = 4`.
- *Additive K-colorable extensions of the rational plane*, Discrete Math. 82
  (1990) 181–195. As summarised by Payne (arXiv 0707.1177), it proves
  `χ(ℚ(√N)²) ≤ 4` for `N ≡ 3 (mod 8)`.

So the entry "An open question closed" above claims too much: the question was
open in the later literature, but Fischer had answered it. What remains ours,
as far as we know, is `χ(ℚ(√2, √3)²) = 4` (Fischer's hypotheses exclude it) and
the short criterion behind both proofs. We have read only the zbMATH summaries
of Fischer's papers. The note, the READMEs and `notes/literature.md` now credit
Fischer, and the two mathematicians who received the first version of the note
are being told.

## The plane over `ℚ(√47)`: first experiments (25 September)

Moorhouse's classes `d ≡ 47, 143, 167 (mod 168)` are the real quadratic fields
where no bound below 5 is known; the smallest is `ℚ(√47)`, with
`3 ≤ χ ≤ 5` (Fischer 1990 or Moorhouse Thm 8.6 below, reduction at 11 above).
At 2 the place splits in `ℚ(i, √47)`, since `ℚ₂(√47) = ℚ₂(i)`, so no 2-adic
colouring exists.
- **Odd cycles.** Unit vectors of height at most 244 give an odd cycle of
  length 13 (CP-SAT; lower bound 11). With only the axes and the
  `ℚ(√−47)` vectors `(±23 ± i√47)/24`, the shortest odd cycle is longer than 37.
- **Balls are 3-colourable.** With the 48 unit vectors of height at most 40,
  the ball of radius 5 in the Cayley graph (2 248 121 points, 10 439 472 edges)
  is 3-colourable; with the 28 vectors of heights 1, 5 and 40, so is the ball of
  radius 6 (751 073 points).
- **No periodic 3-colouring.** For the 48-vector module `M ≅ ℤ⁴`, no coset
  3-colouring exists, and `Cay(M/pM, U)` is not 3-colourable for
  `p = 8, 9, 11, 16, 17, 18, 19, 22, 23, 24, 25, 27`; for the other `p ≤ 27`
  some unit vector lies in `pM`. So a 3-colouring of the module, if one exists,
  has no period lattice of exponent at most 27.
- **Cores wrap around.** The smallest non-3-colourable core of `Cay(M/8M, U)`
  has 7 vertices, but 4 of its 12 edges do not lift to `M`.

The prototype scripts for these experiments are not yet in the repository.
Whether `χ(ℚ(√47)²)` is 3 or 4 (or 5) remains open.

## Source check of the note (25 September)

Before the correction was sent, every statement in the note about other work
was checked against the source itself: the zbMATH entries for Fischer (1990,
1994) and Johnson (1987); the papers of Payne, Moorhouse, Madore,
Exoo–Ismailescu, Heule, Dúcz and Voronov–Neopryatnaya–Dergachev; and the
Polymath16 comments of Speyer, Gibbs and Voronov. The prime decompositions
were recomputed with PARI/GP, and the tests of both theorems pass. Changes:
- Fischer (1994) *proves the existence* of an additive 4-colouring; the note
  said "constructed".
- Madore's ¶6.6 already states his Prop. 3.2 for any quadratic form, so
  Theorem 2 is a case of it. The note now says so, and describes its own
  contribution as the choice of coordinates.
- Speyer coloured the Moser ring by projecting `R/2R ≅ 𝔽₄ × 𝔽₄`; the
  coordinates `α + βω` are ours. Gibbs reported Hubai's "analysis and computer
  search". Voronov is now quoted word for word.
- Two related works are now cited: Fischer on the connected components of
  `ℚ(√N₁, …, √N_d)²` (Congr. Numer. 72, 1990) and Johnson's status report
  (Geombinatorics 9, 2000). We have not seen the report, nor the full texts
  of Fischer's papers; the abstract of Fischer (1990) is not available online.

## Fischer (1990) read in full (25 September)

The full text of *Additive K-colorable extensions of the rational plane*
(Discrete Math. 82 (1990) 181–195), which is in Elsevier's open archive, was
read. It treats `ℚᵈ` and the quadratic fields `ℚ(√N)`; it says nothing about
`ℚ(√2, √3)` or other biquadratic fields.
- **Thm 1.** The components of `Fᵈ` are the translates of the group `C₀`
  generated by the unit vectors, so it suffices to colour `C₀`. This is the
  coset step, before Moorhouse and Madore. Note v4 credits it.
- **Thm 8.** `ℚ(√N)²` is 2-colourable exactly when `N ≢ 3 (mod 4)`. Johnson
  proved the "if" part independently.
- **Thms 9, 10(i).** An additive 3-colouring for `N ≢ 2 (mod 3)`, and an
  additive 4-colouring for `N ≡ 3 (mod 8)`: a linear functional into `ℤ₍₂₎`
  followed by reduction modulo 4, so already a reduction at 2. Note v4 says so.
- **Thm 10(ii).** For `N ≡ −1 (mod 24)` there is no additive `k`-colouring
  with `k ≤ 6`: with `a = (N + 1)/2`, `b = (N − 1)/2`, the unit vector
  `(b + i√N)/a` gives `1/2, 1/3 ∈ C₀`. This covers `N = 47`, so any
  3-colouring of `ℚ(√47)²` is non-additive, which fits the experiments above.
  For `N = 167` (Example 4) there is none with `k < 11`.

## Level 2 at 11; two square roots; the plane over `ℚ(√3, √5)` (25 September)

- **Level 2 at 11.** Reduction modulo 11 gives `χ(ℚ(√47)²) ≤ χ(G₁₁) = 5`,
  because 11 splits in `ℚ(√47)`. Reduction modulo 121 gives the bound
  `χ(Cay((ℤ/121)², U₂))`, where `U₂` is the 132 unit vectors modulo 121. That
  graph (14 641 vertices, 966 306 edges) has no proper 4-colouring: kissat,
  CaDiCaL and Glucose all return UNSAT with a unit triangle pinned. Level 1
  checks as 5- but not 4-colourable. So the second level gives nothing below 5.
- **Voronov's question** (Polymath16, thread 17, comment 29283): can a field
  generated by two square roots of primes carry a 5-chromatic graph? For
  `ℚ(√3, √q)` the answer is no unless `q ≡ 5, 23 (mod 24)`: reduction at 3
  gives `χ = 3` when `q ≡ 1 (mod 3)`, and Theorem 2 gives `χ ≤ 4` when `q = 2`
  or `q ≡ 1, 3 (mod 8)` (`notes/local_colourings.md` §11).
- **`ℚ(√3, √5)`, the smallest open case.** `τ = (2 + i√5)/3` makes `1/3` a sum
  of unit vectors, and a chain of six unit rhombi (`data/chain35.json`,
  19 vertices) has no proper 3-colouring. Reduction at 11 gives `χ ≤ 5`. So
  `4 ≤ χ(ℚ(√3, √5)²) ≤ 5`.
- **False evidence first.** Before the chain, Minkowski balls and growths at
  three and four colours up to 83 000 points were all 3-colourable. They stayed
  near the origin; the chain reaches radius about 3.5.
- A four-colour growth from the twelve turns of the chain about the origin is
  running.

## The planes over `ℚ(√3, √q)` (25 September)

The lemma behind the `ℚ(√3, √5)` chain works for every squarefree `q ≡ 2 (mod 3)`.
There 3 splits in `ℚ(√−q) ⊂ ℚ(√3, √q)(i)`. A generator `α` of the `h`-th
power of a prime above 3 has norm `3^h`, and `Tr(α²)` is prime to 3, so
`α/ᾱ + ᾱ/α` has exact denominator `3^h`. Hence `1/3` is a sum of unit vectors,
and the rhombus chain shows `χ ≥ 4`. With reduction at 3 and Theorem 2:
- `χ(ℚ(√3, √q)²) = 3` for `q ≡ 1 (mod 3)`;
- `= 4` for `q` even or `q ≡ 1, 3 (mod 8)`: for primes, `q = 2` or
  `q ≡ 11, 17 (mod 24)`;
- `≥ 4` for `q ≡ 5, 23 (mod 24)`, the only open class.

`tests/test_q3q.py` checks the generators for `q ≤ 113`, and for `q = 17` a
chain of 90 rhombi built from `(8 + i√17)/9`, which has no proper 3-colouring.
For `q = 17` the value 4 is new as far as we know; Fischer's family covers
`q ≡ 11 (mod 32)`, and states the lower bound only for `q = 11`.

The search for a 5-chromatic graph over `ℚ(√3, √5)` goes on. The growths at
four colours stay easy. Their colourings are not periodic at the places over
7, 13 or 17, and no degree-1 place up to 61 four-colours the module of the
chain seed.

## More square roots (25 September)

The dichotomy of Theorem 5 (of `notes/local_colourings.md` §11) holds for every multiquadratic `L ∋ √3`: `χ(L²) = 3`
when every other generator, with any factor 3 removed, is `≡ 1 (mod 3)`, and
`χ(L²) ≥ 4` as soon as one is `≡ 2 (mod 3)` (the subfield `ℚ(√−q)` of `L(i)`
supplies `1/3`). PARI/GP gives residue degree 1 above 3 for `ℚ(√3, √7)`,
`ℚ(√3, √7, √13)` and `ℚ(√3, √7, √13, √19)`, and 2 for six fields with a
generator `≡ 2 (mod 3)`. `notes/local_colourings.md` §11.

## The prime 2 was the wall (25 September)

The growths over `ℚ(√3, √5)` stayed 4-colourable up to 112 000 points, and the
module gate found no periodic colouring at the odd primes it checked. The gate
never looks at 2, and 2 is the answer. It splits in `ℚ(i, √3, √5)` with residue field `𝔽₄`.
The unit vectors `ζ^a τ^b` are all units there. Reduction `z ↦ z mod w₂ ∈ 𝔽₄`
sends a unit step `u` to `ρ(u) ≠ 0`, so it is a proper 4-colouring of every
graph those units can build. The two runs targeting pairs at `d² = 7/3` and
`d = 2` were stopped: they could not succeed.

In general (`notes/local_colourings.md` §12, Proposition C), a split place with
residue field `𝔽_q` colours the graph of the edges that are units there with
`χ(H_q)` colours. Here `χ(H₄) = 4` and `χ(H₉) = 3` (`tests/test_split_places.py`).
So over `ℚ(√3, √5)` a 4-chromatic graph needs an edge vector that is not a unit
above 3, and a 5-chromatic one also needs one that is not a unit above 2.
The same holds in de Grey's, Voronov–Neopryatnaya–Dergachev's and
Exoo–Ismailescu's fields, where 2 splits with residue field `𝔽₄`. That is why
their spindle rotations carry a power of 2 in the denominator.

Other findings of the day:
- `1/√3` is a sum of four unit vectors of `ℚ(i, √3, √5)`, which gives a
  13-vertex chain of four rhombi with no proper 3-colouring
  (`tests/test_q35.py`). Three would mean a point on a genus-1 curve; none has
  small height.
- A split-place gate with deep units (`scripts/experiments/split_gate_q35.py`)
  scales by a power of the uniformiser and reduces at `w^A w̄^B`. Results:
  - `ζ, τ, (1 + i√15)/4`: 4-colourable at 2 at the second level.
  - `ζ, o₁, (1 + i√15)/4`: no 4-colouring at 2 up to 1 048 576 points (kissat).
    Here `o₁` has odd valuation above 2.
  - With `τ` added (252 units), none at 3 (59 049 points) or at 5 (up to
    390 625 points) either.
- The Voronov–Neopryatnaya–Dergachev construction `M₃ ∪ ψM₃` does not transfer
  directly. From `ζ, τ` the set `M₃` is 3-colourable. With
  `(1 + i√15)/4` added, `M₂` (5 833 points) is 4-chromatic, containing a
  19-vertex 4-chromatic subgraph of radius 2. But `M₂ ∪ (7 + i√15)/8 · M₂` is
  still 4-colourable.

A colouring-guided growth at four colours inside the 252-unit module is running.
Unlike every earlier run, its colourings are tight: about a thousand rainbow
candidates at each step, and tabu gives way to kissat from 13 000 points on.

## Two square roots without `√3` (25 September)

Proposition A, run over every pair of primes `a < b < 60`, settles most of
Voronov's second question locally (`notes/local_colourings.md` §13). The input
is the non-split places of `K = ℚ(i, √a, √b)` below 200.
- **A ramified place above 2** makes the plane bipartite. The residues `±1`
  coincide in characteristic 2.
- **An unramified place with residue field `𝔽_q`** gives `χ(G_q)`.

Of the 120 fields without `√3`:
- 28 are bipartite, 14 are at most 3-chromatic and 33 at most 4-chromatic;
- 28 are bounded by 5, through `G₁₁` or through `G₁₉`;
- 17 have no non-split place with `q ≤ 19`.

New finite-plane values:
- **`χ(G₁₉) = 5`.** The plane has no triangle and a linear 5-colouring. So
  `ℚ(√5, √7)`, bounded only above 19, is locally consistent with a
  triangle-free 5-chromatic graph.
- **`χ(G₄) = 4`.**

With `√3`, the table reproduces Theorem 5 (§11). The one exception is `ℚ(√3, √29)`,
whose first non-split place is above 23.

Also recorded today:
- The four-rhombus units are square roots of `τ` up to roots of unity:
  `c₁ c₂ = −τ`, `c₂ = ζ⁹ c₁`, and `c₃ = conj(ζτ)`. They have odd valuation
  above 3, as `τ` has even valuation there.
- A pool of 2 954 small unit vectors was sorted by depth above 2 and 3.
  - The 15 204 units of depth at most 1 above 2 have no periodic 4-colouring at
    the levels `(2, 2)`, `(3, 3)` or `(4, 4)` above 2. At `(3, 3)` they meet
    every residue class. The 36 units `ζ^a o₁^{±1}` alone are 4-colourable at
    `(2, 2)`, so richness matters.
  - The 15 300 units of depth at most 2 above 2 and at most 1 above 3 have no
    periodic 4-colouring at 2 up to `(4, 4)`, at 3 at `(2, 2)`, or at 5 up to
    `(2, 1)`.
- The Voronov–Neopryatnaya–Dergachev shape `M₃ ∪ (7 + i√15)/8 · M₃`, built
  from `ζ`, `c₁` and `τ`, has 42 913 points and 588 new edges. It is still
  4-colourable, but CaDiCaL needs 172 s to find the colouring.
- The four-colour growth in the 252-unit module reached 13 708 points. From
  then on each step needs kissat for up to half an hour.

## Finite planes: Moorhouse's table continued (25 September)

For `q ≡ 1 (mod 4)` the substitution `(a, b) = (x + iy, x − iy)` makes the
hyperbola graph `H_q` of Proposition C the unit-distance graph of `𝔽_q²`, the
graph of Moorhouse's Table 6.1 (`notes/local_colourings.md` §12). This gives:
- **`χ(𝔽₁₃²) = 5`**, settling Moorhouse's "5 or 6". The value `χ(H₁₃) = 5` was
  already certified.
- **`χ(𝔽₁₇²) ≤ 6`**, from Moorhouse's "5, 6 or 7". The colouring is
  `(x, y) ↦ ⌊((3x + 6y) mod 17)/3⌋`: the line `αx + βy = r` meets the unit
  circle exactly when `α² + β² − r²` is a square, and `45 − r²` is not a
  square for `r = 0, 1, 2`.

The same *interval colourings*, with `m` consecutive parallel lines per colour,
extend Vinh's pairs of lines. They are optimal for `q = 7, 13, 19`, and give
the upper bounds of the continued table.

| `q` | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
|---|---|---|---|---|---|---|---|
| `χ(𝔽_q²)` | 5 | 5–8 | 5–6 | 5–8 | 5–8 | 5–7 | 5–8 |

The lower bounds come from SAT: no 4-colouring for `q = 23, 29, 31, 37, 41, 43`.

**Five colours for `𝔽₁₇²`.** Still open. What was tried:
- **Local search.** Tabu search and a hybrid evolutionary search stop at exactly
  two monochromatic edges, and never reach one in 5·10⁷ moves. In all 28 such
  colourings the two edges touch a common isotropic line `x ± 4y = c`, which
  two random edges do with probability about 0.4.
- **Plain CDCL.** kissat, kissat on a satsuma lex-leader version and CP-SAT
  ran for 25 to 50 minutes each without an answer; kissat on the hard cube
  below is still running.
- **Cube and conquer.** With one edge pinned and value precedence on the other
  colours, `march_cu` returns 8614 cubes, and `drat-trim` checks that they cover
  every case. They are almost all degenerate cases, where a colour first
  appears late, and the first 7000 fall within 13 seconds each. The last cube
  only asks that colour 2 appear among the first 23 vertices, and it is as hard
  as the whole problem.
  Cubing the plain CNF on colour variables alone gives cubes of depth 12 or 16
  that each take minutes.
- **Symmetric colourings.** No 5-colouring satisfies `c(gx) = π(c(x))` for a
  rotation `g` of order 4, 8 or 16 about a point and any colour permutation `π`
  whose order divides that of `g` (SAT). For the half-turn only `π = 1` was
  checked, and reflections were not checked.

**The plane over `ℚ(√5, √7)`.** It contains a unit 5-cycle,
`1 + s̄ − t + s − t̄ = 0` with `t = (2 + i√5)/3` and `s = (1 + i√35)/6`, so
`3 ≤ χ ≤ 5`. The Minkowski ball `U + U + U` of radius 1.5 for the 130 unit
vectors of height at most 60 (187 003 points) is 3-colourable. A three-colour
growth from `U + U` (8 581 points) found candidate points seeing all three
colours only in its first 19 steps; after that none, up to 50 000 points. So a
4-chromatic graph over this field, if there is one, is not found by local
growth from these units.

**Also.** `H_q` is the hyperbola graph `HG(𝔽_q)` of Bardestani and
Mallahi-Karai (arXiv 1507.05300), which appears inside the graph of every
isotropic quadratic form.

## Finite planes: six colours from the three-point bound (25 September)

The lower bounds of Moorhouse's table (continued in
`notes/local_colourings.md` §12) stopped at 5 because SAT cannot see the
obstruction. Local search suggests that `α(𝔽_q²)` is far below `q²/5` (it finds
`0.13–0.17 q²` for `23 ≤ q ≤ 43`), but a small `α` is a counting fact that
CDCL does not find. kissat ran 40 minutes on `𝔽₂₉²` with five colours, and 30 on `𝔽₂₃²`
with six, without an answer. Both kinds of plane are vertex-transitive, so
`χ ≥ q²/α`, and it is enough to bound `α` (`notes/local_colourings.md` §14).

**What did not work.**
- **Hoffman and Delsarte.** The ratio bound is `0.21–0.25 q²` for
  `23 ≤ q ≤ 53`; the Delsarte linear programme on the circles gives the same
  numbers.
- **Local configurations.** Conditional constraints from unit triangles,
  Moser spindles and rigid templates, over all placements, lower the linear
  programme by at most 2%.
- **Local ratios.** A subgraph needs about half the plane before its
  independence ratio falls below 1/5.
- **CP-SAT on `α`.** Its bound for `𝔽₂₃²` stayed at 163 after 10 minutes.

**The three-point bound.** Schrijver's semidefinite bound for codes, moved to
the plane: the pair frequencies `g` and triangle frequencies `z` of a
recentred, randomly turned independent set satisfy two positive semidefinite
conditions (`M₁ = [z(a, b)]`, `M₀ = [g(b − a) − z(a, b)]`) and linear ones. The
rotation group splits both matrices into real blocks of size about `q`
(`scripts/threepoint.py`).
- **Close for small `q`.** `𝔽₇²`: 14.09 (`α = 14`); `𝔽₁₁²`: 29.18 (28);
  `𝔽₁₃²`: 42.32 (39 found).
- **Solvers.** Clarabel needs about 2 GB per thousand triangle classes and was
  killed at `q = 37`; cvxopt's interior point takes five minutes an iteration
  there; SCS stalls. DSDP solves the primal of `q = 37` in six minutes, but its
  dual is poor: at `q = 23` a certificate built from DSDP's dual gave 726, and
  276 after an LP polish, where the optimum is 108.55. The dual is therefore
  recomputed on the null spaces of the blocks at DSDP's optimum, a small
  semidefinite programme. `scripts/threepoint_verify.py` checks it in interval
  arithmetic, with an exact rational `LDLᵀ` for each dual matrix.

**Results, certified** (`data/threepoint/`):

| plane | `n` | `α ≤` | `n/5` |
|---|---|---|---|
| `G₂₉` | 841 | 163.25 | 168.2 |
| `𝔽₃₇²` | 1 369 | 259.90 | 273.8 |
| `G₃₇` | 1 369 | 263.64 | 273.8 |
| `𝔽₄₁²` | 1 681 | 327.68 | 336.2 |
| `G₄₁` | 1 681 | 300.73 | 336.2 |
| `𝔽₄₃² = G₄₃` | 1 849 | 347.79 | 369.8 |
| `𝔽₄₇² = G₄₇` | 2 209 | 371.42 | 441.8 |

So `χ ≥ 6` for all seven. With the colourings of `notes/local_colourings.md`
§12, `χ(𝔽₄₁²) ∈ {6, 7}`. With Proposition B, `χ(G_q) ≥ 6` for every prime
`q ≥ 29` except 31. The same programme gives `α(G₁₃) ≤ 42.64`, so
`χ(G₁₃) ≥ 5` (`data/threepoint/inert13.npz`).

> **Corrected later** (pre-release audit, 26 September). The step needs
> integrality: `169/42.64 < 4`, but `α` is an integer, so `α(G₁₃) ≤ 42` and
> `χ(G₁₃) ≥ 169/42 > 4`. The test floors the bound in the same way.

**Spectral bounds for large `q`** (interval arithmetic,
`scripts/finite_hoffman.py`): `χ(𝔽_q²) ≥ 6` for `q = 59` and, by Weil, every
prime `q ≥ 67`; `χ ≥ 7` for `q = 71, 97, 101` and every prime `q ≥ 103`.

**A correction.** The proof of Proposition B (`χ(G_q) ≥ 6` for `q ≥ 53`) said
that Weil's bound covers `q > 62` and that 53 and 59 were computed. It did not
mention 61, which Weil's bound does not reach. The computation gives
`χ_f(G₆₁) ≥ 5.20`, so the statement stands; the proofs in both notes now name
61, and `tests/test_finite_planes.py` checks the three cases.

**Not yet** (solver values, not certified).
- **`𝔽₂₃² = G₂₃`.** The plain bound is 108.55. With the localizing matrices of
  a unit edge and a unit triangle, the conditional triangle inequalities and
  the target 106 in the corners, it is 107.04; proving `χ ≥ 6` this way needs
  `α ≤ 105`.
  Almost all of the gain comes from the triangle's localizing matrix (107.07
  with it alone; 108.39 with the edge's, 108.47 with the triangle
  inequalities, 108.55 with the target alone). A unit pentagon adds nothing.
- **`G₃₁ = 𝔽₃₁²`:** 200.89 against 192.2. **`𝔽₂₉²`:** 184.77 against 168.2.
- **`G₁₇`:** 64.56, and 63.33 with the conditional triangles and the
  localizing matrices, against 57.8.
- **Colourings instead of independent sets.** A 5-colouring says more than
  that each class is independent: recentred at a random point, with the other
  four colours in random order, the colour indicators have a component on the
  standard representation of `S₄`, which gives
  `[3 p_ABB(a, b) − p_ABC(a, b)] ⪰ 0`, and `p_ABC ≥ 0` on every triangle.
  (`p_ABB`: the other two vertices share a second colour; `p_ABC`: three
  colours.) A real 5-colouring of `𝔽₁₃²` satisfies every constraint, but the
  bound moves by about 0.01 at most (`G₁₇`: 64.556 to 64.545; `𝔽₂₃²`: 108.547
  to 108.542). Outside the
  rotation-invariant block the new condition reads `4 M₀ ⪰ M₁`, which the
  optimum of the plain programme already nearly satisfies
  (`scripts/experiments/threepoint_colouring.py`).
- **Larger `q`.** `𝔽₅₃²` has 12 827 triangle variables and `𝔽₆₁²` 19 501,
  against 8 625 at `q = 47`, where DSDP used 6 GB; not run. Seven
  colours for `𝔽₄₇²` would need the bound below 369 instead of 371.42; the
  triangle's localizing matrix might give it, but its blocks there have sizes
  up to about 740.

**`G₁₃`.** The anisotropic plane over `𝔽₁₃` was missing from the table of
`notes/local_colourings.md` §4. It is 6-colourable (SAT; the best linear
colouring needs 7), and the three-point certificate gives `χ(G₁₃) ≥ 5`. No
5-colouring turned up: tabu search stops at six monochromatic edges in six runs
of `2·10⁷` moves, and CaDiCaL ran 25 minutes without an answer. Independent
sets of 36 points exist, above `169/5`, so the bound on `α` cannot decide it.
kissat is running on the 5-colouring problem, with one edge pinned and value
precedence on the other colours. Its DRAT proof grew by 1.4 MB a second, too
fast to keep, so it runs without one; an answer of "unsatisfiable" would have
to be re-proved by cube and conquer with a checked proof per cube. Split into
1 024 cubes by `march_cu`, each cube takes kissat more than two minutes.

## Pre-release audit (26 September)

Before release 1.0.0, five independent reviews read the repository: the two
theorems and the note; the local colourings and `notes/rigidity.md`; the
certificates and data; the literature and attributions; and the consistency of
the documentation. A sixth check went back to the primary sources for every
new bibliographic claim. The mathematics of both theorems of the note stands.
What had to change:
- **Two results were not new.** `χ(ℚ(√−3, √−11)) = 4` follows from Fischer's
  theorem, since that field lies in `ℚ(i, √3, √11)`, the plane over
  `ℚ(√3, √11)`. For `χ(ℚ(√−3, √−11, √−247)) = 5`, Madore's argument for
  `ℚ(√3, √11)` (his Cor. 3.4, Lemma 4.5, Prop. 4.6) gives the upper bound at a
  place over 11, and Exoo–Ismailescu's graph lies in the field after a quarter
  turn. The entries "The gate in relation space, and the local criterion" and
  "The field of `five_rho7` is 5-colourable: reduction at 11" state them
  without that earlier work; the notes and READMEs now credit it.
- **The note overstated its checks.** It said the colourings had been tested
  on graphs of up to about 12 000 vertices, and that a second check had used
  PARI/GP; the repository recorded neither. The largest tested graph has 3 134
  vertices. `scripts/decompositions.gp` now recomputes the prime decompositions
  with PARI/GP, checked by `tests/test_decompositions.py`, and
  `tests/test_q311.py` now also checks `data/ei_rho7.json` and the inertness
  of the places over 2 in `ℚ(i, √3, √11)`.
- **A certificate.** The 803-vertex graph `five_247_c` has no 4-colouring:
  kissat's DRAT proof of the plain formula, 4 006 246 lemmas, checked by
  drat-trim in 463 s (`certificates/five_247_c_no4coloring.json`).
- **Local colourings.** In `notes/local_colourings.md` and `notes/rigidity.md`:
  - the Question of §6 was asked for a minimum of at most `χ(ℝ²)`, which made
    the step from `ℚ(√−3, √−7, √−11)` to `χ(ℝ²) ≥ 6` circular; it is now asked
    for a minimum of at most 6;
  - the three-point bounds at 17 and 41 are for the finite planes, level 1;
    the deeper levels map onto them and may need fewer colours, so "no local
    obstruction" holds only at level 1 there;
  - next to `√−3`, a non-split place can also have residue field `𝔽₃` (for
    example in `ℚ(i, √3)`) or `𝔽₈`; the claim "`q ≡ 5 mod 6` or `q = 2`" in
    "The field of `five_rho7` is 5-colourable: reduction at 11" was wrong;
  - the eigenvalues of `G_q` are Kloosterman sums, `λ(n) = −Kl(1, n; q)`, not
    Salié-type sums;
  - in "Blind edges, the level-2 plane at 17, and fields with no local
    obstruction", the list of `d` misses 21, 33 and 77 (the same field as
    `d = 1`); `Γ(F8) □ Γ(F8)` has the chromatic number of `Γ(F8)`, not 5; the
    places of `L16` over 3 and 7 have residue fields `𝔽₉` and `𝔽₄₉`, so the
    targets there are `H₉` and `H₄₉`; and the `L16` seed has 918 unit vectors,
    459 directions, not 471;
  - `chain35.json` reaches distance `√7` from its start, not 3.5;
  - `χ(G₁₇)` is 5 or 6: SAT finds no 4-colouring and a 6-colouring.
- **Attributions.** Reduction modulo a prime goes back to Woodall (1973) and
  Fischer (1990), and extension by cosets to Fischer's Theorem 1 (1990). The
  spectra of the finite planes are Medrano–Myers–Stark–Terras's (1996), and
  Galois symmetry between distances is Tao's (Polymath16, thread 7). Ismailescu
  had a 103-vertex graph with edges at 1 and `2/√3` and no 4-colouring (thread
  3). The reduction to two-distance witnesses is Exoo–Ismailescu's (arXiv
  1805.00157 and 1805.06055), not stated in 1909.13177. Cranston–Rabern had
  also asked about `ℚ(√3, √11)`.
- **Code.** `scripts/tabucol.py` and `scripts/quotient.py` ran their whole
  experiment when imported; they no longer do. The time limits of
  `hn.coloring` never stopped CaDiCaL, the default solver: pysat's interrupt
  does not reach it (Glucose and MapleChrono stop at once). A limited solve now
  runs CaDiCaL in a separate process. The full test suite had waited hours on
  de Grey's graph under a nominal 30-minute limit; its slow tests now have a
  four-hour limit that holds. The maintained scripts write to
  `HN_OUT` (default `/tmp/hn`) instead of this session's scratch directory, and
  the exploratory scripts name `/tmp/hn` too (a rewrite checked to change no
  other part of their syntax trees).
- **Certificates and data.** The review of the certificates and data rebuilt
  them from their coordinates, and every mathematical claim it could check
  holds. What was wrong was in the descriptions:
  - `certificates/five_247_c_no4coloring.json` and the certificates README
    called `five_247_c` a subgraph of `five_247.json`; only 317 of its 803
    points lie there;
  - `five_247_b.json` was built from `Sa` cut to 327 vertices, not 340: that
    peel, glued to its 120° image with overlap 178, reproduces the file's
    476-point union exactly;
  - the sixteen further points of `pressure3_witness_47.json` lie at
    (9 ± √33)/6, not (3 ± √33)/6 ("The machine that makes pressure, taken
    apart");
  - `five_tuned_16` and `five_rho7` share 404 points, 402 of them the carrier
    they both contain, so their union has 4 081 + 2 403 − 404 = 6 080 points;
  - the kissat run on de Grey's formula without the pinned triangle was
    stopped after 110 minutes without a verdict (Status);
  - `α(G₁₃) ≤ 42.64` gives `χ(G₁₃) ≥ 5` only because `α` is an integer:
    `α ≤ 42 < 169/4`;
  - in nine data files the coordinates generate a smaller field than
    `field_generators` names; the index now gives the field of the
    coordinates.
- **Checks added.** `scripts/check_no4.py` rebuilds each of the 27 graphs in
  `data/` said to have no proper 4-colouring, pins one triangle, and has
  drat-trim check kissat's DRAT proof: all 27 are verified
  (`certificates/data_no4_checks.txt`). drat-trim logs are now stored for the
  Moser spindle, the 19-vertex graph, the two pressure certificates, the two
  multi-distance witnesses and the forced pair of Exoo–Ismailescu's graph H.
- **Code, found by the full test run.** `scripts/six.py` shadowed the `six`
  package whenever a test put `scripts/` first on the import path, and two
  tests failed; it is now `scripts/degrey_forced_pair.py`. `hn.cli verify`
  printed `VERIFIED` even when drat-trim was missing and the proof had not
  been checked; it now prints the solver's `UNSAT` and exits with status 2.
  `hn.cli demo` and `degrey`, `scripts/ei_rebuild.py` and
  `scripts/orbit_witness_test.py` no longer write into the repository.

## The four-colour growth in `ℚ(√3, √5)`: its hard step is 4-colourable (26 September)

The growth in `ℚ(√3, √5)` of 25 September (`MODE=same`) tries to make the
vertices 0 and 10 585, at squared distance 3, alike in every 4-colouring. At
14 148 points and 51 727 edges its solver calls ran out of time, and the step
was left undecided. kissat, given eight hours on it, found in about five hours
a 4-colouring in which the two vertices differ. The colouring is proper on all
51 727 edges, recomputed exactly. The pair is not forced at this size, and
`4 ≤ χ(ℚ(√3, √5)²) ≤ 5` stands.

## `α(G₁₇)`: the best sets are rosettes (26 September)

`χ(G₁₇) ≥ 6` would follow from `α(G₁₇) ≤ 57`, since `5 · 57 < 289`. The place
above 17 would then no longer exclude six colours for `ℚ(√−3, √−7, √−11)` at
level 1 (`notes/local_colourings.md` §5). Plain SAT on "58 independent points"
did not finish, so the search was cut down first.

- **Symmetry.** `Aut(G₁₇)` has order 10 404 = 289 · 36, the affine orthogonal
  group, counted by individualisation and refinement. Lex-leader clauses help:
  - there is one for every group element, on the first positions of a vertex
    order that starts with 0 and one whole circle;
  - they keep a formula satisfiable if it was;
  - on `G₁₁` they bring the proof of `α ≤ 28` from 91 s to about 8 s; longer
    prefixes are slower.
- **Symmetric sets.** No independent set of 58 points is invariant under a
  rotation of order 3. The involution `z ↦ −z` and the two kinds of reflection
  were left undecided after ten minutes each.
- **The rosette.** Two kissat runs found 57-point sets, one with symmetry
  breaking and one without. Both lie in the same orbit (stabiliser of order 2).
  The set is a *rosette*:
  - a point `u`;
  - the whole circle `N(z − u) = 12`, 18 points with no unit distance among
    them;
  - 38 points on the circles `N = 3, 4, 5, 6, 9, 14` around `u`.

  It is `ROSETTE` in `scripts/g17_alpha.py`.
- **Part A.** An independent set that contains a point and its whole circle
  `N = c` has at most 57 points.
  - Only the seven circles with no unit distance inside can occur
    (`c = 4, 5, 9, 11, 12, 14, 15`).
  - The other points then lie among the 90 or 108 vertices adjacent to none of
    the 19, and at most 38 of those are independent.
  - kissat and drat-trim verify all seven cases
    (`certificates/g17_part_a_checks.txt`). For `c = 12` the 38 are reached, by
    the rosette.
- **Part B, open.** What remains is an independent set of 58 points containing
  no point together with its whole circle. `g17_alpha.py --write-b` writes it
  with vertex 0 pinned: 2 206 variables, 21 363 clauses.
  - march_cu cut it into 280 cubes, and 270 of them are refuted in under a
    minute each.
  - The other ten gave no answer in 60 s. Cut again, recursively, they produced
    hard cubes faster than they closed: 2 056 refuted leaves and 2 151 open
    cubes when the run was paused.
  - kissat has run on the whole formula for two and a half hours without an
    answer. Runs on each of the ten hard cubes have started.
- **What the rosette suggests.** Suppose every independent set of 57 points
  contained a point together with its whole circle. A 58-point set contains
  57-point ones, so part A would give `α(G₁₇) = 57` at once. The smaller planes
  behave differently: the largest sets found in `G₁₁` (28 points) and `G₁₃` (36
  points) contain no whole circle.

`χ(G₁₇)` stays 5 or 6.

## Formal proofs of the two theorems (27 September)

- `lean/` proves both theorems of the note in Lean 4 with Mathlib (Lean and
  Mathlib v4.34.1): `Q23.chromaticNumber_eq_four` and
  `Q311.chromaticNumber_eq_four`. `#print axioms` lists only `propext`,
  `Classical.choice` and `Quot.sound` (`lean/axioms.expected`).
- The formal proof follows the note except at one step. The note reads the
  residue field of a prime above 2 from the decomposition of 2 in the ring of
  integers. Lean takes a valuation subring with 2 in its maximal ideal, from
  Chevalley's extension theorem, and computes its residue field directly:
  - over `ℚ(√2, √3)`, `π = (√2 + √6)/2 − 1` is a root of the 2-Eisenstein
    polynomial `X⁴ + 4X³ + 2X² − 4X − 2`, and `2 = π⁴ε` with `ε` a unit;
  - over `ℚ(√3, √11)`, 2 has two places, and the proof works at either: with
    `π = √3 − 1`, and a Hensel-type argument for `(1 + √33)/2`.
- The lower bounds are the 10-vertex chain of `data/chain23.json`, whose 16
  edges are proved to have length exactly 1, and the Moser spindle.
- CI (`.github/workflows/lean.yml`) builds the proofs, compares their axioms
  with `lean/axioms.expected`, and replays each file in Lean's kernel with
  `leanchecker`.
- The note is now typeset in LaTeX (version 5, five pages), with a figure of
  the 10-vertex graph and the citation audit's fixes.
