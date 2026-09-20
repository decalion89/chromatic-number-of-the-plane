# Hadwiger–Nelson: a machine-checkable attack

**How many colours does the plane need, so that no two points at distance exactly 1
share a colour?**

Posed around 1950. Still open. The answer, written χ(ℝ²), is known only to lie in
**{5, 6, 7}**.

| bound | value | who, when | how |
|---|---|---|---|
| lower | ≥ 4 | Nelson, 1950 | the 7-vertex Moser spindle |
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

- ❌ χ(ℝ²) ≥ 6 — the actual goal. Not found.

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

| orbit | `|N|` | hyperplanes missing | blocks |
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

## Honest odds

Polymath16 worked on this for years. The chance that this finds a 6-chromatic
unit-distance graph is small. What it does provide is a correct, fast, fully
certifying search whose negative results are recorded precisely enough to be worth
something on their own — and against which any future claim, from anyone, can be
checked in one command.
