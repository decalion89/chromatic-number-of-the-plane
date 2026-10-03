# The search code

These scripts found and certified the graphs of `data/quadratic_planes/`
(`notes/quadratic_planes.md` §3–§4). They are the code as it was run on
1–2 October 2026, kept as a record; the maintained checker is
`scripts/verify_quadratic_planes.py`, which needs none of them.

**Layout.** A working directory `SC` holds kissat at `SC/kissat/build/kissat`
and drat-trim at `SC/drat-trim/drat-trim` (the versions that
`scripts/worker_setup.sh` builds), and these scripts sit in a subdirectory of
it (it was `SC/q47/`). The shell scripts and `min3.py` and `certify_q.py` take
`SC` from the environment, or else use the parent directory. Compile the tabu
search in the same subdirectory: `gcc -O2 -o tabucol tabucol.c`. Temporary
files go to `/dev/shm`.

| file | role |
|---|---|
| `kd.py` | arithmetic in `ℚ(√d)(i)`, with `d` from the environment variable `QD` |
| `k47.py` | the same for `d = 47` only, used by the first runs |
| `units_fast.py` | every unit vector `((a + b√d)/D, (c + e√d)/D)` with integer coordinates |
| `tabu.py`, `tabucol.c` | tabu search for `k`-colourings |
| `grow3.py`, `grow3d.py` | the first growth runs over `ℚ(√47)`; `grow3d.py` with `D = 240` reached the 9 139-point graph |
| `grow3q.py` | colouring-guided growth for any `d` (`QD=d`); it found the graphs for `d = 59, 71, 119, 191` (and the first ones for `d = 11, 23`). Points are looked up by exact keys (four 16-bit fields); the first version used a hash, and two runs with `D = 1170` stopped on a collision of keys (the code detects them, so no edge was wrong) |
| `grow3r.py` | the same, but when no candidate is blocked it 3-colours the graph again from scratch, and then adds candidates whose neighbours see two colours; it found the graphs for `d = 23, 35, 95, 131, 155, 179, 239, 251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959` |
| `exactcheck.py` | an exact check of a grown graph and of its formula |
| `min3.py` | shrinking: rounds of drat-trim cores, then deletion of vertices one at a time, to a vertex-critical graph (two copies with fixed temporary names were run for `d = 11` and `d = 23`) |
| `min3fast.py` | the same shrinking, faster on large graphs: several random seeds per core round, the two vertices of the fixed edge always kept in the core (`min3.py` drops them, and then stops too early), and tabu search before kissat in the deletion phase; used for `d = 431` |
| `min3fast2.py` | `min3fast.py`, but it saves the vertex set of every core round (`OUT_core{N}_pts.npy`, `N` its size), and in the deletion phase it replaces the vertex set by the drat-trim core of each refutation; a fifth argument starts it from a saved core. With `min3multi.py` started from one of these cores it gave the graphs for `d = 191` and `455` (2 October) |
| `min3inc.py` | the same shrinking with one incremental solver (CaDiCaL 1.5.3 through pysat): a selector literal per vertex, each test a solve under assumptions, and the core of failed assumptions of every refutation as the new vertex set; seconds instead of minutes |
| `min3multi.py`, `reshrink.sh` | `min3inc.py` repeated with up to 400 random deletion orders, keeping the smallest vertex-critical graph; `reshrink.sh` ran it on every grown graph and certified the results; the published graphs come from it (for `d = 191`, `239`, `455` and `935`, `python3 min3multi.py STATE.json CORE_pts.npy OUT 1000 1200` (`3000` for `d = 455`), with `CORE_pts.npy` a core of `min3fast2.py` of 122, 449, 78 and 315 vertices) |
| `certify_q.py` | the certificate of a critical graph: exact edges, unlisted unit pairs, triangles, two encodings with kissat and drat-trim, a 4-colouring, and 3-colourings of `G − v`. Since the certificate of `d = 155` it finds the colourings of `G − v` with one incremental CaDiCaL solver (a selector per vertex) and by rotation (a colouring of `G − v` in which exactly one neighbour `w` of `v` has colour `k` gives, with `v` coloured `k`, one of `G − w`); the earlier certificates called kissat once per vertex |
| `export_q.py` | writes `q{d}.json`, `q{d}.cnf` and `q{d}.logs/` from a certificate directory |
| `pipeline_q.sh` | `exactcheck.py`, kissat with a DRAT proof, drat-trim, `min3.py` and `certify_q.py` for one grown graph |
| `scan_q.sh`, `scan_r.sh`, `scan_r2.sh` | the growth over a list of pairs `d:D` (with `grow3q.py`, and with `grow3r.py`; `scan_r2.sh` takes the size limit from `MAXN`), and the pipeline on the first success for each `d` |
| `scanD.py` | the number of unit vectors for each denominator `D` |
| `gates.py`, `gates2.py` | which directions are integral at the places above 2, 3 and 5 (the gates of §3 of the note) |
| `gate2adic.py` | for `d ≡ 3 (mod 8)`: is the Cayley graph of the directions modulo `2^m` 3-colourable? (it was not, for `m ≤ 3`, in every case we tried except `d = 203`, `D = 1105`) |
| `cohen_check.py` | a numerical check of the remark on Cohen's conjecture (§1 of the note) |
| `oddmitm.py`, `oddcount.py` | the shortest odd closed walk among the directions `U_D` (meeting in the middle), and how many points at 3 steps are also at 4 steps |
| `colstruct.py` | is the 3-colouring saved by a stopped growth run a function of the point modulo `mM`? |
| `prescreen.py` | for each `d` and `D`: the number of directions, the gates of `gates2.py` and the shortest odd closed walk; it chose the denominators of the scan of 1–2 October |
| `scan_r5.sh`, `shrink_worker.sh` | the growth over a list of pairs `d:D` with `grow3r.py`, queueing every graph with no 3-colouring for the worker, which shrinks it with `min3inc.py` and certifies it with `certify_q.py`, one at a time |
| `scan5.py` | for each `d`, the denominators `D` (up to a bound, with small prime factors) whose directions pass every gate and close a 5-cycle, five unit vectors with sum 0; it chose the denominators for `d = 263, 299, 407, 599, 935, 959` |
| `odd_published.py` | for each published graph: does `U_D` close a 5-cycle, and the shortest odd cycle of the graph (5 for `d = 11, 23, 35, 71, 119, 131, 191, 251, 263, 455, 599, 935, 959`, 7 for the others) |
| `gatescan.py` | for each `d`: the denominators `D` (multiples of 6, small prime factors) with 40 to 400 directions and every gate of `gates2.py` open, by the number of directions; the best one for each open field went to the scan of 2 October |
| `scan_r8.sh` | that scan: `grow3r.py` over a list `d:D:RAD:SOFT` (radius of the seed and number of soft rounds per field), two copies sharing a list of finished fields, every graph with no 3-colouring queued for `shrink_worker.sh`; it found the graph for `d = 443` (`D = 2652`) |
| `targets.py` | the vectors `m` of the module spanned by `U_D` (sums of at most three directions) at a spindle distance: `4|m|² − 1 = t²` with `t` in `ℚ(√d)` |
| `growforce.py` | `grow3r.py` aimed at a forced pair: the graph contains `0` and `m` (`FORCE`), every colouring it uses gives them different colours, and it stops when that is impossible; it then solves the graph again without the pair, and reports a forced pair only if the graph itself is 3-colourable (the first runs did not check, and reported pairs over `ℚ(√11)`, `ℚ(√119)` and `ℚ(√455)` in graphs with no 3-colouring) |
| `spin.py` | the spindle of a forced pair: the graph and its copy rotated about `0` by `u = (t + i)/(t − i)`, which moves `m` by exactly 1, checked exactly, confirmed not 3-colourable, shrunk over `RUNS` random orders and written in the format of `certify_q.py`, and rotated back when every point lies in the copy. It refuses a graph with no 3-colouring unless `PLAIN=1`; `PLAIN=1 RUNS=300` on the `d = 11` graph of `growforce.py` (`D = 30`, `FORCE=-15,0,24,-9`) gives the published graph for `d = 11` (rotated back) |
| `growforce4.py`, `spin4.py` | the same at K colours (`KCOL`, default 4), for `ℚ(√47)`: a search for a pair forced to the same colour in every 4-colouring at the spindle distance `|m|² = 403/100 − (99/200)√47` (`m = (−371, 36, −180, 5)/240`), whose spindle would need five colours (§6 of the note). Several colourings per round (`MULTI`), soft rounds, 64-bit hash keys checked exactly; `SEED`, `KISSAT` and `OUTDIR` come from the environment, so that copies with other seeds can run on other machines |
| `mkcnf_state.py` | rebuilds the forcing formula of a saved `growforce4.py` state (a proper K-colouring plus `c(0) ≠ c(m)`, that pair fixed to colours 0 and 1), to decide a formula on which the growth timed out with a longer run elsewhere |
| `rotunion.py`, `overlap.py` | unions of a graph with copies rotated about `0` (by any `t`), and the rotations `w v̄` (`v, w ∈ U_D`) ranked by how many points they keep in the graph |
| `pentagons.py` | unit pentagons over `ℚ(√d)` on three directions of `U_D`: `Q = u₁ + u₂ − v₁` closes with two unit vectors of any denominator when `|Q|²(4 − |Q|²)` is a square |
| `periodicq.py` | is there a 3-colouring periodic modulo `mM`, where `M` is the lattice spanned by the directions `U_D`? (for `d = 83` there is none for the `m` and `D` of §6 of the note) |
| `periodic4.py` | the same at K colours, with kissat: for `d = 47`, `D = 240` and four colours there is no periodic colouring modulo `mM` for `m = 7, 9, 11` (for the other `m ≤ 10` a direction lies in `mM`) |
| `padic_planes.py` | bounds for the chromatic number of the `p`-adic plane `ℚ_p²` (note §5): checks the colourings of `𝔽₃²`, `𝔽₇²`, `𝔽₁₁²`, the 5-cycle over `ℚ(√7)` and which certified fields embed in `ℚ_p`, and prints the table up to `p = 83` (`χ(ℚ₂²) = 2`, `χ(ℚ₃²) = 3`, `χ(ℚ₇²) = 4`) |
| `padic_reach.c`, `padic_reach.py` | which primes `p ≡ 3 (mod 4)` the 27 certified fields reach (note §5): the C program (a segmented sieve with Jacobi symbols) finds the first prime after 3 where no `d` is a nonzero square, `2 129 503 819`, so `χ(ℚ_p²) ≥ 4` for every prime `p ≡ 3 (mod 4)` from 7 below it; the Python script checks the C program's list of `d` against the data, that prime (Miller–Rabin, Euler's criterion), and every prime below `10⁷` again. `padic_reach 30000000000 100` lists the 17 primes below `3·10¹⁰` that no `d` reaches (about six minutes); the script checks them too, and that the `d` have rank 25 modulo squares, so that these primes have density `2⁻²⁵` among the primes `≡ 3 (mod 4)` |
| `padic_measurable.py` | lower bounds for the measurable chromatic number of `ℚ_p²`, `p ≡ 3 (mod 4)` (note §5): the least eigenvalue of `𝔽_p²` in interval arithmetic, one character per norm class, and `χ_m ≥ 1 + 1/|m|` (5 from `p = 23`, 6 at `p = 59` and from 67, 7 at 71 and from 103) |
| `level3.py` | the 11-adic plane at level `k`, `Cay((ℤ/11^k)², T_k)`, modulo a group of rotations and translations that fixes no edge: a 4-colouring of it would 4-colour every plane over a field that embeds in `ℚ₁₁`, such as `ℚ(√47)` and `ℚ(√3, √5)`; at level 3 both groups tried give no 4-colouring (§6 of the note) |
| `grow4r.py` | `grow3r.py` at four colours: colouring-guided growth towards a graph with no 4-colouring (over `ℚ(√47)`, for §6 of the note) |
| `levelp.py` | the same for any prime `p ≡ 3 (mod 4)`: the `p`-adic plane at level `k` modulo the rotations `t ≡ 1 (mod p^s)`, written as a 4-colouring CNF (at `p = 19`, level 2, the rotations `t ≡ 1 (mod 19)` leave no 4-colouring) |
| `levelt.py` | the `p`-adic plane at level `k` modulo the translations by `p^j·(1, 0)`, as a 4-colouring CNF (at `p = 19`, level 2, `j = 1`: no 4-colouring) |
| `quotlevel.py` | the same for the translations by `p^(k−1)·(v_x, v_y)` in any direction, as a quotient graph (at `p = 19`, level 2, direction `(1, 1)`, whose norm is not a square mod 19: kissat UNSAT, drat-trim VERIFIED) |
| `two_colour_criterion.gp` | Theorem A of `notes/local_global.md`: `χ(F²) = 2` iff a prime of `F` above 2 ramifies in `F(i)`, computed from the relative discriminant of `F(i)/F` (PARI/GP); checks it against Johnson and Fischer for the 242 squarefree `d ≤ 400` and against Moorhouse's Theorem 7.1 for nine fields of odd degree |
| `odd_walks.gp` | the new half of Theorem A on quartic fields: exact integer relations among unit vectors `z/z̄` with odd coefficient sum (closed walks of odd length) where no prime above 2 ramifies in `F(i)`, and only even relations where one does |
| `admissible.py` | the real multiquadratic fields with no local 4-colouring and no unit triangle (`notes/local_global.md` §4): `√3 ∉ F` and `i` in every completion above 2, 3, 7, 11 and 19; the first quadratic one is `ℚ(√167)`; also the real-place example `d = 186 023` of `papers/two-colour-planes` |
| `ramified_levels.py` | the levels of the plane over a ramified quadratic extension `ℚ_p(√(pm))` and their colourability (levels 1 to 3 of `ℚ₇(√7)` have no 3-colouring) |
| `hyperbola_plane.py` | `H_q = Cay(𝔽_q², {(t, 1/t)})` for `q = p²` and its colourability (`H₂₅` has no 4-colouring) |
| `liftline.py` | the line encoding of level `k + 1` over a set of level-`k` points: a point of level `k + 1` is a level-`k` block with a next digit `x ∈ 𝔽_p²`, and between adjacent blocks the edges join each line `u·x = s` of one block completely to the line `u·y = s − e/2` of the other, and nothing else; line variables keep the 4-colouring formula linear in size (the whole of level 2 at `p = 19`: 795 644 variables, 5 617 523 clauses). Writes the CNF and runs kissat |
| `checkline.py` | checks that line description against the definition of the graph: at `p = 11` every block edge of level 1 with all `11⁴` pairs of digits, at `p = 19` 400 random block edges with all `19⁴` pairs, plus 20 000 random pairs of points; no mismatch |
| `liftcore.py` | obstruction lifting: the preimage at level `k` of a set of level-`k − 1` points with no proper 4-colouring, refuted by kissat (drat-trim checks the proof and gives a core, which is shrunk while it shrinks); at 11 it took a 69-point set of level 1 to levels 2 and 3 (`data/quadratic_planes/padic11.json`) |
| `check_lift.py` | an independent check of one lifting step (its own enumeration of the unit vectors and of the points, the colour-major encoding, another pinned triangle; kissat and drat-trim, with logs) |
| `crit4.py` | shrinks a set of points of a level with no proper 4-colouring to a vertex-critical one, with one incremental CaDiCaL solver (the 244 points of level 2 in `padic11.json`) |
| `spectrum.py` | the eigenvalues of the `p`-adic plane level by level, as fractions of the degree: the new ones at level `ℓ ≥ 2` are at most `2/(p + 1)`, so Hoffman's bound is the same at every level |
| `hoffman_padic.py` | Hoffman's ratio for the level-1 planes `𝔽_p²`, `p ≡ 3 (mod 4)`: below `1/4` (so no level is 4-colourable) for every `p ≥ 23` (checked below 700; the bound `2√p` on the eigenvalues covers `p ≥ 37`), and not for `p = 3, 7, 11, 19` |

A typical run, for `d = 191`:

```sh
QD=191 python3 grow3q.py q191d240 16000 900 240 2 200
mkdir unsat_q191d240
cp grow3_q191d240_state.json grow3_q191d240_unsat_pts.npy unsat_q191d240/
cp /dev/shm/g3_q191d240.cnf unsat_q191d240/q191d240.cnf
bash pipeline_q.sh 191 q191d240
python3 export_q.py 191 unsat_q191d240/cert191 OUTDIR
```
