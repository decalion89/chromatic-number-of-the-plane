# Notes

Self-contained technical notes. Each states what is proved and how it is checked, and credits
prior work.

| note | contents |
|---|---|
| [`local_colourings.md`](local_colourings.md) | Colouring the plane over a number field through one prime. Finite planes; which fields can hold a 6-chromatic graph; a short proof of Fischer's theorem χ(ℚ(√3, √11)²) = 4 (§8), and χ(ℚ(√2, √3)²) = 4 (§10). |
| [`g13.md`](g13.md) | `α(G₁₃) = 36`: the anisotropic plane over 𝔽₁₃ has no 37 independent points. A computer proof (SAT with checked proofs): the argument, what was computed, how to check it and what has to be trusted. |
| [`g13_chi.md`](g13_chi.md) | `χ(G₁₃) = 6`: the anisotropic plane over 𝔽₁₃ has no proper 5-colouring. A computer proof (a case split on the size of the largest colour class, then cube and conquer: 136 548 leaves, each refuted by kissat with a DRAT proof checked by drat-trim): the argument, what was computed, how to check it and what has to be trusted. |
| [`flat852.md`](flat852.md) | A Moser-spindle-free 5-chromatic unit-distance graph with 852 vertices, in ℚ(ζ₂₁): why Haugland's directions need their mirror images (reduction at the places above 2), how the graph was found, its certificate (DRAT proofs checked by drat-trim, twice), and a comparison with the other spindle-free graphs. |
| [`rigidity.md`](rigidity.md) | Coset colourings of unit-distance modules: the rotation κ that every coset colouring is blind to, coarse rigidity, and propagation. |
| [`literature.md`](literature.md) | What the literature already had, and what we have not found elsewhere, with sources. |
| [`worker_jobs.md`](worker_jobs.md) | Working note (September 2026): the parallel search jobs for χ(ℝ²) ≥ 6, via the Exoo–Ismailescu route at repulsive distances and Galois orbits. A status record, not a result. |

The full chronological account, including approaches that failed, is in [`../docs/research-log.md`](../docs/research-log.md).
