# Three-point certificates

Each file is a dual solution of the three-point semidefinite programme of
[`scripts/threepoint.py`](../../scripts/threepoint.py) for one finite plane. It proves an upper bound
on the independence number α of the plane, and so a lower bound on its chromatic number, since the
plane is vertex-transitive and χ ≥ n/α ([`notes/local_colourings.md`](../../notes/local_colourings.md)
§14).

| file | plane | n = q² | bound on α | proves |
|---|---|---|---|---|
| `inert13.npz` | G₁₃ | 169 | 42.64, so α ≤ 42 < 169/4 | χ ≥ 5 |
| `inert29.npz` | G₂₉ | 841 | 163.25 < 841/5 = 168.2 | χ ≥ 6 |
| `std37.npz` | 𝔽₃₇² | 1 369 | 259.90 < 273.8 | χ ≥ 6 |
| `inert37.npz` | G₃₇ | 1 369 | 263.64 < 273.8 | χ ≥ 6 |
| `std41.npz` | 𝔽₄₁² | 1 681 | 327.68 < 336.2 | χ ≥ 6 |
| `inert41.npz` | G₄₁ | 1 681 | 300.73 < 336.2 | χ ≥ 6 |
| `std43.npz` | 𝔽₄₃² = G₄₃ | 1 849 | 347.79 < 369.8 | χ ≥ 6 |
| `std47.npz` | 𝔽₄₇² = G₄₇ | 2 209 | 371.42 < 441.8 | χ ≥ 6 |

Here `std` is the plane 𝔽_q² with the form x² + y², and `inert` the anisotropic plane G_q; the
bounds are rounded up. Each archive holds `meta`, `z`, `value` and `status`, described in
[`../README.md`](../README.md). `SHA256SUMS` fixes the contents: run `sha256sum -c SHA256SUMS` in
this folder.

**Checking.** From `research/hadwiger-nelson`,
`python3 scripts/threepoint_verify.py data/threepoint/std47.npz` rebuilds the blocks in interval
arithmetic, proves the dual matrices positive definite by an exact rational LDLᵀ factorisation and
prints a rigorous bound on α (3–5 minutes). `python3 -m pytest -q tests/test_threepoint_certificates.py`
checks all eight against the table above (about 15 minutes).
