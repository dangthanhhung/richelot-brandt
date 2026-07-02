# The Richelot isogeny graph as a Brandt matrix — computational artifact
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20798967.svg)](https://doi.org/10.5281/zenodo.20798967)

This repository contains the scripts, the assembled matrices, and the
machine-checkable certificates accompanying the paper

> Hung T. Dang,
> *The Richelot isogeny graph as a Brandt matrix: computation, spectral
> catalogue, and type numbers.*

Everything reported in the paper is regenerated from these files: the
eight Brandt matrices `B_2(2)` for `11 <= p <= 37`, the larger run at
`p = 61`, the spectral catalogue, the sign law, the trace identity, the
type numbers, and the closed-form predictions. No numerical value in the
paper was entered by hand.

## Requirements

Two independent toolchains are used, by design, so that the geometric
assembly and the arithmetic certification do not share an
implementation.

| Layer | Tool | Used for |
|-------|------|----------|
| Geometric assembly | **SageMath ≥ 9.5** | vertex enumeration, Richelot edges, matrix assembly, Atkin–Lehner provenance |
| Arithmetic certification | **Python ≥ 3.9** with `sympy`, `numpy`, `galois` | the self-contained re-derivation, the dimension/trace identities, the positivity bound |
| Cross-check (optional) | **Magma** | the per-level `magma/verify_h2_*.m` scripts (self-contained; quaternion-order certification of `h2(p)`) |

Install the Python dependencies with

```
pip install sympy numpy galois
```

The Python certification layer (`verify_three_pillars.py`,
`closed_identities_audit.py`, `prove_trace_identity.py`,
`certificate_positivity.py`, `verify_IK94_typenumber.py`,
`check_outofsample.py`, `gen_appendix_A.py`) runs **without Sage**: it
re-derives the published matrices from the stored data block and checks
every identity in exact arithmetic. Sage is required only to *re-assemble*
the matrices from the geometry.

## Quick start

```
# 1. exact-arithmetic re-derivation of all three pillars (no Sage needed)
python scripts/verify_three_pillars.py

# 2. trace identity, all primes 7 <= p <= 2500
python scripts/closed_identities_audit.py
python scripts/prove_trace_identity.py

# 3. positivity of the block degrees (Lemma)
python scripts/certificate_positivity.py

# 4. type numbers, independent IK94 route
python scripts/verify_IK94_typenumber.py

# 5. regenerate Appendix A verbatim
python scripts/gen_appendix_A.py

# 6. independently verify the p = 61 result (no Sage needed):
#    the (x+7) general-type factor, the block structure, and agreement
#    with the elliptic newform data
python scripts/verify_p61.py
```

To re-assemble the matrices from the geometry (requires Sage):

```
sage scripts/assemble_brandt.sage          # matrices for 11..37 and p=61 + battery (i)-(vi)
sage scripts/step_b3_signs_conjugation.sage # sigma-action, odd-part polynomials
sage scripts/al_provenance.sage             # Atkin-Lehner signs from a_p
```

## Contents

### `scripts/` — the scripts of Appendix B, plus helpers

| Script | Proves / produces |
|--------|-------------------|
| `assemble_brandt.sage` | enumerates vertices, assembles the eight matrices by the four-piece construction, runs the acceptance battery (i)–(vi), prints the elliptic newform `a_2`-data |
| `step_b3_signs_conjugation.sage` | the action of `σ` on vertices and the odd-part characteristic polynomials |
| `al_provenance.sage` | derives every Atkin–Lehner sign from `a_p` and cross-checks against Sage's routine |
| `verify_three_pillars.py` | self-contained pure-Python re-derivation: structure, characteristic polynomials, catalogue (with `Mys_p` as a residue), "cleanly absent" via gcd, trace identity |
| `verify_general_prop.py` | all-`p` ingredients on the eight matrices: regularity, strong connectivity, loop bound, primitivity, simplicity of 15 with Sturm count, the kernel and polarization identities, and the dimension-identity reduction |
| `verify_signlaw.py` | certifies the sign law: independent exhaustive re-derivation of `σ`, fixed/2-cycle counts, the signed characteristic polynomial on the odd part, and the gcd sign-separations |
| `verify_IK94_typenumber.py` | reimplements the analytic `2T−H` of IK94, validates against the published table for `7 ≤ p ≤ 53`, feeds the type-number corollary |
| `check_outofsample.py` | the out-of-sample tests at `p = 41,43,47,53` |
| `closed_identities_audit.py` | exact-arithmetic closed weight-3 dimension formulas; verifies the trace identity for every prime `7 ≤ p ≤ 2500`; generates the predictions table |
| `prove_trace_identity.py` | symbolic proof of the trace identity per residue class (8 mod 24, 32 mod 120), with `B_{2,χ}`, `h(−p)`, `h(−2p)`, `h(−3p)` as independent symbols |
| `certificate_positivity.py` | proves the positivity lemma: polynomial lower bound per class mod 120, finite check, the residual `d(p) ≥ 0` for `p ≥ 61` to `p = 2500` |
| `check_lmfdb_readings.py` | cross-checks dimension, Fricke sign and `a_2`-trace of all 38 trivial-character newform orbits against LMFDB |
| `gen_appendix_A.py` | regenerates Appendix A, re-asserting row sums, Mestre symmetry, both mass identities, recomputing supersingular `j`-invariants independently |
| `worked_edge_p11.sage` | helper: produces one certified Richelot edge of G_11 (Hasse-Witt superspecial test + (2,2)-splitting + codomain), the edge displayed in the worked example of Sec. 3.3 |

### `matrices/` — the assembled Brandt matrices

Plain-text and JSON forms of `B_2(2)` for `p ∈ {11,13,17,19,23,29,31,37}`
and the `p = 61` run, each with its vertex labels and automorphism
weights.

### `data/` — supporting data

The pinned newform `a_2`-data, the LMFDB download records used by
`check_lmfdb_readings.py`, and the residue-class symbolic forms.

### `magma/` — independent quaternion-order certification (Magma)

Self-contained Magma scripts `verify_h2_p.m` for `p \in {7,17,31,37,61}`
that certify the class number `h_2(p)` directly from a maximal order of
the quaternion algebra `B_{p,\infty}`, by a route entirely independent of
the geometric assembly. Each is self-contained (the structure constants
are embedded) and runs in Magma with no external input; every check is an
`assert`. They provide a third, lattice-theoretic cross-check of the
vertex counts.

## What is verified, and to what extent

- **The eight matrices** are reproduced verbatim and certified by the
  full battery; at `p = 11` the matrix agrees with Jordan–Zaytman.
- **The catalogue, sign law and trace identity** are checked as exact
  identities, the trace identity symbolically (all primes) and
  numerically through `p = 2500`.
- **The type numbers** are confirmed by three independent routes.
- The newform inputs are cross-checked against the **LMFDB**.

## Extending to other primes

The scripts differ in how far they run new primes without code changes.

- **`assemble_brandt.sage` (geometric engine).** Edit the `PRIMES` list and
  re-run; the matrix, weights, and factored characteristic polynomial are built
  from scratch per prime. The implementation assumes the six Weierstrass points
  of every genus-2 vertex are rational over `GF(p^2)`; this holds for all primes
  tested through `p = 61` and first fails at `p = 73` (the `assert len(bp) == 6`
  in `neighbors()` aborts), where some vertices need `GF(p^4)`. Going further
  requires generalizing `branch()` and the pairing/splitting routines to carry
  Weierstrass points over `GF(p^4)`; the rest is prime-independent.
- **`verify_three_pillars.py` (Python, `p <= 37`).** Does not recompute geometry;
  it re-derives the three pillars from the integer matrices stored in its `DATA`
  dictionary. To add a prime, run `assemble_brandt.sage` for it, paste the printed
  matrix/weights/charpoly/newform data into a new `DATA[p]` entry, and append the
  prime. It verifies pre-computed matrices, it does not generate them.
- **`verify_p61.py` (Python).** Specific to `p = 61`; the factors and newform
  minimal polynomials are inline. It is a template for other primes but as
  distributed checks only `p = 61`.
- **`closed_identities_audit.py` (Python, `7 <= p <= 2500`).** The range is a loop
  bound limited only by running time, not by any structural assumption; raising it
  (e.g. `range(5, 100000)`) verifies the identities and the Section 7 growth
  statements further out. This script never touches the graph.

In short: the closed-form audit extends freely, the geometric engine extends up
to the `GF(p^2)` Weierstrass bound (`p <= 61` in practice), and the two Python
certificates replay pre-computed geometric data rather than generating it.

## License

The code is released under the MIT License (see `LICENSE`). The paper
itself is not covered by this license.

## How to cite

If you use this artifact, please cite both the paper and the archived
release (see `CITATION.cff` and the Zenodo DOI on the release page).
