# Richelot--Brandt artifact, version 3.0.0

Data and code accompanying the paper *"A structural trace identity and
certified spectra for the Richelot--Brandt graph"*. This record contains
the complete certified computation of the degree-2 Brandt matrices
B_2(2) on superspecial principally polarized abelian surfaces at every
prime 11 <= p <= 149, their full spectral catalogues with Fricke signs,
the two-sided verification of the eigenvalue--sign refinement
(Conjecture 2.2 of the paper) on that range, the exact positivity
certificate for the signed defect d(p) at every prime 7 <= p <= 2500
(the finite input of Theorem 8.5), and an independent Magma
certification of the class numbers h_2(p) on a subrange.

Concept DOI (always resolves to the latest version):
10.5281/zenodo.20798967. This version: 10.5281/zenodo.21927983.

## Contents

Production engines (as in v2.0.2, unchanged):
- `richelot.py`   -- graph engine (adaptive extension fields); one prime per run argument
- `w2side.py`     -- weight-2 eigenvalue systems + Fricke signs (elliptic Brandt)
- `msym.py`       -- weight-4 eigenvalue systems + Fricke signs (Manin symbols, Merel--Heilbronn)
- `assembly.py`   -- two-sided catalogue matching: exact division + trace-decoded sign certificates modulo 2^61 - 1

Frozen data:
- `h2_result_p.json` (31 files, 11 <= p <= 149) -- per-prime record:
  matrix `M`, weights `wts`, involution `sigma`, and the summary `res`
  (h2, vP, vJ, TrB = trace of B_2(2), TrR = trace of R(pi) = 2 T_1 - h_2,
  factored characteristic polynomial, block degrees). Note TrB and TrR
  are traces of two different operators (at p = 61: 126 and 38).
- `w2data.json`, `w4data_part1.json` -- weight-2/4 orbit polynomials with
  signs, all 31 primes (`part1` is a historical name; there is no part 2).
- `assembly_results.json` -- final ledger: block degrees sk/yos/mys and
  the general-type polynomial N_p at all 31 primes.
- `appendix_A_matrices.tex` -- the eight matrices of Appendix A.

Verifiers:
- `verify_all.py` -- independent structural + spectral verifier, all 31
  primes (python-flint). Per prime: row sums 15; Mestre symmetry
  e_j M_ij = e_i M_ji; the Hashimoto--Ibukiyama mass as an exact
  rational; sigma an involution commuting with M and preserving e, with
  Fix(sigma) = TrR; trace(M) = TrB = sum of the roots of the stored
  factorization; the EXACT identity charpoly(M) = product of the stored
  factors; the Eisenstein factor (x - 15) simple; the recomputed
  square-excess degree; the ledger identity
  1 + sk + yos + 2*mys + deg(N_p) = h_2(p). Globally: the Fix counts of
  Lemma 8.1 (5, 4, 8, 8, 14, 18, 18, 11 at p = 11..37 and 38 at p = 61),
  max_p 2 h_2(p) = 2866 at p = 149 (the modulus bound of Section 6), and
  N_61 = x + 7, N_73 = x + 6, N_79 = x + 5. Run with `--verify-manifest`
  to recompute every SHA-256 in the manifest first.
- `verify_positivity_2500.py` -- the finite input of Theorem 8.5: exact
  evaluation of the Proposition 8.2 class-number formula for d(p) at
  every prime 7 <= p <= 2500 (standard library only); checks that d(p)
  is a non-negative integer throughout, that the crude bound of the
  proof is positive for every prime p > 673, and that 673 is the largest
  prime where the crude bound fails.
- `verify_reduction.py` -- minimal standalone verifier of hypotheses
  (a)--(d) of Proposition 6.2 at p = 11; mirrors Lemma 6.1 and
  Proposition 6.2 line by line and shares no code with the engines.
- `verify_p61.py` -- the worked example of Section 6.9 replayed.
- `verify_three_pillars.py`, `verify_IK94_typenumber.py`,
  `verify_signlaw.py`, `verify_general_prop.py`, `closed_identities_audit.py`,
  `certificate_positivity.py` -- the supporting audits of the closed
  dimension identities and of Lemma 4.1 (delta(p) >= 0, d(p) >= 0), as
  in v2.0.2.
- `certify_P4_omf5.py` -- external comparison of every N_p against the
  ALRTV quinary database.
- `gen_appendix_A.py` -- regenerates `appendix_A_matrices.tex`,
  byte-identical to Appendix A of the paper.

Independent certification:
- `magma/` -- Magma certification of h_2(p), the mass, and every |Aut|
  via quaternion Hermitian lattices, self-contained for
  p in {7, 17, 31, 37, 61}, with the underlying lattice-class records
  and a generator script; see `magma/README_magma.md`. This line shares
  no code and no data path with the Python engines.

Reference transcripts and metadata:
- `expected_output_verify_all.txt`, `expected_output_positivity.txt`
  (timings are machine-dependent; every other line should match).
- `MANIFEST.sha256` -- SHA-256 of every file in this record (except itself).
- `LICENSE` (MIT for code; data additionally CC-BY-4.0), `CITATION.cff`,
  `CHANGELOG.md`, `requirements.txt`.

## Requirements

Python >= 3.9. `pip install -r requirements.txt`: python-flint (engines
and `verify_all.py`), sympy (engines and `verify_reduction.py`).
`verify_positivity_2500.py` needs only the standard library. The
`magma/` scripts need Magma; the free online calculator suffices, one
stage at a time.

## Five-minute referee path

    python3 verify_all.py --verify-manifest --primes 11 61
    python3 verify_reduction.py
    python3 verify_positivity_2500.py

Expected: the manifest line, two per-prime `charpoly OK  blocks OK`
lines, the asserts of `verify_reduction.py`, and
`POSITIVITY (Theorem 8.5, finite input): PASSED, 7 <= p <= 2500`.
The p = 11 record can additionally be checked by eye against Appendix A
of the paper and against Jordan--Zaytman, J. Math. Soc. Japan,
Theorem 39 and Section 10.1.

## Full verification and reproduction

    python3 verify_all.py --verify-manifest       # all 31 primes; about 8 minutes (dominated by the exact 1433 x 1433 characteristic polynomial at p = 149)
    python3 verify_positivity_2500.py             # about 5 seconds

Full reproduction from scratch (order matters; every step saves
incrementally; interruption is safe and rerunning resumes):

    python3 richelot.py 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149
    python3 w2side.py
    python3 msym.py
    python3 assembly.py

Approximate engine runtimes on a 2024 laptop: 30--45 min, 3 min,
25--40 min, 1--2 h (unchanged from v2.0.2).

## Claim-by-claim verification map

| Paper statement | Data | Command | Expected line |
|---|---|---|---|
| Appendix A matrices, p = 11..37 | h2_result_p.json | `python3 gen_appendix_A.py` | byte-identical .tex |
| Lemma 8.1 Fix counts; Tr R(pi) = 2T_1 - h_2 | sigma in records | `python3 verify_all.py` | `global: Fix counts of Lemma 8.1 OK` |
| Section 6 modulus bound 2 h_2(p) <= 2866 < 2^61 - 1 | records | `python3 verify_all.py` | `max 2*h2 = 2866 at p=149 OK` |
| Section 6.9 worked example (Tr B_2(2) = 126, blocks (1,15,33,2*39,1), N_61 = x+7) | h2_result_61.json | `python3 verify_p61.py` and `verify_all.py` | asserts pass; `p= 61 ... charpoly OK` |
| Theorem 2.3 / Conjecture 2.2 on 11 <= p <= 149 (factorizations, signs) | records + w2/w4 + ledger | `python3 assembly.py` | `catalogue EXACT, signs=OK` |
| Proposition 6.2 hypotheses (a)--(d) | records | `python3 verify_reduction.py` (p = 11) and `verify_all.py` (all) | asserts pass; `ALL CHECKS PASSED (31/31 primes)` |
| Theorem 8.5, finite input: d(p) >= 0 for 7 <= p <= 2500, cutoff 673 | -- | `python3 verify_positivity_2500.py` | `POSITIVITY ... PASSED, 7 <= p <= 2500` |
| Lemma 4.1 (delta, d >= 0; cubic growth of Delta_III) | -- | `python3 certificate_positivity.py` | `POSITIVITY CERTIFICATE PASSED` |
| External corroboration (Table 7 vs ALRTV) | omf5 data | `python3 certify_P4_omf5.py ...` | `(P4) PASSED ... 31/31` |
| h_2(p), mass, Aut independently (lattices) | magma/ | run `magma/verify_h2_p.m` | `... INDEPENDENTLY CERTIFIED BY MAGMA` |

## Acceptance run (this version)

2026-08-13, clean Linux container, Python 3.12, python-flint 0.9.0:
`verify_all.py` ALL CHECKS PASSED (31/31 primes);
`verify_positivity_2500.py` PASSED (364 primes 7 <= p <= 2500, largest
crude-bound failure at 673 confirmed); `verify_reduction.py` hypotheses
(a)--(d) verified at p = 11. The engine replay figures are those of the
v2.0.2 acceptance run; no data file changed in this version, so the
frozen records are those certified there.

## How to cite

Please cite the paper and this artifact (see `CITATION.cff`).
Concept DOI 10.5281/zenodo.20798967; version DOI of this record:
10.5281/zenodo.21927983.
