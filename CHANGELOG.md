# Changelog

## v3.0.0 (2026-08-13)
Synchronized with paper v4 ("A structural trace identity and certified
spectra for the Richelot--Brandt graph") and hardened for verification.

Added:
- `verify_all.py`: independent structural + spectral verifier for all 31
  primes, including the EXACT characteristic-polynomial identity
  charpoly(M) = product of the stored factors (python-flint), the Mestre
  and mass checks, the sigma/R(pi) checks with Fix = Tr R(pi), and the
  headline numbers of the paper (Fix counts of Lemma 8.1, max 2 h_2 = 2866,
  N_61 = x+7, N_73 = x+6, N_79 = x+5).
- `verify_positivity_2500.py`: the finite input of Theorem 8.5 -- exact
  evaluation of the Proposition 8.2 class-number formula for d(p) at every
  prime 7 <= p <= 2500, with the crude-bound cutoff pinned at 673.
- `expected_output_verify_all.txt`, `expected_output_positivity.txt`:
  reference transcripts of the two verifiers.
- `magma/`: independent Magma certification of h_2(p), mass and |Aut|
  (self-contained scripts for p in {7, 17, 31, 37, 61}; lattice-class
  records; generator script; own README).
- `LICENSE`, `CITATION.cff`, `CHANGELOG.md`, `requirements.txt`.

Changed:
- `README.md` rewritten: paper title and all statement numbers now follow
  the v4 manuscript (Theorem 2.1, Conjecture 2.2, Theorem 2.3, Lemma 4.1,
  Proposition 6.2, Lemma 8.1, Proposition 8.2, Theorem 8.5).
- `MANIFEST.sha256` regenerated over the full v3.0.0 file set.

Removed / clarified:
- The side file `p61_charpoly.txt` of some working copies is not part of
  this record; the factored characteristic polynomial at p = 61 lives in
  `h2_result_61.json` (single source of truth). Note: `res.TrB = 126` is
  the trace of B_2(2); `res.TrR = 38 = 2*T_1(61) - h_2(61)` is the trace
  of the involution R(pi); the two are different operators.
- `w4data_part1.json` keeps its historical name; it covers all 31 primes
  (there is no part 2).

## v2.0.2 (2026-07-23)
As published: production engines, 31 per-prime records, appendix
generator, reduction verifier, ALRTV comparison.
