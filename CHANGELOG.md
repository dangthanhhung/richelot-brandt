# Changelog

## v2.0.2 (2026-07-23)

- Add `scripts/certify_P4_omf5.py` -- certificate (P4): at every prime
  11 <= p <= 149 the general-type factor N_p is compared against the
  type-(G) paramodular eigensystems of the quinary database of
  Assaf-Ladd-Rama-Tornaria-Voight (arXiv:2308.09824; data commit
  1907b8812c3ceb45fd09d47e24f96593116e816a). All 31 levels agree,
  including N_p = 1 at the seventeen primes with delta(p) = 0.
- Add `scripts/verify_reduction.py` -- standalone reference verifier for
  hypotheses (a)-(d) of the reduction proposition (Lemma 5.3 /
  Proposition 5.4), with an embedded self-test at p = 11.
- `scripts/gen_appendix_A.py`: emit the F_{p^2} generator as \omega,
  wrap the per-prime headings in \texorpdfstring, and close the file with
  \bigskip. The regenerated appendix is byte-identical to Appendix A of
  the submitted manuscript.
- `scripts/assembly.py`: corrected stale docstring (the script runs all
  thirty-one primes 11 <= p <= 149; the previous "T3-to-101" header
  predated the extension).
- `scripts/MANIFEST.sha256`: now covers every script and data file in
  `scripts/` (previously only the 34 data files), so the verifying code
  itself is checksummed.
- Housekeeping: removed `__pycache__/`, the `.orig` backup, the duplicate
  `gen_appendix_A_patched.py`, the redundant `backup_run1/` tree, and the
  nested `result-11-149.zip`.
- No numerical output of any certificate changed in this version.

## v2.0.1 (2026-07-02)

- First archived release accompanying the submitted manuscript.
