# Richelot--Brandt artifact, version 2.0.2

Data and code accompanying the paper *"The Richelot isogeny graph as a
Brandt matrix: a global eigenvalue correspondence and a structural trace
formula"*. This record contains the complete certified computation of the
degree-2 Brandt matrices B_2(2) on superspecial principally polarized
abelian surfaces, their full spectral catalogues, and the two-sided
verification of the refined Ibukiyama conjecture (C1) at every prime
11 <= p <= 149.

## Contents

- `richelot.py`   -- graph engine (adaptive extension fields); one prime per run argument
- `w2side.py`     -- weight-2 eigenvalue systems + Fricke signs (elliptic Brandt)
- `msym.py`       -- weight-4 eigenvalue systems + Fricke signs (Manin symbols, Merel--Heilbronn)
- `assembly.py`   -- two-sided catalogue matching: exact division + trace-decoded sign certificates
- `h2_result_p.json` (31 files) -- matrix M, weights, Frobenius permutation sigma, factored charpoly
- `w2data.json`, `w4data_part1.json` -- weight-2/4 orbit polynomials with signs
- `assembly_results.json` -- final ledger: block degrees and N_p at all 31 primes
- `certify_P4_omf5.py` -- certificate (P4): external comparison of every N_p against the ALRTV quinary database
- `verify_reduction.py` -- reference verifier for hypotheses (a)-(d) of the reduction proposition (Lemma 5.3 / Prop. 5.4)
- `MANIFEST.sha256` -- SHA-256 checksums of every data file and script in this record

## Requirements and reproduction

Python >= 3.9 and one library: `pip install python-flint`.
Full reproduction from scratch (order matters; every step saves incrementally):

    python3 richelot.py 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149
    python3 w2side.py
    python3 msym.py
    python3 assembly.py

Approximate runtimes on a 2024 laptop: 30--45 min, 3 min, 25--40 min,
1--2 h. Any interruption is safe; rerunning resumes.

## What is certified

Every prime passes, in exact arithmetic: row sums 15; Mestre symmetry
e_j M_ij = e_i M_ji; Eichler and Hashimoto--Ibukiyama mass formulas as
exact rationals; sigma a permutation with Tr R(pi) equal to the closed
formula of Ibukiyama--Katsura (Compositio 1994) and to L_0(p)+d(p)
(the paper's Theorem T2); dimension formulas for the weight-2/4 sign
splits; exact divisibility of the characteristic polynomial by the
Eisenstein, Saito--Kurokawa and Yoshida blocks built from independent
weight-2/4 data; the residual factor an exact square times N_p of degree
delta(p); and the R(pi)-sign of every irreducible block via projector
traces (integer targets bounded by block degree, so a single 61-bit
modulus is a proof). External anchors: the p=11 matrix of
Jordan--Zaytman Section 10.1 (entrywise), and the paramodular non-lift
eigenvalues of Poor--Yuen and Poor--Schmidt--Yuen at p=61,73,79.
A failed certificate aborts the run; no silent output exists.

## Ledger of general-type factors

|p|N_p|
|--|--|
|61|x + 7|
|73|x + 6|
|79|x + 5|
|89|x + 4|
|97|x^2 + 9x + 19|
|101|x^2 + 7x + 11|
|103|x^2 + 9x + 19|
|109|x^3 + 10x^2 + 31x + 29|
|113|x + 3|
|127|x^3 + 9x^2 + 24x + 19|
|131|x^2 + 6x + 7|
|137|x^2 + 4x + 2|
|139|x^4 + 14x^3 + 63x^2 + 110x + 61|
|149|x^4 + 6x^3 + 10x^2 + 3x + (-1)|

(N_p = 1 at the other seventeen primes.)

## Determinism note

All quantities are deterministic. Vertex labelings may permute between
runs that group the primes differently on the command line (the RNG state
at each prime depends on the primes processed before it in the same
invocation); characteristic polynomials, weights multisets, traces and
all certified invariants are the comparison standard and are invariant.


## Certificate (P4): external corroboration

`certify_P4_omf5.py` checks that at every prime 11 <= p <= 149 the
general-type factor N_p equals the characteristic polynomial of T(2,1) on
the type-(G) paramodular newforms of the public database of
Assaf-Ladd-Rama-Tornaria-Voight (arXiv:2308.09824), including N_p = 1 at
the seventeen primes with delta(p) = 0. The data are not redistributed
here; fetch them and run:

    git clone --depth 1 https://github.com/assaferan/omf5_data
    # pinned commit: 1907b8812c3ceb45fd09d47e24f96593116e816a
    python certify_P4_omf5.py omf5_data/hecke_evs_3_0/data_nl_200

Expected final line: `(P4) PASSED: N_p agrees with [ALRTV] at all 31/31
primes 11 <= p <= 149`.

## Reference verifier for the reduction proposition

`verify_reduction.py` is a standalone implementation (pure Python +
sympy, under 200 lines) of hypotheses (a)-(d) of the reduction
proposition of the paper, with an embedded self-test at p = 11 taken from
Appendix A. It mirrors the printed Lemma and Proposition line by line and
shares no code with the generating engines; the production replay at all
thirty-one primes is the flint routine in `assembly.py` (same
mathematics, same Mersenne modulus 2^61 - 1).

    python verify_reduction.py

## Acceptance run (this version)

All four scripts were re-run from a clean environment for v2.0.2:
`verify_three_pillars.py` (all asserts pass, eight primes),
`gen_appendix_A.py` (all asserts pass; the emitted `appendix_A_matrices.tex`
is byte-identical to Appendix A of the submitted manuscript),
`verify_reduction.py` (hypotheses (a)-(d) verified at p = 11), and
`assembly.py` (catalogue EXACT with `signs=OK`; `N_p = x + 7` at p = 61).
