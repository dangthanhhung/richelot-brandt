# Independent Magma certification of h_2(p)

This directory is an independent verification line for the class numbers
h_2(p), the Hashimoto--Ibukiyama mass, and every |Aut| value, built from
quaternion Hermitian lattices. It shares no code and no data path with
the Python graph engines in the parent directory: the two computations
meet only in the integers they certify.

Contents:
- `verify_h2_p.m` for p in {7, 17, 31, 37, 61}: self-contained Magma
  scripts (structure constants and class representatives embedded).
  EVERY check is an assert; silence means a failure was found; the final
  line prints the verdict, e.g.
  `VERDICT: h2(61) = 128 INDEPENDENTLY CERTIFIED BY MAGMA`.
  Stages can be toggled at the top of each file; on the free Magma
  online calculator run them one at a time if you hit the CPU limit.
- `h2_lattice_p.json` for p in {11, 13, 17, 23, 29, 31, 37, 61}:
  the underlying lattice-class records (class representatives, |Aut|
  multiset, exact mass). These use a different schema from the graph
  records `../h2_result_p.json`; the two families agree on h_2(p) and on
  the |Aut| multiset wherever both exist.
- `gen_magma_verify.py`: emits a `verify_h2_p.m` script from any
  `h2_lattice_p.json`, so the certification can be extended to the
  remaining primes of the lattice family.

What is certified per script: the embedded structure constants define
the stated maximal order; det(T4) = p^2 and det(2H) = p^4; a positive
control (an explicit GL2(O) transform lands in the same class) and a
negative control (two distinct classes are non-isometric); |Aut| of
every class equals the recorded value; the Eichler--Hashimoto--Ibukiyama
mass closes EXACTLY as a rational number; and all classes are pairwise
non-isometric (theta-series collisions resolved by explicit isometry
tests).
