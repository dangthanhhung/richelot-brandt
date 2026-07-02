# Assembled Brandt matrices

`B2_p11.txt` is included verbatim as a self-contained sample (the 5x5
case, with vertex labels, weights, regularity, weighted symmetry, and the
factorization of its characteristic polynomial).

The full set for p in {11,13,17,19,23,29,31,37} and the p=61 run is
emitted by `scripts/gen_appendix_A.py` (which regenerates Appendix A of
the paper from the verified data block) and by `scripts/assemble_brandt.sage`
(which re-assembles them from the geometry). Either script writes the full set, in text and JSON form, into this directory.
