# Manifest

```
richelot-brandt/
├── README.md                    overview, requirements, quick start
├── CITATION.cff                 citation metadata (Zenodo DOI added after upload)
├── LICENSE                      MIT
├── .zenodo.json                 Zenodo deposition metadata
├── MANIFEST.md                  this file
├── scripts/                     the certification and assembly scripts
│   ├── verify_three_pillars.py        pure-Python re-derivation of the three pillars
│   ├── verify_p61.py                   pure-Python check of the p=61 result (x+7 factor, blocks)
│   ├── closed_identities_audit.py     dimension/trace identities, 7 <= p <= 2500
│   ├── prove_trace_identity.py        symbolic proof of the trace identity
│   ├── certificate_positivity.py      positivity certificate for the block degrees
│   ├── verify_IK94_typenumber.py      IK94 2T-H type numbers, independent route
│   ├── check_outofsample.py           out-of-sample tests p = 41,43,47,53
│   ├── verify_general_prop.py         all-p ingredients on the eight matrices
│   ├── verify_signlaw.py              sign-law certification
│   ├── check_lmfdb_readings.py        LMFDB cross-check of the newform data
│   ├── gen_appendix_A.py              regenerates Appendix A
│   ├── assemble_brandt.sage           matrix assembly from the geometry (Sage)
│   ├── step_b3_signs_conjugation.sage sigma-action, odd-part polynomials (Sage)
│   ├── al_provenance.sage             Atkin-Lehner provenance (Sage)
│   └── worked_edge_p11.sage           one certified Richelot edge of G_11 (Sage)
├── magma/                       independent quaternion-order certification (Magma)
│   ├── verify_h2_7.m
│   ├── verify_h2_17.m
│   ├── verify_h2_31.m
│   ├── verify_h2_37.m
│   └── verify_h2_61.m
├── matrices/                    assembled Brandt matrices
│   ├── B2_p11.txt                     the 5x5 sample
│   └── README.md
└── data/                        output location for regenerated data
    └── README.md
```

## Toolchains

- The pure-Python certification layer (`verify_three_pillars.py`,
  `closed_identities_audit.py`, `prove_trace_identity.py`,
  `certificate_positivity.py`, `verify_IK94_typenumber.py`,
  `check_outofsample.py`, `verify_general_prop.py`, `verify_signlaw.py`,
  `gen_appendix_A.py`) runs under Python >= 3.9 with `sympy` and `numpy`,
  and reproduces every number in the paper without Sage.
- The geometric assembly (`assemble_brandt.sage`,
  `step_b3_signs_conjugation.sage`, `al_provenance.sage`,
  `worked_edge_p11.sage`) and `check_lmfdb_readings.py` run under
  SageMath >= 9.5.
- The Magma scripts in `magma/` are self-contained and certify the class
  numbers `h_2(p)` from a maximal quaternion order, independently of the
  geometric and arithmetic layers.

All theorem references in the script headers match the paper's numbering
(Theorem 4.10 catalogue, Theorem 5.1 sign law, Theorem 5.4 trace formula).
