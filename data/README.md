# Supporting data

The certification scripts are self-contained: the pinned newform
`a_2`-data and the LMFDB readings used by `check_lmfdb_readings.py` are
embedded verbatim in the scripts themselves, so the Python layer runs
without querying any external service.

This directory is the output location for the regenerated artifacts:
running `scripts/assemble_brandt.sage` writes the per-level newform
`a_2`-data here, and `scripts/closed_identities_audit.py` writes the
residue-class symbolic forms here.

`p61_charpoly.txt` is the factored characteristic polynomial of `B_2(2)`
at `p = 61`, exactly as printed by `assemble_brandt.sage`, annotated with
the block structure; it is the input checked by `scripts/verify_p61.py`.
