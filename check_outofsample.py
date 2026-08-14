#!/usr/bin/env python3
# Out-of-sample test of the trace identity (Theorem 5.4 /
# Remark 5.x) at p = 41, 43, 47, 53 -- primes with NO graph computation.
#
# Inputs and their provenance:
#  - ORBITS: (dim, eps) per newform orbit, verbatim from the user's
#    al_provenance.sage (Sage 10.7).  eps derived from
#    a_p = -eps * p^(k/2-1) and cross-checked by Sage's independent
#    atkin_lehner_eigenvalue at every orbit (38/38 agreement, 11<=p<=53).
#  - External anchor #2 (LMFDB newspace 41/2/a, user screenshot
#    2026-06-10): single orbit 41.2.a.a, dim 3, A-L sign '-' at 41,
#    Tr a2 = -1, matching minpoly x^3+x^2-5x-1.
#  - IK94: (2T-H, H) from verify_IK94_typenumber.py (analytic formula,
#    validated against the published table on all 13 primes 7<=p<=53).
ORBITS = {
    41: {2: [(3, -1)],            4: [(3, -1), (7, +1)]},
    43: {2: [(1, +1), (2, -1)],   4: [(4, -1), (6, +1)]},
    47: {2: [(4, -1)],            4: [(3, -1), (8, +1)]},
    53: {2: [(1, +1), (3, -1)],   4: [(1, -1), (4, -1), (8, +1)]},
}
IK94 = {41: (32, 50), 43: (19, 55), 47: (44, 72), 53: (33, 93)}
GENUS_X0 = {41: 3, 43: 3, 47: 4, 53: 4}     # dim S2^new = genus X0(p)

total_orbits = 0
for p in sorted(ORBITS):
    S2, S4 = ORBITS[p][2], ORBITS[p][4]
    total_orbits += len(S2) + len(S4)
    S2p = sum(d for d, s in S2 if s > 0); S2m = sum(d for d, s in S2 if s < 0)
    S4p = sum(d for d, s in S4 if s > 0); S4m = sum(d for d, s in S4 if s < 0)
    assert S2p + S2m == GENUS_X0[p], (p, "genus mismatch")
    rhs = 1 + (S4p + S4m) + S2m*S4p - S2p*S4m
    tgt, H = IK94[p]
    assert rhs == tgt, (p, rhs, tgt)
    assert (H + tgt) % 2 == 0
    print(f"p={p}: (S2+,S2-;S4+,S4-)=({S2p},{S2m};{S4p},{S4m})   "
          f"1+{S4p+S4m}+{S2m}*{S4p}-{S2p}*{S4m} = {rhs} = 2T-H  OK   "
          f"(T = {(H+tgt)//2})")
assert total_orbits == 15            # 3+4+3+5 orbits at 41,43,47,53
print("OUT-OF-SAMPLE: the trace identity holds at p=41,43,47,53.")
print("Together with 11<=p<=37 (graphs): TWELVE primes, 11 <= p <= 53.")
print("Subtracted channel S2+*S4- nonzero at p=37 (=4), 43 (=4), 53 (=5).")
