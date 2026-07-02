#!/usr/bin/env python3
# check_lmfdb_readings.py -- cross-check ALL newform data of the paper
# against the LMFDB, for every prime 11 <= p <= 53 and weights k in {2,4}.
#
# EXPECTED side: (p, k, a2-minimal-polynomial, eps_p) exactly as produced by
#   al_provenance.sage (Sage 10.7, 2026-06-10); Tr a2 is computed here with
#   sympy as -[coeff of x^(d-1)].
# LMFDB side: rows transcribed VERBATIM from three download files fetched on
#   10 June 2026 (format: [Label, Dim, A, Field, CM, (RM,) Traces, Fricke
#   sign, q-expansion], sorted by analytic conductor):
#   https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/?download=1&query={'weight': 4}&weight=4
#   https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/?download=1&query={'level': {'$gte': 11, '$lte': 100}}&level=11-100
#   (the level=11-100 file fetched twice with increasing truncation limits).
# Matching convention: orbits are matched by (dim, Tr a2), NEVER by label
#   order (protocol rule).  All comparisons are asserted.

import sympy as sp

x = sp.symbols('x')

# ---- expected: 38 orbits from al_provenance.sage --------------------------
EXPECTED = [
 (11,2,'x+2',-1), (11,4,'x**2-2*x-2',+1),
 (13,4,'x+5',-1), (13,4,'x**2-x-4',+1),
 (17,2,'x+1',-1), (17,4,'x+3',-1), (17,4,'x**3-x**2-24*x+32',+1),
 (19,2,'x',-1),   (19,4,'x+3',-1), (19,4,'x**3-3*x**2-18*x+38',+1),
 (23,2,'x**2+x-1',-1), (23,4,'x+2',-1),
 (23,4,'x**4-2*x**3-24*x**2+61*x+2',+1),
 (29,2,'x**2+2*x-1',-1), (29,4,'x**2+2*x-1',-1),
 (29,4,'x**5-33*x**3+28*x**2+192*x-256',+1),
 (31,2,'x**2-x-1',-1), (31,4,'x**2+5*x+2',-1),
 (31,4,'x**5-3*x**4-30*x**3+79*x**2+167*x-386',+1),
 (37,2,'x+2',+1), (37,2,'x',-1),
 (37,4,'x**4+6*x**3-x**2-16*x+6',-1),
 (37,4,'x**5-4*x**4-21*x**3+74*x**2+102*x-296',+1),
 (41,2,'x**3+x**2-5*x-1',-1), (41,4,'x**3+3*x**2-5*x-3',-1),
 (41,4,'x**7-x**6-49*x**5+33*x**4+720*x**3-320*x**2-3200*x+512',+1),
 (43,2,'x+2',+1), (43,2,'x**2-2',-1),
 (43,4,'x**4+4*x**3-9*x**2-14*x+2',-1),
 (43,4,'x**6-6*x**5-17*x**4+124*x**3+26*x**2-608*x+540',+1),
 (47,2,'x**4-x**3-5*x**2+5*x-1',-1), (47,4,'x**3+5*x**2-2*x-12',-1),
 (47,4,'x**8-3*x**7-50*x**6+124*x**5+844*x**4-1549*x**3-5393*x**2+5418*x+10316',+1),
 (53,2,'x+1',+1), (53,2,'x**3+x**2-3*x-1',-1),
 (53,4,'x',-1), (53,4,'x**4+4*x**3-16*x**2-42*x+49',-1),
 (53,4,'x**8-2*x**7-52*x**6+94*x**5+823*x**4-1136*x**3-4480*x**2+3296*x+6784',+1),
]
assert len(EXPECTED) == 38

# ---- LMFDB readings: label -> (dim, fricke, first trace = Tr a2) ----------
LMFDB = {
 '11.2.a.a':(1,-1,-2), '11.4.a.a':(2,+1, 2),
 '13.4.a.a':(1,-1,-5), '13.4.a.b':(2,+1, 1),
 '17.2.a.a':(1,-1,-1), '17.4.a.a':(1,-1,-3), '17.4.a.b':(3,+1, 1),
 '19.2.a.a':(1,-1, 0), '19.4.a.a':(1,-1,-3), '19.4.a.b':(3,+1, 3),
 '23.2.a.a':(2,-1,-1), '23.4.a.a':(1,-1,-2), '23.4.a.b':(4,+1, 2),
 '29.2.a.a':(2,-1,-2), '29.4.a.a':(2,-1,-2), '29.4.a.b':(5,+1, 0),
 '31.2.a.a':(2,-1, 1), '31.4.a.a':(2,-1,-5), '31.4.a.b':(5,+1, 3),
 '37.2.a.a':(1,+1,-2), '37.2.a.b':(1,-1, 0),
 '37.4.a.a':(4,-1,-6), '37.4.a.b':(5,+1, 4),
 '41.2.a.a':(3,-1,-1),
 '41.4.a.a':(3,-1,-3), '41.4.a.b':(7,+1, 1),
 '43.2.a.a':(1,+1,-2), '43.2.a.b':(2,-1, 0),
 '43.4.a.a':(4,-1,-4), '43.4.a.b':(6,+1, 6),
 '47.2.a.a':(4,-1, 1),
 '47.4.a.a':(3,-1,-5), '47.4.a.b':(8,+1, 3),
 '53.2.a.a':(1,+1,-1), '53.2.a.b':(3,-1,-1),
 '53.4.a.a':(1,-1, 0), '53.4.a.b':(4,-1,-4), '53.4.a.c':(8,+1, 2),
}
assert len(LMFDB) == 38

PRIMES = [11,13,17,19,23,29,31,37,41,43,47,53]

def tr_a2(mp):
    P = sp.Poly(sp.sympify(mp), x); d = P.degree()
    return (P.degree(), int(-P.coeff_monomial(x**(d-1))) if d >= 1 else 0)

# group both sides by (p,k); compare as sorted multisets of (dim, sign, Tr a2)
exp_by, lm_by = {}, {}
for (p,k,mp,e) in EXPECTED:
    d,t = tr_a2(mp); exp_by.setdefault((p,k),[]).append((d,e,t))
for lab,(d,f,t) in LMFDB.items():
    p,k = int(lab.split('.')[0]), int(lab.split('.')[1])
    lm_by.setdefault((p,k),[]).append((d,f,t,lab))

print('=== Per-orbit comparison (matched by (dim, Tr a2)) ===')
print(f"{'p':>3} {'k':>2}  {'label':<10} {'LMFDB (d,sign,Tra2)':>20}   {'expected':>14}  verdict")
n_ok = 0
for (p,k) in sorted(set(exp_by) | set(lm_by)):
    E = sorted(exp_by.get((p,k),[]))
    L = sorted((d,f,t) for (d,f,t,_) in lm_by.get((p,k),[]))
    assert E == L, f'MISMATCH at (p,k)=({p},{k}): expected {E}, LMFDB {L}'
    for (d,f,t,lab) in sorted(lm_by[(p,k)]):
        s = '+' if f > 0 else '-'
        print(f'{p:>3} {k:>2}  {lab:<10} {f"({d},{s},{t:+d})":>20}   {"same":>14}  OK')
        n_ok += 1
# weight-2 space at p=13 must be EMPTY on both sides
assert (13,2) not in exp_by and (13,2) not in lm_by
print(f'... and S_2^new(Gamma_0(13)) = 0 on both sides (no orbit listed).')

print()
print('=== Atkin-Lehner cells (dim+, dim-) reconstructed from LMFDB rows ===')
print(f"{'p':>3}  {'k=2 (+,-)':>10}  {'k=4 (+,-)':>10}")
for p in PRIMES:
    row = []
    for k in (2,4):
        c = [0,0]
        for (d,f,t,_) in lm_by.get((p,k),[]):
            c[0 if f > 0 else 1] += d
        ce = [0,0]
        for (d,e,t) in exp_by.get((p,k),[]):
            ce[0 if e > 0 else 1] += d
        assert c == ce, f'cell mismatch at ({p},{k})'
        row.append(tuple(c))
    print(f'{p:>3}  {str(row[0]):>10}  {str(row[1]):>10}')

print()
print(f'ALL {n_ok}/38 ORBITS AND ALL 24 AL-CELLS MATCH THE LMFDB,'
      f' for every prime 11 <= p <= 53.')
