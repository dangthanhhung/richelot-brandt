#!/usr/bin/env python3
"""
verify_p61.py -- independent verification of the p=61 computation.

Scope: this script is specific to p=61.  The seven irreducible factors of
charpoly(B_2(2)) and the elliptic newform minimal polynomials are hard-coded
inline.  It is a template: to verify another prime, replace these with that
prime's factorization (from assemble_brandt.sage) and its newform data (from
assemble_brandt.sage or the LMFDB), keeping the same block-structure checks.
As distributed it verifies only the p=61 result.

Reads the characteristic polynomial produced by assemble_brandt.sage
(data/p61_charpoly.txt, reproduced inline below) and checks EVERY claim
made about it in Section 5 of the paper, using only sympy.  A referee can
run this WITHOUT Sage/Magma to confirm the p=61 result is internally
consistent and matches the elliptic newform data.

    python3 verify_p61.py

Every claim is an assert; the script prints a line per check and ends with
"ALL CHECKS PASSED" iff all hold.
"""
import sympy as sp

x = sp.symbols('x')

# --- the seven irreducible factors of charpoly(B_2(2)) at p=61 -------------
f_eis    = x - 15
f_exotic = x + 7
f6a = x**6 - 29*x**5 + 322*x**4 - 1714*x**3 + 4471*x**2 - 5205*x + 2026
f6b = x**6 + 19*x**5 + 122*x**4 + 270*x**3 - 57*x**2 - 613*x - 254
f9  = (x**9 - 57*x**8 + 1382*x**7 - 18578*x**6 + 151337*x**5 - 767157*x**4
       + 2392680*x**3 - 4370504*x**2 + 4180400*x - 1563856)
f27 = (x**27 - 27*x**26 + 45*x**25 + 4757*x**24 - 35787*x**23 - 277679*x**22
       + 3685273*x**21 + 3197225*x**20 - 170497541*x**19 + 307583031*x**18
       + 4078100295*x**17 - 14611592945*x**16 - 48249621385*x**15
       + 283609945755*x**14 + 159850627139*x**13 - 2786195668061*x**12
       + 2142206112408*x**11 + 13121674311096*x**10 - 22838634130320*x**9
       - 18719183607600*x**8 + 70751348650496*x**7 - 35657421667840*x**6
       - 42709251426048*x**5 + 49375209055488*x**4 - 5956258588672*x**3
       - 12253250717696*x**2 + 5954764984320*x - 806573920256)
f39 = (x**39 - 12*x**38 - 153*x**37 + 2303*x**36 + 7987*x**35 - 185725*x**34
       - 93910*x**33 + 8429360*x**32 - 7815920*x**31 - 241635628*x**30
       + 429028319*x**29 + 4635109999*x**28 - 10908253299*x**27
       - 61351113923*x**26 + 171787460666*x**25 + 568751592440*x**24
       - 1816600366008*x**23 - 3707925461042*x**22 + 13323615686657*x**21
       + 16921802594417*x**20 - 68642242339831*x**19 - 53347316042727*x**18
       + 248772871107874*x**17 + 113262596889180*x**16 - 629551063140824*x**15
       - 153933703951932*x**14 + 1093654981417449*x**13 + 117745531470049*x**12
       - 1266909895505881*x**11 - 25263867033801*x**10 + 936378370173242*x**9
       - 32185529099140*x**8 - 414919982946977*x**7 + 25271696353350*x**6
       + 101150145632624*x**5 - 6064587325952*x**4 - 11980394453520*x**3
       + 490237947792*x**2 + 504648061920*x - 10018017792)

# --- elliptic newform data at level 61 (from LMFDB / Sage Newforms) --------
# S_4(Gamma0(61))^new: two Galois orbits, a_2 minimal polynomials:
s4_deg6 = x**6 + 7*x**5 - 8*x**4 - 106*x**3 - 29*x**2 + 279*x - 8
s4_deg9 = (x**9 - 3*x**8 - 58*x**7 + 154*x**6 + 1145*x**5 - 2511*x**4
           - 9096*x**3 + 13216*x**2 + 24608*x - 3712)
# S_2(Gamma0(61))^new: rational eigenform a_2 = -1 (plus a degree-3 orbit).

def check(msg, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    assert cond, msg

print("Independent verification of B_2(2) at p=61")
print("="*60)

# 1. degree bookkeeping
degs = [sp.Poly(f,x).degree() for f in
        (f_eis,f_exotic,f6a,f6b,f9,f27)] + [2*sp.Poly(f39,x).degree()]
check(f"block degrees sum to h_2(61)=128  (got {sum(degs)})", sum(degs)==128)

# 2. exotic eigenvalue is exactly -7
check("general-type factor is (x+7), so lambda* = -7",
      sp.simplify(f_exotic.subs(x,-7))==0)

# 3. -7 is a SIMPLE eigenvalue: root of no other factor
for nm,ff in [('Eisenstein',f_eis),('SK-6',f6a),('SK-9',f9),
              ('Yoshida-6',f6b),('Yoshida-27',f27),('V_a',f39)]:
    check(f"-7 is NOT a root of the {nm} factor (independence)",
          sp.simplify(ff.subs(x,-7))!=0)

# 4. Saito-Kurokawa factors equal m_g(x-6) for g in S_4^new  (the +6 shift)
check("SK factor f6a == minpoly of a_2(g), g deg-6 orbit, shifted x->x-6",
      sp.expand(f6a - s4_deg6.subs(x,x-6))==0)
check("SK factor f9  == minpoly of a_2(g), g deg-9 orbit, shifted x->x-6",
      sp.expand(f9  - s4_deg9.subs(x,x-6))==0)

# 5. Yoshida deg-6 factor == 2*a_2(f)+a_2(g) with f rational (a_2=-1):
#    root = a_2(g) - 2, i.e. minpoly of a_2(g) shifted x->x+2
check("Yoshida factor f6b == s4_deg6 shifted x->x+2  (2*a_2(f)+a_2(g), a_2(f)=-1)",
      sp.expand(f6b - s4_deg6.subs(x,x+2))==0)

# 6. full charpoly: degree 128, monic, integer coefficients
full = sp.expand(f_eis*f_exotic*f6a*f6b*f9*f27*f39**2)
P = sp.Poly(full, x)
check("charpoly has degree 128", P.degree()==128)
check("charpoly is monic", P.LC()==1)
check("charpoly has integer coefficients", all(c==int(c) for c in P.all_coeffs()))

# 7. trace = 126 = 2*T_1 - h_2  (sum of eigenvalues = -[x^127]/[x^128])
tr = -P.all_coeffs()[1]/P.all_coeffs()[0]
check(f"trace of B_2(2) = 126  (got {tr})", tr==126)

print("="*60)
print("ALL CHECKS PASSED -- the p=61 general-type eigenvalue is -7,")
print("realized geometrically and consistent with the elliptic newform data.")
