# prove_trace_identity.py -- SYMBOLIC proofs of the two closed identities:
#   (A)  H_npg(p) = 1 + dim S3(K(p))          [Ibu22 Thm 5.1] vs [Ibu07dim 2.1]
#   (D)  tr R(pi)|_pr = c0(p) + d(p)          the trace identity (Conj 5.3)
#        equivalently  trRpr(p) - trRnpg(p) = nu(p) * (p - (-3|p)) / 12
# Method: every constituent is a closed expression in p whose coefficients
# depend only on p mod 120 (for A) or p mod 24 (for D), linearly in the
# arithmetic blocks {B_{2,chi}, h(-p), h(-2p), h(-3p)}.  For each residue
# class we build both sides symbolically (sympy) and assert the difference
# is the zero polynomial.  Transcriptions are first sanity-checked against
# the pinned numeric functions of closed_identities_audit.py on many primes.
import re, sympy as sp
from fractions import Fraction
from math import gcd

# ---------------------------------------------------------------------------
# Step 0: import the PINNED numeric functions verbatim from the audit script
# (no manual retyping -- extract the function sources by name).
import os
_here = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(_here, "closed_identities_audit.py")).read()

def grab(name, nxt):
    i = SRC.index(f"def {name}(")
    j = SRC.index(f"def {nxt}(")
    return SRC[i:j]

code = "".join([
    grab("jacobi", "is_prime"), grab("is_prime", "squarefree_part"),
    grab("squarefree_part", "fundamental_discriminant"),
    grab("fundamental_discriminant", "h_disc"), grab("h_disc", "h_field"),
    grab("h_field", "B2_fast").split("for D, v in")[0],
    grab("B2_fast", "dim_new").split("assert B2_fast")[0],
    grab("dim_new", "nu_of"), grab("nu_of", "trace_w"),
    grab("trace_w", "al_split"), grab("al_split", "tr_R_pi"),
    grab("tr_R_pi", "closed_rhs").split("TABLE_2TmH")[0],
    grab("closed_rhs", "d_of"), grab("d_of", "S3_K").split("assert all")[0],
    grab("S3_K", "S3_G0"), grab("trR_npg", "H_npg"),
    grab("H_npg", "plus3"),
])
ns = {"Fraction": Fraction, "gcd": gcd}
exec(code, ns)
(jacobi, is_prime, h_field, B2_fast, dim_new, nu_of, trace_w, al_split,
 tr_R_pi, closed_rhs, S3_K, trR_npg, H_npg) = (
    ns["jacobi"], ns["is_prime"], ns["h_field"], ns["B2_fast"],
    ns["dim_new"], ns["nu_of"], ns["trace_w"], ns["al_split"],
    ns["tr_R_pi"], ns["closed_rhs"], ns["S3_K"], ns["trR_npg"], ns["H_npg"])
print("Step 0: pinned numeric functions imported from the audit script.")

# quick re-anchor (the audit's own anchors, repeated here):
assert tr_R_pi(61) == 38 and trR_npg(61) == 8 and H_npg(61) == 8
assert S3_K(61) == 7 and al_split(61) == (1, 3, 9, 6)
print("        re-anchors at p=61 OK (trRpr=38, trRnpg=8, Hnpg=8, S3K=7).")

# ---------------------------------------------------------------------------
# Symbolic environment: all Jacobi/Kronecker symbols are functions of the
# residue class r alone (p > 5 prime):
#   e1=(-1|p): r mod 4;  e3=(-3|p)=(p|3): r mod 3;  e2=(2|p): r mod 8;
#   em2=(-2|p)=e1*e2;    j3=(3|p): r mod 12;        p5=(p|5): r mod 5.
p = sp.symbols('p')
B2s, hp, h2p, h3p = sp.symbols('B2 hp h2p h3p')

def env(r):
    e1 = 1 if r % 4 == 1 else -1
    e3 = 1 if r % 3 == 1 else -1
    e2 = 1 if r % 8 in (1, 7) else -1
    j3 = 1 if r % 12 in (1, 11) else -1
    p5 = 1 if r % 5 in (1, 4) else -1
    return dict(e1=e1, e3=e3, e2=e2, em2=e1 * e2, j3=j3, p5=p5)

# Symbolic mirrors of the closed formulas (same branch structure as the
# numeric transcriptions; sanity-checked below against them).
def sym_S3K(r):
    E = env(r)
    v = (sp.Integer(-1) + (p*p - 1)/sp.Integer(2880)
         + sp.Rational(1, 64)*(p + 1)*(1 - E['e1'])
         + sp.Rational(5, 192)*(p - 1)*(1 + E['e1'])
         + sp.Rational(1, 72)*(p + 1)*(1 - E['e3'])
         + sp.Rational(1, 36)*(p - 1)*(1 + E['e3']))
    if r % 5 in (2, 3):
        v += sp.Rational(2, 5)
    v += sp.Rational(1, 8)*(1 - E['e2'])
    if r % 12 == 5:
        v += sp.Rational(1, 6)
    return sp.expand(v)

def sym_Hnpg(r):
    E = env(r)
    return sp.expand((p*p - 1)/sp.Integer(2880)
        + (p - E['e1'])/sp.Integer(24) + (p*E['e1'] - 1)/sp.Integer(96)
        + (p - E['e3'])/sp.Integer(24) + (p*E['e3'] - 1)/sp.Integer(72)
        + sp.Rational(1, 5)*(1 - E['p5'])
        + sp.Rational(1, 8)*(1 - E['e2'])
        + sp.Rational(1, 24)*(1 - E['j3'] + E['e1'] - E['e3']))

def sym_trRnpg(r):
    E = env(r)
    if r % 4 == 1:
        return sp.expand(sp.Rational(1, 96)*(9 - 2*E['e2'])*B2s
            + hp/sp.Integer(16) + h2p/sp.Integer(8)
            + sp.Rational(1, 12)*(3 + E['e2'])*h3p)
    return sp.expand(sp.Rational(1, 96)*B2s
        + sp.Rational(1, 16)*(1 - E['e2'])*hp + h2p/sp.Integer(8)
        + h3p/sp.Integer(12))

def sym_trRpr(r):
    E = env(r)
    if r % 4 == 3:
        return sp.expand(sp.Rational(1, 96)*B2s + h2p/sp.Integer(8)
            + h3p/sp.Integer(12)
            + ((p - 1)/sp.Integer(48)*(9 - 4*E['e2'])
               + (p - E['e2'])/sp.Integer(16)
               + sp.Rational(1, 12)*(1 - E['e3'])*(3 - E['e2']))*hp)
    return sp.expand(sp.Rational(1, 96)*(9 - 2*E['e2'])*B2s
        + (4*p - 1)/sp.Integer(48)*hp + h2p/sp.Integer(8)
        + sp.Rational(1, 12)*(3 + E['em2'])*h3p
        + sp.Rational(1, 12)*(1 - E['e3'])*hp)

def nu_coef(r):
    # conductor relation: nu(p) = a * h(-p)_field, a depending on p mod 8
    if r % 4 == 1:
        return 1
    return 2 if r % 8 == 7 else 4

def sym_d4(r):
    E = env(r)
    return sp.expand((p - 1)*sp.Rational(3, 12)
                     - sp.Rational(1, 4)*(1 - E['e1']))  # dim_new(4,p)

def sym_S4m(r):
    # S4^- = (d4 - trW4)/2,  trW4 = +nu/2 = nu_coef*hp/2
    return sp.expand((sym_d4(r) - nu_coef(r)*hp/2)/2)

def sym_c0(r):
    E = env(r)
    nu = nu_coef(r)*hp
    return sp.expand(1 + (p - 2 + E['e1'])/sp.Integer(8)
                     + nu*(p - 3 - E['e3'])/sp.Integer(12))  # closed_rhs

# ---------------------------------------------------------------------------
# Step 1 (two-layer sanity): symbolic mirrors == pinned numerics, many primes
prs = [q for q in range(7, 1300) if is_prime(q)]
for q in prs:
    r120, r24 = q % 120, q % 24
    sub = {p: q, B2s: (lambda F: sp.Rational(F.numerator, F.denominator))(B2_fast(q)), hp: h_field(-q),
           h2p: h_field(-2*q), h3p: h_field(-3*q)}
    assert sym_S3K(r120).subs(sub) == S3_K(q), ("S3K", q)
    assert sym_Hnpg(r120).subs(sub) == H_npg(q), ("Hnpg", q)
    assert sym_trRnpg(r24).subs(sub) == trR_npg(q), ("trRnpg", q)
    assert sym_trRpr(r24).subs(sub) == tr_R_pi(q), ("trRpr", q)
    assert sym_d4(r24).subs(sub) == dim_new(4, q), ("d4", q)
    assert sym_S4m(r24).subs(sub) == al_split(q)[3], ("S4m", q)
    assert sym_c0(r24).subs(sub) == closed_rhs(q), ("c0", q)
    assert nu_coef(q % 8) * h_field(-q) == nu_of(q), ("nu", q)
print(f"Step 1: symbolic mirrors == pinned numerics at {len(prs)} primes;")
print("        and nu(p) = a(p mod 8) * h(-p), a in {1,2,4}, verified there.")

# ---------------------------------------------------------------------------
# Step 2 = PROOF of (A): H_npg(p) - 1 - S3K(p) == 0 on every class mod 120
classes120 = [r for r in range(1, 120) if gcd(r, 120) == 1]
for r in classes120:
    diff = sp.simplify(sym_Hnpg(r) - 1 - sym_S3K(r))
    assert diff == 0, ("A fails on class", r, diff)
print(f"Step 2 (PROOF A): H_npg - 1 - S3K == 0 identically on all "
      f"{len(classes120)} classes mod 120.")

# ---------------------------------------------------------------------------
# Step 3 = PROOF of (D): on every class mod 24,
#   (i)  RHS reduction: c0 - 1 - S4m... wait, the identity is
#        trRpr = c0 + d  with  d = trRnpg - 1 - S4m  (granted (A));
#        equivalently  trRpr - trRnpg = c0 - 1 - S4m.
#   (ii) we also confirm the clean form: c0 - 1 - S4m == nu*(p - e3)/12.
classes24 = [r for r in range(1, 24) if gcd(r, 24) == 1]
for r in classes24:
    E = env(r)
    nu = nu_coef(r)*hp
    rhs = sp.expand(sym_c0(r) - 1 - sym_S4m(r))
    clean = sp.expand(nu*(p - E['e3'])/sp.Integer(12))
    assert sp.simplify(rhs - clean) == 0, ("clean form fails", r)
    lhs = sp.expand(sym_trRpr(r) - sym_trRnpg(r))
    assert sp.simplify(lhs - rhs) == 0, ("D fails on class", r,
                                         sp.simplify(lhs - rhs))
print(f"Step 3 (PROOF D): trRpr - trRnpg == c0 - 1 - S4^- == nu*(p-(-3|p))/12")
print(f"        identically on all {len(classes24)} classes mod 24;")
print(f"        all B2, h(-2p), h(-3p) terms cancel exactly.")

# ---------------------------------------------------------------------------
# Step 4: independent end-to-end numeric confirmation (sanity, not proof)
bad = []
for q in prs:
    s4m = al_split(q)[3]
    lhs = tr_R_pi(q)
    rhs = closed_rhs(q) + (trR_npg(q) - 1 - s4m)
    if lhs != rhs:
        bad.append(q)
assert bad == [], bad[:5]
print(f"Step 4: end-to-end numeric identity re-confirmed at {len(prs)} primes.")

print("\nALL ASSERTS PASSED -- prove_trace_identity.py complete.")
print("CONCLUSION: granted the published closed formulas")
print("  [IK94 Rem.3 + HI80], [Ibu22 Thms 1.1, 5.1], [Ibu07dim Thm 2.1],")
print("  [Yamauchi via Ibu22 (4),(5)] and the conductor relation for h(-4p),")
print("the trace identity   tr R(pi)|pr = c0(p) + d(p)   is a THEOREM,")
print("and so is the dual form  H_npg(p) = 1 + dim S3(K(p)).")
