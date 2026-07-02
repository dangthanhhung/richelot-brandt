#!/usr/bin/env python3
# ===========================================================================
# verify_general_prop.py -- backs every CHECKABLE ingredient of the general
# (all-p) results: the all-prime catalogue (Theorem 4.10) and the
# accompanying remarks.
#
# The Proposition itself is a pen-and-paper theorem; this script verifies
# (a) its computational ingredients on all eight matrices, and (b) the exact
# algebraic identities quoted in its proof and in the Remark.
#
#   G0  import verify_three_pillars  -> re-runs the full certification.
#   G1  row sums = 15 = (2+1)(2^2+1)  (Lagrangian count in Sp(4,F_2)).
#   G2  irreducibility (strong connectivity) of every matrix  [general: JZ20].
#   G3  square product vertices P{a,a}: weight e = 2|Aut|^2 with |Aut| in
#       {2,4,6}, and diagonal entry  M_vv >= |Aut|/2 >= 1  (the loop bound
#       proved in the Proposition).
#   G4  primitivity: some boolean power of M is everywhere positive
#       (irreducible + positive diagonal => aperiodic).
#   G5  15 is a SIMPLE root of charpoly; the cofactor C_p satisfies
#       C_p(15) != 0 and C_p(-15) != 0; Sturm count: all deg(C_p) roots of
#       C_p lie in the open interval (-15, 15).
#   G6  the two exact identities used in the loop construction:
#       (i)  deg(alpha - 1) = 2 - tr(alpha) for an automorphism alpha,
#            and 2 - tr is not divisible by 4 for tr in {-1, 0, 1}
#            (=> ker(Aut E -> Aut E[2]) = {+-1});
#       (ii) for a norm-1 quaternion a, the matrix F = [[1, a^{-1}],
#            [1, -a^{-1}]] satisfies F^dagger F = 2*Id  (=> the quotient by
#            the graph kernel carries the product principal polarization).
#   G7  dimension identities (II), (III), (VI) on the verified range, and
#       the SYMBOLIC equivalence:  given (VI) [= Ibu18RMS Thm 3.1], the
#       degree consistency of Theorem 4.10 is equivalent to (III).
#   G8  Eisenstein numerology: 15 = 2*sigma(2) + sigma_3(2) = 2*3 + 9
#       (the 'Yoshida of the two Eisenstein series' reading of
#       the Eisenstein eigenvalue in Sec. 4.1).
# ===========================================================================
import re
import sys, os

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_three_pillars as V3P  # noqa: E402  (G0: full re-certification)

x = sp.symbols('x')
PRIMES = [11, 13, 17, 19, 23, 29, 31, 37]

LABEL_RE = re.compile(r'(P\{[^}]*\}|J\d+):(\d+)')


def parse_labels_weights(evec_str):
    out = LABEL_RE.findall(evec_str)
    return [lab for (lab, _) in out], [int(w) for (_, w) in out]


def is_square_product(label):
    if not label.startswith('P{'):
        return False
    inner = label[2:-1]
    parts = re.findall(r'\([^)]*\)', inner)
    assert len(parts) == 2, label
    return parts[0] == parts[1]


def bool_mat_mul(A, B):
    n = len(A)
    return [[any(A[i][t] and B[t][j] for t in range(n)) for j in range(n)]
            for i in range(n)]


print("=" * 78)
print("GENERAL PROPOSITION -- checkable ingredients on all eight matrices")
print("=" * 78)

# ---- G1 (the constant 15) ------------------------------------------------
assert 15 == (2 + 1) * (2 ** 2 + 1)          # Lagrangians in Sp(4, F_2)

for p in PRIMES:
    d = V3P.DATA[p]
    M = V3P.parse_matrix(d['rows'])
    labels, weights = parse_labels_weights(d['evec'])
    n = len(M)
    assert len(labels) == n == len(weights), (p, "label parse")
    assert all(sum(r) == 15 for r in M), (p, "row sums")          # G1

    # ---- G2: irreducibility (strong connectivity) -----------------------
    A = [[M[i][j] > 0 for j in range(n)] for i in range(n)]
    R = [[A[i][j] or i == j for j in range(n)] for i in range(n)]
    for _ in range(n.bit_length() + 1):                 # transitive closure
        R = bool_mat_mul(R, R)
    assert all(all(row) for row in R), (p, "not strongly connected")

    # ---- G3: loops at square product vertices ----------------------------
    n_squares = 0
    for i, lab in enumerate(labels):
        if is_square_product(lab):
            n_squares += 1
            e = weights[i]
            assert e % 2 == 0, (p, lab, e)
            aut2 = e // 2                                # = |Aut|^2
            aut = sp.sqrt(aut2)
            assert aut.is_integer and int(aut) in (2, 4, 6), (p, lab, e)
            assert M[i][i] >= int(aut) // 2 >= 1, \
                (p, lab, "diagonal", M[i][i], int(aut) // 2)
    assert n_squares >= 1, (p, "no square vertex found")

    # ---- G4: primitivity (=> unique eigenvalue of max modulus) ----------
    assert any(M[i][i] > 0 for i in range(n)), (p, "no loop at all")
    P_ = A
    k = 1
    while not all(all(row) for row in P_):
        P_ = bool_mat_mul(P_, P_)
        k *= 2
        assert k <= 4 * n * n, (p, "not primitive within bound")
    # (irreducible + positive diagonal already implies this; double-checked)

    # ---- G5: simplicity of 15, strictness, Sturm count -------------------
    cp = sp.Matrix(M).charpoly(x)
    cp_poly = sp.Poly(cp.as_expr(), x)
    assert cp_poly.eval(15) == 0, (p, "15 not a root")
    C, rem = sp.div(cp_poly, sp.Poly(x - 15, x))
    assert rem.is_zero
    C = sp.Poly(C, x)
    assert C.eval(15) != 0, (p, "15 not simple")
    assert C.eval(-15) != 0, (p, "-15 is a root")
    # Sturm counts DISTINCT real roots, so test on the squarefree part;
    # together with C(+-15) != 0 this places ALL roots of C in (-15, 15).
    Csf = sp.Poly(sp.prod(f.as_expr() for (f, _) in C.factor_list()[1]), x)
    n_in = Csf.count_roots(-15, 15)          # exact Sturm count
    assert n_in == Csf.degree(), (p, "roots escape (-15,15)",
                                  n_in, Csf.degree())

    print(f"  p={p:2d}: rowsum=15  irreducible  squares={n_squares} "
          f"(loop>=|Aut|/2 OK)  primitive  15 simple  "
          f"all {C.degree()} other roots (counted with mult.) in (-15,15)")

# ---- G6(i): kernel of Aut E -> Aut E[2] is {+-1} --------------------------
# alpha an automorphism: deg(alpha) = Nm(alpha) = 1, so
# deg(alpha - 1) = Nm(alpha - 1) = Nm(alpha) - tr(alpha) + 1 = 2 - tr(alpha).
# alpha trivial on E[2]  =>  alpha - 1 = 2*eta  =>  deg(alpha-1) in 4Z.
# Possible traces of alpha != +-1 (orders 3, 4, 6): tr in {-1, 0, 1}.
for tr in (-1, 0, 1):
    assert (2 - tr) % 4 != 0, ("kernel argument fails", tr)
assert (2 - (-2)) % 4 == 0                    # alpha = -1: consistent

# ---- G6(ii): F^dagger F = 2 Id for norm-1 quaternions ---------------------
from sympy.algebras.quaternion import Quaternion  # noqa: E402

one = Quaternion(1, 0, 0, 0)
test_alphas = [
    Quaternion(0, 1, 0, 0),                                   # i, order 4
    Quaternion(sp.Rational(1, 2), sp.Rational(1, 2),
               sp.Rational(1, 2), sp.Rational(1, 2)),         # order 6
    Quaternion(sp.Rational(1, 2), sp.sqrt(3) / 2, 0, 0),      # zeta_6
]
for a in test_alphas:
    na = sp.simplify(a.norm())
    assert sp.simplify(na - 1) == 0, ("not norm 1", a)
    ainv = a.inverse()
    Fm = [[one, ainv], [one, -ainv]]
    Fd = [[Fm[j][i].conjugate() for j in range(2)] for i in range(2)]  # F^dag
    prod = [[Fd[i][0] * Fm[0][j] + Fd[i][1] * Fm[1][j] for j in range(2)]
            for i in range(2)]
    for i in range(2):
        for j in range(2):
            target = Quaternion(2 if i == j else 0, 0, 0, 0)
            diff = prod[i][j] - target
            assert all(sp.simplify(c) == 0 for c in
                       (diff.a, diff.b, diff.c, diff.d)), \
                ("F^dag F != 2I", a, i, j)

# ---- G7: dimension identities and the symbolic reduction ------------------
# (h2, G'0=Klingen, K=paramodular, G0=Siegel, S2+, S2-, S4+, S4-, g_mys)
# source: verify_v06_numbers.py / MEMORY 2026-06-10 (project, re-asserted).
DIMS = {
    11: (5, 0, 0, 0, 0, 1, 2, 0, 0),
    13: (4, 1, 1, 0, 0, 0, 2, 1, 0),
    17: (8, 1, 1, 1, 0, 1, 3, 1, 0),
    19: (10, 2, 1, 1, 0, 1, 3, 1, 1),
    23: (16, 2, 1, 2, 0, 2, 4, 1, 1),
    29: (24, 5, 2, 4, 0, 2, 5, 2, 3),
    31: (26, 6, 2, 4, 0, 2, 5, 2, 4),
    37: (37, 13, 4, 9, 1, 1, 5, 4, 9),
}
for p in PRIMES:
    h2, Gp, K, G0, S2p, S2m, S4p, S4m, g = DIMS[p]
    S2, S4 = S2p + S2m, S4p + S4m
    assert g == Gp - K, (p, "(II)")
    assert S2p * S4p + S2m * S4m == G0, (p, "(III)")
    assert 2 * Gp - G0 - 2 * K == h2 - 1 - (S2 + 1) * S4, (p, "(VI)")
    # degree consistency of the all-p theorem, instance check:
    assert h2 == 1 + S4 + (S2m * S4p + S2p * S4m) + 2 * g, (p, "deg")
    # deg Mys from the matrices agrees with g = Klingen - paramodular:
    assert sp.Poly(sp.sympify(V3P.DATA[p]['Mys']), x).degree() == g, \
        (p, "deg Mys != Klingen - paramodular")

# symbolic:  D = T + E, i.e. given (VI) [Thm 3.1], deg-consistency <=> (III)
h2s, Gps, Ks, G0s, a, b, c, e = sp.symbols(
    'h2 Gp K G0 S2p S2m S4p S4m')           # a=S2+, b=S2-, c=S4+, e=S4-
S2s, S4s = a + b, c + e
D = h2s - 1 - S4s - (b * c + a * e) - 2 * (Gps - Ks)
T = (h2s - 1 - (S2s + 1) * S4s) - (2 * Gps - G0s - 2 * Ks)
E = (a * c + b * e) - G0s
assert sp.expand(D - (T + E)) == 0, "symbolic reduction D = T + E fails"

# ---- G8: Eisenstein numerology and spin-sum dictionary --------------------
sigma1_2 = 1 + 2          # sigma(2)
sigma3_2 = 1 + 2 ** 3     # sigma_3(2)
assert sigma1_2 == 3 and sigma3_2 == 9
assert 2 * sigma1_2 + sigma3_2 == 15
assert 1 + 2 + 4 + 8 == 15                      # Eisenstein spin parameters
ag, bg, af, bf = sp.symbols('alpha_g beta_g alpha_f beta_f')
assert sp.expand((ag + bg + 2 + 4) - ((ag + bg) + 6)) == 0          # SK
assert sp.expand((ag + bg + 2 * af + 2 * bf)
                 - ((ag + bg) + 2 * (af + bf))) == 0                # Yoshida

print("=" * 78)
print("ALL GENERAL-PROPOSITION INGREDIENTS AND REMARK IDENTITIES PASSED.")
print("  G6(i): ker(Aut E -> Aut E[2]) = {+-1} via deg(alpha-1)=2-tr(alpha)")
print("  G6(ii): F^dag F = 2I verified for three norm-1 quaternions")
print("  G7: (II),(III),(VI) 8/8; deg Mys = dim Klingen - dim paramodular 8/8;")
print("      symbolically: given Thm 3.1, [deg-consistency] <=> [(III)]")
print("  G8: 15 = 2*sigma(2) + sigma_3(2) = 2*3 + 9")
print("=" * 78)
