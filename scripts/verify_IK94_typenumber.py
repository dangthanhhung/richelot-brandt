#!/usr/bin/env python3
# verify_IK94_typenumber.py
# ---------------------------------------------------------------------------
# Independent verification of Ibukiyama--Katsura, Compositio Math. 91 (1994),
# 37-46 ("IK94"), Theorem 2 + Remark 3, PRINCIPAL genus, n = 2, prime level p:
#
#     #{ppas (A,Theta) with a model over F_p} = tr(R(pi)) = 2T - H.
#
# Three independent sources are cross-checked:
#   (a) the analytic formula of Remark 3, recomputed here FROM SCRATCH
#       (generalized Bernoulli number B_{2,chi}; class numbers h(-D) by
#        counting reduced binary quadratic forms -- no tables imported);
#   (b) the published table of IK94, p. 42 (rows H and 2T-H, p = 7..53);
#   (c) our Brandt-graph computation (steps B1-B3): h_2 = #vertices,
#       fixed points and 2-cycles of the Galois-conjugation involution.
# Final target:  T(p) = (H + tr R(pi))/2  ==  predicted T1(p), p = 11..37.
#
# Protocol: NO mental arithmetic -- every claimed number is an assert.
# A red assert is data, not something to patch by hand.
# ---------------------------------------------------------------------------

from fractions import Fraction
from math import gcd
from sympy import isprime, jacobi_symbol


# --------------------------- elementary helpers ----------------------------

def squarefree_part(n):
    """Squarefree d (carrying the sign of n) with n = d * (perfect square)."""
    assert n != 0
    sign = -1 if n < 0 else 1
    m = abs(n)
    d, f = 1, 2
    while f * f <= m:
        if m % f == 0:
            e = 0
            while m % f == 0:
                m //= f
                e += 1
            if e % 2 == 1:
                d *= f
        f += 1
    return sign * d * m          # leftover m is prime (exponent 1) or 1


def fundamental_discriminant(m):
    """Fundamental discriminant of Q(sqrt(m)), m a nonsquare integer."""
    d = squarefree_part(m)
    return d if d % 4 == 1 else 4 * d        # Python: (-7) % 4 == 1


def h_imag(D):
    """Class number of the imaginary quadratic field of fundamental
    discriminant D < 0, by counting reduced primitive forms (a,b,c):
    b^2 - 4ac = D, |b| <= a <= c, with b >= 0 whenever |b| == a or a == c."""
    assert D < 0 and D % 4 in (0, 1)
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a, a + 1):
            num = b * b - D
            if num % (4 * a) == 0:
                c = num // (4 * a)
                if c >= a:
                    if b < 0 and (abs(b) == a or a == c):
                        continue                      # avoid double count
                    if gcd(gcd(a, abs(b)), c) == 1:   # primitive (auto if D fund.)
                        h += 1
        a += 1
    return h


# sanity: classical class numbers (rock-solid values)
for D, hval in [(-3, 1), (-4, 1), (-7, 1), (-8, 1), (-11, 1), (-19, 1),
                (-43, 1), (-67, 1), (-163, 1), (-15, 2), (-20, 2),
                (-23, 3), (-24, 2), (-47, 5), (-71, 7)]:
    assert h_imag(D) == hval, (D, h_imag(D), hval)


def B2_chi(p):
    """Second generalized Bernoulli number B_{2,chi} for the quadratic
    character chi attached to the real quadratic field Q(sqrt(p)).
    Conductor f = p (p = 1 mod 4) or f = 4p (p = 3 mod 4).
    B_{2,chi} = f * sum_{a=1}^{f} chi(a) * B2(a/f),  B2(x) = x^2 - x + 1/6."""
    if p % 4 == 1:
        f = p
        chi = lambda a: int(jacobi_symbol(a, p))
    else:
        f = 4 * p

        def chi(a):                       # Kronecker(4p|a) = (p|a) for odd a
            if a % 2 == 0:
                return 0
            return int(jacobi_symbol(p, a))
    s = Fraction(0)
    for a in range(1, f + 1):
        c = chi(a)
        if c:
            x = Fraction(a, f)
            s += c * (x * x - x + Fraction(1, 6))
    return f * s


# sanity: zeta_{Q(sqrt 5)}(-1) = 1/30 = B_{2,chi}/24  ==>  B_{2,chi} = 4/5
assert B2_chi(5) == Fraction(4, 5), B2_chi(5)


# ------------------- IK94 Remark 3: the analytic formula -------------------

def tr_R_pi(p):
    """tr(R(pi)) = 2T - H for the principal genus, n = 2 (IK94 Remark 3,
    valid for primes p >= 7).  Transcribed verbatim from the paper."""
    assert isprime(p) and p >= 7
    leg2   = int(jacobi_symbol(2, p))      # (2/p)
    legm2  = int(jacobi_symbol(-2, p))     # (-2/p)
    p_on_3 = int(jacobi_symbol(p, 3))      # (p/3)
    h_p  = h_imag(fundamental_discriminant(-p))       # h(Q(sqrt(-p)))
    h_2p = h_imag(fundamental_discriminant(-2 * p))   # h(Q(sqrt(-2p)))
    h_3p = h_imag(fundamental_discriminant(-3 * p))   # h(Q(sqrt(-3p)))
    B2 = B2_chi(p)
    if p % 4 == 3:
        val = (Fraction(1, 96) * B2
               + Fraction(1, 8) * h_2p
               + Fraction(1, 12) * h_3p
               + (Fraction(p - 1, 48) * (9 - 4 * leg2)
                  + Fraction(1, 16) * (p - leg2)
                  + Fraction(1, 12) * (1 - p_on_3) * (3 - leg2)) * h_p)
    else:                                  # p = 1 mod 4 ; chi(2) = (2/p)
        val = (Fraction(1, 96) * (9 - 2 * leg2) * B2
               + Fraction(4 * p - 1, 48) * h_p
               + Fraction(1, 8) * h_2p
               + Fraction(1, 12) * (3 + legm2) * h_3p
               + Fraction(1, 12) * (1 - p_on_3) * h_p)
    assert val.denominator == 1, ("non-integral tr(R(pi))", p, val)
    return int(val)


# ----------------- (b) published table, IK94 p. 42 -------------------------
TABLE_H = {2: 1, 3: 1, 5: 2, 7: 2, 11: 5, 13: 4, 17: 8, 19: 10, 23: 16,
           29: 24, 31: 26, 37: 37, 41: 50, 43: 55, 47: 72, 53: 93}
TABLE_2TmH = {2: 1, 3: 1, 5: 2, 7: 2, 11: 5, 13: 4, 17: 8, 19: 8, 23: 14,
              29: 18, 31: 18, 37: 11, 41: 32, 43: 19, 47: 44, 53: 33}

# ----------------- (c) our Brandt-graph data (B1-B3 output) -----------------
PRIMES = [11, 13, 17, 19, 23, 29, 31, 37]
H2_GRAPH = {11: 5, 13: 4, 17: 8, 19: 10, 23: 16, 29: 24, 31: 26, 37: 37}
FIXED = {11: 5, 13: 4, 17: 8, 19: 8, 23: 14, 29: 18, 31: 18, 37: 11}
TWOCYC = {11: 0, 13: 0, 17: 0, 19: 1, 23: 1, 29: 3, 31: 4, 37: 13}
T1_PRED = {11: 5, 13: 4, 17: 8, 19: 9, 23: 15, 29: 21, 31: 22, 37: 24}

# secondary data (MEMORY section 3) for the closed-form T1 cross-check
GAMMA_PRIME = {11: 0, 13: 1, 17: 1, 19: 2, 23: 2, 29: 5, 31: 6, 37: 13}   # dim S3(Gamma'_0(p))
K_PARA      = {11: 0, 13: 1, 17: 1, 19: 1, 23: 1, 29: 2, 31: 2, 37: 4}    # dim S3(K(p))
S2_PLUS     = {11: 0, 13: 0, 17: 0, 19: 0, 23: 0, 29: 0, 31: 0, 37: 1}
S4_MINUS    = {11: 0, 13: 1, 17: 1, 19: 1, 23: 1, 29: 2, 31: 2, 37: 4}


# ------------------------------- checks -------------------------------------

print("=" * 78)
print("STEP 1: formula (Remark 3) vs published table, ALL primes 7 <= p <= 53")
print("=" * 78)
ok_formula = True
for p in sorted(TABLE_2TmH):
    if p < 7:
        continue                      # formula stated only for p >= 7
    t = tr_R_pi(p)
    tag = "OK " if t == TABLE_2TmH[p] else "!! MISMATCH"
    if t != TABLE_2TmH[p]:
        ok_formula = False
    print(f"  p={p:>2}:  formula tr(R(pi)) = {t:>3}   table 2T-H = "
          f"{TABLE_2TmH[p]:>3}   {tag}")
assert ok_formula, "Remark-3 formula does not reproduce the published table"
print("  --> formula transcription/implementation validated on 13 primes.\n")

print("=" * 78)
print("STEP 2: three-way comparison and the type number, p = 11..37")
print("=" * 78)
hdr = (f"  {'p':>2} | {'H(tab)':>6} {'h2(graph)':>9} | {'tr(formula)':>11} "
       f"{'2T-H(tab)':>9} {'fixed(B3)':>9} | {'T=(H+tr)/2':>10} "
       f"{'T1 pred':>7} {'fix+2cyc':>8}")
print(hdr)
print("  " + "-" * (len(hdr) - 2))
all_ok = True
for p in PRIMES:
    H = TABLE_H[p]
    t = tr_R_pi(p)
    # (i) class number: published H = our vertex count h_2
    assert H == H2_GRAPH[p], ("H vs h2", p, H, H2_GRAPH[p])
    # (ii) tr(R(pi)): formula = table = sigma-fixed vertices (Osaka Thm 4.3)
    assert t == TABLE_2TmH[p] == FIXED[p], ("trace triple", p)
    # (iii) involution bookkeeping: fixed + 2 * (#2-cycles) = h_2
    assert FIXED[p] + 2 * TWOCYC[p] == H2_GRAPH[p], ("orbit count", p)
    # (iv) type number: T = (H + tr)/2, integral, equals prediction
    T2 = Fraction(H + t, 2)
    assert T2.denominator == 1, ("T not integral", p)
    T = int(T2)
    assert T == T1_PRED[p], ("T vs prediction", p, T, T1_PRED[p])
    assert T == FIXED[p] + TWOCYC[p], ("T vs sigma-orbits", p)
    # (v) closed form  T1 = h2 - (Gamma'_0 - K) - S2plus * S4minus
    assert T == (H2_GRAPH[p] - (GAMMA_PRIME[p] - K_PARA[p])
                 - S2_PLUS[p] * S4_MINUS[p]), ("closed form T1", p)
    print(f"  {p:>2} | {H:>6} {H2_GRAPH[p]:>9} | {t:>11} {TABLE_2TmH[p]:>9} "
          f"{FIXED[p]:>9} | {T:>10} {T1_PRED[p]:>7} "
          f"{FIXED[p] + TWOCYC[p]:>8}")
print()
print("ALL ASSERTS PASSED: 8/8 match, three independent sources agree.")
print("  T(p) = 5, 4, 8, 9, 15, 21, 22, 24  ==  predicted T1(p).  CLOSED.")

# bonus: type numbers at the remaining table primes (for the record)
print()
print("For the record, T at the other table primes (>= 7):")
for p in [7, 41, 43, 47, 53]:
    T = (TABLE_H[p] + tr_R_pi(p)) // 2
    print(f"  p={p:>2}: H={TABLE_H[p]:>3}, tr(R(pi))={tr_R_pi(p):>3}, T={T:>3}")
