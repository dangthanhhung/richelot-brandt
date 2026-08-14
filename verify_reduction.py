#!/usr/bin/env python3
# ============================================================================
# verify_reduction.py -- reference verifier for Proposition [prop:reduction]
# of "The Richelot Isogeny Graph as a Brandt Matrix".
#
# Checks, from finite data (M, R, e, block polynomials with signs), the four
# hypotheses of the proposition, in exact arithmetic:
#   (a) row sums 15 and Mestre symmetry  e_j M_ij = e_i M_ji;
#   (b) R^2 = I, RM = MR, e_{sigma(i)} = e_i;
#   (c) chi_M = (x-15) * P_SK * P_Y * Q^2 * N exactly in Z[x],
#       with Q monic and deg N = delta(p);
#   (d) for every irreducible f | chi_M:
#       Tr(R * h_f(M)) = s_f * deg f, certified modulo the Mersenne prime
#       q = 2^61 - 1, sufficient by Lemma [lem:onemod] since 2*h_2(p) < q.
#
# This file is the *reference* implementation (sympy + python ints), meant to
# mirror the printed Lemma/Proposition line by line; it is practical for
# p <= 61. The production replay at all 31 primes is the flint routine of
# assembly.py in the archived artifact (same mathematics, same modulus).
#
# Usage: python3 verify_reduction.py            # runs the embedded p=11 test
#        or import verify() and feed archived (M, R, e, blocks).
# ============================================================================
import sympy as sp

x = sp.symbols('x')
Q61 = (1 << 61) - 1                      # Mersenne prime 2^61 - 1


def verify(p, M, R, e, sk_blocks, y_blocks, Np, delta, verbose=True):
    """sk_blocks / y_blocks: lists of (poly_in_x, sign in {+1,-1}) per
    Galois orbit; Np: the general-type factor (sp.Integer(1) if delta=0)
    with predicted sign +1 for p < 167; returns True or raises."""
    n = len(M)
    E = [row[:] for row in M]

    # ---- (a) ----
    assert all(sum(row) == 15 for row in M), 'row sums'
    assert all(e[j]*M[i][j] == e[i]*M[j][i] for i in range(n)
               for j in range(n)), 'Mestre symmetry'

    # ---- (b) ----
    R2 = [[sum(R[i][k]*R[k][j] for k in range(n)) for j in range(n)]
          for i in range(n)]
    assert R2 == [[int(i == j) for j in range(n)] for i in range(n)], 'R^2=I'
    RM = [[sum(R[i][k]*M[k][j] for k in range(n)) for j in range(n)]
          for i in range(n)]
    MR = [[sum(M[i][k]*R[k][j] for k in range(n)) for j in range(n)]
          for i in range(n)]
    assert RM == MR, 'RM = MR'
    sigma = [next(j for j in range(n) if R[i][j] == 1) for i in range(n)]
    assert all(e[sigma[i]] == e[i] for i in range(n)), 'e sigma-invariant'

    # ---- (c): exact factorization over Z[x] ----
    chi = sp.Matrix(M).charpoly(x).as_expr()
    pred = sp.expand((x - 15)
                     * sp.prod([f for f, _ in sk_blocks] or [1])
                     * sp.prod([f for f, _ in y_blocks] or [1]))
    quo, rem = sp.div(chi, pred, x)
    assert sp.expand(rem) == 0, 'exact division by Eis*SK*Y'
    quoN, remN = sp.div(quo, Np, x)
    assert sp.expand(remN) == 0 and sp.degree(Np, x) == delta, 'N_p factor'
    # quoN must be an exact square Q^2
    fl = sp.factor_list(quoN)
    assert all(m % 2 == 0 for _, m in fl[1]), 'quotient is an exact square'
    Qpoly = sp.expand(fl[0]**sp.Rational(1, 2)
                      * sp.prod([f**(m//2) for f, m in fl[1]] or [1])) \
        if fl[1] or fl[0] == 1 else sp.Integer(1)
    assert sp.expand(Qpoly**2 - quoN) == 0, 'Q^2 reconstruction'

    # ---- predicted eigenvalue-sign multiset, per irreducible factor ----
    signed = {}                          # irreducible poly -> sum of eps

    def acc(poly, sign, mult=1):
        for f, m in sp.factor_list(sp.expand(poly))[1]:
            key = sp.Poly(f, x)
            signed[key] = signed.get(key, 0) + sign*m*mult*key.degree()
    acc(x - 15, +1)
    for f, s in sk_blocks: acc(f, s)
    for f, s in y_blocks: acc(f, s)
    # Q^2: each root twice with signs +1 and -1 -> signed contribution 0
    acc(Qpoly, +1); acc(Qpoly, -1)
    if delta:
        acc(Np, +1)                      # general-type sign +1 for p < 167

    # ---- (d): projector traces modulo Q61 (Lemma lem:onemod) ----
    assert 2*n < Q61
    dom = sp.GF(Q61)
    chiP = sp.Poly(chi, x)
    for f, m in sp.factor_list(chi)[1]:
        fP = sp.Poly(f, x)
        cof = sp.Poly(sp.div(chi, f**m, x)[0], x, domain=dom)
        fm = sp.Poly(f**m, x, domain=dom)
        inv = sp.invert(cof, fm)         # cof^{-1} mod f^m over GF(Q61)
        h = (sp.Poly(cof*inv, x, domain=dom) % sp.Poly(chiP, x, domain=dom))
        # evaluate h(M) mod Q61 by Horner, then trace of R*h(M)
        cs = [int(c) % Q61 for c in h.all_coeffs()]   # descending
        H = [[0]*n for _ in range(n)]
        for c in cs:
            H = [[(sum(H[i][k]*M[k][j] for k in range(n))
                   + (c if i == j else 0)) % Q61
                  for j in range(n)] for i in range(n)]
        tr = sum(R[i][k]*H[k][i] for i in range(n) for k in range(n)) % Q61
        target = signed.get(sp.Poly(f, x), 0) % Q61
        assert tr == target, ('projector trace', p, str(f), tr, target)
        if verbose:
            print(f'  (d) f = {f}:  Tr(R h_f(M)) = s_f deg f = '
                  f'{signed.get(sp.Poly(f, x), 0)}  OK')
    if verbose:
        print(f'p={p}: hypotheses (a)-(d) of Proposition [prop:reduction] '
              f'verified; conclusion follows.')
    return True


# ---------------------------------------------------------------------------
# Embedded self-test: p = 11 (matrix and weights printed in Appendix A;
# R(pi) = identity at p = 11 since all five classes are sigma-fixed).
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    M11 = [[3, 0, 9, 3, 0],
           [0, 6, 3, 0, 6],
           [4, 4, 3, 4, 0],
           [1, 0, 3, 3, 8],
           [0, 3, 0, 4, 8]]
    e11 = [72, 24, 32, 24, 12]
    R11 = [[int(i == j) for j in range(5)] for i in range(5)]
    verify(11, M11, R11, e11,
           sk_blocks=[(x**2 - 14*x + 46, +1)],
           y_blocks=[(x**2 + 6*x + 6, +1)],
           Np=sp.Integer(1), delta=0)
