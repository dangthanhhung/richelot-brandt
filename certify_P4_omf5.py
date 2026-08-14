#!/usr/bin/env python3
# ============================================================================
# certify_P4_omf5.py -- Certificate (P4): external corroboration of the
# general-type factors N_p of Theorem 5.1 ("The Richelot Isogeny Graph as a
# Brandt Matrix") against the public paramodular database of
#
#   E. Assaf, W. Ladd, G. Rama, G. Tornaria, J. Voight,
#   "A database of paramodular forms from quinary orthogonal modular forms",
#   LuCaNT proceedings, Contemp. Math. 796 (2024); arXiv:2308.09824.
#   Data: https://github.com/assaferan/omf5_data
#   (pinned commit used for this certificate:
#    1907b8812c3ceb45fd09d47e24f96593116e816a)
#
# For each prime 11 <= p <= 149 the script reads the nonlift (type-G)
# eigensystems of weight (3,0) at level p from
# hecke_evs_3_0/data_nl_200/hecke_ev_3_0_nl_200_{p}.dat, reconstructs the
# eigenvalue a_2 of T(2,1) in the Hecke field (converting through the
# supplied Hecke-ring basis: numerators/denominators need NOT be the power
# basis -- e.g. at level 139), forms the characteristic polynomial of a_2
# over each Galois orbit, multiplies over orbits, and asserts equality with
# the general-type factor N_p of Theorem 5.1 (with N_p = 1, i.e. an empty
# nonlift file, at the seventeen primes with delta(p) = 0).
#
# The Atkin--Lehner strings of the database are REPORTED but not asserted:
# the sign convention of [ALRTV] relative to the paramodular Fricke sign
# used in the paper has not been pinned here.
#
# Usage:  python3 certify_P4_omf5.py [path/to/omf5_data/hecke_evs_3_0/data_nl_200]
# Requires: sympy.
# ============================================================================
import ast
import sys

import sympy as sp

x = sp.symbols('x')

# General-type factors N_p of Theorem 5.1, as coefficient lists
# [a_0, ..., a_n] with N_p(x) = a_0 + a_1 x + ... + a_n x^n.
PAPER_NP = {
    61: [7, 1], 73: [6, 1], 79: [5, 1], 89: [4, 1],
    97: [19, 9, 1], 101: [11, 7, 1], 103: [19, 9, 1],
    109: [29, 31, 10, 1], 113: [3, 1], 127: [19, 24, 9, 1],
    131: [7, 6, 1], 137: [2, 4, 1],
    139: [61, 110, 63, 14, 1], 149: [-1, 3, 10, 6, 1],
}
PRIMES = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
          73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149]


def a2_charpoly(form):
    """Characteristic polynomial of a_2 = lambda_{T(2,1)} acting on the
    Hecke field of one Galois orbit, from the [ALRTV] record."""
    fp = form['field_poly']                       # a_0 + ... + a_n x^n, monic
    d = len(fp) - 1
    K = sp.Matrix.zeros(d, d)                     # companion matrix
    for i in range(d - 1):
        K[i + 1, i] = 1
    for i in range(d):
        K[i, d - 1] = -sp.Rational(fp[i], fp[-1])
    nums = form['hecke_ring_numerators']
    dens = form['hecke_ring_denominators']
    gens = [sp.Matrix([sp.Rational(nij, dens[i]) for nij in nums[i]])
            for i in range(d)]                    # Hecke-ring basis, power coords
    a2 = sp.Matrix.zeros(d, 1)
    for i, ci in enumerate(form['lambda_p'][0]):  # first entry of lambda_p is p = 2
        a2 += sp.Integer(ci) * gens[i]
    A = sp.Matrix.zeros(d, d)                     # multiplication by a_2
    for j in range(d):
        A += a2[j] * (K ** j)
    return sp.expand(A.charpoly(x).as_expr())


def main(data_dir):
    ok = 0
    for p in PRIMES:
        forms = ast.literal_eval(
            open(f'{data_dir}/hecke_ev_3_0_nl_200_{p}.dat').read())
        gforms = [f for f in forms if f['aut_rep_type'] == 'G']
        np_db = sp.Integer(1)
        for f in gforms:
            np_db = sp.expand(np_db * a2_charpoly(f))
        cs = PAPER_NP.get(p)
        np_paper = (sp.expand(sum(c * x ** i for i, c in enumerate(cs)))
                    if cs else sp.Integer(1))
        assert sp.expand(np_db - np_paper) == 0, (p, np_db, np_paper)
        ok += 1
        al = ','.join(f['atkin_lehner_string'] for f in gforms) or '-'
        print(f'  p={p:>3}: #G={len(gforms)}  N_p = {sp.factor(np_db)}  '
              f'[AL(db): {al}]  OK')
    print(f'\n(P4) PASSED: N_p agrees with [ALRTV] at all {ok}/31 primes '
          f'11 <= p <= 149, including N_p = 1 at the 17 primes with delta = 0.')


if __name__ == '__main__':
    d = sys.argv[1] if len(sys.argv) > 1 else 'omf5_data/hecke_evs_3_0/data_nl_200'
    main(d)
