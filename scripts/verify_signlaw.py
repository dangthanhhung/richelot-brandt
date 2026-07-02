#!/usr/bin/env python3
# ===========================================================================
# verify_signlaw.py -- complete certification of the sign law (Theorem 5.3),
# upgrading it from tabulated sign data to a fully asserted proof.
#
#   T0  import verify_three_pillars -> re-runs the full certification
#       (matrices, catalogue, Mys^2, cleanly-absent, trace).
#   T1  sigma re-derived INDEPENDENTLY of the original Sage run: exhaustive
#       search over all involutions of the vertex set that are
#       (P1) graph automorphisms  M[s(a)][s(b)] = M[a][b],
#       (P2) weight-preserving    e[s(a)] = e[a],
#       (P3) equal to j-conjugation on product vertices (forced by labels:
#            rational j fixed; the unique non-rational pair swapped).
#       The number of candidates is reported; the fixed-point data of
#       Sec. 5.1 selects sigma if needed.
#   T2  fixed = 5,4,8,8,14,18,18,11 and 2-cycles = 0,0,0,1,1,3,4,13
#       re-derived;  fixed + 2*cyc = h2;  Tr(P) = fixed re-proves the LHS
#       of the trace identity directly from the involution; degree bookkeeping
#       #2cyc = deg Mys + dim S2+ * dim S4-  (identity (IV)).
#   T3  SIGNED (ODD-PART) IDENTITY:  charpoly(B|V^-) =
#         prod_{(f,g): eps(f)=+1, eps(g)=-1} Res_y(t_f, m_g(x-y)) * Mys_p,
#       computed exactly in the basis {e_i - e_{sigma(i)}}, and equal to
#       the entries of Table tab:odd (hardcoded verbatim from the source).
#   T4  EVEN PART:  charpoly/odd = (x-15) * prod m_g(x-6)
#                                * prod_{eps(f)=-1 pairs} Res * Mys_p.
#   T5  SEPARATIONS (these place each eigenvector in the right sign space):
#       gcd(prod m_g(x-6), odd) = 1;  gcd(prod_{f^-} Res, odd) = 1;
#       gcd(Res of each f^+ pair, even) = 1;  Mys squarefree and coprime
#       to every lift factor (=> each irreducible mystery factor has
#       multiplicity exactly 1 in each of odd and even: the pair splits).
#   T6  EISENSTEIN: P*ones = ones (sign +1; unconditional for every p),
#       odd(15) != 0.
# ===========================================================================
import re
import sys, os

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_three_pillars as V3P  # noqa: E402

x, y = sp.symbols('x y')
PRIMES = [11, 13, 17, 19, 23, 29, 31, 37]
LABEL_RE = re.compile(r'(P\{[^}]*\}|J\d+):(\d+)')

# Sec. 5.1 reference data (fixed-point and 2-cycle counts), re-derived below.
FIXED_STORED = dict(zip(PRIMES, [5, 4, 8, 8, 14, 18, 18, 11]))
CYC_STORED = dict(zip(PRIMES, [0, 0, 0, 1, 1, 3, 4, 13]))

# Table tab:odd, verbatim from the source (trivial = 1 for p <= 17).
TAB_ODD = {
    11: "1", 13: "1", 17: "1",
    19: "x+2",
    23: "x",
    29: "x**3-5*x**2+3*x+5",
    31: "x**4+4*x**3-10*x**2-23*x-8",
    37: "(x**4+22*x**3+167*x**2+520*x+566)"
        "*(x**9-6*x**8-23*x**7+162*x**6+30*x**5-772*x**4+268*x**3"
        "+520*x**2-200*x-16)",
}


def parse_labels(evec_str):
    return [lab for (lab, _) in LABEL_RE.findall(evec_str)]


def product_components(label):
    """P{(a, b),(c, d)} -> [(a,b),(c,d)] as int tuples."""
    parts = re.findall(r'\(([^)]*)\)', label)
    return [tuple(int(t) for t in q.split(',')) for q in parts]


def monic_resultant(tf, mg):
    R = sp.expand(sp.resultant(tf, mg.subs(x, x - y), y))
    P_ = sp.Poly(R, x)
    if P_.LC() < 0:
        P_ = sp.Poly(-R, x)
    assert P_.LC() == 1
    return P_.as_expr()


def find_involutions(M, e, prod_map, jac_idx, cap=400):
    """All involutions s with (P1) M[s a][s b]=M[a][b], (P2) e-preserving,
    (P3) s|products = prod_map.  Returns list of dicts a->s(a)."""
    n = len(M)
    rows_ms = [tuple(sorted(M[a])) for a in range(n)]
    cols_ms = [tuple(sorted(M[b][a] for b in range(n))) for a in range(n)]

    def compatible(i, j, s):
        if e[i] != e[j] or M[i][i] != M[j][j] or M[i][j] != M[j][i]:
            return False
        if rows_ms[i] != rows_ms[j] or cols_ms[i] != cols_ms[j]:
            return False
        for a, sa in s.items():
            if M[j][sa] != M[i][a] or M[sa][j] != M[a][i]:
                return False
        return True

    sols = []

    def extend(s, remaining):
        if len(sols) >= cap:
            return
        if not remaining:
            full = dict(s)
            ok = all(M[full[a]][full[b]] == M[a][b]
                     for a in range(n) for b in range(n))
            if ok:
                sols.append(full)
            return
        i = remaining[0]
        # try fixing i
        if compatible(i, i, s):
            s2 = dict(s); s2[i] = i
            extend(s2, remaining[1:])
        # try pairing i with a later unassigned j
        for j in remaining[1:]:
            if compatible(i, j, s):
                s2 = dict(s); s2[i] = j; s2[j] = i
                rem2 = [r for r in remaining[1:] if r != j]
                extend(s2, rem2)

    base = dict(prod_map)
    extend(base, jac_idx)
    return sols


print("=" * 78)
print("SIGN LAW (THEOREM 5.3) -- full certification")
print("=" * 78)

for p in PRIMES:
    d = V3P.DATA[p]
    M = V3P.parse_matrix(d['rows'])
    e = V3P.parse_evec(d['evec'])
    labels = parse_labels(d['evec'])
    n = len(M)
    assert len(labels) == n == len(e)

    # ---- T1: forced action on product vertices ---------------------------
    prod_idx = [i for i, L in enumerate(labels) if L.startswith('P{')]
    jac_idx = [i for i, L in enumerate(labels) if L.startswith('J')]
    comps = {i: product_components(labels[i]) for i in prod_idx}
    nonrat = sorted({c for i in prod_idx for c in comps[i] if c[1] != 0})
    if p <= 31:
        assert nonrat == [], (p, "unexpected non-rational j below 37")
        conj = {}
    else:
        assert len(nonrat) == 2, (p, nonrat)   # forced pairing (order 2)
        conj = {nonrat[0]: nonrat[1], nonrat[1]: nonrat[0]}

    def conj_j(c):
        return conj.get(c, c)

    key = {frozenset_or_pair(comps[i]): i for i in prod_idx} if False else {}
    # build lookup by sorted component multiset
    lookup = {tuple(sorted(comps[i])): i for i in prod_idx}
    prod_map = {}
    for i in prod_idx:
        target = tuple(sorted(conj_j(c) for c in comps[i]))
        assert target in lookup, (p, labels[i], "conjugate vertex missing")
        prod_map[i] = lookup[target]
    assert all(prod_map[prod_map[i]] == i for i in prod_idx)

    # ---- T1: exhaustive involution search --------------------------------
    sols = find_involutions(M, e, prod_map, jac_idx)
    assert sols, (p, "no involution found")
    matching = [s for s in sols
                if sum(1 for a in range(n) if s[a] == a) == FIXED_STORED[p]]
    assert matching, (p, "no involution with the stored fixed count")
    sigma = matching[0]
    ident = {a: a for a in range(n)}
    for s2 in sols:
        if s2 != sigma:
            assert s2 == ident, (p, "unexpected extra automorphism")
    if p <= 17:
        assert sigma == ident, (p, "sigma should be trivial here")
    if p == 37:
        assert ident not in sols, (p, "identity should violate (P3) at 37")
    # if several candidates, they must all give the same signed split
    # (checked below after computing the odd part for sigma).

    # basic certified properties
    assert all(sigma[sigma[a]] == a for a in range(n))                 # P^2=I
    assert all(e[sigma[a]] == e[a] for a in range(n))                  # P2
    assert all(M[sigma[a]][sigma[b]] == M[a][b]
               for a in range(n) for b in range(n))                    # P1

    # ---- T2: counts, trace identity, degree bookkeeping ------------------
    fixed = sum(1 for a in range(n) if sigma[a] == a)
    cycles = [(a, sigma[a]) for a in range(n) if a < sigma[a]]
    c = len(cycles)
    assert fixed == FIXED_STORED[p] and c == CYC_STORED[p], (p, fixed, c)
    assert fixed + 2 * c == n
    # trace identity LHS from sigma itself (Theorem 5.4, left side):
    def deg_of(q):
        return sp.Poly(q, y).degree() if sp.sympify(q).has(y) else \
            sp.Poly(sp.sympify(q), y).degree()
    S2p = sum(sp.Poly(q, y).degree() for (q, s) in d['S2'] if s == +1)
    S2m = sum(sp.Poly(q, y).degree() for (q, s) in d['S2'] if s == -1)
    S4p = sum(sp.Poly(q, y).degree() for (q, s) in d['S4'] if s == +1)
    S4m = sum(sp.Poly(q, y).degree() for (q, s) in d['S4'] if s == -1)
    S4 = S4p + S4m
    assert fixed == 1 + S4 + S2m * S4p - S2p * S4m, (p, "trace identity")
    Mys = sp.expand(sp.sympify(d['Mys']))
    gmys = sp.Poly(Mys, x).degree() if Mys != 1 else 0
    assert c == gmys + S2p * S4m, (p, "identity (IV): #2cyc")

    # ---- T3: odd-part characteristic polynomial --------------------------
    if c == 0:
        odd = sp.Integer(1)
    else:
        ii = [a for (a, b) in cycles]
        jj = [b for (a, b) in cycles]
        N = sp.Matrix(c, c, lambda k, l: M[ii[k]][ii[l]] - M[ii[k]][jj[l]])
        odd = sp.expand(N.charpoly(x).as_expr())

    # predicted odd part from newform data
    odd_lift = sp.Integer(1)
    even_lift = sp.Integer(1)
    for (qf, sf) in d['S2']:
        tf = V3P.minpoly_of_double(qf)
        for (mg, sg) in d['S4']:
            if sf * sg == -1:
                R = monic_resultant(tf, sp.sympify(mg).subs(y, x))
                if sf == +1:
                    odd_lift = sp.expand(odd_lift * R)
                else:
                    even_lift = sp.expand(even_lift * R)
    odd_pred = sp.expand(Mys * odd_lift)
    assert sp.expand(odd - odd_pred) == 0, (p, "odd identity fails")
    assert sp.expand(odd - sp.expand(sp.sympify(TAB_ODD[p]))) == 0, \
        (p, "Table tab:odd mismatch")

    # all stored-count candidates give the same odd part
    for s2 in matching[1:]:
        cyc2 = [(a, s2[a]) for a in range(n) if a < s2[a]]
        ii2 = [a for (a, b) in cyc2]; jj2 = [b for (a, b) in cyc2]
        N2 = sp.Matrix(len(cyc2), len(cyc2),
                       lambda k, l: M[ii2[k]][ii2[l]] - M[ii2[k]][jj2[l]])
        odd2 = sp.expand(N2.charpoly(x).as_expr()) if cyc2 else sp.Integer(1)
        assert sp.expand(odd2 - odd) == 0, (p, "candidate split differs")

    # ---- T4: even part ----------------------------------------------------
    cp = sp.expand(sp.Matrix(M).charpoly(x).as_expr())
    even, rem = sp.div(sp.Poly(cp, x), sp.Poly(odd, x))
    assert rem.is_zero
    even = sp.expand(sp.Poly(even, x).as_expr())
    sk = sp.Integer(1)
    for (mg, sg) in d['S4']:
        sk = sp.expand(sk * sp.sympify(mg).subs(y, x - 6))
    even_pred = sp.expand((x - 15) * sk * even_lift * Mys)
    assert sp.expand(even - even_pred) == 0, (p, "even identity fails")

    # ---- T5: separations ---------------------------------------------------
    Podd, Peven = sp.Poly(odd, x), sp.Poly(even, x)
    assert sp.gcd(sp.Poly(sk, x), Podd).as_expr() == 1, (p, "SK leaks to odd")
    assert sp.gcd(sp.Poly(even_lift, x), Podd).as_expr() == 1, \
        (p, "f^- Yoshida leaks to odd")
    for (qf, sf) in d['S2']:
        if sf != +1:
            continue
        tf = V3P.minpoly_of_double(qf)
        for (mg, sg) in d['S4']:
            if sg == -1:
                R = monic_resultant(tf, sp.sympify(mg).subs(y, x))
                assert sp.gcd(sp.Poly(R, x), Peven).as_expr() == 1, \
                    (p, "f^+ Yoshida leaks to even")
    if Mys != 1:
        PM = sp.Poly(Mys, x)
        assert sp.gcd(PM, sp.Poly(sp.diff(Mys, x), x)).as_expr() == 1, \
            (p, "Mys not squarefree")
        lifts = sp.expand((x - 15) * sk * even_lift * odd_lift)
        assert sp.gcd(PM, sp.Poly(sp.expand(lifts), x)).as_expr() == 1, \
            (p, "Mys shares a root with a lift factor")
        # => each irreducible factor of Mys has multiplicity exactly 1
        #    in odd and exactly 1 in even: the mystery pair splits (+1,-1).

    # ---- T6: Eisenstein ----------------------------------------------------
    assert all(sum(1 for b in range(n) if sigma[b] == a) == 1
               for a in range(n))               # permutation => P*ones = ones
    assert sp.Poly(odd, x).eval(15) != 0, (p, "15 in odd part")

    print(f"  p={p:2d}: involutions(P1-P3)={len(sols):3d} "
          f"(fixed-count match={len(matching)})  fixed={fixed:2d} 2cyc={c:2d} "
          f" odd==pred==tab:odd OK  even OK  separations OK  trace OK")

print("=" * 78)
print("SIGN LAW (THEOREM 5.3) FULLY CERTIFIED: sigma re-derived, counts re-derived,")
print("odd/even identities, all sign separations, Eisenstein clause,")
print("and the trace-identity left-hand side re-proved from the involution itself.")
print("=" * 78)
