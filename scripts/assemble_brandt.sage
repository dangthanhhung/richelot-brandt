# ============================================================
#  ASSEMBLE THE BRANDT MATRIX  B_2(2)  for the eight primes 11 <= p <= 37
#  of the paper (and, with the same code, 41 <= p <= 53).
#  B_2(2) is the adjacency operator of the superspecial Richelot graph
#  G_p (= the (2,2)-isogeny graph; Jordan-Zaytman [JZ20], Def. 2 + Thm. 36).
#
#  WHY: diagonal = loop counts (ours at Jacobians + products),
#       trace and SPECTRUM are the objects matched against the
#       paramodular correspondence ([RW21], [DPRT21]) -> structured
#       explanation of the loop totals and the spectral catalogue.
#
#  HOW THE MISSING ROWS (product vertices) ARE BUILT -- four pieces:
#   (1) Jacobian rows: previous Richelot scan, UPGRADED so that edges
#       into the product locus record WHICH product (pair of j-invariants).
#       Split-codomain identification: when delta(F1,F2,F3)=0 the three
#       quadratics lie in a pencil; the pencil defines an involution of P^1
#       whose fixed points are the roots of the Jacobian covariant of two
#       members; a Moebius map sending those fixed points to {0,oo} makes
#       all three forms EVEN (checked by runtime assertion, instance by
#       instance), and then  y^2 = prod(a_i x^2 + b_i)  covers the two
#       elliptic factors   E1: y^2 = prod(a_i u + b_i),
#                          E2: w^2 = prod(b_i v + a_i)   (u=x^2, v=1/x^2).
#       j-invariants are twist-invariant, so scalings drop out.
#   (2) The 9 "product-type" kernels C1 x C2 of a product vertex map to
#       (E_a/C1) x (E_b/C2):  pure elliptic 2-isogenies.
#   (3) The 6 "graph-type" kernels going to Jacobians are recovered FREE
#       from the Mestre/Brandt symmetry  e_j*B_ij = e_i*B_ji
#       (Jordan-Zaytman Thm 18(c)) applied to the Jacobian rows of (1).
#   (4) Graph-type kernels with SPLIT codomain: occur (at these p) only at
#       square vertices E_a x E_a, where phi extends to an automorphism of
#       E_a and the quotient is the SAME vertex (a loop); expected count
#       |Aut E_a|/2  (=3 for j=0, 2 for j=1728, 1 otherwise).
#       This is NOT assumed silently: the script computes
#       6 - sum_J a(P->J) and ASSERTS it equals the expected count
#       (0 at mixed vertices).  Any mismatch aborts validity for that p.
#       (Cross-checked against the Katsura-Takashima count [KT20], C4.)
#
#  FIVE AUTOMATIC CHECKS per prime:
#   C1 row sums = 15            (JZ Thm 18(a), N_2(2) = 3*5)
#   C2 Mestre symmetry on the WHOLE matrix with weights e_i
#   C3 mass: products sum to M1^2/2, total = (p-1)(p^2+1)/5760
#      (Eichler; Hashimoto-Ibukiyama / Ekedahl, JZ Thm 18(d))
#   C4 #(edges J -> products) = Katsura-Takashima long-involution count [KT20]
#   C5 p=11 must reproduce VERBATIM the 5x5 matrix computed independently
#      by Jordan-Zaytman ([JZ20], Sec. 10.1)  -- calibration.
#
#  |Aut| facts used:  supersingular elliptic |Aut| = 6 (j=0), 4 (j=1728),
#  2 otherwise (Silverman III.10);  |Aut(E x E)| = 2*|Aut E|^2
#  (XOZ Lemma 3.2 gives 32 for E_1728, calibrated);  |Aut(E_a x E_b)| =
#  |Aut E_a|*|Aut E_b| for a != b.  Jacobian |Aut| derived from balance +
#  mass exactly as in the previous validated scan.
#
#  Run:   sage assemble_brandt.sage
# ============================================================
import itertools
from collections import Counter
from sage.all import (GF, PolynomialRing, Matrix, HyperellipticCurve, QQ, ZZ,
                      EllipticCurve, RDF, CDF, Newforms, Gamma0)

# Primes computed here.  To run others, edit this list.  The engine assumes the
# six Weierstrass points of every genus-2 vertex are rational over GF(p^2); this
# holds for all primes tested through p = 61 and first fails at p = 73 (the
# assert len(bp)==6 in neighbors() aborts), where some vertices need GF(p^4).
# The eight primes 11..37 are the Appendix-A matrices; p = 61 is the 128x128
# target.  Everything except branch()/pairing is prime-independent.
PRIMES = [11, 13, 17, 19, 23, 29, 31, 37, 61]

# ---------- small helpers ----------------------------------------------------

def jkey(j):
    """Canonical hashable key for an element of GF(p^2)."""
    pol = j.polynomial().list()
    pol = pol + [0]*(2-len(pol))
    return (int(pol[0]), int(pol[1]))

def aut_e(j, Fq):
    """|Aut| of a supersingular elliptic curve with this j (char p > 3)."""
    if j == Fq(0):
        return 6
    if j == Fq(1728):
        return 4
    return 2

def ss_j_list(Fq):
    """All supersingular j-invariants in GF(p^2)."""
    out = []
    for j in Fq:
        E = EllipticCurve(Fq, j=j)
        if E.is_supersingular(proof=False):
            out.append(j)
    return out

def two_isog_targets(j, Fq):
    """Multiset (list, length 3) of j-invariants of the three 2-isogeny
    codomains of a supersingular E with invariant j.  All 2-torsion x's lie
    in GF(p^2) for supersingular curves (Frobenius acts as +-p = +-1 on E[2])."""
    E = EllipticCurve(Fq, j=j).short_weierstrass_model()
    psi = E.division_polynomial(2)
    xs = [r for r, m in psi.roots() for _ in range(m)]
    assert len(xs) == 3, "2-torsion not rational -- unexpected for supersingular"
    out = []
    for x0 in xs:
        P = E([x0, 0])
        out.append(E.isogeny(P).codomain().j_invariant())
    return out

# ---------- split-codomain identification (delta = 0 case) -------------------

def split_factors(Fs, emb, Fq4, back):
    """Given the 3 (linearly dependent) quadratics of a degenerate Richelot
    splitting, return the sorted pair of j-invariant keys of the two elliptic
    factors of the split codomain, or (None, reason) on any degeneracy.
    emb  : the embedding GF(p^2) -> GF(p^4) (representation-independent);
    back : lookup {emb(j): jkey(j)} over the supersingular j's, used to map
           factor j-invariants back to GF(p^2) without touching internals.
    Every algebraic step is verified by an in-place assertion."""
    RT = PolynomialRing(Fq4, 'T')
    T = RT.gen()

    def bq(F):
        # binary quadratic (A,B,C) for A x^2 + B xz + C z^2 ;
        # a linear factor x-a stands for roots {a, oo}: (c1 x + c0 z) z
        cs = F.list()
        if F.degree() == 2:
            return (emb(cs[2]), emb(cs[1]), emb(cs[0]))
        return (Fq4(0), emb(cs[1]), emb(cs[0]))

    Qs = [bq(F) for F in Fs]

    # two independent pencil members
    pair = None
    for i, k in itertools.combinations(range(3), 2):
        A1, B1, C1 = Qs[i]; A2, B2, C2 = Qs[k]
        if (A1*B2-A2*B1, A1*C2-A2*C1, B1*C2-B2*C1) != (Fq4(0), Fq4(0), Fq4(0)):
            pair = (Qs[i], Qs[k]); break
    if pair is None:
        return None, "all three forms proportional"
    (A1, B1, C1), (A2, B2, C2) = pair

    # Jacobian covariant of two binary quadratics: its roots are the fixed
    # points of the involution defined by the pencil.
    jA = 2*(A1*B2 - A2*B1)
    jB = 4*(A1*C2 - A2*C1)
    jC = 2*(B1*C2 - B2*C1)

    if jA != 0:
        rts = [r for r, m in (jA*T**2 + jB*T + jC).roots() for _ in range(m)]
        if len(rts) != 2 or rts[0] == rts[1]:
            return None, "degenerate fixed-point form"
        r1, r2 = rts
        def tau(x):
            if x == 'inf':
                return Fq4(1)
            if x == r2:
                return 'inf'
            return (x - r1) / (x - r2)
    else:
        if jB == 0:
            return None, "fixed-point form vanished"
        r1 = -jC/jB
        def tau(x):
            if x == 'inf':
                return 'inf'
            return x - r1

    Cs = []
    for (A, B, C) in Qs:
        if A != 0:
            rs = [r for r, m in (A*T**2 + B*T + C).roots() for _ in range(m)]
            if len(rs) != 2:
                return None, "branch pair not rational over GF(p^4)"
            roots = rs
        else:
            roots = [-C/B, 'inf']
        s = [tau(x) for x in roots]
        if 'inf' in s:
            return None, "branch point hit a fixed point"
        if s[0] + s[1] != 0:
            return None, "EVENNESS CHECK FAILED"   # the in-place verification
        Cs.append(s[0]*s[1])

    e1 = Cs[0]+Cs[1]+Cs[2]
    e2 = Cs[0]*Cs[1]+Cs[0]*Cs[2]+Cs[1]*Cs[2]
    e3 = Cs[0]*Cs[1]*Cs[2]
    if e3 == 0:
        return None, "branch point hit the other fixed point"

    def j_of_cubic(a3, a2, a1, a0):
        # y^2 = a3 x^3 + a2 x^2 + a1 x + a0  ->  monic via X=a3*x, Y=a3*y
        if a3 == 0:
            return None
        try:
            E = EllipticCurve(Fq4, [0, a2, 0, a1*a3, a0*a3**2])
        except (ArithmeticError, ValueError):
            return None
        return E.j_invariant()

    j1 = j_of_cubic(Fq4(1), e1, e2, e3)   # E1: y^2 = prod(u + C_i)
    j2 = j_of_cubic(e3, e2, e1, Fq4(1))   # E2: w^2 = prod(C_i v + 1)
    if j1 is None or j2 is None:
        return None, "factor cubic degenerate"

    out = []
    for j4 in (j1, j2):
        k = back.get(j4)                  # in GF(p^2) AND supersingular?
        if k is None:
            return None, "factor j not a supersingular j of GF(p^2) -- inconsistent"
        out.append(k)
    return tuple(sorted(out)), None

# ---------- Jacobian side (upgraded previous scan) ----------------------------

def jacobian_side(p, Fq, ssj):
    R = PolynomialRing(Fq, 'x'); x = R.gen()
    # GF(p^4) plumbing, robust to how Sage represents the extension:
    Fq4 = Fq.extension(2, 'v')
    emb = Fq4.coerce_map_from(Fq)
    if emb is None:                          # fallback: explicit hom via a root
        rts = Fq.modulus().roots(Fq4)
        assert rts, "cannot embed GF(p^2) into GF(p^4)"
        emb = Fq.hom([rts[0][0]], Fq4)
    back = {emb(j): jkey(j) for j in ssj}    # inverse lookup on supersingular j's
    m_exp = (p-1)//2

    def is_ss(f):
        h = f**m_exp; cs = h.list()
        def c(k): return cs[k] if 0 <= k < len(cs) else Fq(0)
        return c(p-1) == 0 and c(p-2) == 0 and c(2*p-1) == 0 and c(2*p-2) == 0

    def igusa_vec(f):
        C = HyperellipticCurve(f); I2, I4, I6, I10 = C.igusa_clebsch_invariants()
        J2 = I2/8; J4 = (4*J2**2-I4)/96; J6 = (8*J2**3-160*J2*J4-I6)/576
        J8 = (J2*J6-J4**2)/4; J10 = I10/4096
        return (J2, J4, J6, J8, J10)

    weights = [2, 4, 6, 8, 10]
    def vkey(v):
        i = None
        for idx, c in enumerate(v):
            if c != 0:
                i = idx; break
        if i is None:
            return ('zero',)
        Ji = v[i]; out = [i]
        for j in range(5):
            if j == i:
                continue
            out.append(v[j]**weights[i] / Ji**weights[j])
        return tuple(out)

    def branch(f):
        rts = []
        for r, m in f.roots():
            rts += [r]*m
        if f.degree() == 5:
            rts.append('inf')
        return rts

    def quad(pr):
        a, b = pr
        if a == 'inf':
            return R([-b, 1])
        if b == 'inf':
            return R([-a, 1])
        return (x-a)*(x-b)

    def co3(G):
        cs = G.list(); return [cs[j] if j < len(cs) else Fq(0) for j in range(3)]

    def delta(F1, F2, F3):
        return Matrix(Fq, 3, 3, [co3(F1), co3(F2), co3(F3)]).determinant()

    def pair_parts(e):
        if not e:
            yield []; return
        f0 = e[0]
        for i in range(1, len(e)):
            for rem in pair_parts(e[1:i]+e[i+1:]):
                yield [(f0, e[i])]+rem

    def codomain(F1, F2, F3, d):
        def Gi(Fj, Fk): return (Fj.derivative()*Fk-Fk.derivative()*Fj)/d
        return Gi(F2, F3)*Gi(F3, F1)*Gi(F1, F2)

    anomalies = []

    def neighbors(f):
        bp = branch(f)
        # Requires all six Weierstrass points rational over GF(p^2).
        # Holds through p=61; fails first at p=73 (some vertices need GF(p^4)).
        # To go further, generalize branch()/pairing to carry points over GF(p^4).
        assert len(bp) == 6
        out = []
        for prt in pair_parts(bp):
            F1, F2, F3 = [quad(prt[k]) for k in range(3)]
            d = delta(F1, F2, F3)
            if d == 0:
                tgt, reason = split_factors([F1, F2, F3], emb, Fq4, back)
                if tgt is None:
                    anomalies.append(('split-id', reason))
                    out.append(('prod?', None))
                else:
                    out.append(('prod', tgt))
            else:
                Cp = codomain(F1, F2, F3, d)
                if Cp.degree() in (5, 6) and Cp.discriminant() != 0:
                    out.append(('jac', Cp))
                else:
                    anomalies.append(('nonzero-delta-degenerate', None))
                    out.append(('prod?', None))
        return out

    def find_start():
        for f in [x**6-1, x*(x**4-1), x**5-x, x**6-x]:
            if f.discriminant() != 0 and is_ss(f):
                return f
        for a in Fq:
            for b in Fq:
                f = x**6+a*x**3+b
                if f.discriminant() != 0 and is_ss(f):
                    return f

    s = find_start()
    curves = [s]; seen = {vkey(igusa_vec(s)): 0}; frontier = [0]
    while frontier:
        nf = []
        for i in frontier:
            for typ, Cp in neighbors(curves[i]):
                if typ == 'jac':
                    k = vkey(igusa_vec(Cp))
                    if k not in seen:
                        seen[k] = len(curves); curves.append(Cp); nf.append(seen[k])
        frontier = nf
    V = len(curves)

    adjJJ = [[0]*V for _ in range(V)]
    adjJP = [Counter() for _ in range(V)]
    for i, f in enumerate(curves):
        for typ, t in neighbors(f):
            if typ == 'jac':
                adjJJ[i][seen[vkey(igusa_vec(t))]] += 1
            elif typ == 'prod':
                adjJP[i][t] += 1
            else:
                adjJP[i]['UNKNOWN'] += 1

    # |Aut| via Mestre balance among Jacobians + mass pinning (validated before)
    w = [None]*V; w[0] = QQ(1)
    changed = True
    while changed:
        changed = False
        for i in range(V):
            if w[i] is None:
                continue
            for j in range(V):
                if adjJJ[i][j] > 0 and adjJJ[j][i] > 0 and w[j] is None:
                    w[j] = w[i]*QQ(int(adjJJ[i][j]))/QQ(int(adjJJ[j][i]))
                    changed = True
    assert all(wi is not None for wi in w), "Jacobian balance graph not connected"
    M2 = QQ((p-1)*(p**2+1))/QQ(5760); M1 = QQ(p-1)/QQ(24)
    jac_mass = M2 - M1**2/2                       # product mass is exactly M1^2/2
    scale = jac_mass / sum(w)
    eJ = [QQ(1)/(wi*scale) for wi in w]
    assert all(e.is_integer() for e in eJ), "non-integral |Aut| at a Jacobian"
    eJ = [ZZ(e) for e in eJ]

    return curves, seen, vkey, igusa_vec, adjJJ, adjJP, eJ, anomalies, R

# ---------- full assembly for one prime ---------------------------------------

KT_LONG_INVOLUTIONS = {2: 0, 4: 1, 8: 2, 10: 0, 12: 3, 24: 4, 48: 6}

def assemble(p, verbose=True):
    Fq = GF(p**2, 'w')
    ssj = ss_j_list(Fq)
    ss_keys = set(jkey(j) for j in ssj)
    eE = {jkey(j): aut_e(j, Fq) for j in ssj}

    curves, seen, vkey, igusa_vec, adjJJ, adjJP, eJ, anomalies, R = \
        jacobian_side(p, Fq, ssj)
    nJ = len(curves)

    # product vertices = unordered pairs of supersingular j's (with repeats)
    keys = sorted(ss_keys)
    prods = [tuple(sorted((a, b))) for idx, a in enumerate(keys)
             for b in keys[idx:]]
    nP = len(prods)
    pidx = {t: i for i, t in enumerate(prods)}
    eP = []
    for (a, b) in prods:
        eP.append(2*eE[a]**2 if a == b else eE[a]*eE[b])

    V = nP + nJ
    e = eP + eJ
    M = [[0]*V for _ in range(V)]

    ok = True
    msgs = []

    if anomalies:
        ok = False
        msgs.append("ANOMALIES in Jacobian scan: %s" % str(Counter(
            a for a, _ in anomalies)))

    # Jacobian rows
    for i in range(nJ):
        for j in range(nJ):
            M[nP+i][nP+j] = adjJJ[i][j]
        for t, c in adjJP[i].items():
            if t == 'UNKNOWN':
                ok = False; continue
            M[nP+i][pidx[t]] = c
        # C4: KT long involutions
        kt = KT_LONG_INVOLUTIONS.get(int(eJ[i]))
        if kt is None or sum(adjJP[i].values()) != kt:
            ok = False
            msgs.append("C4 FAIL at Jacobian %d: ->prod=%d, KT expects %s for |Aut|=%s"
                        % (i, sum(adjJP[i].values()), kt, eJ[i]))

    # product rows, piece (2): the 9 product-type kernels
    tg = {jkey(j): [jkey(t) for t in two_isog_targets(j, Fq)] for j in ssj}
    for P, (a, b) in enumerate(prods):
        for ta in tg[a]:
            for tb in tg[b]:
                M[P][pidx[tuple(sorted((ta, tb)))]] += 1

    # piece (3): product -> Jacobian via Mestre symmetry
    for P in range(nP):
        for i in range(nJ):
            val = QQ(e[P]*M[nP+i][P]) / QQ(eJ[i])
            if not val.is_integer():
                ok = False
                msgs.append("Mestre fill non-integral at P=%d, J=%d" % (P, i))
                val = QQ(0)
            M[P][nP+i] = ZZ(val)

    # piece (4): split graph-type kernels (loops at square vertices)
    for P, (a, b) in enumerate(prods):
        g = 6 - sum(M[P][nP+i] for i in range(nJ))
        expected = eE[a]//2 if a == b else 0
        if g != expected:
            ok = False
            msgs.append("GRAPH-SPLIT COUNT at P=%d (pair %s): got %d, expected %d"
                        " -- investigate (Kani-type split? see KT 2003.00633)"
                        % (P, str((a, b)), g, expected))
        else:
            M[P][P] += g

    # ----- checks C1-C3 -----
    for i in range(V):
        if sum(M[i]) != 15:
            ok = False; msgs.append("C1 FAIL row %d sums to %d" % (i, sum(M[i])))
    for i in range(V):
        for j in range(V):
            if e[j]*M[i][j] != e[i]*M[j][i]:
                ok = False
                msgs.append("C2 Mestre FAIL at (%d,%d)" % (i, j))
    mass_prod = sum(QQ(1)/QQ(x) for x in eP)
    mass_all = sum(QQ(1)/QQ(x) for x in e)
    M1 = QQ(p-1)/QQ(24); M2 = QQ((p-1)*(p**2+1))/QQ(5760)
    if mass_prod != M1**2/2:
        ok = False; msgs.append("C3 FAIL product mass %s != %s" % (mass_prod, M1**2/2))
    if mass_all != M2:
        ok = False; msgs.append("C3 FAIL total mass %s != %s" % (mass_all, M2))

    # ----- report -----
    Mz = Matrix(ZZ, V, V, M)
    labels = ["P{%s,%s}" % (a, b) for (a, b) in prods] + \
             ["J%d" % i for i in range(nJ)]
    if verbose:
        print("="*72)
        print("p = %d : %d product vertices + %d Jacobians = %d (h_2)" % (p, nP, nJ, V))
        print("vertex | e=|Aut| :",
              ", ".join("%s:%s" % (labels[i], e[i]) for i in range(V)))
        print("matrix rows (out-edges):")
        for i in range(V):
            print("   %-14s %s" % (labels[i], M[i]))
        loops = [M[i][i] for i in range(V)]
        print("loops:", loops, "  trace =", sum(loops),
              " (Jacobian part %d + product part %d)"
              % (sum(loops[nP:]), sum(loops[:nP])))
        cp = Mz.charpoly()
        print("charpoly factors over QQ:", cp.factor())
        # numeric spectrum via the weighted symmetrization (real by JZ 18(c,f))
        import math
        D = [math.sqrt(float(x)) for x in e]
        S = Matrix(RDF, V, V, [[float(M[i][j])*D[j]/D[i] for j in range(V)]
                               for i in range(V)])
        asym = max(abs(S[i][j]-S[j][i]) for i in range(V) for j in range(V))
        raw = S.change_ring(CDF).eigenvalues()
        max_im = max(abs(float(t.imag())) for t in raw)
        ev = sorted(float(t.real()) for t in raw)
        print("numeric spectrum:", [round(t, 6) for t in ev],
              " (sym defect %.1e, max|imag| %.1e)" % (asym, max_im))
        print("checks C1-C4:", "ALL PASS" if ok and not msgs else "PROBLEMS")
        for m in msgs:
            print("   !!", m)
    return Mz, labels, e, prods, nP, ok, (curves, seen, vkey, igusa_vec, R, Fq)

# ---------- C5: calibration certificate against Jordan-Zaytman, p = 11 --------

def jz_certificate():
    p = 11
    Mz, labels, e, prods, nP, ok, (curves, seen, vkey, igusa_vec, R, Fq) = \
        assemble(p, verbose=True)
    x = R.gen()
    # JZ Sec. 10.1 vertex order: E1xE1 (j=0), E2xE2 (j=1728), E1xE2, J(C1), J(C2)
    JZ = Matrix(ZZ, 5, 5, [[3,9,0,3,0],
                           [4,3,4,4,0],
                           [0,3,6,0,6],
                           [1,3,0,3,8],
                           [0,0,3,4,8]])
    k0 = jkey(Fq(0)); k1728 = jkey(Fq(1728))
    order = [ (k0, k0), (k1728, k1728), tuple(sorted((k0, k1728))) ]
    perm = [prods.index(t) for t in order]
    # map JZ's named curves to our BFS vertices via the Igusa key
    for f, name in [(x**6+1, "C1: y^2=x^6+1"), (x**6+3*x**3+1, "C2: y^2=x^6+3x^3+1")]:
        k = vkey(igusa_vec(f))
        assert k in seen, "JZ curve %s not among our vertices!" % name
        perm.append(nP + seen[k])
        print("identified %-22s = our %s" % (name, labels[nP+seen[k]]))
    Mperm = Matrix(ZZ, 5, 5, [[Mz[perm[i]][perm[j]] for j in range(5)]
                              for i in range(5)])
    same = (Mperm == JZ)
    print("C5 CALIBRATION vs Jordan-Zaytman 5x5 (their order):",
          "PASS -- matrices identical" if same else "FAIL")
    if not same:
        print("ours (permuted):"); print(Mperm)
        print("JZ:"); print(JZ)
    return same and ok

# ---------- elliptic newform data feeding the eigenvalue catalogue ------------

def newform_data(p):
    print("-"*72)
    print("elliptic newform data at level %d (for the eigenvalue catalogue):" % p)
    for k in (2, 4):
        try:
            nfs = Newforms(Gamma0(p), k, names='a')
            print("  S_%d(Gamma0(%d))^new : %d Galois orbit(s)" % (k, p, len(nfs)))
            for f in nfs:
                a2 = f[2]
                try:
                    deg = a2.parent().degree()
                except AttributeError:
                    deg = 1
                if deg == 1:
                    print("     a2 = %s   (rational orbit)" % a2)
                else:
                    print("     a2 minpoly = %s   (orbit degree %d)"
                          % (a2.minpoly(), deg))
        except Exception as ex:
            print("  S_%d computation failed: %s" % (k, ex))

# ---------- main ---------------------------------------------------------------

print("#"*72)
print("# assemble B_2(2): calibrate at p=11 against [JZ20], then 13 <= p <= 53")
print("#"*72)

cal = jz_certificate()
newform_data(11)
if not cal:
    print("\n*** CALIBRATION FAILED at p=11 -- aborting (matrix does not match [JZ20]).")
else:
    for p in PRIMES[1:]:
        assemble(p, verbose=True)
        newform_data(p)
    print("\nDone.")
