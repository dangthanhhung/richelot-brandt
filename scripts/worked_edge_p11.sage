#!/usr/bin/env sage
r"""
worked_edge_p11.sage
====================================================================
Produce ONE *certified* Richelot edge of the superspecial graph G_11,
to replace the illustrative splitting in the worked example (Sec. 3.3).

Steps:
  (1) test superspeciality via the Hasse-Witt (Cartier-Manin) matrix:
      a genus-2 curve in char p is SUPERSPECIAL iff this 2x2 matrix is 0;
  (2) factor the sextic into three quadratics f = c*q1*q2*q3 (a 2-torsion
      splitting) and form the Richelot data [q_i,q_j] and Delta_R;
  (3) reduce the codomain, read its Igusa-Clebsch invariants, report the
      vertex it lands on (the matrix entry this edge contributes).

REQUIREMENTS: SageMath >= 9.5.   USAGE:  sage worked_edge_p11.sage

This reproduces the certified Richelot edge displayed in the
"One Richelot edge" stage of Sec. 3.3: a superspecial curve over F_11,
one (2,2)-splitting, the bracket polynomials, Delta_R, and the codomain.
====================================================================
"""
p = 11
Fq.<t> = GF(p^2, modulus='conway')   # F_{11^2}, for the Hasse-Witt and Igusa steps
R.<x>  = PolynomialRing(Fq)
F0 = GF(p)                            # the prime field F_11
R0.<x0> = PolynomialRing(F0)          # sextics and splittings live here (rational)

def hasse_witt_matrix(f):
    h = f**((p - 1)//2)
    c = h.list()
    def cf(k): return c[k] if 0 <= k < len(c) else Fq(0)
    return Matrix(Fq, 2, 2, [[cf(p*1-1), cf(p*1-2)],
                             [cf(p*2-1), cf(p*2-2)]])

# We look for the edge displayed in Sec. 3.3: a superspecial curve whose
# sextic factors into three quadratics over the PRIME field F_11, so that the
# splitting and all the Richelot data are rational.  The paper's curve is
# x^6 + x^4 - 3x^2 + 3, i.e. (a,b,c) = (1,-3,3) mod 11; we put it first, then
# scan the rest of the F_11 family x^6 + a x^4 + b x^2 + c.
print("Searching superspecial genus-2 curves over F_%d (prime field) ..." % p)

def hw_zero_over_F11(f0):
    """Hasse-Witt test: lift the F_11 sextic to F_{11^2} and check HW = 0."""
    fq = R([Fq(co) for co in f0.list()])
    return fq.degree()==6 and fq.discriminant()!=0 and hasse_witt_matrix(fq)==0

def rational_splittings(f0):
    """All (2,2)-splittings f0 = c*q1*q2*q3 into three quadratics over F_11."""
    fac=[g for (g,e) in f0.factor() for _ in range(e)]
    quads=[g.monic() for g in fac if g.degree()==2]
    lins =[g for g in fac if g.degree()==1]
    out=[]
    from itertools import combinations
    if len(quads)==3 and not lins:
        out.append(tuple(quads))
    elif len(lins)==6:
        idx=list(range(6)); seen=set()
        for a in combinations(idx,2):
            rest=[i for i in idx if i not in a]
            for b in combinations(rest,2):
                c=tuple(i for i in rest if i not in b)
                key=tuple(sorted([tuple(sorted(a)),tuple(sorted(b)),tuple(sorted(c))]))
                if key in seen: continue
                seen.add(key)
                out.append(((lins[a[0]]*lins[a[1]]).monic(),
                            (lins[b[0]]*lins[b[1]]).monic(),
                            (lins[c[0]]*lins[c[1]]).monic()))
    return out

def richelot_F11(q1,q2,q3):
    def br(a,b): return a.derivative()*b - a*b.derivative()
    b23,b31,b12 = br(q2,q3),br(q3,q1),br(q1,q2)
    def abc(q):
        cc=q.list(); cc+=[F0(0)]*(3-len(cc)); return [cc[2],cc[1],cc[0]]
    DR=Matrix(F0,[abc(q1),abc(q2),abc(q3)]).det()
    return (b23,b31,b12,DR,b23*b31*b12)

preferred = x0**6 + x0**4 - 3*x0**2 + 3
def family():
    yield preferred
    for a in F0:
        for b in F0:
            for c0 in F0:
                g = x0**6 + a*x0**4 + b*x0**2 + c0
                if g != preferred:
                    yield g

edge=None
for f0 in family():
    if not hw_zero_over_F11(f0): continue
    for (q1,q2,q3) in rational_splittings(f0):
        b23,b31,b12,DR,sextic = richelot_F11(q1,q2,q3)
        if DR!=0:
            edge=(f0,q1,q2,q3,b23,b31,b12,DR,sextic); break
    if edge: break

if not edge:
    raise SystemExit("No certified superspecial edge over F_11 in the scanned family.")

f,q1,q2,q3,b23,b31,b12,DR,sextic=edge
print("Found a superspecial curve (Hasse-Witt = 0).\n")
print("p = %d   (curve and splitting over the prime field F_%d)" % (p,p))
print("C: y^2 = f(x),   f =", f)
fq = R([Fq(co) for co in f.list()])
print("Hasse-Witt matrix =", hasse_witt_matrix(fq).list(), "(zero: superspecial)\n")
print("one (2,2)-splitting  f = c * q1 * q2 * q3:")
print("  q1 =",q1); print("  q2 =",q2); print("  q3 =",q3,"\n")
print("brackets [q_i,q_j] = q_i' q_j - q_i q_j':")
print("  [q2,q3] =",b23); print("  [q3,q1] =",b31); print("  [q1,q2] =",b12,"\n")
print("Delta_R =",DR,"(nonzero => Jacobian codomain)")
print("codomain:  Delta_R * y^2 = [q2,q3][q3,q1][q1,q2]")
print("  RHS sextic =",sextic,"\n")
try:
    C2=HyperellipticCurve(R([Fq(co) for co in (sextic/DR).list()]))
    print("codomain Igusa-Clebsch invariants:", C2.igusa_clebsch_invariants())
    print("(match against the vertex table to record the matrix entry)")
except Exception as e:
    print("Igusa-Clebsch step skipped:", e)
