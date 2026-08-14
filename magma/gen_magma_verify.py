#!/usr/bin/env python3
# gen_magma_verify.py -- emit a self-contained Magma script that INDEPENDENTLY
# certifies our h2(p) result using Magma's Plesken-Souvignier isometry engine.
#
#   usage:  py gen_magma_verify.py <p>     (needs h2_result_<p>.json + h2_pure.py here)
#   output: verify_h2_<p>.m  -- paste into Magma (free calculator works; staged)
#
# What the .m file checks, with assert at every step:
#   STAGE 1  structure constants of our O-basis vs native quaternion arithmetic,
#            det(T4) = p^2, det(2H) = p^4, a positive control (an explicit
#            GL2(O) transform must be detected as isometric) and a negative
#            control (two classes with different |Aut| must not be).
#   STAGE 2  for every class: #AutomorphismGroup([2H] + O-linearity forms)
#            equals our |Aut|; Eichler mass closes exactly.
#   STAGE 3  classes are pairwise NON-isometric (theta-bucketed).
# O-linearity encoding: a Z-isometry U preserves G and F_c := G*Transpose(A_c)
# (A_c = right multiplication by basis e_c) iff U commutes with every right
# multiplication, i.e. U is left multiplication by an element of GL2(O).
import sys, os, json
from fractions import Fraction
import h2_pure

p = int(sys.argv[1])
R = json.load(open(os.path.join(os.getcwd(), f"h2_result_{p}.json")))
assert R["p"] == p
M = h2_pure.h2_of(p, inspect=True)
Q, O, den, gram8 = M["Q"], M["O"], M["den"], M["gram8"]
q = int(-Q.A); assert Q.B == -p
n_cls = len(R["classes"])

# ---- exact 4x4 linear algebra over Fractions (solve coords in O-basis) ----
def solve4(Bm, y):
    A = [row[:] + [yy] for row, yy in zip([list(r) for r in Bm], y)]
    for c in range(4):
        piv = next(r for r in range(c, 4) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        inv = Fraction(1) / A[c][c]
        A[c] = [v * inv for v in A[c]]
        for r in range(4):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [v - f * w for v, w in zip(A[r], A[c])]
    return [A[r][4] for r in range(4)]

Bm = [list(map(Fraction, e)) for e in O]          # rows: e_a over (1,i,j,k)
Bm_T = [[Bm[r][c] for r in range(4)] for c in range(4)]
def coords(x):                                    # x (4-tuple over 1,i,j,k)
    return solve4(Bm_T, list(map(Fraction, x)))

# LEFT-mult structure constants L_c[a][b]: e_c * e_a = sum_b L_c[a][b] e_b
# (automorphisms of h(v) = v g v^dagger are v |-> v*u, which commute with
#  LEFT scalar multiplication; the aux forms must therefore encode L, not R)
Ac, Rc = [], []
for c in range(4):
    rl, rr = [], []
    for a in range(4):
        col = coords(Q.mul(O[c], O[a]))          # LEFT:  e_c * e_a
        cor = coords(Q.mul(O[a], O[c]))          # RIGHT: e_a * e_c
        assert all(x.denominator == 1 for x in col + cor), (a, c)
        rl.append([int(x) for x in col]); rr.append([int(x) for x in cor])
    Ac.append(rl); Rc.append(rr)
# conjugation closure sanity (maximal orders are conj-stable)
for a in range(4):
    co = coords(Q.conj(O[a]))
    assert all(x.denominator == 1 for x in co), ("conj", a, co)
# T4 and dets
T4 = [[Q.trd(Q.mul(O[a], Q.conj(O[b]))) for b in range(4)] for a in range(4)]
assert all(v.denominator == 1 for row in T4 for v in row)
T4 = [[int(v) for v in row] for row in T4]
def detn(Min):
    A = [list(map(Fraction, r)) for r in Min]; n = len(A); d = Fraction(1)
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        if piv != c: A[c], A[piv] = A[piv], A[c]; d = -d
        d *= A[c][c]; inv = Fraction(1)/A[c][c]
        for r in range(c+1, n):
            f = A[r][c]*inv
            if f: A[r] = [v - f*w for v, w in zip(A[r], A[c])]
    return d
assert detn(T4) == p*p, detn(T4)
# cross-check our 2H block construction against gram8 on two classes
def twoH_blocks(s, t, rc):
    r = (Fraction(0),)*4
    for cco, e in zip(rc, O):
        r = Q.add(r, Q.smul(Fraction(cco), e))
    G = [[None]*8 for _ in range(8)]
    for a in range(4):
        for b in range(4):
            G[a][b]     = Fraction(s) * T4[a][b]
            G[4+a][4+b] = Fraction(t) * T4[a][b]
            v = Q.trd(Q.mul(Q.mul(O[a], r), Q.conj(O[b])))
            G[a][4+b] = Fraction(v); G[4+b][a] = Fraction(v)
    return G
for k in (0, n_cls // 2):
    s, t, rc, aut = R["classes"][k]
    G = twoH_blocks(s, t, rc)
    H = gram8(s, t, tuple(rc))
    assert all(G[i][j] == 2*H[i][j] for i in range(8) for j in range(8)), k
    assert all(v.denominator == 1 for row in G for v in row)
    if k == 0: assert detn(G) == p**4, detn(G)

# ---- convention proof (exact): u = [[1,0],[alpha,1]], alpha = e_2 --------
def mat8_mul(X, Y):
    return [[sum(X[i][k]*Y[k][j] for k in range(8)) for j in range(8)]
            for i in range(8)]
def mat8_T(X):
    return [[X[i][j] for i in range(8)] for j in range(8)]
def lift44(m):
    Z = [[Fraction(0)]*8 for _ in range(8)]
    for i in range(4):
        for j in range(4):
            Z[i][j] = Fraction(m[i][j]); Z[4+i][4+j] = Fraction(m[i][j])
    return Z
def forms_of(G, mats44):
    return [G] + [mat8_mul(G, mat8_T(lift44(m))) for m in mats44]
s1, t1, rc1, aut1 = R["classes"][0]
al = O[1]
r0 = (Fraction(0),)*4
for cco, e in zip(rc1, O):
    r0 = Q.add(r0, Q.smul(Fraction(cco), e))
r2q = Q.add(r0, Q.smul(Fraction(s1), Q.conj(al)))
t2 = Fraction(t1) + Q.trd(Q.mul(al, r0)) + Fraction(s1) * Q.nrd(al)
assert t2.denominator == 1 and t2 > 0
rc2 = coords(r2q)
assert all(x.denominator == 1 for x in rc2)
G1 = twoH_blocks(s1, t1, rc1)
G2 = twoH_blocks(s1, int(t2), [int(x) for x in rc2])
# witness matrix U8 for phi(v) = v*u : (x, y) -> (x + y*alpha, y); the
# coordinate matrix of y -> y*alpha is the RIGHT-mult matrix of alpha = Rc[1]
U8 = [[Fraction(1 if i == j else 0) for j in range(8)] for i in range(8)]
for i in range(4):
    for j in range(4):
        U8[4+i][j] = Fraction(Rc[1][i][j])
def transports(U, Fs1, Fs2):
    Ut = mat8_T(U)
    return all(mat8_mul(mat8_mul(U, F2), Ut) == F1
               for F1, F2 in zip(Fs1, Fs2))
# direction: h_{g'}(v) = h_g(v u)  =>  U8 (forms of g) U8^t = (forms of g')
okL = transports(U8, forms_of(G2, Ac), forms_of(G1, Ac))
okR = transports(U8, forms_of(G2, Rc), forms_of(G1, Rc))
assert okL, "LEFT-aux tuple does not transport -- convention still wrong!"
assert not okR, "RIGHT-aux also transports?! ambiguity -- investigate"
print("convention proof: LEFT-mult aux forms transport under the GL2(O) "
      "witness; RIGHT-mult do not  (as required)")

# --------------------------------------------------------------- emit .m ----
def frac_m(x):
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else f"{x.numerator}"
L = []
A = L.append
A(f"// verify_h2_{p}.m -- INDEPENDENT Magma certification of h2({p}) = {n_cls}")
A(f"// generated by gen_magma_verify.py from h2_result_{p}.json  (q = {q})")
A("// Stages can be toggled; on the free calculator run them one at a time")
A("// if you hit the CPU limit.  EVERY check is an assert: silence = failure")
A("// found nothing; the final lines print the verdict.")
A("STAGE1 := true;  STAGE2 := true;  STAGE3 := true;")
A("RLO := 1;  RHI := %d;  // STAGE2 class range (chunk if needed)" % n_cls)
A("")
A(f"p := {p};  EXPECTED := {n_cls};")
A(f"B<i,j,k> := QuaternionAlgebra< Rationals() | {-q}, {-p} >;")
A("e := [")
for a in range(4):
    A("  B![%s, %s, %s, %s]%s" % tuple([frac_m(c) for c in O[a]] +
                                       ["," if a < 3 else ""]))
A("];")
A("// LEFT-multiplication structure constants: e[c]*e[a] = &+ A[c][a][b]*e[b]")
A("Amats := [")
for c in range(4):
    flat = ", ".join(str(Ac[c][a][b]) for a in range(4) for b in range(4))
    A(f"  Matrix(Integers(), 4,4, [{flat}]){',' if c < 3 else ''}")
A("];")
A("classes := [")
for kk, (s, t, rc, aut) in enumerate(R["classes"]):
    A("  <%d, %d, [%d,%d,%d,%d], %d>%s" %
      (s, t, rc[0], rc[1], rc[2], rc[3], aut, "," if kk < n_cls-1 else ""))
A("];")
A("""
// ------------------------------------------------------------ machinery ----
T4 := Matrix(Rationals(), 4,4, [Trace(e[a]*Conjugate(e[b])) : b in [1..4],
                                                              a in [1..4]]);
assert forall{ x : x in Eltseq(T4) | IsIntegral(x) };
assert Determinant(T4) eq p^2;
Vbas := Matrix(Rationals(), 4,4, &cat[ Eltseq(e[a]) : a in [1..4] ]);
function OCoords(x)  // coordinates of x in our O-basis (must be integral)
  sol := Solution(Vbas, Vector(Rationals(), Eltseq(x)));
  assert forall{ c : c in Eltseq(sol) | IsIntegral(c) };
  return [ Integers()!c : c in Eltseq(sol) ];
end function;
function QuatOf(rc)
  return &+[ rc[a]*e[a] : a in [1..4] ];
end function;
function Gram2H(s, t, rc)  // integral 8x8 Gram, v*G*v^t = 2*h(v)
  r := QuatOf(rc);
  G := ZeroMatrix(Rationals(), 8, 8);
  for a in [1..4] do for b in [1..4] do
    G[a,b]     := s*T4[a,b];
    G[4+a,4+b] := t*T4[a,b];
    v := Trace(e[a]*r*Conjugate(e[b]));
    G[a,4+b] := v;  G[4+b,a] := v;
  end for; end for;
  assert forall{ x : x in Eltseq(G) | IsIntegral(x) };
  return Matrix(Integers(), 8, 8, [Integers()!x : x in Eltseq(G)]);
end function;
RA := [ DiagonalJoin(Amats[c], Amats[c]) : c in [1..4] ];
function FormsOf(G)  // PD form + O-linearity enforcers
  return [G] cat [ G*Transpose(RA[c]) : c in [1..4] ];
end function;
function AutOrder(s, t, rc)
  fs := FormsOf(Gram2H(s, t, rc));
  L := LatticeWithGram(fs[1]);
  return #AutomorphismGroup(L, [fs[n] : n in [2..5]]);
end function;
function SameClass(c1, c2)
  f1 := FormsOf(Gram2H(c1[1], c1[2], c1[3]));
  f2 := FormsOf(Gram2H(c2[1], c2[2], c2[3]));
  L1 := LatticeWithGram(f1[1]);  L2 := LatticeWithGram(f2[1]);
  ok, _ := IsIsometric(L1, [f1[n] : n in [2..5]], L2, [f2[n] : n in [2..5]]);
  return ok;
end function;

// --------------------------------------------------------------- STAGE 1 ---
if STAGE1 then
  for a in [1..4] do for c in [1..4] do
    assert e[c]*e[a] eq &+[ Amats[c][a,b]*e[b] : b in [1..4] ];
  end for; end for;
  print "STAGE1: structure constants OK";
  c1 := classes[1];
  G1 := Gram2H(c1[1], c1[2], c1[3]);
  assert Determinant(G1) eq p^4;
  print "STAGE1: det(T4) = p^2, det(2H) = p^4 OK";
  // positive control: explicit GL2(O) transform u = [[1, e2],[0,1]]
  al := e[2];
  s2 := c1[1];                       // u g u* for g = [[s, r],[rbar, t]]:
  r0 := QuatOf(c1[3]);
  r2 := r0 + c1[1]*Conjugate(al);    // [[s, r + s*albar],[..., ...]]
  t2 := c1[2] + Trace(al*r0) + c1[1]*Norm(al);
  assert IsIntegral(t2) and t2 gt 0;
  cT := < c1[1], Integers()!t2, OCoords(r2), c1[4] >;
  assert SameClass(c1, cT);
  print "STAGE1: positive control (unipotent GL2(O) transform) OK";
  i2 := 2;  while classes[i2][4] eq c1[4] do i2 +:= 1; end while;
  assert not SameClass(c1, classes[i2]);
  print "STAGE1: negative control (distinct classes) OK";
end if;

// --------------------------------------------------------------- STAGE 2 ---
if STAGE2 then
  mass := Rationals()!0;
  for idx in [RLO..RHI] do
    c := classes[idx];
    a := AutOrder(c[1], c[2], c[3]);
    assert a eq c[4];
    mass +:= 1/a;
  end for;
  printf "STAGE2: |Aut| verified for classes %o..%o\\n", RLO, RHI;
  if RLO eq 1 and RHI eq EXPECTED then
    assert mass eq (p-1)*(p^2+1)/5760;
    print "STAGE2: Eichler mass closes EXACTLY:", mass;
  else
    print "STAGE2 partial mass:", mass;
  end if;
end if;

// --------------------------------------------------------------- STAGE 3 ---
if STAGE3 then
  keys := [];
  for idx in [1..EXPECTED] do
    c := classes[idx];
    L := LatticeWithGram(Gram2H(c[1], c[2], c[3]));
    th := ThetaSeries(L, 12);
    Append(~keys, <[Coefficient(th, n) : n in [1..12]], c[4]>);
  end for;
  npairs := 0;
  for x in [1..EXPECTED] do
    for y in [x+1..EXPECTED] do
      if keys[x] eq keys[y] then
        npairs +:= 1;
        assert not SameClass(classes[x], classes[y]);
      end if;
    end for;
  end for;
  printf "STAGE3: pairwise distinct (theta-collisions tested: %o) OK\\n", npairs;
end if;
printf "VERDICT: h2(%o) = %o INDEPENDENTLY CERTIFIED BY MAGMA\\n", p, EXPECTED;
""")
out = os.path.join(os.getcwd(), f"verify_h2_{p}.m")
open(out, "w").write("\n".join(L))
print(f"wrote {out}  ({n_cls} classes, q={q}, den={den})")
