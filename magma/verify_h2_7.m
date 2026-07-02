// verify_h2_7.m -- INDEPENDENT Magma certification of h2(7) = 2
// self-contained certification; the embedded structure constants use q = 1
// Stages can be toggled; on the free calculator run them one at a time
// if you hit the CPU limit.  EVERY check is an assert: silence = failure
// found nothing; the final lines print the verdict.
STAGE1 := true;  STAGE2 := true;  STAGE3 := true;
RLO := 1;  RHI := 2;  // STAGE2 class range (chunk if needed)

p := 7;  EXPECTED := 2;
B<i,j,k> := QuaternionAlgebra< Rationals() | -1, -7 >;
e := [
  B![1/2, 0, 1/2, 0],
  B![0, 1/2, 0, 1/2],
  B![0, 0, 1, 0],
  B![0, 0, 0, 1]
];
// LEFT-multiplication structure constants: e[c]*e[a] = &+ A[c][a][b]*e[b]
Amats := [
  Matrix(Integers(), 4,4, [-3, 0, 2, 0, 0, 4, 0, -2, -7, 0, 4, 0, 0, 7, 0, -3]),
  Matrix(Integers(), 4,4, [0, -3, 0, 2, -4, 0, 2, 0, 0, -7, 0, 4, -7, 0, 3, 0]),
  Matrix(Integers(), 4,4, [-7, 0, 4, 0, 0, 7, 0, -4, -14, 0, 7, 0, 0, 14, 0, -7]),
  Matrix(Integers(), 4,4, [0, -7, 0, 4, -7, 0, 4, 0, 0, -14, 0, 7, -14, 0, 7, 0])
];
classes := [
  <1, 1, [0,0,0,0], 32>,
  <2, 2, [-2,-1,1,0], 48>
];

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
  printf "STAGE2: |Aut| verified for classes %o..%o\n", RLO, RHI;
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
  printf "STAGE3: pairwise distinct (theta-collisions tested: %o) OK\n", npairs;
end if;
printf "VERDICT: h2(%o) = %o INDEPENDENTLY CERTIFIED BY MAGMA\n", p, EXPECTED;
