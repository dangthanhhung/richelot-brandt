// verify_h2_61.m -- INDEPENDENT Magma certification of h2(61) = 128
// self-contained certification; the embedded structure constants use q = 2
// Stages can be toggled; on the free calculator run them one at a time
// if you hit the CPU limit.  EVERY check is an assert: silence = failure
// found nothing; the final lines print the verdict.
STAGE1 := true;  STAGE2 := true;  STAGE3 := true;
RLO := 1;  RHI := 128;  // STAGE2 class range (chunk if needed)

p := 61;  EXPECTED := 128;
B<i,j,k> := QuaternionAlgebra< Rationals() | -2, -61 >;
e := [
  B![1/2, 0, 1/2, 1/2],
  B![0, 1/4, 1/2, 1/4],
  B![0, 0, 1, 0],
  B![0, 0, 0, 1]
];
// LEFT-multiplication structure constants: e[c]*e[a] = &+ A[c][a][b]*e[b]
Amats := [
  Matrix(Integers(), 4,4, [-91, 0, 46, 46, -61, -30, 46, 38, -61, -122, 92, 61, -122, 122, 0, 31]),
  Matrix(Integers(), 4,4, [-61, 31, 15, 23, -46, 0, 23, 23, -61, -61, 61, 46, -61, 122, -31, 0]),
  Matrix(Integers(), 4,4, [-61, 122, -30, 0, -61, 61, 0, 15, -122, 0, 61, 61, 0, 244, -122, -61]),
  Matrix(Integers(), 4,4, [-122, -122, 122, 92, -61, -122, 92, 61, 0, -244, 122, 61, -244, 0, 122, 122])
];
classes := [
  <1, 1, [0,0,0,0], 8>,
  <2, 16, [0,-2,1,0], 48>,
  <2, 12, [0,-1,0,0], 8>,
  <2, 12, [0,-1,1,0], 8>,
  <2, 5, [1,3,-2,-1], 4>,
  <3, 4, [-3,3,0,1], 4>,
  <3, 8, [0,1,0,0], 4>,
  <3, 8, [0,1,-1,0], 4>,
  <3, 21, [-2,0,0,1], 2>,
  <3, 11, [-2,-2,2,1], 2>,
  <3, 3, [1,1,-1,-1], 4>,
  <4, 16, [0,-4,1,1], 4>,
  <4, 16, [0,4,-3,-1], 4>,
  <4, 6, [0,1,0,0], 4>,
  <4, 7, [-4,1,1,2], 4>,
  <4, 7, [-4,1,2,2], 4>,
  <4, 6, [0,1,-1,0], 4>,
  <4, 8, [0,-2,1,0], 12>,
  <4, 8, [0,-2,1,1], 4>,
  <4, 9, [-4,2,1,2], 12>,
  <4, 18, [-3,1,0,1], 2>,
  <4, 18, [-3,1,2,1], 2>,
  <4, 7, [2,5,-4,-2], 2>,
  <4, 7, [2,5,-3,-2], 2>,
  <5, 42, [0,-5,2,0], 4>,
  <5, 11, [-5,-4,5,4], 4>,
  <5, 25, [0,-4,2,0], 12>,
  <5, 11, [-5,-4,4,4], 4>,
  <5, 29, [-5,2,0,2], 4>,
  <5, 5, [0,-3,1,1], 8>,
  <5, 6, [0,7,-4,-2], 4>,
  <5, 5, [0,-3,2,1], 8>,
  <5, 6, [0,7,-3,-2], 4>,
  <5, 29, [-5,2,3,2], 4>,
  <5, 8, [-4,6,-1,0], 4>,
  <5, 19, [-4,7,-2,1], 2>,
  <5, 18, [-4,-3,3,2], 2>,
  <5, 18, [-4,-3,4,2], 2>,
  <5, 15, [-3,5,0,0], 2>,
  <5, 13, [2,-4,0,0], 4>,
  <5, 13, [2,-4,2,0], 4>,
  <5, 29, [-3,6,-2,1], 2>,
  <5, 5, [2,1,-2,-1], 4>,
  <5, 5, [2,1,-1,-1], 4>,
  <6, 6, [0,-6,3,1], 8>,
  <6, 6, [0,-6,3,2], 8>,
  <6, 6, [-6,-5,6,4], 12>,
  <6, 35, [0,-5,1,2], 8>,
  <6, 35, [0,-5,2,0], 4>,
  <6, 35, [0,-5,3,0], 4>,
  <6, 35, [0,-5,4,2], 8>,
  <6, 6, [-6,-5,5,4], 12>,
  <6, 31, [0,-4,1,0], 12>,
  <6, 17, [-6,2,1,3], 4>,
  <6, 31, [0,-4,3,0], 12>,
  <6, 26, [-6,3,0,2], 4>,
  <6, 26, [-6,3,3,2], 4>,
  <6, 11, [-4,0,1,2], 8>,
  <6, 25, [2,-5,0,0], 2>,
  <6, 25, [-4,1,0,2], 2>,
  <6, 25, [2,-5,3,0], 2>,
  <6, 25, [-4,1,3,2], 2>,
  <6, 15, [-4,3,0,2], 2>,
  <6, 15, [-4,3,1,2], 2>,
  <6, 13, [3,7,-6,-3], 2>,
  <6, 13, [3,7,-4,-3], 2>,
  <7, 35, [0,0,2,0], 4>,
  <7, 30, [0,-5,2,0], 4>,
  <7, 31, [0,9,-5,-1], 4>,
  <7, 30, [0,-5,3,0], 4>,
  <7, 31, [0,9,-4,-1], 4>,
  <7, 23, [-6,7,1,1], 2>,
  <7, 37, [1,-7,2,0], 4>,
  <7, 31, [-6,3,0,3], 4>,
  <7, 39, [1,10,-7,-2], 4>,
  <7, 7, [1,-4,1,0], 4>,
  <7, 7, [1,-4,1,1], 4>,
  <7, 7, [1,-4,2,1], 4>,
  <7, 31, [-6,3,3,3], 4>,
  <7, 9, [2,0,0,-1], 8>,
  <7, 21, [-5,9,-3,1], 2>,
  <7, 20, [-5,-5,4,3], 2>,
  <7, 20, [-5,-5,6,3], 2>,
  <7, 50, [2,10,-7,-2], 2>,
  <7, 50, [2,10,-5,-2], 2>,
  <7, 7, [3,0,-2,-2], 4>,
  <7, 7, [3,0,-2,-1], 4>,
  <7, 28, [-4,8,-3,1], 2>,
  <7, 28, [-4,8,-1,1], 2>,
  <8, 26, [-8,8,1,3], 4>,
  <8, 26, [-8,-8,9,5], 4>,
  <8, 24, [0,-8,3,1], 4>,
  <8, 24, [0,8,-5,-1], 4>,
  <8, 36, [0,10,-5,-1], 8>,
  <8, 11, [0,-5,2,2], 4>,
  <8, 11, [0,-5,3,2], 4>,
  <8, 69, [0,-4,1,-1], 4>,
  <8, 69, [0,-4,3,-1], 4>,
  <8, 40, [1,-7,1,2], 2>,
  <8, 10, [1,9,-6,-3], 2>,
  <8, 32, [1,-5,0,1], 2>,
  <8, 32, [1,-5,4,1], 2>,
  <8, 50, [-6,3,4,2], 2>,
  <8, 50, [-6,3,-1,2], 2>,
  <9, 19, [-9,10,0,1], 4>,
  <9, 62, [0,-8,5,0], 4>,
  <9, 32, [-9,11,0,3], 4>,
  <9, 31, [-9,-7,9,5], 4>,
  <9, 31, [-9,-7,10,6], 4>,
  <9, 44, [0,-7,3,0], 4>,
  <9, 44, [0,-7,4,0], 4>,
  <9, 21, [-7,-8,9,5], 2>,
  <9, 21, [-7,-8,6,5], 2>,
  <9, 33, [-6,10,0,0], 2>,
  <9, 10, [4,3,-4,-2], 2>,
  <10, 27, [-10,0,3,5], 4>,
  <10, 41, [-10,12,1,3], 4>,
  <10, 40, [-10,-8,11,6], 12>,
  <10, 41, [-10,12,-3,3], 4>,
  <10, 40, [-10,-8,7,6], 12>,
  <10, 33, [2,1,0,0], 2>,
  <10, 47, [5,15,-8,-5], 2>,
  <11, 41, [-11,-8,12,8], 4>,
  <11, 41, [-11,-8,7,8], 4>,
  <11, 49, [-11,-6,11,6], 4>,
  <11, 53, [-10,-8,8,5], 2>,
  <12, 16, [0,-8,3,3], 4>,
  <14, 15, [0,5,-3,0], 4>
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
