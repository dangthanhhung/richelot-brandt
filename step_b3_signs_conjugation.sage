# ===========================================================================
#  step_b3_signs_conjugation.sage -- the Galois-conjugation involution and
#  the signed spectrum (the geometric side of the sign law, Theorem 5.3).
#
#  Runs on top of assemble_brandt.sage and produces two things:
#
#  (A) The Atkin-Lehner signs eps_p(f) of all newforms in S_2(Gamma0(p))
#      and S_4(Gamma0(p)).  These realize the Yoshida selection rule of the
#      catalogue: the pair (f,g) contributes an opposite-sign Yoshida factor
#      exactly when eps_p(f) eps_p(g) = -1 (e.g. the would-be eigenvalues
#      -5 at p=17, -3 at p=19, -3+-sqrt5 at p=23 are absent, as the rule
#      predicts).
#
#  (B) The Galois-conjugation involution sigma on G_p.  Frobenius x -> x^p
#      of GF(p^2)/GF(p) permutes the vertices and commutes with B_2(2)
#      (asserted at runtime); the script splits the module into the (+1)-
#      and (-1)-eigenspaces of sigma and factors charpoly(B_2(2)) on each.
#      This computes the odd part of the catalogue directly from the
#      geometry and exhibits the sign on each block (Theorem 5.3): for
#      instance the anti-invariant eigenvectors at p=19 (J5-J6) and p=23
#      (J3-J5) place the corresponding eigenvalues in the (-1)-space.
#
#  Run:   sage step_b3_signs_conjugation.sage
#  (it loads assemble_brandt.sage first, so that output prints first).
# ===========================================================================
load("assemble_brandt.sage")


print()
print("#"*72)
print("# (A) Atkin-Lehner eigenvalues eps_p of elliptic newforms")
print("#"*72)
for p in PRIMES:
    for k in (2, 4):
        try:
            nfs = Newforms(Gamma0(p), k, names='a')
        except Exception as ex:
            print("p=%d k=%d: Newforms failed: %s" % (p, k, ex)); continue
        for f in nfs:
            a2 = f[2]
            try:
                w = f.atkin_lehner_eigenvalue(p)
            except Exception:
                try:
                    # fallback: for prime level, a_p = -eps_p * p^(k/2 - 1)
                    ap = f[p]
                    w = -ap / p**((k-2)//2) if k == 4 else -ap
                except Exception:
                    w = '?'
            print("  p=%2d  S_%d  a2=%-22s  AL_p = %s" % (p, k, str(a2), w))

print()
print("#"*72)
print("# (B) spectrum split by the Galois-conjugation involution sigma")
print("#"*72)
for p in PRIMES:
    Mz, labels, e, prods, nP, ok, (curves, seen, vkey, igusa_vec, R, Fq) = \
        assemble(p, verbose=False)
    V = Mz.nrows()
    w = Fq.gen()

    def conj_jkey(k):
        j = Fq(k[0]) + Fq(k[1])*w
        return jkey(j**p)

    perm = []
    for (a, b) in prods:
        perm.append(prods.index(tuple(sorted((conj_jkey(a), conj_jkey(b))))))
    for f in curves:
        fc = R([c**p for c in f.list()])          # Frobenius-conjugate curve
        perm.append(nP + seen[vkey(igusa_vec(fc))])

    P = Matrix(ZZ, V, V, lambda i, j: 1 if perm[i] == j else 0)
    assert P*P == 1, "conjugation not an involution?!"
    assert all(e[perm[i]] == e[i] for i in range(V)), "|Aut| not conj-invariant?!"
    assert P*Mz == Mz*P, "B_2(2) does not commute with conjugation -- investigate"

    nfix = sum(1 for i in range(V) if perm[i] == i)
    print("p=%d: conjugation fixes %d of %d vertices; 2-cycles: %s"
          % (p, nfix, V, [(i, perm[i]) for i in range(V) if perm[i] > i]))

    MQ = Mz.change_ring(QQ)
    for s, name in [(1, 'even (+1)'), (-1, 'odd  (-1)')]:
        Kb = (P.change_ring(QQ) - s).right_kernel().basis_matrix().transpose()
        if Kb.ncols() == 0:
            print("   %s part: dimension 0" % name); continue
        A = (Kb.transpose()*Kb).inverse()*Kb.transpose()*MQ*Kb
        print("   %s part: dim %d, charpoly = %s"
              % (name, Kb.ncols(), A.charpoly().factor()))
print()
print("Done.")
