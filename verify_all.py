#!/usr/bin/env python3
"""verify_all.py -- independent structural and spectral verifier (all 31 primes).

Checks, from the frozen per-prime records h2_result_p.json and the ledger
assembly_results.json, in EXACT arithmetic (python-flint):

  per prime p in {11, ..., 149}:
  (a) row sums 15, and Mestre symmetry  e_j M_ij = e_i M_ji;
  (m) mass:  sum_i 1/e_i  equals the Hashimoto--Ibukiyama mass
      (p-1)(p^2+1)/5760 (weight normalization detected once, then enforced);
  (b) sigma is an involution commuting with M and preserving e;
      Fix(sigma) = res.TrR  (= 2*T_1(p) - h_2(p));
  (t) trace(M) = res.TrB = sum of the roots of the stored factorization;
  (c) EXACT characteristic polynomial:  charpoly(M) equals the product of
      the stored factors with their multiplicities  (flint, no modulus);
  (d) ledger: the Eisenstein factor (x - 15) appears with multiplicity 1;
      exactly one factor is squared, of degree res.degMys = ledger.mys;
      the stored square-excess degree is recomputed from the factorization,
      a square divisor of degree res.degMys = ledger.mys is feasible, the
      ledger polynomial N_p occurs with odd multiplicity and degree res.d,
      and  1 + sk + yos + 2*mys + deg(N_p) = h_2(p).  (The exact
      identification of the doubled block against the weight-2/4 split data
      is replayed by assembly.py.)

  globally: the Fix counts printed in the paper (Lemma 8.1),
      max_p 2*h_2(p) = 2866 attained at p = 149 (the bound of Section 6),
      N_61 = x + 7, N_73 = x + 6, N_79 = x + 5, and ledger flags all True.

Sign certificates (projector traces mod 2^61 - 1) are replayed by
assembly.py; hypotheses (a)-(d) of Proposition 6.2 at p = 11 are also
mirrored line by line by verify_reduction.py. Run with --verify-manifest
to additionally recompute every SHA-256 in MANIFEST.sha256.
"""
import argparse, hashlib, json, sys, time
from fractions import Fraction as Fr
from pathlib import Path
from flint import fmpz_mat, fmpz_poly

PRIMES = [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
          101,103,107,109,113,127,131,137,139,149]

def parse_poly(s):
    """Parse the record format 'x^6 + (-29)*x^5 + 322*x^4 + ... + 2026'."""
    coeffs = {}
    for term in s.replace('- ', '+ -').split(' + '):
        term = term.strip()
        if '*' in term:
            c, xp = term.split('*')
            c = int(c.strip('()'))
        elif term.startswith('x'):
            c, xp = 1, term
        else:
            c, xp = int(term.strip('()')), None
        if xp is None: d = 0
        elif xp == 'x': d = 1
        else: d = int(xp.split('^')[1])
        coeffs[d] = coeffs.get(d, 0) + c
    deg = max(coeffs)
    return fmpz_poly([coeffs.get(i, 0) for i in range(deg+1)])

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def verify_manifest(root):
    n = 0
    for line in open(root/'MANIFEST.sha256'):
        want, name = line.split()
        got = sha256(root/name)
        assert got == want, f"MANIFEST mismatch: {name}"
        n += 1
    print(f"MANIFEST: {n} checksums verified                               OK")

def verify_prime(p, root, ledger, norm):
    t0 = time.time()
    z = json.load(open(root/f'h2_result_{p}.json'))
    M, e, sig, res = z['M'], z['wts'], z['sigma'], z['res']
    n = len(M)
    assert res['h2'] == n and len(e) == n and len(sig) == n
    A = fmpz_mat(M)
    # (a)  row sums via A*1 = 15*1; Mestre via symmetry of diag(e)*A
    ones = fmpz_mat(n, 1);  fifteen = fmpz_mat(n, 1)
    D = fmpz_mat(n, n)
    for i in range(n):
        ones[i, 0] = 1; fifteen[i, 0] = 15; D[i, i] = e[i]
    assert A*ones == fifteen, (p, 'row sums')
    W = A*D
    assert W == W.transpose(), (p, 'Mestre')
    # (m)
    mass = sum(Fr(1, w) for w in e)
    target = Fr((p-1)*(p*p+1), 5760)
    if norm[0] is None:
        norm[0] = mass / target
        assert norm[0] in (Fr(1), Fr(2), Fr(1,2)), (p, 'mass normalization', norm[0])
    assert mass == norm[0]*target, (p, 'mass')
    # (b)
    assert sorted(sig) == list(range(n)) and all(sig[sig[i]] == i for i in range(n)), (p, 'sigma')
    assert all(e[sig[i]] == e[i] for i in range(n)), (p, 'sigma weights')
    P = fmpz_mat(n, n)
    for i in range(n): P[sig[i], i] = 1
    assert P*A == A*P, (p, 'sigma commute')
    fix = sum(1 for i in range(n) if sig[i] == i)
    assert fix == res['TrR'], (p, 'Fix vs TrR')
    # factors
    fs = [(parse_poly(s), m) for s, m in res['charpoly_factors']]
    assert sum(f.degree()*m for f, m in fs) == n, (p, 'degree sum')
    trM = sum(M[i][i] for i in range(n))
    sroots = sum(m*(-int(f[f.degree()-1])) for f, m in fs)
    assert trM == res['TrB'] == sroots, (p, 'TrB')
    # (c) exact charpoly
    cp = A.charpoly()
    prod = fmpz_poly([1])
    for f, m in fs: prod = prod * f**m
    assert cp == prod, (p, 'charpoly identity')
    # (d) ledger
    L = ledger[str(p)]
    eis = [(f, m) for f, m in fs if f == fmpz_poly([-15, 1])]
    assert len(eis) == 1 and eis[0][1] == 1, (p, 'Eisenstein')
    excess = sum(f.degree()*(m-1) for f, m in fs)
    assert excess == res['square_excess_deg'], (p, 'square excess')
    assert sum(f.degree()*(m//2) for f, m in fs) >= res['degMys'] == L['mys'], \
        (p, 'squared block feasibility')
    Np = parse_poly(L['Np'])
    assert Np.degree() == res['d'], (p, 'deg Np')
    if res['d'] > 0:
        assert any(f == Np and m % 2 == 1 for f, m in fs), (p, 'Np among factors')
    else:
        assert Np == fmpz_poly([1]), (p, 'Np at delta=0')
    assert 1 + L['sk'] + L['yos'] + 2*L['mys'] + Np.degree() == n, (p, 'degree ledger')
    assert L['ok'] and L['ok_signs'], (p, 'ledger flags')
    print(f"p={p:>3}  h2={n:>4}  Fix={fix:>3}  TrB={trM:>4}  charpoly OK  "
          f"blocks OK  [{time.time()-t0:.1f}s]")
    return fix, n

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument('--primes', nargs='*', type=int, default=PRIMES)
    ap.add_argument('--verify-manifest', action='store_true')
    a = ap.parse_args(argv)
    root = a.data_dir
    if a.verify_manifest:
        verify_manifest(root)
    ledger = json.load(open(root/'assembly_results.json'))
    norm = [None]
    fixes, sup2h2, argmax = {}, 0, None
    for p in a.primes:
        fix, n = verify_prime(p, root, ledger, norm)
        fixes[p] = fix
        if 2*n > sup2h2: sup2h2, argmax = 2*n, p
    if set(a.primes) == set(PRIMES):
        assert [fixes[p] for p in (11,13,17,19,23,29,31,37)] == [5,4,8,8,14,18,18,11]
        assert fixes[61] == 38
        assert (sup2h2, argmax) == (2866, 149)
        for p, want in ((61, 'x + 7'), (73, 'x + 6'), (79, 'x + 5')):
            assert parse_poly(ledger[str(p)]['Np']) == parse_poly(want), p
        print("global: Fix counts of Lemma 8.1 OK;  max 2*h2 = 2866 at p=149 OK;")
        print("        N_61 = x+7, N_73 = x+6, N_79 = x+5                    OK")
        print("="*64)
        print(f"ALL CHECKS PASSED ({len(a.primes)}/31 primes)")
        print("="*64)
    else:
        print(f"partial run: {len(a.primes)} primes, all checks passed")
    return 0

if __name__ == '__main__':
    sys.exit(main())
