#!/usr/bin/env python3
"""verify_positivity_2500.py -- finite input of Theorem 8.5 (positivity of the bias).

Evaluates, in exact rational arithmetic and with no dependencies beyond the
standard library, the closed class-number formula of Proposition 8.2 for the
signed defect d(p), at every prime 7 <= p <= 2500, and checks:

  (1) d(p) is an integer and d(p) >= 0                      [the finite input
      quoted in the proof of Theorem 8.5];
  (2) the crude lower bound obtained by dropping the (positive) class-number
      terms and using B_{2,chi} >= D^{3/2}/15 is positive for every prime
      p > 673, and 673 is the largest prime where it fails  [the cutoff used
      in the proof];
  (3) regression pins: d(61) = d(73) = d(79) = 1 (the certified range, where
      d(p) = delta(p)).

Conventions. h(sqrt(-m)) is the class number of the field Q(sqrt(-m)),
computed as the number of reduced primitive binary quadratic forms of the
field discriminant (-m if -m = 1 mod 4, else -4m). B_{2,chi} is the
generalized Bernoulli number for the quadratic character chi of conductor
D0 = disc(Q(sqrt(p))) (= p for p = 1 mod 4, = 4p otherwise), via
B_{2,chi} = D0 * sum_{a=1}^{D0} chi(a) B_2(a/D0), B_2(x) = x^2 - x + 1/6.
"""
from fractions import Fraction as Fr
import math, sys

PMAX = 2500

def sieve(n):
    s = bytearray([1])*(n+1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = b'\x00'*len(s[i*i::i])
    return [i for i in range(2, n+1) if s[i]]

def kronecker(a, n):
    if n == 0: return 1 if a in (1, -1) else 0
    if n < 0:
        return kronecker(a, -n) * (-1 if a < 0 else 1)
    r = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0: return 0
        if a % 8 in (3, 5): r = -r
    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def h_field(m):
    """Class number of Q(sqrt(-m)), m > 4 squarefree-free-form not required:
    uses the field discriminant of Q(sqrt(-m))."""
    D = -m if (-m) % 4 == 1 else -4*m
    cnt = 0
    amax = math.isqrt(-D//3) + 1
    for a in range(1, amax+1):
        for b in range(abs(D) % 2, a+1, 2):
            if (b*b - D) % (4*a): continue
            c = (b*b - D)//(4*a)
            if c < a: continue
            if math.gcd(math.gcd(a, b), c) != 1: continue
            cnt += 1 if (b == 0 or a == b or a == c) else 2
    return cnt

def B2chi(p):
    D0 = p if p % 4 == 1 else 4*p
    s = Fr(0)
    for a in range(1, D0+1):
        ch = kronecker(D0, a)
        if ch:
            x = Fr(a, D0)
            s += ch * (x*x - x + Fr(1, 6))
    return D0 * s

def d_exact(p):
    s = kronecker(2, p)
    B = B2chi(p)
    h1, h2m, h3 = h_field(p), h_field(2*p), h_field(3*p)
    if p % 4 == 1:
        d = Fr(9-2*s, 96)*B + Fr(5,16)*h1 + Fr(1,8)*h2m + Fr(3+s,12)*h3 - Fr(p+7, 8)
    else:
        d = Fr(1,96)*B + Fr(13-5*s,16)*h1 + Fr(1,8)*h2m + Fr(1,12)*h3 - Fr(p+5, 8)
    return d, B, s

def crude_positive(p, s):
    """True iff the dropped-terms bound  coef * D^{3/2}/15 - linear  is > 0,
    with D = p (p=1 mod 4) or 4p, tested exactly via squaring."""
    if p % 4 == 1:
        coef, lin = Fr(9-2*s, 96*15), Fr(p+7, 8)   # coef * p^{3/2} > lin ?
        # p^{3/2} > (lin/coef)  <=>  p^3 > (lin/coef)^2
        rhs = lin / coef
        return p**3 > rhs*rhs
    else:
        coef, lin = Fr(8, 96*15), Fr(p+5, 8)       # D^{3/2} = 8 p^{3/2}
        rhs = lin / coef
        return p**3 > rhs*rhs

def main():
    primes = [q for q in sieve(PMAX) if q >= 7]
    worst_fail = None; nneg = 0
    for p in primes:
        d, B, s = d_exact(p)
        assert B > 0, ("B2chi <= 0", p)
        assert d.denominator == 1, ("d not integral", p, d)
        assert d >= 0, ("d negative", p, d)
        if not crude_positive(p, s):
            worst_fail = p
        if p in (61, 73, 79):
            assert d == 1, ("pin", p, d)
    assert worst_fail == 673, ("largest crude failure", worst_fail)
    print(f"checked {len(primes)} primes 7 <= p <= {PMAX}")
    print("(1) d(p) integral and d(p) >= 0 at EVERY prime            OK")
    print("(2) crude bound positive for every prime p > 673;")
    print(f"    largest prime where the crude bound fails: {worst_fail}      OK")
    print("(3) pins d(61) = d(73) = d(79) = 1                        OK")
    print("="*64)
    print(f"POSITIVITY (Theorem 8.5, finite input): PASSED, 7 <= p <= {PMAX}")
    print("="*64)

if __name__ == "__main__":
    if len(sys.argv) > 1: PMAX = int(sys.argv[1])
    main()
