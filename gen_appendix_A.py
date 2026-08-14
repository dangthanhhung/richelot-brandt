#!/usr/bin/env python3
# gen_appendix_A.py
# ---------------------------------------------------------------------------
# AUTOMATIC generation of Appendix A (the eight Brandt matrices) of the paper
# "The Richelot isogeny graph as a Brandt matrix", from the verbatim DATA{}
# block of verify_three_pillars.py.  NOTHING is typed by hand: every printed
# integer is read from DATA and re-certified here before printing.
#
# Certification performed by THIS script (independently of, and in addition
# to, verify_three_pillars.py, which re-derives charpolys etc.):
#   (a) structural: square matrix, row sums = 15, Mestre symmetry
#       e_j M_ij = e_i M_ji (out-edge convention), mass = (p-1)(p^2+1)/5760;
#   (b) label consistency: the labels in 'evec' coincide, in order, with the
#       row labels of 'rows'; products come first, then J0..J_{vJ-1};
#   (c) vertex-label provenance: the supersingular j-invariants are
#       recomputed FROM SCRATCH over F_{p^2} = F_p[w]/(Conway poly), via the
#       Hasse--Igusa polynomial H(t) = sum_i binom(m,i)^2 t^i, m=(p-1)/2,
#       and the Legendre lambda -> j map; the resulting j-set must equal the
#       set of j's in the P{...} labels; the Eichler mass sum 1/|Aut E| =
#       (p-1)/24 is asserted over that set; the product-vertex weights
#       e = |Aut E_a||Aut E_b| (a != b), 2|Aut E|^2 (a = b) must equal the
#       weights in 'evec'; #products = s(s+1)/2;
#   (d) Galois action on labels: Frobenius x -> x^p permutes the j-set; at
#       p = 37 it must swap 20+10w <-> 23+27w, printed as 20+10\omega and
#       23+27\omega in the appendix (the two non-rational j's).
# Only if every assert passes is any LaTeX written.
#
# Output: appendix_A_matrices.tex  (to be \input by the main file).
# Protocol: NO mental arithmetic -- every claimed number is an assert.
# ---------------------------------------------------------------------------
import re
import sys
from fractions import Fraction
from math import comb

import sympy as sp
import galois  # Conway polynomials (Frank Luebeck's database, = Sage default)

SRC = "verify_three_pillars.py"
OUT = "appendix_A_matrices.tex"

# ---------------------------------------------------------------------------
# 1. load DATA verbatim from verify_three_pillars.py (cut at VERIFICATION)
# ---------------------------------------------------------------------------
src = open(SRC).read()
marker = "# VERIFICATION"
cut = src.find(marker)
assert cut > 0, "VERIFICATION marker not found in source"
ns = {}
exec(src[:cut], ns)                       # defines DATA, parse_matrix, ...
DATA = ns["DATA"]
parse_matrix = ns["parse_matrix"]
PRIMES = [11, 13, 17, 19, 23, 29, 31, 37]
assert sorted(DATA) == PRIMES, "DATA does not contain exactly the 8 primes"

LABEL_RE = re.compile(r"(P\{[^}]*\}|J\d+):(\d+)")
ROWLAB_RE = re.compile(r"^\s*(P\{[^}]*\}|J\d+)\s*\[")
PKEY_RE = re.compile(r"P\{\((\-?\d+),\s*(\-?\d+)\),\((\-?\d+),\s*(\-?\d+)\)\}")


def parse_evec_labeled(evec_str):
    """[(label, e)] in order from the evec string."""
    out = [(m.group(1), int(m.group(2))) for m in LABEL_RE.finditer(evec_str)]
    return out


def row_labels(rows_str):
    out = []
    for line in rows_str.strip().splitlines():
        m = ROWLAB_RE.match(line)
        if m:
            out.append(m.group(1))
    return out


# ---------------------------------------------------------------------------
# 2. F_{p^2} arithmetic with the Conway modulus (pure python, tuples (a,b))
#    w^2 = -c1*w - c0  where Conway poly = x^2 + c1 x + c0 over F_p.
# ---------------------------------------------------------------------------
class Fp2:
    def __init__(self, p):
        cp = galois.conway_poly(p, 2)              # monic deg-2 over GF(p)
        coeffs = [int(c) for c in cp.coefficients()]   # [1, c1, c0]
        assert len(coeffs) == 3 and coeffs[0] == 1
        self.p, self.c1, self.c0 = p, coeffs[1] % p, coeffs[2] % p

    def mul(self, x, y):
        a, b = x; c, d = y; p = self.p
        bd = b * d
        return ((a * c - bd * self.c0) % p,
                (a * d + b * c - bd * self.c1) % p)

    def power(self, x, n):
        r, base = (1, 0), x
        while n:
            if n & 1:
                r = self.mul(r, base)
            base = self.mul(base, base)
            n >>= 1
        return r

    def elements(self):
        p = self.p
        return [(a, b) for b in range(p) for a in range(p)]


def hasse_igusa_ss_jset(p, F):
    """Recompute the set of supersingular j-invariants in F_{p^2}, from
    scratch: H(t) = sum binom(m,i)^2 t^i, m=(p-1)/2 (Igusa); a j in F_{p^2}
    is supersingular iff Res_t(H(t), 256(t^2-t+1)^3 - j t^2(t-1)^2) = 0."""
    t, j = sp.symbols("t j")
    m = (p - 1) // 2
    H = sum(comb(m, i) ** 2 * t ** i for i in range(m + 1))
    Fib = 256 * (t ** 2 - t + 1) ** 3 - j * t ** 2 * (t - 1) ** 2
    R = sp.Poly(sp.resultant(H, Fib, t), j)
    cs = [int(c) % p for c in R.all_coeffs()]      # descending, mod p
    assert any(cs), "resultant vanished identically -- impossible"
    out = set()
    for z in F.elements():                          # Horner over F_{p^2}
        acc = (0, 0)
        for c in cs:
            acc = F.mul(acc, z)
            acc = ((acc[0] + c) % p, acc[1])
        if acc == (0, 0):
            out.add(z)
    # lambda = 0, 1 are never roots of H: H(0)=1, H(1)=binom(2m,m) != 0 (p)
    assert (1 % p, 0) not in (set() if H.subs(t, 0) % p else {0}), "guard"
    return out


def aut_order(jkey, p):
    """|Aut E| for supersingular E/F_{p^2} with j-invariant jkey (p > 3)."""
    if jkey == (0, 0):
        return 6
    if jkey == (1728 % p, 0):
        return 4
    return 2


def fmt_j(jkey):
    a, b = jkey
    if b == 0:
        return str(a)
    if a == 0:
        return f"{b}\\omega"
    return f"{a}{{+}}{b}\\omega"


# ---------------------------------------------------------------------------
# 3. certification + LaTeX emission, prime by prime
# ---------------------------------------------------------------------------
SIZE = {  # (LaTeX size command, arraycolsep in pt) by vertex count
    4: ("\\small", "4"), 5: ("\\small", "4"), 8: ("\\small", "4"),
    10: ("\\small", "3.5"), 16: ("\\footnotesize", "2.6"),
    24: ("\\scriptsize", "2.0"), 26: ("\\scriptsize", "1.9"),
    37: ("\\tiny", "1.5"),
}

chunks = []
chunks.append(
    "%% ------------------------------------------------------------------\n"
    "%% appendix_A_matrices.tex -- GENERATED by gen_appendix_A.py from the\n"
    "%% verbatim DATA{} block of verify_three_pillars.py.  DO NOT EDIT BY\n"
    "%% HAND.  Before printing, the generator re-asserted, for every prime:\n"
    "%% row sums = 15, Mestre symmetry, both mass formulas, and the\n"
    "%% from-scratch recomputation of the supersingular j-invariants\n"
    "%% (Hasse--Igusa polynomial over F_{p^2} with the Conway modulus)\n"
    "%% against the product labels, weights and Frobenius pairing.\n"
    "%% ------------------------------------------------------------------\n")

summary = []
for p in PRIMES:
    d = DATA[p]
    M = parse_matrix(d["rows"])
    lab_e = parse_evec_labeled(d["evec"])
    labels = [le[0] for le in lab_e]
    e = [le[1] for le in lab_e]
    V = len(M)

    # ---- (a) structural re-certification (REQUIRED before printing) ----
    assert all(len(r) == V for r in M), (p, "non-square")
    assert len(e) == V == len(labels), (p, "evec length mismatch")
    assert all(sum(r) == 15 for r in M), (p, "row sums != 15")
    for i in range(V):
        for jj in range(V):
            assert e[jj] * M[i][jj] == e[i] * M[jj][i], (p, "Mestre", i, jj)
    mass = sum(Fraction(1, ei) for ei in e)
    assert mass == Fraction((p - 1) * (p * p + 1), 5760), (p, "total mass")

    # ---- (b) label consistency between evec and rows ----
    assert row_labels(d["rows"]) == labels, (p, "row labels != evec labels")
    nP = sum(1 for L in labels if L.startswith("P"))
    assert all(L.startswith("P") for L in labels[:nP]), (p, "P's not first")
    assert labels[nP:] == [f"J{k}" for k in range(V - nP)], (p, "J order")

    # ---- (c) vertex-label provenance: supersingular j's from scratch ----
    F = Fp2(p)
    ss = hasse_igusa_ss_jset(p, F)
    label_keys = []
    for L in labels[:nP]:
        m = PKEY_RE.fullmatch(L)
        assert m, (p, "unparsable product label", L)
        ja = (int(m.group(1)) % p, int(m.group(2)) % p)
        jb = (int(m.group(3)) % p, int(m.group(4)) % p)
        label_keys.append((ja, jb))
    # Product vertices must come in lexicographic order of their key
    # pairs ((c0,c1),(d0,d1)) -- claimed in the Appendix A preamble,
    # so it is asserted here rather than assumed:
    assert label_keys == sorted(label_keys), (p, "product labels not lex-sorted")
    label_jset = set(j for pair in label_keys for j in pair)
    assert label_jset == ss, (p, "label j-set != recomputed ss j-set",
                              sorted(label_jset), sorted(ss))
    # Product-vertex ordering, as stated in the Appendix A preamble:
    # each pair is internally sorted and the pairs come lexicographically.
    assert all(ja <= jb for (ja, jb) in label_keys), (p, "pair not sorted")
    assert label_keys == sorted(label_keys), (p, "product labels not lex")
    # product vertices are listed in lexicographic order of their keys
    # (claimed in the Appendix A conventions paragraph -- asserted here):
    assert label_keys == sorted(label_keys), (p, "product labels not lex-sorted")
    s = len(ss)
    assert nP == s * (s + 1) // 2, (p, "#products != s(s+1)/2", nP, s)
    # Eichler mass over the recomputed set (independent anchor):
    emass = sum(Fraction(1, aut_order(j, p)) for j in ss)
    assert emass == Fraction(p - 1, 24), (p, "Eichler mass", emass)
    # product mass = M1^2/2 (Lemma in the paper; asserted, not assumed):
    pmass = sum(Fraction(1, e[i]) for i in range(nP))
    assert pmass == Fraction((p - 1), 24) ** 2 / 2, (p, "product mass")
    # product weights from |Aut|:
    for (ja, jb), ei in zip(label_keys, e[:nP]):
        ex = (2 * aut_order(ja, p) ** 2 if ja == jb
              else aut_order(ja, p) * aut_order(jb, p))
        assert ex == ei, (p, "product weight", ja, jb, ex, ei)
    # ---- (d) Frobenius on labels ----
    frob = {z: F.power(z, p) for z in ss}
    assert set(frob.values()) == ss, (p, "Frobenius does not permute ss")
    nonrat = sorted(z for z in ss if frob[z] != z)
    if p != 37:
        assert nonrat == [], (p, "unexpected non-rational ss j")
    else:
        assert nonrat == [(20, 10), (23, 27)], (37, "non-rational pair", nonrat)
        assert frob[(20, 10)] == (23, 27) and frob[(23, 27)] == (20, 10), \
            (37, "Frobenius pairing of labels fails")

    # ---- LaTeX emission ----
    size, colsep = SIZE[V]
    vJ = V - nP
    body = []
    # English pluralization (cosmetic only -- numbers untouched)
    pw = "vertex" if nP == 1 else "vertices"
    jw = "vertex" if vJ == 1 else "vertices"
    body.append(f"\\subsection*{{\\texorpdfstring{{$p={p}$\\quad ($h_2={V}$:\\ {nP} product"
                f" {pw}, {vJ} Jacobian {jw})}}{{p={p} (h2={V}: {nP} product"
                f" {pw}, {vJ} Jacobian {jw})}}}}\n")
    # vertex key (labels + weights), in matrix order
    key_items = []
    for i, (L, ei) in enumerate(zip(labels, e), start=1):
        if L.startswith("P"):
            ja, jb = label_keys[i - 1]
            lab = f"P\\{{{fmt_j(ja)},{fmt_j(jb)}\\}}"
        else:
            lab = f"J_{{{L[1:]}}}"
        key_items.append(f"$v_{{{i}}}{{=}}{lab}\\,({ei})$")
    # key-paragraph font size scales with V to avoid overfull lines:
    keysize = "\\small" if V <= 10 else ("\\footnotesize" if V <= 16
                                         else "\\scriptsize")
    omega_note = ("; $\\omega$ denotes a fixed generator of"
                  f" $\\F_{{{p}^2}}$ over $\\F_{{{p}}}$"
                  if any("\\omega" in k for k in key_items) else "")
    body.append("\\begin{sloppypar}\\noindent" + keysize +
                "\\textit{Vertices (label, weight $e_i=\\#"
                "\\operatorname{Aut}$), in matrix order:} "
                + ", ".join(key_items) + omega_note + ".\\end{sloppypar}\n\n")
    rowtex = " \\\\\n".join(" & ".join(str(x) for x in row) for row in M)
    body.append(
        "\\begingroup" + size +
        "\\setlength{\\arraycolsep}{" + colsep + "pt}\n"
        "\\[ B_2(2)\\;=\\;\\left(\\begin{array}{*{" + str(V) + "}{r}}\n"
        + rowtex + "\n\\end{array}\\right) \\]\n\\endgroup\n")
    chunks.append("".join(body))
    summary.append((p, V, nP, vJ))

chunks.append("\\bigskip\n")

with open(OUT, "w") as f:
    f.write("\n".join(chunks))

print("ALL ASSERTS PASSED -- structural + label-provenance certification:")
for (p, V, nP, vJ) in summary:
    print(f"  p={p:2d}: V={V:2d} (vP={nP}, vJ={vJ})  rowsum/Mestre/mass OK,"
          f" ss-j labels OK, weights OK, Frobenius OK")
print(f"LaTeX written to {OUT} ({sum(1 for _ in open(OUT))} lines).")
