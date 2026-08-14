"""Richelot (2,2)-isogeny graph engine, from scratch.
Fields F_{p^D} as F_p[x]/(mod) via flint nmod_poly; adaptive extensions."""
import random, itertools, json
from fractions import Fraction
from math import gcd as _g
from flint import nmod_poly, fmpz_mat, fmpz_poly
random.seed(20260718)

def lcm(a,b): return a*b//_g(a,b)

class Fld:
    def __init__(s, p, D):
        s.p, s.D = p, D
        s.zero = nmod_poly([0], p); s.one = nmod_poly([1], p)
        s.x = nmod_poly([0,1], p)
        s.mod = s._find_irred(D)
        s.emb_cache = {}
    def _find_irred(s, D):
        p = s.p
        if D == 1: return nmod_poly([0,1], p)
        while True:
            c = [random.randrange(p) for _ in range(D)] + [1]
            m = nmod_poly(c, p)
            X = s.x
            t = X
            for _ in range(D): t = t.pow_mod(p, m)
            if t != X % m: continue
            ok = True
            for r in set(_pf(D)):
                t2 = X
                for _ in range(D//r): t2 = t2.pow_mod(p, m)
                if (t2 - X).gcd(m).degree() != 0: ok = False; break
            if ok: return m
    def red(s,a): return a % s.mod
    def c(s,n): return nmod_poly([n % s.p], s.p)
    def add(s,a,b): return a+b
    def sub(s,a,b): return a-b
    def mul(s,a,b): return (a*b) % s.mod
    def neg(s,a): return -a
    def is0(s,a): return a.degree() < 0
    def inv(s,a):
        g, u, v = a.xgcd(s.mod)
        assert g.degree() == 0
        gi = pow(int(g[0]), s.p-2, s.p)
        return (u * gi) % s.mod
    def pw(s,a,e):
        r, b = s.one, a % s.mod
        while e:
            if e & 1: r = s.mul(r,b)
            b = s.mul(b,b); e >>= 1
        return r
    def rnd(s): return nmod_poly([random.randrange(s.p) for _ in range(s.D)], s.p)
    def frob_ord(s,a):
        t = s.pw(a, s.p); k = 1
        while t != a: t = s.pw(t, s.p); k += 1
        return k

def _pf(n):
    f=[]; d=2
    while d*d<=n:
        while n%d==0: f.append(d); n//=d
        d+=1
    if n>1: f.append(n)
    return f

# ---------- polynomials over Fq (dense lists, index=degree) ----------
def pnorm(F,f):
    while f and F.is0(f[-1]): f.pop()
    return f
def pmul(F,f,g):
    if not f or not g: return []
    r=[F.zero]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        if F.is0(a): continue
        for j,b in enumerate(g): r[i+j]=r[i+j]+F.mul(a,b)
    return pnorm(F,r)
def psub(F,f,g):
    n=max(len(f),len(g)); r=[]
    for i in range(n):
        a=f[i] if i<len(f) else F.zero; b=g[i] if i<len(g) else F.zero
        r.append(a-b)
    return pnorm(F,r)
def pdivmod(F,f,g):
    g=pnorm(F,list(g)); f=pnorm(F,list(f))
    assert g
    li=F.inv(g[-1]); gm=[F.mul(c,li) for c in g]
    q=[F.zero]*max(1,len(f)-len(gm)+1); r=list(f)
    while len(r)>=len(gm) and r:
        c=r[-1]; d=len(r)-len(gm)
        q[d]=q[d]+c
        for i,cg in enumerate(gm):
            r[d+i]=r[d+i]-F.mul(c,cg)
        r=pnorm(F,r)
    return pnorm(F,[F.mul(c,li) for c in q]), r
def pmod(F,f,g): return pdivmod(F,f,g)[1]
def pmonic(F,f):
    f=pnorm(F,list(f))
    if not f: return f
    li=F.inv(f[-1]); return [F.mul(c,li) for c in f]
def pgcd(F,f,g):
    f,g=pnorm(F,list(f)),pnorm(F,list(g))
    while g: f,g=g,pmod(F,f,g)
    return pmonic(F,f)
def ppow_mod(F,b,e,m):
    r=[F.one]; b=pmod(F,b,m)
    while e:
        if e&1: r=pmod(F,pmul(F,r,b),m)
        b=pmod(F,pmul(F,b,b),m); e>>=1
    return r
def pderiv(F,f):
    return pnorm(F,[F.mul(F.c(i),f[i]) for i in range(1,len(f))])

def roots_in(F,f):
    """distinct roots of f lying in F"""
    f=pmonic(F,f)
    if len(f)<=1: return []
    q=F.p**F.D
    X=[F.zero,F.one]
    g=pgcd(F,f,psub(F,ppow_mod(F,X,q,f),X))
    return _cz(F,g,q)
def _cz(F,g,q):
    if len(g)<=1: return []
    if len(g)==2: return [F.neg(g[0])]
    while True:
        r=F.rnd()
        h=ppow_mod(F,[r,F.one],(q-1)//2,g)
        h=psub(F,h,[F.one])
        d=pgcd(F,g,h)
        if 0<len(d)-1<len(g)-1:
            g2,_=pdivmod(F,g,d)
            return _cz(F,d,q)+_cz(F,g2,q)
def ddf_lcm(F,f):
    """lcm of degrees of irreducible factors of squarefree monic f over F"""
    f=pmonic(F,f); X=[F.zero,F.one]; q=F.p**F.D
    t=list(X); L=1; ff=list(f); d=0
    while len(ff)>2:
        d+=1
        t=ppow_mod(F,t,q,f)
        g=pgcd(F,ff,psub(F,t,X))
        if len(g)>1:
            L=lcm(L,d); ff,_=pdivmod(F,ff,g)
    if len(ff)==2: L=lcm(L,1)
    return L

# ---------- field registry, embeddings, downmaps ----------
class Tower:
    def __init__(s,p):
        s.p=p; s.flds={}
    def F(s,D):
        if D not in s.flds: s.flds[D]=Fld(s.p,D)
        return s.flds[D]

class Emb:
    """embedding Fs -> Fb (Fs.D | Fb.D), with exact down-map"""
    def __init__(s,Fs,Fb):
        s.Fs,s.Fb=Fs,Fb; d,p=Fs.D,Fs.p
        lift=[Fb.c(int(Fs.mod[i])) for i in range(d+1)]
        rts=roots_in(Fb,lift); assert rts, "no embedding root"
        s.g=rts[0]
        s.pows=[Fb.one]
        for _ in range(d-1): s.pows.append(Fb.mul(s.pows[-1],s.g))
        # gaussian solver: matrix (Fb.D x d) over Fp of coeff vectors of pows
        A=[[int(s.pows[j][i]) for j in range(d)] for i in range(Fb.D)]
        s.solver=_gauss_prepare(A,p,d)
    def up(s,a):
        r=s.Fb.zero
        for i in range(a.degree()+1):
            r=r+s.pows[i]*int(a[i])
        return r % s.Fb.mod
    def down(s,z):
        vec=[int(z[i]) for i in range(s.Fb.D)]
        co=_gauss_solve(s.solver,vec,s.Fs.p)
        assert co is not None, "element not in subfield"
        return nmod_poly(co,s.Fs.p)

def _gauss_prepare(A,p,d):
    rows=len(A); M=[row[:] for row in A]
    piv=[]; r=0
    for c in range(d):
        pr=None
        for i in range(r,rows):
            if M[i][c]%p: pr=i; break
        assert pr is not None
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],p-2,p)
        M[r]=[(v*inv)%p for v in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(d)]
        piv.append(c); r+=1
    return (A,p,d)
def _gauss_solve(prep,vec,p):
    A,_,d=prep; rows=len(A)
    M=[A[i][:]+[vec[i]%p] for i in range(rows)]
    r=0
    for c in range(d):
        pr=None
        for i in range(r,rows):
            if M[i][c]%p: pr=i; break
        if pr is None: return None
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(d+1)]
        r+=1
    sol=[0]*d
    for i in range(r): 
        c=next(j for j in range(d) if M[i][j]); sol[c]=M[i][d]
    for i in range(r,rows):
        if M[i][d]%p: return None
    return sol

EMBS={}
def emb(T,Fs,Fb):
    k=(Fs.D,Fb.D)
    k=(Fs.p,)+k
    if k not in EMBS: EMBS[k]=Emb(Fs,Fb)
    return EMBS[k]

# ---------- combinatorics ----------
IDX=list(range(6))
def _pairings():
    out=[]
    def rec(rem,cur):
        if not rem: out.append(tuple(cur)); return
        a=rem[0]
        for b in rem[1:]:
            rec([x for x in rem if x not in (a,b)],cur+[(a,b)])
    rec(IDX,[])
    return out
PAIR15=_pairings()
TRI10=[((0,)+t, tuple(x for x in range(1,6) if x not in t)) for t in itertools.combinations(range(1,6),2)]
PERM3=list(itertools.permutations(range(3)))

# ---------- Igusa-type root invariants; key over abstract F_{p^2} ----------
def inv_key(T,F,rs):
    d={} 
    for i in range(6):
        for j in range(6):
            if i!=j: d[(i,j)]=rs[i]-rs[j]
    sq=lambda a: F.mul(a,a)
    Ddisc=F.one
    for i in range(6):
        for j in range(i+1,6): Ddisc=F.mul(Ddisc,sq(d[(i,j)]))
    A=F.zero
    for P in PAIR15:
        t=F.one
        for (i,j) in P: t=F.mul(t,sq(d[(i,j)]))
        A=A+t
    def qtri(t):
        (i,j,k)=t
        return sq(F.mul(F.mul(d[(i,j)],d[(j,k)]),d[(k,i)]))
    B=F.zero; C=F.zero
    for (T1,T2) in TRI10:
        qq=F.mul(qtri(T1),qtri(T2))
        B=B+qq
        m=F.zero
        for s in PERM3:
            t=F.one
            for l in range(3): t=F.mul(t,sq(d[(T1[l],T2[s[l]])]))
            m=m+t
        C=C+F.mul(qq,m)
    Di=F.inv(Ddisc)
    A2=F.mul(A,A); A3=F.mul(A2,A); A5=F.mul(A2,A3)
    B2=F.mul(B,B); B5=F.mul(F.mul(B2,B2),B)
    C2=F.mul(C,C); C5=F.mul(F.mul(C2,C2),C)
    r=[F.mul(A5,Di), F.mul(F.mul(A3,B),Di), F.mul(F.mul(A2,C),Di),
       F.mul(F.mul(B,C),Di), F.mul(B5,F.mul(Di,Di)),
       F.mul(C5,F.mul(Di,F.mul(Di,Di)))]
    assert F.D%2==0
    e=emb(T,T.F(2),F) if F.D>2 else None
    out=[]
    for z in r:
        if F.D==2: w=z
        else: w=e.down(z)
        a=int(w[0]) if w.degree()>=0 else 0
        b=int(w[1]) if w.degree()>=1 else 0
        out.append((a,b))
    return ('J',)+tuple(out)

# ---------- Moebius ----------
def mob_triple(F,tr):
    a,b,c=tr
    return [[b-c, F.neg(F.mul(a,b-c))],[b-a, F.neg(F.mul(c,b-a))]]
def mob_inv(F,M):
    return [[M[1][1],F.neg(M[0][1])],[F.neg(M[1][0]),M[0][0]]]
def mob_mul(F,M,N):
    return [[F.mul(M[0][0],N[0][0])+F.mul(M[0][1],N[1][0]), F.mul(M[0][0],N[0][1])+F.mul(M[0][1],N[1][1])],
            [F.mul(M[1][0],N[0][0])+F.mul(M[1][1],N[1][0]), F.mul(M[1][0],N[0][1])+F.mul(M[1][1],N[1][1])]]
def mob_ap(F,M,z):
    if z is None:
        den=M[1][0]
        return None if F.is0(den) else F.mul(M[0][0],F.inv(den))
    num=F.mul(M[0][0],z)+M[0][1]; den=F.mul(M[1][0],z)+M[1][1]
    if F.is0(den): return None
    return F.mul(num,F.inv(den))

def norm6(F,rs,ninf):
    """Given <=6 finite roots + ninf points at infinity, Moebius-normalize
    to 6 finite distinct roots (z -> 1/(z-c))."""
    assert len(rs)+ninf==6
    if ninf==0:
        assert len(set(map(str,rs)))==6
        return rs
    c=0
    while True:
        cc=F.c(c)
        if all(not F.is0(r-cc) for r in rs): break
        c+=1
    out=[F.inv(r-cc) for r in rs]+[F.zero]*ninf
    assert len(set(map(str,out)))==6
    return out

def dn2(T,F,z):
    if F.D==2: w=z
    else: w=emb(T,T.F(2),F).down(z)
    a=int(w[0]) if w.degree()>=0 else 0
    b=int(w[1]) if w.degree()>=1 else 0
    return (a,b)

def cr4(F,P):
    def s(x,y):
        if x is None or y is None: return F.one
        return x-y
    num=F.mul(s(P[0],P[2]),s(P[1],P[3])); den=F.mul(s(P[0],P[3]),s(P[1],P[2]))
    return F.mul(num,F.inv(den))
def j_from_lam(F,l):
    l2=F.mul(l,l); t=l2-l+F.one
    num=F.mul(F.c(256),F.pw(t,3))
    d0=F.mul(l,l-F.one); den=F.mul(d0,d0)
    return F.mul(num,F.inv(den))
def j_of_branch4(F,P): return j_from_lam(F,cr4(F,P))

def split_quads(T,F,quads):
    quads=[pnorm(F,list(q)) for q in quads]
    need=False
    for q in quads:
        if len(q)==3 and len(roots_in(F,q))<2: need=True; break
    if need:
        Fb=T.F(2*F.D); e0=emb(T,F,Fb)
        quads=[[e0.up(c) for c in q] for q in quads]; F=Fb
    roots=[]; ninf=0
    for q in quads:
        if len(q)==3:
            rr=roots_in(F,q); assert len(rr)==2, "quad not split"
            roots+=rr
        elif len(q)==2:
            roots.append(F.mul(F.neg(q[0]),F.inv(q[1]))); ninf+=1
        else: assert False, "H constant"
    return roots,ninf,F

def det3(F,M):
    return (F.mul(M[0][0],F.mul(M[1][1],M[2][2])-F.mul(M[1][2],M[2][1]))
          - F.mul(M[0][1],F.mul(M[1][0],M[2][2])-F.mul(M[1][2],M[2][0]))
          + F.mul(M[0][2],F.mul(M[1][0],M[2][1])-F.mul(M[1][1],M[2][0])))

def jac_edge(T,F,rs,P):
    G=[]
    for (i,j) in P:
        G.append([F.mul(rs[i],rs[j]), F.neg(rs[i]+rs[j]), F.one])
    if F.is0(det3(F,[[g[2],g[1],g[0]] for g in G])):
        return ('prod', degen_split(T,F,rs,P))
    H=[]
    for i in range(3):
        j,k=(i+1)%3,(i+2)%3
        H.append(psub(F,pmul(F,pderiv(F,G[j]),G[k]),pmul(F,G[j],pderiv(F,G[k]))))
    roots,ninf,Fc=split_quads(T,F,H)
    return ('jac', norm6(Fc,roots,ninf), Fc)

def degen_split(T,F,rs,P):
    (a,b),(c,d),(e,f)=[(rs[i],rs[j]) for (i,j) in P]
    A=mob_triple(F,(a,c,e)); B=mob_triple(F,(b,d,f))
    S=mob_mul(F,mob_inv(F,B),A)
    assert mob_ap(F,S,b)==a and mob_ap(F,S,d)==c and mob_ap(F,S,f)==e, "not involution"
    assert F.is0(S[0][0]+S[1][1]), "involution must have trace 0"
    quad=[F.neg(S[0][1]), S[1][1]-S[0][0], S[1][0]]
    quad=pnorm(F,quad)
    if len(quad)==3:
        fps,Ff=(roots_in(F,quad),F)
        if len(fps)<2:
            Fb=T.F(2*F.D); e0=emb(T,F,Fb)
            fps=roots_in(Fb,[e0.up(cq) for cq in quad]); Ff=Fb
            a,b,c,d,e,f=[e0.up(z) for z in (a,b,c,d,e,f)]
        q1,q2=fps
        iv=Ff.inv
        def phi(z): return Ff.mul(z-q1, iv(z-q2))
    else:
        # affine involution s(z) = -z + beta: fixed points beta/2 and infinity
        assert len(quad)==2, "parabolic involution (impossible in odd char)"
        Ff=F; q1=F.mul(F.neg(quad[0]),F.inv(quad[1]))
        def phi(z): return z-q1
    u=[]
    for (x,y) in ((a,b),(c,d),(e,f)):
        px,py=phi(x),phi(y)
        assert px == Ff.neg(py), "phi anti-symmetry fails"
        u.append(Ff.mul(px,px))
    j1=j_of_branch4(Ff,[u[0],u[1],u[2],None])
    j2=j_of_branch4(Ff,[Ff.zero,u[0],u[1],u[2]])
    k1,k2=dn2(T,Ff,j1),dn2(T,Ff,j2)
    return tuple(sorted([k1,k2]))

def velu2_j(F,et,x0):
    s=et[0]+et[1]+et[2]; sh=F.mul(s,F.inv(F.c(3)))
    E=[z-sh for z in et]; X0=x0-sh
    A=F.mul(E[0],E[1])+F.mul(E[0],E[2])+F.mul(E[1],E[2])
    Bb=F.neg(F.mul(F.mul(E[0],E[1]),E[2]))
    t=F.mul(F.c(3),F.mul(X0,X0))+A
    a2=A-F.mul(F.c(5),t); b2=Bb-F.mul(F.c(7),F.mul(X0,t))
    a3=F.mul(F.mul(a2,a2),a2)
    den=F.mul(F.c(4),a3)+F.mul(F.c(27),F.mul(b2,b2))
    return F.mul(F.mul(F.c(1728),F.mul(F.c(4),a3)),F.inv(den))

def glue(T,F,Ea,Eb,sg):
    e=Ea; f=[Eb[sg[i]] for i in range(3)]
    a11=F.neg(e[0]-e[1]); a12=f[0]-f[1]
    a21=F.neg(e[0]-e[2]); a22=f[0]-f[2]
    det=F.mul(a11,a22)-F.mul(a12,a21)
    if F.is0(det): return None
    r1=F.neg(F.mul(e[0],f[0])-F.mul(e[1],f[1]))
    r2=F.neg(F.mul(e[0],f[0])-F.mul(e[2],f[2]))
    di=F.inv(det)
    tt=F.mul(di,F.mul(a22,r1)-F.mul(a12,r2))
    bb=F.mul(di,F.mul(a11,r2)-F.mul(a21,r1))
    u=[z+bb for z in e]
    if any(F.is0(z) for z in u): return None
    k=F.mul(u[0],f[0]-tt)
    if F.is0(k): return None
    for i in (1,2): assert F.mul(u[i],f[i]-tt)==k, "glue inconsistency"
    roots,ninf,Fc=split_quads(T,F,[[F.neg(z),F.zero,F.one] for z in u])
    assert ninf==0
    return ('jac', norm6(Fc,roots,0), Fc)

def wgt(k,p):
    if k==(0,0): return 6
    if k==(1728%p,0): return 4
    return 2

def ss_data(T):
    from math import comb
    p=T.p; F2=T.F(2); m=(p-1)//2
    H=[F2.c(comb(m,i)**2 % p) for i in range(m+1)]
    lams=roots_in(F2,H)
    assert len(lams)==m, (len(lams),m)
    byj={}
    for l in lams:
        byj.setdefault(dn2(T,F2,j_from_lam(F2,l)), l)
    assert sum(Fraction(1,wgt(k,p)) for k in byj)==Fraction(p-1,24), "Eichler mass"
    return byj

def jac_weight(F,rs):
    def kz(z): return tuple(int(z[i]) if z.degree()>=i else 0 for i in range(F.D))
    rset={kz(r) for r in rs}
    A=mob_triple(F,(rs[0],rs[1],rs[2])); cnt=0
    for (y1,y2,y3) in itertools.permutations(rs,3):
        phi=mob_mul(F,mob_inv(F,mob_triple(F,(y1,y2,y3))),A)
        ok=True
        for k in range(3,6):
            im=mob_ap(F,phi,rs[k])
            if im is None or kz(im) not in rset: ok=False; break
        if ok: cnt+=1
    return 2*cnt

# ---------- arithmetic side (Ibukiyama closed forms; inlined, self-contained) ----------
#!/usr/bin/env python3
# certify_section4.py -- self-contained, cached, fast.
# [P] TrR-TrRnpg = nu(p)(p-e3)/12 ; [A][C] integrality of AL splits from
# TrW2=1-nu/2, TrW4=nu/2 ; [B] d2+d4+1=(p-e3)/3 ; [S3] L0-1-b4 = nu(p-e3)/12 ;
# [X] cross-family ; anchors p=11,13,37,61 and both printed tables.
from fractions import Fraction as F
from functools import lru_cache
from math import gcd

def jac(a, n):                     # Jacobi symbol, n odd > 0
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def kron(D, n):                    # Kronecker (D/n), D>0, n>=1
    r = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0: return 0
        if D % 8 in (3, 5): r = -r
    return r * jac(D % n, n) if n > 1 else r

@lru_cache(maxsize=None)
def h_field(m):                    # class number of Q(sqrt(-m)), m>0 squarefree
    D = -m if (-m) % 4 == 1 else -4*m
    h = 0; a = 1
    while 3*a*a <= -D:
        for b in range(-a, a+1):
            num = b*b - D
            if num % (4*a) == 0:
                c = num // (4*a)
                if c >= a:
                    if b < 0 and (abs(b) == a or a == c): continue
                    if gcd(gcd(a, abs(b)), c) == 1: h += 1
        a += 1
    return h

@lru_cache(maxsize=None)
def B2chi(p):
    D0 = p if p % 4 == 1 else 4*p
    s2 = s1 = 0
    for a in range(1, D0+1):
        k = kron(D0, a); s2 += k*a*a; s1 += k*a
    return F(s2, D0) - s1

def nu(p):
    return h_field(p) if p % 4 == 1 else (4 if p % 8 == 3 else 2)*h_field(p)

def TrR_pr(p):
    s = jac(2, p); t = jac(p % 3, 3) if p % 3 else 0
    t = jac(p, 3)
    hp, h2p, h3p = h_field(p), h_field(2*p), h_field(3*p); B = B2chi(p)
    if p % 4 == 3:
        c = F(1,48)*(p-1)*(9-4*s)+F(1,16)*(p-s)+F(1,12)*(1-t)*(3-s)
        return F(1,96)*B+F(1,8)*h2p+F(1,12)*h3p+c*hp
    return (F(9-2*s,96)*B+F(4*p-1,48)*hp+F(1,8)*h2p
            +F(1,12)*(3+jac(-2 % p, p)* (1 if True else 1))*h3p+F(1,12)*(1-t)*hp) if False else \
           (F(9-2*s,96)*B+F(4*p-1,48)*hp+F(1,8)*h2p
            +F(1,12)*(3+jac((-2) % p, p))*h3p+F(1,12)*(1-t)*hp)

def TrR_npg(p):
    s = jac(2, p); hp,h2p,h3p = h_field(p),h_field(2*p),h_field(3*p); B=B2chi(p)
    if p % 4 == 1:
        return F(9-2*s,96)*B+F(1,16)*hp+F(1,8)*h2p+F(1,12)*(3+s)*h3p
    return F(1,96)*B+F(1,16)*(1-s)*hp+F(1,8)*h2p+F(1,12)*h3p

def sieve(n):
    P=[]; is_c=[False]*(n+1)
    for i in range(2,n+1):
        if not is_c[i]:
            P.append(i)
            for j in range(i*i,n+1,i): is_c[j]=True
    return P



def arith_pack(p):
    F=Fraction
    e1=jac(-1%p,p); e3=jac(-3%p,p); N=nu(p)
    d4=F(p-2+e1,4)
    nu2=1+jac(-1%p,p); nu3=1+jac(-3%p,p)
    d2=F(p+1,12)-F(nu2,4)-F(nu3,3)
    a2=(d2+1-F(N,2))/2; b2=(d2-1+F(N,2))/2
    a4=(d4+F(N,2))/2;  b4=(d4-F(N,2))/2
    for z in (d2,d4,a2,b2,a4,b4): assert z.denominator==1 and z>=0,(p,z)
    L0=1+d4+b2*a4-a2*b4
    d=TrR_pr(p)-(1+d4+b2*a4-a2*b4)  # placeholder; real d below
    dd=TrR_npg(p)-1-b4
    cross=b2*a4+a2*b4
    return dict(d2=int(d2),d4=int(d4),a2=int(a2),b2=int(b2),a4=int(a4),b4=int(b4),
                L0=int(L0),d=int(dd),cross=int(cross),TrRpr=TrR_pr(p))

# ---------- main build ----------
from collections import deque
def build(p, verbose=True):
    T=Tower(p); F2=T.F(2)
    byj=ss_data(T); jl=sorted(byj)
    verts={}; order=[]; models=[]; wts=[]
    def add(key,mdl,w):
        verts[key]=len(order); order.append(key); models.append(mdl); wts.append(w)
        return verts[key]
    for i in range(len(jl)):
        for j2 in range(i,len(jl)):
            ka,kb=jl[i],jl[j2]
            w=2*wgt(ka,p)**2 if i==j2 else wgt(ka,p)*wgt(kb,p)
            add(('P',ka,kb),('P',byj[ka],byj[kb],ka,kb),w)
    Q=deque(range(len(order))); rows={}
    maxdeg=[2]
    def newjac(rs,Fc):
        L=1
        for r in rs: L=lcm(L,Fc.frob_ord(r))
        L=lcm(L,2)
        if L<Fc.D:
            Ts=T.F(L); em0=emb(T,Ts,Fc)
            rs=[em0.down(r) for r in rs]; Fc=Ts
        maxdeg[0]=max(maxdeg[0],Fc.D)
        key=inv_key(T,Fc,rs)
        if key not in verts:
            Q.append(add(key,('J',rs,Fc),jac_weight(Fc,rs)))
        return key
    nproc=0
    while Q:
        i=Q.popleft(); mdl=models[i]; row={}
        def bump(k): row[k]=row.get(k,0)+1
        if mdl[0]=='P':
            la,lb,ka,kb=mdl[1],mdl[2],mdl[3],mdl[4]
            Ea=(F2.zero,F2.one,la); Eb=(F2.zero,F2.one,lb)
            for x0 in Ea:
                ja=dn2(T,F2,velu2_j(F2,Ea,x0))
                for y0 in Eb:
                    jb=dn2(T,F2,velu2_j(F2,Eb,y0))
                    tk=('P',)+tuple(sorted([ja,jb]))
                    assert tk in verts, ("velu target missing",tk)
                    bump(tk)
            for sg in PERM3:
                res=glue(T,F2,Ea,Eb,sg)
                if res is None: bump(('P',ka,kb))
                else: bump(newjac(res[1],res[2]))
        else:
            rs,Fc=mdl[1],mdl[2]
            for P in PAIR15:
                res=jac_edge(T,Fc,rs,P)
                if res[0]=='prod':
                    tk=('P',)+res[1]
                    assert tk in verts, ("degen target missing",tk)
                    bump(tk)
                else:
                    bump(newjac(res[1],res[2]))
        assert sum(row.values())==15,(i,row)
        rows[i]=row; nproc+=1
        if verbose and nproc%25==0: print(f"  processed {nproc}, verts {len(order)}")
    n=len(order)
    M=[[0]*n for _ in range(n)]
    for i,row in rows.items():
        for k,c in row.items(): M[i][verts[k]]=c
    # certificates
    for i in range(n): assert sum(M[i])==15
    for i in range(n):
        for j2 in range(n):
            assert wts[j2]*M[i][j2]==wts[i]*M[j2][i], ("Mestre",i,j2)
    mass=sum(Fraction(1,w) for w in wts)
    assert mass==Fraction((p-1)*(p*p+1),5760), ("mass",mass)
    # R(pi): frobenius on keys
    g2=F2.mod; gp=F2.pw(F2.x,p)
    c0=int(gp[0]) if gp.degree()>=0 else 0
    c1=int(gp[1]) if gp.degree()>=1 else 0
    def fr(ab): a,b=ab; return ((a+b*c0)%p,(b*c1)%p)
    def frkey(k):
        if k[0]=='P': return ('P',)+tuple(sorted([fr(k[1]),fr(k[2])]))
        return ('J',)+tuple(fr(x) for x in k[1:])
    sigma=[verts[frkey(k)] for k in order]
    assert sorted(sigma)==list(range(n))
    TrR=sum(1 for i in range(n) if sigma[i]==i)
    ar=arith_pack(p)
    assert TrR==ar['TrRpr']==ar['L0']+ar['d'], ("TrR",TrR,ar['TrRpr'],ar['L0']+ar['d'])
    vP=sum(1 for k in order if k[0]=='P'); vJ=n-vP
    delta=ar['d']  # p<167
    degM=Fraction(n-1-ar['d4']-ar['cross']-delta,2)
    assert degM.denominator==1 and degM>=0
    cp=fmpz_mat(M).charpoly(); fac=cp.factor()
    sq=[(f,e) for (f,e) in fac[1] if e>=2]
    sqdeg=sum(f.degree()*e for f,e in [(f,e-1) for f,e in sq])
    res=dict(p=p,h2=n,vP=vP,vJ=vJ,TrB=sum(M[i][i] for i in range(n)),
             TrR=TrR,d=ar['d'],degMys=int(degM),maxdeg=maxdeg[0],
             charpoly_factors=[(str(f),e) for f,e in fac[1]],
             square_excess_deg=sqdeg, weights=sorted(wts))
    return res,M,order,wts,sigma

if __name__ == "__main__":
    import sys as _s, time as _t
    args=[a for a in _s.argv[1:]]
    ps=[int(a) for a in args if a.isdigit()] or [11]
    for p in ps:
        t0=_t.time()
        res,M,order,wts,sigma=build(p,verbose=("-v" in args))
        print(f"p={p}: h2={res['h2']} vP={res['vP']} vJ={res['vJ']} "
              f"TrB2={res['TrB']} TrR={res['TrR']} d={res['d']} "
              f"degMys={res['degMys']} sq={res['square_excess_deg']} "
              f"maxdeg={res['maxdeg']}  ({_t.time()-t0:.1f}s)")
        print("  linear rational factors:",
              [f for f,e in res['charpoly_factors'] if 'x^' not in f])
        json.dump({'p':p,'res':res,'M':M,'wts':wts,'sigma':sigma},
                  open(f"h2_result_{p}.json","w"))
        print(f"  certificates: ALL PASSED  ->  h2_result_{p}.json")
