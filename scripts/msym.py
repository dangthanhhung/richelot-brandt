"""Modular symbols for Gamma_0(p), weights 2 and 4: Manin symbols with
Merel's Heilbronn set; T_2 and U_p; cuspidal charpoly split by Fricke sign.
Convention auto-calibrated against the elliptic-graph weight-2 data."""
import json, itertools
from fractions import Fraction
from flint import fmpz_mat, fmpz_poly

def merel(l):
    out=[]
    for a in range(1,l+1):
        for b in range(0,a):
            # ad - bc = l, d>c>=0
            for c in range(0,l+1):
                num=l+b*c
                if num%a: continue
                d=num//a
                if d>c>=0 and a*d-b*c==l: out.append((a,b,c,d))
    return out

def p1(p):
    L=[(0,1)]+[(1,d) for d in range(p)]
    idx={}
    for i,(u,v) in enumerate(L): idx[(u,v)]=i
    def norm(u,v):
        u%=p; v%=p
        if u==0:
            assert v%p!=0
            return (0,1)
        ui=pow(u,p-2,p)
        return (1,(v*ui)%p)
    return L,norm

def act_poly(P,g,var):
    # P: coeff list of poly in (X,Y), homogeneous deg 2: [c_X2, c_XY, c_Y2]
    # variant 0: (X,Y) -> (aX+bY, cX+dY); variant 1: (X,Y) -> (dX-bY, -cX+aY)
    a,b,c,d=g
    if var==1: a,b,c,d = d,-b,-c,a
    # X' = aX+bY, Y' = cX+dY ; P(X',Y')
    cX2,cXY,cY2=P
    # (aX+bY)^2 = a2X2+2abXY+b2Y2 ; (aX+bY)(cX+dY)=acX2+(ad+bc)XY+bdY2 ; (cX+dY)^2
    return (cX2*a*a+cXY*a*c+cY2*c*c,
            cX2*2*a*b+cXY*(a*d+b*c)+cY2*2*c*d,
            cX2*b*b+cXY*b*d+cY2*d*d)

def build(p,k,var):
    L,norm=p1(p); n1=len(L)
    dP = 1 if k==2 else 3
    N = n1*dP
    def sid(i,j): return i*dP+j
    def act_sym(vec_updates, i, j, g, coef):
        u,v=L[i]; a,b,c,d=g
        uu,vv=(u*a+v*c)%p,(u*b+v*d)%p
        if uu==0 and vv==0: return
        ii=None
        uu,vv=norm(uu,vv)
        ii={(0,1):0}.get((uu,vv), 1+vv if uu==1 else None)
        if k==2:
            vec_updates.append((sid(ii,0),coef))
        else:
            P=[0,0,0]; P[j]=1
            Q=act_poly(tuple(P),g,var)
            for jj in range(3):
                if Q[jj]: vec_updates.append((sid(ii,jj),coef*Q[jj]))
    # relations: x + x.sigma = 0 ; x + x.tau + x.tau^2 = 0
    sig=(0,-1,1,0); tau=(0,-1,1,-1); tau2=(-1,1,-1,0)
    rows=[]
    for i in range(n1):
        for j in range(dP):
            r={} 
            def add(t):
                for (s,c) in t: r[s]=r.get(s,0)+c
            add([(sid(i,j),1)]); t=[]; act_sym(t,i,j,sig,1); add(t)
            rows.append(r)
            r={}; add([(sid(i,j),1)])
            t=[]; act_sym(t,i,j,tau,1); add(t)
            t=[]; act_sym(t,i,j,tau2,1); add(t)
            rows.append(r)
    A=[[0]*N for _ in range(len(rows))]
    for ri,r in enumerate(rows):
        for s,c in r.items(): A[ri][s]=c
    M=fmpz_mat(A)
    K=M.nullspace()  # columns spanning kernel? flint: returns (mat, nullity)?
    return L,norm,N,dP,sid,act_sym,K

def hecke_on_quotient(p,k,l,var):
    L,norm=p1(p); n1=len(L); dP=1 if k==2 else 3; N=n1*dP
    def sid(i,j): return i*dP+j
    # build relation matrix, get quotient via rref pivots
    sig=(0,-1,1,0); tau=(0,-1,1,-1); tau2=(-1,1,-1,0)
    def sym_act(i,j,g):
        u,v=L[i]; a,b,c,d=g
        uu,vv=(u*a+v*c)%p,(u*b+v*d)%p
        if uu==0 and vv==0: return []
        uu,vv=norm(uu,vv)
        ii=0 if (uu,vv)==(0,1) else 1+vv
        if k==2: return [(sid(ii,0),1)]
        P=[0,0,0]; P[j]=1
        Q=act_poly(tuple(P),g,var)
        return [(sid(ii,jj),Q[jj]) for jj in range(3) if Q[jj]]
    rows=[]
    for i in range(n1):
        for j in range(dP):
            r={sid(i,j):1}
            for s,c in sym_act(i,j,sig): r[s]=r.get(s,0)+c
            rows.append(r)
            r={sid(i,j):1}
            for s,c in sym_act(i,j,tau): r[s]=r.get(s,0)+c
            for s,c in sym_act(i,j,tau2): r[s]=r.get(s,0)+c
            rows.append(r)
    A=[[0]*N for _ in range(len(rows))]
    for ri,r in enumerate(rows):
        for s,c in r.items(): A[ri][s]=c
    from fractions import Fraction as Fr
    # rational rref to find pivots/free variables
    from flint import fmpq_mat
    RRm,rank=fmpq_mat(A).rref()
    piv=[]
    for r in range(rank):
        c=next(cc for cc in range(N) if RRm[r,cc]!=0)
        piv.append(c)
    pivset=set(piv)
    free=[c for c in range(N) if c not in pivset]
    dim=len(free)
    expr=[[Fr(0)]*dim for _ in range(N)]
    for fi,c in enumerate(free): expr[c][fi]=Fr(1)
    for r in range(rank):
        c=piv[r]
        for fi,cc in enumerate(free):
            v=RRm[r,cc]
            if v!=0: expr[c][fi]=-Fr(int(v.p),int(v.q))
    def op(l):
        H=merel(l)
        T=[[Fr(0)]*dim for _ in range(dim)]
        for t,c in enumerate(free):
            i,j=divmod(c,dP)
            acc={}
            for g in H:
                for s,co in sym_act(i,j,g):
                    acc[s]=acc.get(s,0)+co
            for s,co in acc.items():
                for tt in range(dim):
                    if expr[s][tt]: T[tt][t]+=co*expr[s][tt]
        return T
    return dim,op

def int_charpoly(T):
    from math import lcm as _lcm
    den=1
    for row in T:
        for x in row: den=_lcm(den,x.denominator)
    n=len(T)
    Mi=fmpz_mat([[int(x*den) for x in row] for row in T])
    cp=Mi.charpoly()
    # char_{den*T}(den*x) = den^n char_T(x): coefficients of char_T integer if T integral-similar
    co=[int(cp[i]) for i in range(n+1)]
    out=[]
    for i in range(n+1):
        v=Fraction(co[i], den**(n-i))
        assert v.denominator==1,("nonintegral charpoly",i)
        out.append(int(v))
    return fmpz_poly(out)

def w4data(p):
    """cusp orbit polys (in a2 variable) with Fricke signs, weight 4 level p"""
    dim,op=hecke_on_quotient(p,4,2,0)
    T=op(2)
    cp=int_charpoly(T)
    fac=cp.factor()[1]
    orbits=[]; eis=0
    for f,e in fac:
        assert e%2==0,("odd multiplicity",p,str(f))
        if str(f)=="x + (-9)": eis=e//2; continue
        orbits.append((f,e//2))
    # signs via U_p trace on each orbit block: eigenvalue -eps*p with mult 2*deg*m
    U=op(p)
    n=len(T)
    from fractions import Fraction as Fr
    def matmul(A,B):
        m=len(A); q=len(B[0]); r=len(B)
        return [[sum(A[i][t]*B[t][j] for t in range(r)) for j in range(q)] for i in range(m)]
    def polyeval(coeffs,M):
        R=[[Fr(1) if i==j else Fr(0) for j in range(n)] for i in range(n)]
        out=[[Fr(0)]*n for _ in range(n)]
        # Horner
        res=[[Fr(coeffs[-1]) if i==j else Fr(0) for j in range(n)] for i in range(n)]
        for c in reversed(coeffs[:-1]):
            res=matmul(res,T)
            for i in range(n): res[i][i]+=c
        return res
    out=[]
    cpz=cp
    for f,m in orbits:
        d=f.degree()
        # projector: h = cof * (cof^{-1} mod f^m), cof = cp/(f^m)
        fm=f**(2*m)
        cof=fmpz_poly([1])
        for g,e in fac:
            if str(g)!=str(f): cof=cof*g**e
        # inverse of cof mod f^(2m) over Q: xgcd in Q[x]
        import flint
        from flint import fmpq_poly
        A=fmpq_poly(cof); B=fmpq_poly(fm)
        g,s,t=A.xgcd(B)
        assert g.degree()==0
        h=(A*s)/g[0]  # h = cof*s/g0 == 1 mod fm, == 0 mod cof-part
        hh=h % fmpq_poly(cpz)
        coeffs=[Fr(int(hh.numer()[i]),int(hh.denom())) for i in range(hh.degree()+1)]
        P=polyeval(coeffs,T)
        tr=sum(matmul(U,P)[i][i] for i in range(n))
        val=tr/Fr(2*d*m)
        assert val in (Fr(p),Fr(-p)),("sign trace",p,str(f),tr)
        eps = -1 if val==p else 1   # U_p eigenvalue = -eps*p
        out.append((str(f),m,eps))
    return dict(dim=dim,eis=eis,orbits=out)

if __name__=="__main__":
    import sys,json,time
    ps=[int(a) for a in sys.argv[1:] if a.isdigit()] or [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149]
    try: res=json.load(open("w4data_part1.json"))
    except Exception: res={}
    for p in ps:
        t0=time.time(); d=w4data(p); res[str(p)]=d
        json.dump(res,open("w4data_part1.json","w"))
        print(f"p={p}: dim M4 = {d['dim']}, orbits = {[(e, f if len(f)<25 else f[:22]+'...') for f,m,e in d['orbits']]}  ({time.time()-t0:.0f}s)")
        sys.stdout.flush()
