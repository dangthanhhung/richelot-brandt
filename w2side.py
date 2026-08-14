"""Weight-2 side: elliptic Brandt matrix + Fricke involution from the
supersingular graph machinery, with per-sign charpolys, validated against
the closed dimension formulas."""
import sys, json
sys.path.insert(0,'/home/claude/gs')
from richelot import *
from flint import fmpz_mat, fmpz_poly, fmpq
from fractions import Fraction as Fr

def w2data(p):
    T=Tower(p); F2=T.F(2)
    byj=ss_data(T); jl=sorted(byj)
    idx={k:i for i,k in enumerate(jl)}; h=len(jl)
    B=[[0]*h for _ in range(h)]
    for k in jl:
        la=byj[k]; E=(F2.zero,F2.one,la)
        for x0 in E:
            jt=dn2(T,F2,velu2_j(F2,E,x0))
            B[idx[k]][idx[jt]]+=1
    for r in B: assert sum(r)==3
    # weights w(j) and Mestre symmetry (Eichler)
    W=[wgt(k,p) for k in jl]
    for i in range(h):
        for j2 in range(h):
            assert W[j2]*B[i][j2]==W[i]*B[j2][i],("mestre-ell",p,i,j2)
    # Fricke: j -> j^p permutation
    g2=F2.mod; gp=F2.pw(F2.x,p)
    c0=int(gp[0]) if gp.degree()>=0 else 0
    c1=int(gp[1]) if gp.degree()>=1 else 0
    fr=lambda ab: ((ab[0]+ab[1]*c0)%p,(ab[1]*c1)%p)
    sig=[idx[fr(k)] for k in jl]
    assert sorted(sig)==list(range(h))
    # split charpoly by P-eigenvalue: P B = B P check, then charpoly on (1±P)/2 spaces
    P=[[1 if sig[i]==j2 else 0 for j2 in range(h)] for i in range(h)]
    MB=fmpz_mat(B); MP=fmpz_mat(P)
    assert MP*MB==MB*MP, ("commute",p)
    # rational projectors: work with 2B(1±P)/2 = B±BP to stay integral: charpoly of B on ker(P-1):
    # eigen-split via charpoly of (B + t P) trick is messy; instead: basis of (P±1)-eigenspaces over Q
    def eigbasis(sign):
        # rows spanning the +/- eigenspace: for orbits {i,sig[i]}
        rows=[]; seen=set()
        for i in range(h):
            if i in seen: continue
            if sig[i]==i:
                if sign==1: rows.append([1 if j2==i else 0 for j2 in range(h)])
            else:
                seen.add(sig[i])
                v=[0]*h; v[i]=1; v[sig[i]]=sign
                rows.append(v)
            seen.add(i)
        return rows
    import numpy as _n
    out={}
    for sgn in (1,-1):
        rows=eigbasis(sgn)
        if not rows: out[sgn]=fmpz_poly([1]); continue
        # restriction of B to span(rows): solve R*B = C*R  =>  C = R B R^+ (rows are 0/±1, orthogonal-ish)
        # rows pairwise disjoint supports except... they ARE orthogonal (distinct orbit supports)
        R=rows; m=len(R)
        C=[[0]*m for _ in range(m)]
        # image of each basis vector under B^T? act by B on column vectors: v_k = e_i + sgn e_{sig i}
        # B acts on functions by (Bf)(i)=sum_j B[i][j] f(j): column action; represent vectors as columns
        cols=[[R[k][j2] for j2 in range(h)] for k in range(m)]
        for k in range(m):
            w=[sum(B[i][j2]*cols[k][j2] for j2 in range(h)) for i in range(h)]
            # express w in basis: coefficient on col l = w[i_l]/1 at the defining index
            for l in range(m):
                il=next(j2 for j2 in range(h) if cols[l][j2]==1)
                C[l][k]=w[il]
        out[sgn]=fmpz_mat(C).charpoly()
    return dict(h=h,B=B,W=W,sig=sig,cp_plus=str(out[1]),cp_minus=str(out[-1]))

if __name__=='__main__':
    import sys, json as _json, re as _re
    from fractions import Fraction as F
    ps=[int(a) for a in sys.argv[1:] if a.isdigit()] or \
       [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149]
    try: res=_json.load(open('w2data.json'))
    except Exception: res={}
    def dg(s):
        m=_re.findall(r'x\^(\d+)',s)
        return max(int(x) for x in m) if m else (1 if 'x' in s else 0)
    OK=True
    for p in ps:
        d=w2data(p)
        e1=jac(-1%p,p); N=nu(p)
        nu2=1+jac(-1%p,p); nu3=1+jac(-3%p,p)
        d2=F(p+1,12)-F(nu2,4)-F(nu3,3)
        a2=int((d2+1-F(N,2))/2); b2=int((d2-1+F(N,2))/2)
        okp = dg(d['cp_plus'])==1+b2 and dg(d['cp_minus'])==a2
        if not okp:
            OK=False
            print(f"p={p}: DIM MISMATCH plus {dg(d['cp_plus'])} vs {1+b2}, minus {dg(d['cp_minus'])} vs {a2}")
        else:
            print(f"p={p}: h={d['h']} dims OK (eps=-1:{b2}, eps=+1:{a2})")
        res[str(p)]=d
        _json.dump(res,open('w2data.json','w'))
    print("weight-2 side:", "ALL VALIDATED" if OK else "FAILURES")
