"""T3-to-149 assembly: exact catalogue matching at every prime 11<=p<=149.
charpoly(graph) == (x-15) * SK * Yoshida * Mys^2 * N_p with SK,Yoshida built
from independent weight-2/weight-4 data; signs certified via one-modulus
projector traces (targets bounded by block degree)."""
import sys,json,re,itertools
sys.path.insert(0,'/home/claude/gs')
from flint import fmpz_mat, fmpz_poly, fmpq_poly, nmod_mat, nmod_poly
from fractions import Fraction as Fr

def parse_fp(s):
    s=s.strip()
    if s=='1': return fmpz_poly([1])
    co={}
    for term in s.split(' + '):
        term=term.strip()
        if '*x^' in term:
            c,k=term.split('*x^'); co[int(k)]=int(c.strip('()'))
        elif term.startswith('x^'):
            co[int(term[2:])]=1
        elif '*x' in term:
            c=term.split('*x')[0]; co[1]=int(c.strip('()'))
        elif term=='x':
            co[1]=1
        else:
            co[0]=int(term.strip('()'))
    d=max(co)
    return fmpz_poly([co.get(i,0) for i in range(d+1)])

def companion(f):
    d=f.degree(); a=[int(f[i]) for i in range(d+1)]
    assert a[d]==1
    C=[[0]*d for _ in range(d)]
    for i in range(1,d): C[i][i-1]=1
    for i in range(d): C[i][d-1]=-a[i]
    return C

def kron_sum_charpoly(g4,f2):
    dg,df=g4.degree(),f2.degree()
    Cg=companion(g4); Cf=companion(f2)
    n=dg*df
    K=[[0]*n for _ in range(n)]
    for i in range(df):
        for a in range(dg):
            for b in range(dg):
                K[i*dg+a][i*dg+b]+=Cg[a][b]
    for i in range(df):
        for j in range(df):
            for a in range(dg):
                K[i*dg+a][j*dg+a]+=2*Cf[i][j]
    return fmpz_mat(K).charpoly()

L0=(1<<61)-1  # Mersenne prime
def run(p,w2,w4,gr):
    M=gr['M']; n=len(M); sig=gr['sigma']; res=gr['res']
    cp=fmpz_mat(M).charpoly()
    # weight-2 orbits with signs
    plus=parse_fp(w2['cp_plus']); minus=parse_fp(w2['cp_minus'])
    q,r=divmod(fmpq_poly(plus),fmpq_poly(fmpz_poly([-3,1])))
    assert r==0
    plus_c=fmpz_poly([int((q[i]).p) for i in range(q.degree()+1)])
    f2orbs=[]
    for f,e in plus_c.factor()[1]: f2orbs += [(f,-1)]*e
    for f,e in minus.factor()[1]: f2orbs += [(f,+1)]*e
    # weight-4 orbits
    g4orbs=[(parse_fp(f),m,e) for f,m,e in w4['orbits']]
    xm6=fmpz_poly([-6,1])
    SK=fmpz_poly([1])
    for g,m,e in g4orbs: SK=SK*(g(xm6))**m
    Yos=fmpz_poly([1]); ypairs=[]
    for (f,e2) in f2orbs:
        for (g,m,e4) in g4orbs:
            if e2*e4==-1:
                blk=kron_sum_charpoly(g,f)**m
                Yos=Yos*blk; ypairs.append((blk,-e2))
    pred=fmpz_poly([-15,1])*SK*Yos
    Q,R=divmod(fmpq_poly(cp),fmpq_poly(pred))
    assert R==0, ("EXACT DIVISION FAILS",p)
    Qz=fmpz_poly([int(Q[i].p) for i in range(Q.degree()+1)])
    facQ={str(f):(f,e) for f,e in Qz.factor()[1]}
    # ---- one-modulus isotypic traces on the graph ----
    Mm=nmod_mat([[x% L0 for x in row] for row in M],L0)
    X=nmod_mat(n,n,L0)
    for i in range(n): X[i,i]=1
    tr=[]
    for k in range(cp.degree()):
        tr.append(sum(int(X[sig[i],i]) for i in range(n))%L0)
        if k<cp.degree()-1: X=Mm*X
    cpm=nmod_poly([int(cp[i])%L0 for i in range(cp.degree()+1)],L0)
    def iso_trace(f,e):
        fm=nmod_poly([int(f[i])%L0 for i in range(f.degree()+1)],L0)
        fe=fm**e
        cof=cpm//fe
        g0,s0,t0=cof.xgcd(fe)
        assert g0.degree()==0
        h=(cof*s0*int(pow(int(g0[0]),L0-2,L0)))%cpm
        v=sum(int(h[k2])*tr[k2] for k2 in range(h.degree()+1))%L0
        return v-L0 if v>L0//2 else v
    # signed lift-multiplicity per irreducible (Eis + SK + Yoshida)
    lift_signed={}; lift_mult={}
    def acc(poly,sign,mult=1):
        for f,e in poly.factor()[1]:
            k=str(f)
            lift_signed[k]=lift_signed.get(k,0)+sign*e*mult
            lift_mult[k]=lift_mult.get(k,0)+e*mult
    acc(fmpz_poly([-15,1]),+1)
    for g,m,e in g4orbs: acc(g(xm6),+1,m)
    for blk,s in ypairs: acc(blk,s)
    # ---- trace-decoded split of Q into Mys^2 * N ----
    mys=fmpz_poly([1]); Np=fmpz_poly([1]); ok_signs=True
    for key,(f,eq) in facQ.items():
        etot=eq+lift_mult.get(key,0)
        tau=iso_trace(f,etot)
        ms=tau//f.degree()
        assert ms*f.degree()==tau, ("trace not multiple of degree",p,key)
        b=ms-lift_signed.get(key,0)      # N-multiplicity (R-sign +1, eps_F=-1 for p<167)
        assert 0<=b<=eq and (eq-b)%2==0, ("split fails",p,key,eq,b)
        a=(eq-b)//2
        mys=mys*f**a; Np=Np*f**b
    # blocks appearing only in the lift part: verify their signs too
    for key,(sm) in list(lift_signed.items()):
        if key in facQ: continue
        fl=[f for f,e in cp.factor()[1] if str(f)==key]
        if not fl: continue
        f=fl[0]; e=lift_mult[key]
        tau=iso_trace(f,e)
        if tau!=sm*f.degree(): ok_signs=False; print(f"   SIGN MISMATCH p={p} {key}: got {tau} want {sm*f.degree()}")
    assert mys.degree()==res['degMys'], ("Mys deg",p,mys.degree(),res['degMys'])
    assert Np.degree()==res['d'], ("Np deg",p,Np.degree(),res['d'])
    return dict(p=p,Np=str(Np),ok=True,ok_signs=ok_signs,
                sk=SK.degree(),yos=Yos.degree(),mys=mys.degree())

if __name__=='__main__':
    w2=json.load(open('w2data.json')); w4=json.load(open('w4data_part1.json'))
    try: out=json.load(open('assembly_results.json'))
    except Exception: out={}
    import sys as _s
    ps=[int(x) for x in _s.argv[1:] if x.isdigit()] or [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149]
    for p in ps:
        try: gr=json.load(open(f'h2_result_{p}.json'))
        except FileNotFoundError:
            print(f'p={p}: THIEU h2_result_{p}.json — chay buoc 1 truoc:  python3 richelot.py {p}'); continue
        if str(p) not in w2: print(f'p={p}: THIEU du lieu weight-2 — chay:  python3 w2side.py {p}'); continue
        if str(p) not in w4: print(f'p={p}: THIEU du lieu weight-4 — chay:  python3 msym.py {p}'); continue
        r=run(p,w2[str(p)],w4[str(p)],gr)
        print(f"p={p}: catalogue EXACT  SK={r['sk']} Yos={r['yos']} Mys={r['mys']}  N_p={r['Np']}  signs={'OK' if r['ok_signs'] else 'FAIL'}")
        sys.stdout.flush()
        out[str(p)]=r
        json.dump(out,open('assembly_results.json','w'))
    json.dump(out,open('assembly_results.json','w'))
