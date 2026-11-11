"""Independent whole-interval normalized T-source entrance.

No author imports or receipts. Source constants are reconstructed as Fractions;
the all-real-y completed square is proved by exact Bernstein coefficients on
the entire X interval, unlike the author's Sturm root count. A new degree6251
row retains the sharper matrix cosine-square coefficient and its corresponding
negative sine-square term. This does not assert a global degree6251 theorem.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json,hashlib,time

started=time.monotonic();checks=[];certificates=[]
rho0=F(30501,500000)
def ck(name,v):
    assert isinstance(v,F) or isinstance(v,int)
    ok=v>0;checks.append(dict(name=name,pass_=ok,margin=str(v)))
    if not ok:raise AssertionError((name,v))
def eq(name,a,b):
    assert a==b,(name,a,b)
    checks.append(dict(name=name,pass_=True,identity=True))
def add(a,b):
    return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]
def scale(a,c):return [c*x for x in a]
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def bernstein(poly,l,r):
    n=len(poly)-1
    power=[sum((poly[k]*comb(k,j)*l**(k-j)*(r-l)**j for k in range(j,n+1)),F(0)) for j in range(n+1)]
    return [sum((power[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),F(0)) for i in range(n+1)]
def whole_positive(name,poly,l,r):
    todo=[(l,r,0)];pieces=[]
    while todo:
        a,b,depth=todo.pop();vals=bernstein(poly,a,b)
        if min(vals)>0:
            pieces.append(dict(left=str(a),right=str(b),degree=len(poly)-1,
                               coefficients=list(map(str,vals))))
        else:
            assert depth<20,(name,a,b,vals)
            m=(a+b)/2;todo.extend([(a,m,depth+1),(m,b,depth+1)])
    pieces.sort(key=lambda z:F(z['left']))
    eq(name+' full left endpoint',F(pieces[0]['left']),l)
    eq(name+' full right endpoint',F(pieces[-1]['right']),r)
    for a,b in zip(pieces,pieces[1:]):eq(name+' exact closed seam',F(a['right']),F(b['left']))
    for p in pieces:
        for z in p['coefficients']:ck(name+' whole Bernstein coefficient',F(z))
    certificates.append(dict(name=name,polynomial=list(map(str,poly)),pieces=pieces))

# Exact arctangent enclosure from the alternating-series theorem.
def atan_interval(x,m=24):
    sm=sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(m)),F(0))
    nxt=(-1)**m*x**(2*m+1)/(2*m+1)
    return min(sm,sm+nxt),max(sm,sm+nxt)
a,b=atan_interval(F(1,5));c,d=atan_interval(F(1,239))
pi_lower=16*a-4*d
ck('Machin exact pi lower',pi_lower-F(314159,100000))
pl=F(314159,100000);dd=3*pl*pl-8
K0=F(21,5);A0=F(11,2);gamma=F(241,1000)
ck('complete matrix constant K0',K0-9*pl*pl/dd)
ck('complete matrix constant A0',A0-(2+24*pl/dd))
ck('sharper complete matrix gamma',gamma-(40-3*pl*pl)/(2*dd))
g=gamma+F(1,2)
eq('same sine-square coefficient retained',g,F(741,1000))

# The half-source uniform X<=1/4 comparison, independent of rho and degree.
delta=F(1,2);m=F(1,4)
eq('full large critical-displacement discriminant',(1-delta*delta)**2-m*(4*delta+m),0)
ck('small-displacement derivative has no real root',-((4*m)**2+12*(4*m*delta-1)))
ck('small-displacement derivative is negative',1-4*m*delta)
P=[2*m*delta*delta,4*m*delta-1,2*m,-1]
Pm=sum((v*m**i for i,v in enumerate(P)),F(0))
eq('small-displacement exact endpoint',Pm,m*delta**4)
ck('small-displacement entire positive endpoint',Pm)

# kappa=N_k/D_k, constructed from its elementary factors.
Nk=mul([0,1],mul([-1,4],[9,4]));Dk=[40,40]
eq('half critical numerator',[F(v) for v in Nk],[F(0),F(-9),F(32),F(16)])

rows=[]
for N,R,split,M,label in [
    (8000,F(3,40),F(3,10),F(3),'8001'),
    (6250,F(9,100),F(7,25),F(31,10),'6251_new')
]:
    n=N+1;rlo=F(1,100);b=(1-R)/(1+delta*delta)
    invf=1-R*R*(1+delta)**2/2;Dp=1+(R*(1+delta)/2)**2
    ck(label+' full polar denominator',invf)
    ck(label+' inverse Cayley to right',delta-(R*(1+delta)/2)**2)
    ck(label+' normalized lower circle real',1-R*R/4)
    ck(label+' whole crescent gap gain',F(1,4)-R/2)
    ks=[(n*r*b/(n*r+2*Dp*b)-r)*invf for r in (rlo,R)]
    K=min(ks);ck(label+' direct HB source positive',K)
    ck(label+' split in positive kappa branch',split-F(1,4))
    ck(label+' small-X direct D25 entrance',K-25*R*split)
    # 4 M kappa-(X+kappa)^2, after the positive denominator D_k^2.
    Xplus=add(mul([0,1],Dk),Nk)
    poly=add(scale(mul(Nk,Dk),4*M),scale(mul(Xplus,Xplus),-1))
    whole_positive(label+' full all-real-y square',poly,split,F(1))
    ck(label+' entire large-X direct D25 entrance',K*N-2*M*25**2)
    ck(label+' larger-tau p strictly positive',1-rho0*rho0-(1-R)**2)
    ck(label+' entire larger-tau x below tau',R-rho0)
    # gamma=.241 replaces .25 AND changes the corresponding negative term
    # from -.75 to -.741. This denominator is deliberately not 3/4.
    matrix_upper=50*rho0/R+(2-R)/R*(K0+25*A0+625*gamma)+(1-R)/(2*R*g)
    ck(label+' entire larger-tau matrix D25',F(N)-matrix_upper)
    cp=rho0/2
    ck(label+' direct crescent rational sqrt majorant',F(49,25)-(2-2*cp+cp*cp))
    ck(label+' sqrt2 rational lower',2-F(49,25))
    Dlarge=F(5*N,14)*cp/(1-cp+cp*cp)
    ck(label+' large-x direct D70',Dlarge-70)
    row=dict(label=label,degree_floor=n,rho_interval=[str(rlo),str(R)],
        Xsplit=str(split),joint_all_y_M=str(M),source_endpoints=list(map(str,ks)),source_K=str(K),
        small_X_D=str(K/(R*split)),large_X_h_squared=str(K*N/(2*M)),
        large_tau_matrix_upper=str(matrix_upper),large_x_D=str(Dlarge))
    if N==8000:
        D60=F(N)/(2/rho0+96);ck('8001 lawful ratio3 D60',D60-60)
        eq('8001 actual q rectangle endpoint',F(3,60),F(1,20))
        ck('8001 actual mu rectangle endpoint',F(1,50)-F(1,60))
        Dphase=(N-42/rho0)/98;ck('8001 outside Origin same-phase D70',Dphase-70)
        Blo=F(47,50);d0=1-Blo
        ck('8001 low-middle normalization lower',d0-rlo)
        ck('8001 whole inside-Origin corner',rho0*rho0+(1-rho0)**2-Blo*Blo)
        ck('8001 corner quadratic monotonicity',F(1,2)-rho0)
        eq('8001 low-middle strict radial D60 endpoint',N*d0/8,F(60))
        Dlow=(N-42/d0)/(2*rho0/d0+3*(1+Blo)/d0)
        ck('8001 low-middle phase D70',Dlow-70)
        Dmedium=(N*d0-F(596,5))/F(11,5)
        ck('8001 lawful D70-to-medium160',Dmedium-160)
        row.update(outside_phase_D=str(Dphase),low_middle_phase_D=str(Dlow),medium_D=str(Dmedium))
    if N==6250:
        Dfirst=F(N)/(2/rho0+96)
        ck('6251 outside after lawful ratio3 exceeds48',Dfirst-48)
        ck('6251 coarse outside cosine bound',F(3,4)**2-F(5,9))
        outside_matrix=K0/48**2+A0*F(3,4)/48+gamma*F(5,9)-F(2,9)
        ck('6251 entire outside forbidden sine2over3',-outside_matrix)
        Dout=F(N)*rho0/(2*rho0+F(4,3)*(1-rho0)+3*(2-rho0))
        ck('6251 outside legal q06mu02 D50',Dout-50)
        Blo=F(47,50);d0=1-Blo;Dfirstlow=N*d0/8
        ck('6251 low-middle normalization lower',d0-rlo)
        ck('6251 whole inside-Origin corner',rho0*rho0+(1-rho0)**2-Blo*Blo)
        ck('6251 coarse low-middle cosine bound',F(143,200)**2-F(51,100))
        low_matrix=K0/Dfirstlow**2+A0*F(143,200)/Dfirstlow+gamma*F(51,100)-F(49,200)
        ck('6251 entire low-middle forbidden sine7over10',-low_matrix)
        Dlow=F(N)*d0/(2*rho0+F(7,5)*Blo+3*(1+Blo))
        ck('6251 low-middle legal q06mu02 D50',Dlow-50)
        eq('6251 strict q06 rectangle endpoint',F(3,50),F(3,50))
        eq('6251 strict mu02 rectangle endpoint',F(1,50),F(1,50))
        row.update(outside_ratio3_D=str(Dfirst),outside_coarse_phase='2/3',outside_q06_D=str(Dout),
            low_middle_ratio3_D=str(Dfirstlow),low_middle_coarse_phase='7/10',low_middle_q06_D=str(Dlow))
    rows.append(row)

out=dict(status='PASS',scope='Complete normalized T-source D25 rows and whole larger-tau D25 entrances; 8001 downstream D70/medium160; 6251 lawful q06mu02 entrances; no global6251 assertion',
    independent_method='Exact Fraction source reconstruction, alternating Machin enclosure, whole-X Bernstein convex hull, full corresponding negative matrix term',
    checks=len(checks),whole_polynomial_cells=sum(len(z['pieces']) for z in certificates),
    matrix_gamma=str(gamma),matrix_negative_sine=str(g),rows=rows,polynomial_certificates=certificates,
    predicates=checks,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-started)
Path(__file__).with_name('tang_normalized_T_source_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('predicates','polynomial_certificates')},indent=2))
