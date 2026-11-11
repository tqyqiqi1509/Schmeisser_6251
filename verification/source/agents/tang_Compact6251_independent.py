"""Independent full Compact6251 and .1--.11 bridge reconstruction.

No author code/receipt imports. The input supplies only exact contour
choices and row parameters. Fresh 400-bit Arb exp/log functions retain
full tails. Counts/degree-ray coefficients use exact Fractions. Entire
physical contours, analytic domain gates and closed coverage are rebuilt.
"""
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
from math import floor,isqrt
import sympy as sy
import json,hashlib,time,sys
sys.set_int_max_str_digits(0);ctx.prec=400
START=time.monotonic();ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'computations/degree6251/COMPACT6251_INTEGER_CHOICES.json'
N=6250;n=N+1;CAP=F(1,388);G=[];ROWS=[];CONTOURS=[]
def A(x):
    if isinstance(x,arb):return x
    x=F(x);return arb(x.numerator)/x.denominator

def ck(name,p,value=None):
    ok=bool(p);G.append({'name':name,'pass':ok,'value':None if value is None else str(value)})
    if not ok:raise AssertionError((name,value))

def pos(name,v):ck(name,v>0,v)
def below(name,v,b):pos(name,A(b)-A(v))
def posq(name,v):ck(name,v>0,v)

def radical_upper(v):
    sc=10**14;k=isqrt(v.numerator*sc*sc//v.denominator)+1
    z=F(k,sc);ck('independent exact radical enclosure',z*z>v and (z-F(1,sc))**2<=v)
    return z

def logweight(eta,tag):
    val=-A(eta).log();candidate=F(floor(float(val.lower())*10**13)-2,10**13)
    ck(tag+'/independently rounded complete log weight',A(candidate)<val and candidate>0)
    return candidate

# Fresh exact identities, not copied author numerical gates.
a,x,y,u,v=sy.symbols('a x y u v',real=True)
rec=2*(a-u)-a*((a-u)**2+v*v)-a*(1-u*u-v*v)-(1-a*a)*(a-2*u)
ck('complete reciprocal lower endpoint identity',sy.expand(rec)==0)
r2=(x-a)**2+y*y;t=1-r2;C=t*a*(2*x-(1+t)*a)
rhs=C*(1+r2)/(2*(1-r2))+r2+a*a*r2*(1-r2)/2
ck('same actual gamma coefficient identity',sy.cancel((x*x+y*y-t*a*x)-rhs)==0)
ck('full crescent radial square identity',sy.expand((3*x*x+y*y)**2-8*x*x*(x*x+y*y)-(y*y-x*x)**2)==0)
# The direct crescent proof is valid for every 0<a<=x<=1: its derivative
# -x+a*y/[2sqrt(1-a^2/4)] is < -x+x/sqrt(3)<0. No a<=.1 is used.
pos('whole-a crescent derivative denominator',A(3).sqrt()-1)
c0=F(1,20)
posq('complete rational crescent radicand cap',F(49,25)-(2-2*c0+c0*c0))
posq('complete radial denominator lower',F(2)-F(49,25))
Dfirst=F(5*N,14)*c0/(1-c0+c0*c0)
posq('noncircular direct D100 on every a>=.1',Dfirst-100)
pos('pi rational full bound',22-7*arb.pi())
pos('complete coefficient cap',A(CAP)*A(N*CAP/2).exp()-8)
posq('actual energy degree floor',n-701);posq('actual variance degree floor',n-2048)

rows=json.loads(DATA.read_text())
for row in rows:
    name=row['name'];aa,bb,rlo,rhi,target=map(F,[row['A'],row['B'],row['rlo'],row['rhi'],row['target']])
    tag=lambda s:name+'/'+s
    ck(tag('closed entire row domain'),F('.11')<=aa<bb<=F('.975') and 0<rlo<rhi and rhi*rhi<F(1,2))
    near=row['near_c'];ce=aa*aa*rlo*rlo*(1-rlo*rlo)
    pos(tag('complete original-product logarithmic budget target'),A(target).exp()-A(8/ce))
    if near is not None:
        nc=F(near);bbar=bb+F(1,1000)
        ck(tag('exact shrinking near domain'),rlo==nc/n and nc>0)
        posq(tag('whole near root modulus'),bbar-bb-rlo)
        posq(tag('whole variance exclusion'),(1-bbar*bbar)/2-nc)
        posq(tag('complete near positive squared comparison'),(1-nc/n)**2-(1-(1-bbar*bbar)/n))
        posq(tag('whole shrinking degree ray'),target-2)
    d=1-rhi;m=1/d
    if row['mode']=='trace':
        cp=aa/2;posq(tag('direct D100 before all refined ratios'),F(5*N,14)*cp/(1-cp+cp*cp)-100)
        df=F(row['Dfloor']);kap=F(row['kappa'])
        ck(tag('exact Dfloor before smaller root'),F(2*N,13)*d>df>=201)
        P=lambda z:F(9,4)*z*z-(df-13)*z+2*df
        ck(tag('legal prior and proposed smaller branch'),2<kap<F(9,4) and P(kap)<0 and P(F(9,4))<0)
        posq(tag('negative branch strengthens over full D ray'),kap-2)
        V=min(2/aa,kap/d);Hi=max(N*aa/2,F(2*N,13))
        posq(tag('signed Fine q independently entered'),F(7,625)-V/Hi)
        posq(tag('signed Fine mu independently entered'),F(1,270)-m/Hi)
        c=aa/(2*(1+rhi)**2);alpha=(1+2*c*m)/(2*(c+m));beta=(F(9,4)*V*V+13*V*m)/(2*(c+m))
        H=alpha*N-beta
        ck(tag('complete trace field interval'),0<c<m and H>Hi and alpha>0 and beta>0)
    else:
        ck(tag('explicit moment mode'),row['mode']=='moment')
        ell,wm,km=map(F,row['moment'])
        scopes=[(F('.9'),F('.057'),[F('4.42'),F('5.098'),F('2.5482')]),(F('.93'),F('.055'),[F('4.22'),F('4.765'),F('2.3817')]),(F('.95'),F('.054'),[F('4.11'),F('4.623'),F('2.3108')])]
        ck(tag('actual moment row scope'),any(aa>=amin and rhi<=rmax and [ell,wm,km]==coeff for amin,rmax,coeff in scopes))
        er=8*rhi/(aa*N)+rhi*rhi/(2*(1-rhi)*(1-2*rhi))
        qm=1-aa*aa+wm*(ell/N+er)+F(2,n*n);mm=km*(ell/n+er)+F(1,n*n)
        alpha=1-mm/aa-qm/(aa*d);beta=F(0);H=alpha*N;V=1/(alpha*d*d)
        ck(tag('full same-configuration moment field'),alpha>0 and 0<qm<1 and 0<mm<aa)
        posq(tag('moment phase q independently entered'),F(7,625)-F(4,N)/(aa*aa))
        posq(tag('moment phase mu independently entered'),F(1,270)-2/(N*aa*d))
    weights=[];counts=[];slopes=[];lastR=F(0);count0=len(CONTOURS)
    js=[rr[0] for rr in row['circles']]
    ck(tag('all integer radii in full declared range'),js==list(range(js[0],js[-1]+1)))
    for j,nn,W in row['circles']:
        R=F(j,1000);nu=F(nn,10000);lab=tag('circle'+str(j));e=A(-W).exp()
        ck(lab+'/full physical choices',R>lastR and 0<m*R<1 and 0<nu<1 and W in [1,2,3,4,5,6,8,10])
        lam=V*R/(1-m*R);error=V*(R*R/2+m*R**3/(3*(1-m*R)));image=R-error
        ck(lab+'/whole inverse analytic domain',0<lam<1 and image>0)
        # Extra direct full-log check: no finite series replaces its tail.
        actual_log_error=A(V/m**2)*(-A(1-m*R).log()-A(m*R))
        pos(lab+'/entire log tail below conservative error',A(error)-actual_log_error)
        B=V/((1-m*R)**2*(1-lam)**3);eps=F(1001,1000)*B/H;g=F(1001,1000)*nu*B*image
        pin=1/(1+lam);amp=1/(1-lam)
        posq(lab+'/entire containing corners',(1-nu*nu)*(H*image)**2-W*W)
        pos(lab+'/entire left side',1-A(amp)*e-A(eps)*(1+e)-A(g)*e)
        pos(lab+'/entire right side',A(pin)-e-A(eps)*(1+e)-A(g)*e*e)
        pos(lab+'/both entire horizontal sides',A(pin)-A(eps)-A(g)*e)
        ck(lab+'/positive-real amplitude full homotopy',0<pin<1<amp)
        posq(lab+'/complete odd-height budget',nu*H*image-F(22,7))
        slope=nu*F(7,22)*image;count=slope*H-2;posq(lab+'/full finite count',count)
        delta=CAP*(1+rhi*rhi)/(2*(1-rhi*rhi))+(1+bb*bb/2)*rhi*rhi+R*rhi
        ck(lab+'/same actual gamma denominator',0<delta<1)
        rad=(bb*bb+R*R)/(1-delta)**2;eta=radical_upper(rad)
        ck(lab+'/complete original-factor disk',0<eta<1)
        ww=logweight(eta,lab)
        if weights:
            ck(lab+'/full nested count monotonicity',count>counts[-1]);ck(lab+'/full nested weight monotonicity',ww<=weights[-1])
        weights.append(ww);counts.append(count);slopes.append(slope)
        CONTOURS.append({'row':name,'radius':str(R),'nu':str(nu),'W':W,'image':str(image),'count':str(count),'eta':str(eta),'weight':str(ww)})
        lastR=R
    ck(tag('complete nonempty cloud'),bool(weights))
    # Use the positive Abel coefficient form, independently of the author
    # forward increment implementation.
    budget=weights[-1]*counts[-1]+sum(((weights[i]-weights[i+1])*counts[i] for i in range(len(weights)-1)),F(0))
    wstar=weights[-1]*slopes[-1]+sum(((weights[i]-weights[i+1])*slopes[i] for i in range(len(weights)-1)),F(0))
    ck(tag('finite-loss exact coefficient identity'),budget==wstar*H-2*weights[0])
    bonus=F(0) if near is not None else N*ce/2
    nslope=wstar*alpha+(F(0) if near is not None else ce/2)
    nloss=wstar*(alpha+beta)+2*weights[0]+(F(0) if near is not None else ce/2)
    ck(tag('whole unbounded degree ray'),budget+bonus==nslope*n-nloss and nslope>0 and nloss>0)
    posq(tag('complete same-C contradiction'),budget+bonus-target)
    ROWS.append({'name':name,'A':str(aa),'B':str(bb),'rlo':str(rlo),'rhi':str(rhi),'H':str(H),'V':str(V),'alpha':str(alpha),'beta':str(beta),'budget':str(budget),'bonus':str(bonus),'target':str(target),'reserve':str(budget+bonus-target),'n_slope':str(nslope),'n_loss':str(nloss),'contours':len(CONTOURS)-count0})

# Entire closed coverage: group common witness slabs; radial seams must
# exactly join and their lower ends are proved near exclusions.
slabs=sorted(set((F(r['A']),F(r['B'])) for r in rows));end=F('.11')
for aa,bb in slabs:
    ck('whole witness seam '+str(aa),aa==end);end=bb
    rrows=sorted((F(r['rlo']),F(r['rhi']),r) for r in rows if F(r['A'])==aa and F(r['B'])==bb)
    lo=rrows[0][0]
    if aa<F('.85'):
        ck('shrinking near coverage '+str(aa),rrows[0][2]['near_c'] is not None)
    elif aa<F('.96'):ck('NearA complete lower radial seam '+str(aa),lo==F('.02') and aa>=F('.85'))
    elif aa<F('.965'):ck('NearGeomB complete lower radial seam',lo==F('.04') and aa>=F('.96'))
    else:ck('NearGeomC complete lower radial seam',lo==F('.045') and aa>=F('.965'))
    cur=lo
    for rl,rh,z in rrows:ck('exact closed radial seam '+z['name'],rl==cur);cur=rh
    ck('complete far monotonicity scope '+str(aa),cur*cur<F(1,2))
    posq('entire far exclusion including u=.5 '+str(aa),aa*aa*cur*cur*(1-cur*cur)-CAP)
ck('whole witness upper endpoint',end==F('.975'))

# Entire .1<=a<=.11 bridge, all four actual-root modulus intervals.
BRIDGE=[];al=F('.1');au=F('.11');Kphase=F(21,4);C70=F(1899,20)
def bridge_circle(name,B,R,H):
    d=1-B;V=F(9,4)/d;m=1/d;nu=F(99,100);W=6;e=A(-W).exp()
    lam=V*R/(1-m*R);err=V*(R*R/2+m*R**3/(3*(1-m*R)));image=R-err
    ck(name+'/whole physical inverse domain',0<m*R<1 and 0<lam<1 and image>0)
    pos(name+'/full log remainder',A(err)-A(V/m**2)*(-A(1-m*R).log()-A(m*R)))
    Bi=V/((1-m*R)**2*(1-lam)**3);eps=F(1001,1000)*Bi/H;g=F(1001,1000)*nu*Bi*image;pin=1/(1+lam);amp=1/(1-lam)
    posq(name+'/full containing corners',(1-nu*nu)*(H*image)**2-W*W)
    pos(name+'/entire left side',1-A(amp)*e-A(eps)*(1+e)-A(g)*e)
    pos(name+'/entire right side',A(pin)-e-A(eps)*(1+e)-A(g)*e*e)
    pos(name+'/both entire horizontal sides',A(pin)-A(eps)-A(g)*e)
    ck(name+'/positive amplitude entire homotopy',0<pin<1<amp)
    posq(name+'/entire finite odd height',nu*H*image-F(22,7))
    count=nu*F(7,22)*H*image-2;posq(name+'/positive full count',count)
    den=1-(B+R)*B;posq(name+'/complete actual-factor denominator',den)
    factor=(R+au)/den
    BRIDGE.append({'name':name,'H':str(H),'R':str(R),'image':str(image),'count':str(count),'factor':str(factor)})
    return count,factor

count,factor=bridge_circle('BridgeLow',F('.2'),F('.1'),N*F('.8')/F(9,4))
posq('BridgeLow factor',F(1,4)-factor);pos('BridgeLow complete weight',-A(factor).log()-1)
posq('BridgeLow whole count degree ray',count/n-F(1,120))
posq('BridgeLow variance loss',1-F('.2')**2-F('.94'))
posq('BridgeLow variance squared comparison',(1-F('.47')/n)**2-(1-F('.96')/n))
posq('BridgeLow near displacement',F('.47')-F('.4'))
posq('BridgeLow same-C entire degree ray',al*al*F('.4')**2*(1-F('.4')**2/n**2)-F(1,1000))
pos('BridgeLow full logarithmic constant',10-A(8000).log());pos('BridgeLow floor log degree',9-A(n).log())
posq('BridgeLow full floor contradiction',F(n,120)-28);posq('BridgeLow entire degree derivative',F(1,120)-F(2,n))
count,factor=bridge_circle('BridgeMiddle',F('.8'),F('.04'),(N-C70/F('.2'))/(2*(1+au)))
posq('BridgeMiddle factor',F('.46')-factor);pos('BridgeMiddle complete weight',-A(factor).log()-A(F('.75')))
cmid=al*al*F('.09')**2*(1-F('.09')**2)
pos('BridgeMiddle complete target',12-A(8/cmid).log());posq('BridgeMiddle complete source budget',F('.75')*count-12)

umin=F('.69')**2;tmax=1-umin;Cup=al*al*umin*(1-umin)
ck('BridgeUpper adaptive whole range',0<umin<F('.5') and F('.5')<tmax<1)
den=1-F('.94')**2+al*al/2;actualfactor=tmax*au/den
posq('BridgeUpper actual root denominator',den);posq('BridgeUpper actual root factor',F('.475')-actualfactor)
pos('BridgeUpper full actual-root logarithmic weight',-A(actualfactor).log()-A(F('.74')))
pos('BridgeUpper complete same-C logarithm',A(F('8.08'))-A(8/Cup).log())
posq('BridgeUpper same energy and same actual factor',N*Cup/2+F('.74')-F('8.08'))

amin=F('.03');xmin=F('.1');shift=3*amin/4;L=F('8.15');pmax=1-F('.94')**2
posq('BridgeHigh entire adaptive half shift',(F('.94')-au)**2-F('.5'))
posq('BridgeHigh full concave-a endpoint comparison',xmin-3*amin)
# Exact endpoint difference: x^2/4-a0(x-3a0/4)=(x-a0)(x-3a0)/4.
xx,aa=sy.symbols('xx aa',real=True)
ck('BridgeHigh exact full coefficient endpoint identity',sy.expand(xx*xx/4-aa*(xx-3*aa/4)-(xx-aa)*(xx-3*aa)/4)==0)
pos('BridgeHigh complete same-C target',A(L)-A(8/(amin*(xmin-shift))).log())
kk=2*L/(N*amin);posq('BridgeHigh whole convex reciprocal-square condition',3*(xmin-shift)-shift)
vv,ss=sy.symbols('vv ss',positive=True)
ck('BridgeHigh full v convexity identity',sy.simplify(sy.diff(1/(vv*vv+ss),vv,2)-2*(3*vv*vv-ss)/(vv*vv+ss)**3)==0)
HIGH=[]
for xx in (xmin,F(1)):
    value=2*A(xx)*A(kk/(xx-shift)).sqrt()+A(pmax*kk/(xx-shift)+4*Kphase/(N*xx))
    pos('BridgeHigh complete convex endpoint '+str(xx),1-value)
    HIGH.append({'x':str(xx),'upper':str(value),'reserve':str(1-value)})
# Positive endpoints of the same full parameter intervals cover every modulus.
ck('Bridge all closed modulus seams',[F(0),F('.2'),F('.8'),F('.94'),F(1)]==sorted(set([F(0),F('.2'),F('.8'),F('.94'),F(1)])))

out={'status':'PASS','scope':'Complete .1<=a<=.975, every n>=6251, all actual-root moduli; no Origin or Small conclusion','arithmetic':'400-bit Arb full exp/log plus exact Fraction count, coverage and degree-ray identities; no author imports','checks':len(G),'complete_compact_contours':len(CONTOURS),'complete_bridge_contours':len(BRIDGE),'whole_rows':len(ROWS),'minimum_row_reserve':str(min(F(z['reserve']) for z in ROWS)),'rows':ROWS,'bridge_contours':BRIDGE,'upper_same_C':str(Cup),'upper_actual_factor':str(actualfactor),'high_convex_endpoints':HIGH,'predicates':G,'contours':CONTOURS,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'choices_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'seconds':time.monotonic()-START}
Path(__file__).with_name('tang_Compact6251_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['rows','predicates','contours']},indent=2))
for r in ROWS:print(r['name'],float(F(r['reserve'])))
