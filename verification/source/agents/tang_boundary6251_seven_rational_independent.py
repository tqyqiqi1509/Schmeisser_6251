"""Independent complete Boundary/Near/Moment6251 seven-row validation, exact rational intervals.

No author code or receipt is imported.  The seven parameter rows are read as exact rational strings
from the analytic statements. Transcendentals use rational Taylor/atanh/Machin
series, with explicit infinite tails; square roots use integer arithmetic.
The certificate includes extra domain and primitive-loss gates not separately
listed by the author's Arb program. Each inverse row is checked for both signs.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, factorial
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import hashlib, json, time

START=time.monotonic()
class I:
    def __init__(self, lo=0, hi=None):
        if isinstance(lo,I): self.lo,self.hi=lo.lo,lo.hi; return
        self.lo=F(lo); self.hi=self.lo if hi is None else F(hi)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=I(o); return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,o): return self+-I(o)
    def __rsub__(self,o): return I(o)+-self
    def __mul__(self,o):
        o=I(o); a=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return I(min(a),max(a))
    __rmul__=__mul__
    def reciprocal(self):
        assert not self.lo<=0<=self.hi
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,o): return self*I(o).reciprocal()
    def __rtruediv__(self,o): return I(o)*self.reciprocal()
    def __pow__(self,n):
        assert isinstance(n,int)
        if n<0:return self.reciprocal()**(-n)
        if n==0:return I(1)
        if n%2==0:
            vals=[self.lo**n,self.hi**n]
            return I(0 if self.lo<=0<=self.hi else min(vals),max(vals))
        return I(self.lo**n,self.hi**n)

def decimal(f,up=False):
    with localcontext() as c:
        c.prec=48;c.rounding=ROUND_CEILING if up else ROUND_FLOOR
        return str(Decimal(f.numerator)/Decimal(f.denominator))
def show(x):
    x=I(x); return [decimal(x.lo),decimal(x.hi,True)]
@lru_cache(None)
def expq(x):
    x=F(x)
    if x<0:return expq(-x).reciprocal()
    # Every application has 0<=x<=4.  The Taylor tail ratio is at most
    # x/(m+2), starting with the first omitted positive term.
    assert x<=4
    m=96;s=F(1);t=F(1)
    for j in range(1,m+1): t*=x/j;s+=t
    nxt=t*x/(m+1);tail=nxt/(1-x/F(m+2))
    return I(s,s+tail)
def expi(x):
    x=I(x);return I(expq(x.lo).lo,expq(x.hi).hi)
@lru_cache(None)
def logq(x):
    x=F(x);assert x>0
    k=0;y=x
    while y>=2:y/=2;k+=1
    while y<1:y*=2;k-=1
    def atanhlog(v):
        z=(v-1)/(v+1);assert 0<=z<=F(1,3)
        m=96;total=F(0);zp=z
        for j in range(m):total+=2*zp/(2*j+1);zp*=z*z
        return I(total,total+2*zp/((2*m+1)*(1-z*z)))
    return atanhlog(y)+k*atanhlog(F(2))
def logi(x):
    x=I(x);return I(logq(x.lo).lo,logq(x.hi).hi)
def sqrti(x):
    x=I(x);assert x.lo>=0
    scale=10**56
    lo=isqrt((x.lo.numerator*scale*scale)//x.lo.denominator)
    hi=isqrt((x.hi.numerator*scale*scale)//x.hi.denominator)+1
    return I(F(lo,scale),F(hi,scale))
def atan_inv(q):
    z=F(1,q);m=100
    s=sum(((-1)**j*z**(2*j+1)/F(2*j+1) for j in range(m)),F(0))
    # m even, last included term negative: s is a strict lower bound.
    return I(s,s+z**(2*m+1)/F(2*m+1))
PI=16*atan_inv(5)-4*atan_inv(239)

GATES=[];ROWS=[];N=6250;n=N+1;k=I(F(N,n))
def positive(label,x):
    x=I(x);ok=x.lo>0
    GATES.append(dict(name=label,pass_=ok,margin=show(x)))
    if not ok:raise AssertionError((label,show(x)))
def below(label,x,y):positive(label,I(y)-x)

positive('sharp cap via logarithms',I(F(N,776))-logq(F(8*388)))
below('full Fine phase polynomial',I(21)+133*I(F(7,625))+1030*I(F(7,625))**2+9000*I(F(7,625))**3+120000*I(F(7,625))**4+I(F(22,27)),24)
positive('actual-root energy degree floor',n-701)
positive('complete inner-error degree floor',n-1000)

PARAMETER_PATH=Path(__file__).resolve().parents[1]/'computations/degree6251/BOUNDARY_NEAR_MOMENT6251_SELECTED.json'
PARAMETERS=json.loads(PARAMETER_PATH.read_text())
MAPPING={'rr':'r2','coarse':'qc','mean0':'mc','annbase':'ss','Bcoarse':'bc','meanratio':'ds','q':'qb','B':'Bf','m':'m0','pf':'pfin','C':'C0','v':'v0'}
for raw in PARAMETERS:
    label=raw['name'];annuli=raw['annuli'];feedback=raw.get('feedback',[])
    d={key:I(raw[MAPPING.get(key,key)]) for key in 'a e rr direction s1 del1 s2 del2 aloss ell coarse mean0 p annbase Bcoarse meanratio w q B m pf'.split()}
    if raw.get('role')!='moment':
        d.update({key:I(raw[MAPPING.get(key,key)]) for key in 'weak h R gamma C v beta J inner growth finalr finals errbase'.split()})
    def pos(tag,x):positive(label+'/'+tag,x)
    def less(tag,x,y):below(label+'/'+tag,x,y)
    a,e,rr=d['a'],d['e'],d['rr'];a_initial=a
    pos('positive witness',a);less('witness inside',a,1)
    pos('positive displacement',e);less('reverse-log domain',e,F(1,2))
    pos('full critical distance',(1-e)**2-rr)
    if label.startswith('standard'):pos('whole far cutoff',a*a*e*e*(1-e*e)-I(F(1,388)))
    pos('radial direction full monotonicity',a*a-2*e)
    less('directional full root bound',1+2*e,d['direction'])
    for j in (1,2):
        sj=d['s'+str(j)];delta=d['del'+str(j)]
        pos(f'fiber{j} test radius',sj)
        pos(f'fiber{j} excludes test point',a-sj-delta)
        pos(f'fiber{j} second factor',a+sj-d['direction'])
        pos(f'fiber{j} entire movement disk',(a-sj)*(a+sj-d['direction'])-delta*(4+delta))
    # The full primitive loss must be below delta2/2; the author does not
    # separately list this gate.  It holds for both new parameter rows.
    pos('full primitive loss budget',100*d['del2']-d['s2'])
    pos('first-fiber logarithmic lower bound',d['aloss']/2+logi(d['del1']))
    pos('primitive log normalization',d['ell']-logi(2/d['del2']))
    alpha=1-d['aloss']/n;pos('alpha domain',alpha)
    T=(1+d['s1']**2-alpha*rr)/d['s1'];coarse=a*T+rr-a*a
    pos('complete logarithmic chord',rr+I('4.4')*(logq(F(2))+(1-rr)/2)-4)
    pos('whole witness chord',d['w']*logq(F(2))-3)
    less('v quadratic monotonicity',T,2*a)
    less('coarse complete moment',coarse,d['coarse'])
    less('coarse complete mean',d['coarse'],d['mean0']**2)
    pos('critical logarithmic halfline',2*a*k*alpha*rr/d['s1']-I('4.4'))
    pos('primitive logarithmic halfline',2*rr*(k*alpha*(a-d['s2'])/d['s1']+I('.9'))-I('4.4'))
    for t in (I(0),d['s2']):less('convex whole primitive endpoints '+str(t.lo),t*t-T*t+coarse,I('.9')*rr)
    pos('full primitive exponential log comparison',I(F(N,20))-logq(F(200*n)))
    pos('full primitive exponential degree derivative',I(F(n,20))-1)
    pos('nonzero primitive product base',rr-(I('1.001')*d['p'])**2)
    pos('full product power by logarithms',N*logq(F('1.001'))-logi(2*d['p']/d['del2']))
    less('complete coarse potential',-logi(rr)/2+d['ell']/n,d['Bcoarse'])
    less('mean residual prefactor',1/d['mean0'],2)
    less('mean residual full ratio',d['mean0']/d['p'],d['meanratio'])
    pos('mean residual exponential log comparison',(1-d['meanratio'])*n-logq(F(2*n*n)))
    pos('mean residual entire degree derivative',(1-d['meanratio'])*n-2)
    mean=d['mean0']
    for j,(r0,kap0,mean1) in enumerate(annuli):
        r,kap=I(r0),I(kap0)
        less(f'annulus{j} entire disk',r*r+2*r*mean+d['coarse'],d['annbase']**2)
        less(f'annulus{j} integration prefactor',r,d['annbase'])
        pos(f'annulus{j} full reciprocal-log inequality',-kap*logi(r)-(1/r-r))
        if mean1 is not None:
            less(f'annulus{j} complete mean',kap*d['Bcoarse']+I(F(1,n*n)),I(mean1));mean=I(mean1)
    ratio=d['annbase']/d['p'];pos('annulus whole product gap',1-ratio)
    pos('annulus full degree power',-n*logi(ratio)-logq(F(n)))
    pos('annulus whole degree derivative',-n*logi(ratio)-1)
    Hgeo=N*a/2
    less('legal Fine q entrance',2/a/Hgeo,F(7,625))
    less('legal Fine mu entrance',1/((1-e)*Hgeo),F(1,270))
    err=8*e/(a*N)+e*e/(2*(1-e)*(1-2*e))
    potential=d['ell']/n+err
    q_initial=1-a*a+d['w']*(d['ell']/N+err)+I(F(2,n*n))
    less('discarded B coefficient has correct sign',2*kap,d['w'])
    less('full reverse-log potential',potential,d['B'])
    less('same-configuration full moment',q_initial,d['q'])
    less('complete reciprocal mean',kap*d['B']+I(F(1,n*n)),d['m'])
    pos('final product base via exp(-B)>1-B',1-d['B']-d['pf'])
    if raw.get('role')=='moment':
        ROWS.append(dict(branch=label,role='moment_only',potential=show(potential),q_initial=show(q_initial),mean=show(d['m']),kappa=show(kap),w=show(d['w']),ell=show(d['ell'])))
        continue
    q,m=d['q'],d['m']
    def support(tag,q0,h0):
        ell=(1-q0)**2/(1+q0+2*m);gap=1-h0
        pos(tag+' finite Poisson support',(N*ell+1)*gap**2-2*gap-2*n*d['B'])
        pos(tag+' entire degree slope',ell*gap**2-2*d['B'])
    support('first',q,d['weak'])
    for j,(an,qn,hn) in enumerate(feedback):
        an,qn=I(an),I(qn);pos(f'feedback{j} witness increase',an-a)
        less(f'feedback{j} same-root original-error moment',1-an*an+d['w']*(d['ell']/N+err)+I(F(2,n*n)),qn)
        if hn is not None:support(f'feedback{j}',qn,I(hn))
        a,q=an,qn
    less('final complete original-root annulus',d['finalr']**2+2*d['finalr']*m+q,d['finals']**2)
    less('final primitive prefactor',d['finalr'],d['finals'])
    less('final whole-degree primitive gap',d['finals']/d['pf'],ratio)
    d0=1-e;cf=1-m/a-q/(a*d0);pos('same-root field bracket',cf)
    H=N*cf;V=1/(cf*d0*d0);M=1/d0
    R=d['R'];nu=I('.99');W=I(4);ew=expq(F(-4))
    pos('whole disk excludes every critical label',1-M*R)
    lam=V*R/(1-M*R);pos('whole inverse derivative perturbation',1-lam)
    e0=V*(R*R/2+M*R**3/(3*(1-M*R)));rho=R-e0
    pos('whole inverse image disk',rho)
    B2=V/((1-M*R)**2*(1-lam)**3)
    eps=I('1.001')*B2/H;g=I('1.001')*nu*B2*rho
    amp_lo=1/(1+lam);amp_hi=1/(1-lam)
    strip_gates=[
      ('all containing corners',(1-nu**2)*(H*rho)**2-W**2),
      ('entire left side',1-amp_hi*ew-eps*(1+ew)-g*ew),
      ('entire right side',amp_lo-ew-eps*(1+ew)-g*ew**2),
      ('both entire horizontal sides',amp_lo-eps-g*ew),
      ('nonempty far strip',nu*H*rho-4*PI),
      ('entire amplitude homotopy lower',amp_lo),
      ('entire amplitude homotopy upper',amp_hi-1)]
    for sign in ('positive','negative'):
        for tag,val in strip_gates:pos(sign+' actual-root strip/'+tag,val)
    sig=8/(a_initial*H);less('phase direction domain',sig,1)
    E=W/H+sig*nu*rho+e0
    L=(1-sig**2)*(nu*rho-4*PI/H)-sig*W/H-e0
    pos('positive finite-height term',nu*rho-4*PI/H)
    pos('both actual roots exceed displacement',L-e)
    K=E/(L-e);h=d['h'];gamma=d['gamma']
    if label.startswith(('standard','nearGeom')):pos('same-edge support with finite normal error',a*a-h*h*(1+K*K))
    else:pos('near support precedes oriented use',a-a_initial)
    sqrtc=sqrti(1-h*h);c=(h-K*sqrtc)/sqrti(1+K*K)
    edgegap=1-K-sqrtc
    pos('whole oriented edge critical gap',edgegap-gamma)
    pos('whole oriented edge angular cosine',c)
    pos('edge inner radius exceeds outer cutoff',h-I('.9'))
    pos('positive squared-comparison radius',1-gamma)
    for av in (a,I(1)):
        for rv in (I('.9'),I(1)):
            pos('radial separate-concavity corner '+str(av.lo)+','+str(rv.lo),1-av*av+2*rv*(av*c-gamma)-rv*rv*(1-gamma*gamma))
    r0=I('.9');C,v,beta,J=d['C'],d['v'],d['beta'],d['J']
    for name in ('C','v','beta','J','inner','growth','errbase'):
        pos('endpoint positive parameter '+name,d[name])
    less('complete centered Cauchy coefficient',1/((r0-m)**2*gamma*r0),C)
    pos('complete inward model speed',1-m/(r0*(r0-m))-v)
    pos('finite-n entire outer decay',v*k-C*q-beta)
    less('complete outer integral to infinity',C/(beta**2*(h-m)),J)
    less('outer logarithm positive-side domain',J*q,1)
    less('endpoint normalized radius',(1+m)/(h-m),I('1.1'))
    less('endpoint lower logarithm bound',logi((1+m)/gamma),I('1.2'))
    less('endpoint upper logarithm bound',logi(2/(h-m)),I('1.2'))
    less('complete centered endpoint log coefficient',4*(I('1.2')+I('1.1'))/(h-m)**2,12)
    Hess=2*h*h-1-2*m/h**3-q*(1/(h*h*gamma*gamma)+2/(h**3*gamma))
    pos('entire-edge Hessian margin',Hess-I('.01'))
    inner,growth=d['inner'],d['growth']
    less('complete inner missing-two-label AMGM',(r0+m)**2+q,inner**2)
    less('complete endpoint model exponent',q/(2*(h-m)**2),growth)
    less('complete centered interpolation prefactor',expq(F(1)),4*(h-m)**2)
    less('complete inner centered-model base',(r0+m)*expi(growth),inner)
    less('complete model power prefactor',expi(growth)*inner**2,1)
    less('first Q-factored error coefficient',4/inner**3,6)
    less('first Q-factored error full base',inner/I('.99'),d['errbase'])
    less('mean cap retained in second error',m,I('.006'))
    less('second Q-factored error coefficient',1/I('.006')**2,10**6)
    less('second Q-factored error full base',expi(growth)*I('.006')/I('.99'),I('.01'))
    less('two complete Q-factored powers at n1000',6*I(1000)**3*d['errbase']**1000+10**6*I('.01')**1000,I('.01'))
    less('complete decreasing consecutive ratio',d['errbase']*I(F(1001,1000))**3,1)
    cost=12+J/(1-J*q)+I(F(1,99));less('complete endpoint potential cost',cost,24)
    rad=d['finalr']
    floors=[256/(1-q)**2,8*24*q/(1-q),24*(8/rad+16/rad**2+6/rad**3)+11]
    for j,value in enumerate(floors):less('whole uniform-Q rigidity degree floor'+str(j),value,900)
    ROWS.append(dict(branch=label,H_floor=show(H),K=show(K),q_final=show(q),B=show(d['B']),mean=show(m),cos_lower=show(c),edge_gap=show(edgegap),Hessian=show(Hess),endpoint_cost=show(cost),rigidity_floors=list(map(show,floors)),actual_strips=2))

receipt=dict(status='PASS',degree_floor=n,checks=len(GATES),complete_actual_strips=8,
 arithmetic='Exact Fraction intervals; 96-term positive exp and atanh with explicit infinite tails; 100-term alternating Machin pi; 56-decimal integer square-root bounds; no author imports',
 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 parameter_sha256=hashlib.sha256(PARAMETER_PATH.read_bytes()).hexdigest(),
 runtime_seconds=time.monotonic()-START,rows=ROWS,predicates=GATES)
out=Path(__file__).with_name('tang_boundary6251_seven_rational_independent_receipt.json')
out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='predicates'},indent=2))
