"""Independent final-review checks. No historical verifier imported."""
from fractions import Fraction as F
from flint import arb,ctx
import sympy as s
import json
from pathlib import Path
ctx.prec=256
checks={}
def need(name,b):
    assert b,name
    checks[name]=True

def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator

# Seven complete energy subinterval envelopes, independently evaluated.
rows=[(F(0),F(1,4),F(659,1000)),(F(1,4),F(3,8),F(396,1000)),
      (F(3,8),F(1,2),F(346,1000)),(F(1,2),F(5,8),F(192,1000)),
      (F(5,8),F(3,4),F(54,1000)),(F(3,4),F(7,8),F(11,1000)),
      (F(7,8),F(19,20),F(6,1000))]
values=[]
m=F(6250)
for i,(lo,hi,target) in enumerate(rows):
    dm=max(F(1,10),1-hi*hi);dp=1-lo*lo
    Lm=ar((m+2)*dm/2);Lp=ar((m+2)*dp/2)
    D=7*Lp/(ar(m)*Lp.log())
    if lo:D=min(D,ar(28/(m*lo*lo)))
    U=ar(F(5,3)*hi*hi*m)*(-ar((m-1)/(m+2)/dp)*Lm.log()).exp()
    V=(-ar(m/(m+2)/dp)*Lm.log()).exp()
    value=D+U+V
    need('energy interval '+str(i+1),value<ar(target))
    need('energy interval target gap '+str(i+1),target<1-hi)
    values.append(str(value))
need('clipped interval coverage',rows[-1][1]**2>F(9,10))

# Symbolic full matrix envelope; p denotes pi, h denotes eta.
p,h=s.symbols('p h',real=True)
lam=(24*p+(40-3*p*p)*h)/(6*p*(3*p+4*h))
al2=(2/p-lam)**2+(2/(3*p)-lam)**2
Delta=3*p*p-8
rhs=9*p*p/Delta+(2+24*p/Delta)*h+(40-3*p*p)/(2*Delta)*h*h
need('optimized full Schur envelope identity',s.factor(2*h+3*(1+lam*h)**2/(1-3*al2)-rhs)==0)
# Positive actual endpoints for lambda range and quadratic convex maximum.
pi=arb.pi();lzero=4/(3*pi);linf=(40-3*pi*pi)/(24*pi)
need('lambda endpoint range',linf>2/(5*pi) and lzero<2/pi)
for ix,l in enumerate((lzero,linf)):
    need('Schur quadratic concavity endpoint '+str(ix),3*((2/pi-l)**2+(2/(3*pi)-l)**2)<1)

# Full-domain algebra in saddle displacement and convex variance.
X,Y,a,b=s.symbols('X Y a b',real=True)
s2=X*X+Y*Y;r2=a*a+b*b;Den=1+2*Y+s2
B=X/Den;A=s.symbols('A',real=True)
dv=1-2*X*a-2*Y*b+s2*r2
lhs=(A*a-B*r2)*dv-(a-X*r2)
rhs=a*(A*dv-1+2*B*X*r2)+B*r2*(2*Y*(1+b)+s2*(1-r2))
need('saddle affine identity',s.factor(lhs-rhs)==0)
u,c=s.symbols('u c',real=True)
need('convex variance affine identity',s.expand((u+c*c)*(1-2*u*c+c*c)-(1+c*c)*(u-c)**2-u*(1-u)*(1+c)**2)==0)

# All D270 comparisons, independently rebuilt from exact rational formulas.
def C(r):
    v=r+17*r*r+137*r**3+313*r**4
    A=(1+17*r+137*r*r+313*r**3)/(1-v)
    B=(18+154*r+450*r*r+313*r**3)/(1-v)
    return 1+(2+F(100,199))*r*A*A+3*r*B+F(600,199)*r*r*A*A+4*r*A+29*r+436*r*r
D=F(270);Y0=F(3523,270)
def P(y):return y*y*C(y/D)+F(50000,3)*y**5/D**3+F(650,9)*y-D*(y-2)
need('affine replacement for C13 at D270',20+11*D<=13*(1+D))
need('coarse q entrance',13*(D+1)/D**2<F(97,2000))
need('coarse left exact fraction',P(F(3))==-F(343077017913365252920829,10130931628072339534650))
need('coarse right endpoint',P(Y0)<-197)
need('fine q entrance',F(3,270)<F(7,625))
need('fine left endpoint',F(9,4)*F(9,4)**2+13*F(9,4)-D*F(1,4)==-F(1719,64))
need('fine right endpoint',F(9,4)*9+39-D==-F(843,4))

# Actual complete-circle scalar gates; absolute row bounds, not discovery margins.
for nm,q,mu,Z0,M0,E0,model,P0,R5,Rnu,Z1,L,beta in [
 ('coarse',F(97,2000),F(1,270),F(39,5),F(101,200),F(43,25),126000,550,300000,1300,F(41,5),F(143,100),F(27,100)),
 ('fine',F(7,625),F(1,200),F(657,100),F(2,25),F(109,100),86200,175,120000,220,F(659,100),F(11,10),F(1,100))]:
    nu=q*mu; R=R5*q**5+Rnu*nu
    r3=nu*Z0**3/(3*(1-mu*Z0))
    need(nm+' full derivative exponent',ar(M0+r3).exp()<ar(E0))
    need(nm+' full actual primitive path',Z0**4*E0/(12*(1-mu*Z0))<P0)
    need(nm+' actual circle critical free',Z0+R<Z1<1/mu)
    need(nm+' nonlinear circle', (ar(L*R).exp()-1-ar(L*R))/ar(L*R)<ar(beta))
    need(nm+' coefficient reserve q5',E0*model<(1-beta)*R5)
    need(nm+' coefficient reserve qmu',E0*P0<(1-beta)*Rnu)

# Paper formulas for degree derivatives and low/mid interval domains.
mm,dd=s.symbols('m delta',positive=True)
LL=(mm+2)*dd/2
VV=mm*s.log(LL)/((mm+2)*dd)
need('V degree derivative',s.simplify(s.diff(VV,mm)-(mm+2*s.log(LL))/((mm+2)**2*dd))==0)
UU=s.log(mm)-(mm-1)*s.log(LL)/((mm+2)*dd)
need('U degree derivative',s.simplify(s.diff(UU,mm)-((mm+2)**2*dd-mm*(mm-1)-3*mm*s.log(LL))/(mm*(mm+2)**2*dd))==0)
need('D30 q upper enlarged domain',F(20+11*30,30**2)<F(403,900))
need('D100 q upper enlarged domain',F(20+11*100,100**2)<F(1313,10000))
need('D25 phi exclusion',F(21,5)/625+F(181,200)/25-F(1,8)==-F(513,6250))

out={'status':'PASS','count':len(checks),'checks':checks,'energy_whole_cell_arb_enclosures':values,
     'scope':'independent symbolic and exact/256-bit endpoint reconstruction; no old proof code imported'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
