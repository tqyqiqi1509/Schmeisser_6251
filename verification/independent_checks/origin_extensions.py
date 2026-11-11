"""Fresh targeted final-review checks; no author verifier is imported."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json
import sympy as s

OUT = Path(__file__).resolve().parent
V = OUT.parent / 'small_extensions'
checks = []
def ck(name, value):
    ok = bool(value)
    checks.append({'name': name, 'passed': ok})
    if not ok:
        raise AssertionError(name)

rho, z, r, t, al = s.symbols('rho z r t al', real=True)
cp = rho/2
dp = (z-r-cp*z*t)**2 + (t-cp*z*r)**2
plain = (z-r)**2+t*t
rhs = rho*z*z*(rho*((r-al)**2+t*t)/2-t) + rho*rho*z*z*(z-al)*(z+al-2*r)/2 + rho*rho*z*z*(r*r+t*t)/4
ck('full critical polar identity', s.expand(dp-(1-rho*rho*z*z/2)*plain-rhs)==0)
ck('parallel source denominator identity', s.expand(dp-plain+rho*z*z*t-rho*rho*z*z*(r*r+t*t)/4)==0)

A,X,T,de = s.symbols('A X T de', real=True)
available = A-(A-X/2)*(A*A+(T-t)**2)/((A+de)**2+t*t)
printed = X/2+(A-X/2)*(2*de*A+de*de+2*T*t-T*T)/((A+de)**2+t*t)
ck('critical-label negative-Q coefficient identity', s.cancel(available-printed)==0)
c = s.symbols('c', real=True)
I = T-c*(1+T*T)
D = 1-c*(T+I)
ss = (1+I)/D
bracket = 1+2*c*(1-T)+c*c*(T*T-4*T-1)+2*c**3*(1+T*T)
ck('moving-test exact entrance identity', s.cancel(D*D*(ss-1-c*c*ss*ss)-I*bracket)==0)

# Original degree is retained when original and critical labels map to infinity.
u,w = s.symbols('u w')
cc = s.Rational(1,20)
beta = s.Rational(4,5)-s.I
P = (u+s.I/cc)**2*(u-beta)
d = 1+s.I*cc*w
G = s.cancel(d**3*P.subs(u,w/d)).expand()
J = s.expand(d*s.diff(G,w)-3*s.I*cc*G)
Jdirect = s.cancel(d**2*s.diff(P,u).subs(u,w/d))
ck('degree-drop finite polynomial identity', s.cancel(J-Jdirect)==0)
ck('two original labels at infinity', s.Poly(G,w).degree()==1)
ck('one critical label at infinity', s.Poly(J,w).degree()==1)
aa = -3*s.I*cc/d
ck('first polar correction identity', s.cancel(s.diff(aa,w)+3*cc*cc/d**2)==0)

# Verify normalized D25 with independently reconstructed Bernstein coefficients.
x = s.symbols('x')
N,n = s.Integer(6250),s.Integer(6251)
R = s.Rational(9,100)
f = 1/(1-s.Rational(9,8)*R*R)
Dp = 1+(3*R/4)**2
b0 = 4*(1-R)/5
K = min((n*rr*b0/(n*rr+2*Dp*b0)-rr)/f for rr in [s.Rational(1,100),R])
ck('normalized D25 exact source', K==s.Rational(239012578698741,379195048000000))
ck('normalized D25 lower-X field', K/(R*s.Rational(7,25))>25)
ck('normalized D25 joint field', K*N/(2*s.Rational(31,10))>625)
Nk=x*(4*x-1)*(4*x+9)
Dk=40*(x+1)
pp=s.Poly(s.expand(((s.Rational(62,5)*Nk*Dk-(x*Dk+Nk)**2).subs(x,s.Rational(7,25)+s.Rational(18,25)*x))),x)
bc=[sum(pp.nth(j)*s.Rational(comb(k,j),comb(6,j)) for j in range(k+1)) for k in range(7)]
printed_bc=[s.Rational(120508031,244140625),s.Rational(1411226329,1953125),s.Rational(185744659,78125),s.Rational(416205433,78125),s.Rational(30945547,3125),s.Rational(407837,25),s.Integer(24527)]
ck('normalized D25 all seven printed Bernstein coefficients', bc==printed_bc and min(bc)>0)

# Independent Sturm checks, rather than reusing the authors' subdivision code.
lo,hi,x0=s.Rational(30501,500000),s.Rational(31,500),s.Rational(69,100)
dd=s.Rational(3,5)
kappa=x*(x-(1-dd)**2)*(x+(1+dd)**2)/(2*(1+dd*dd)*(x+2*dd))
def positive_sturm(name, expr, left, right):
    num,den=s.cancel(expr).as_numer_denom()
    for part,label in [(num,'numerator'),(den,'denominator')]:
        pp=s.Poly(part,x)
        ck(name+'/'+label+' endpoints', pp.eval(left)>0 and pp.eval(right)>0)
        ck(name+'/'+label+' no roots', pp.degree()<=0 or pp.count_roots(left,right)==0)
positive_sturm('highX kappa whole interval',kappa-s.Rational(12,25)*x*x,x0,1)
b0=(1-hi)/(1+dd*dd)
f=1/(1-hi*hi*(1+dd)**2/2)
Dp=1+(hi*(1+dd)/2)**2
K=min((n*rr*b0/(n*rr+2*Dp*b0)-rr)/f for rr in [lo,hi])
ck('highX exact source',K==s.Rational(1793507148275409343,2882630131640625000))
cc=s.Rational(12,25)
b=s.Rational(7,10)+s.Rational(3,10)*hi
iv=1-lo*(1+x0*x0)/2
positive_sturm('highX q05 whole interval',N*K-800*x*x*(1-cc*x)*(b/x0+iv*x/(x0*x0)),0,1/cc)
RR,mm=s.Rational(17,98),s.Rational(19,490)
B=2*Dp*(95*f+1)/n
upper=(f*2*s.Rational(25,48)*95*102/(N-190)+hi)/(1-B)
Is=2*(1-hi)*(1-mm)-hi*RR*RR
Ds=2*(1+dd*dd)+4*s.Rational(7,6)*mm+2*RR*RR
lower=(1-hi)/(1+dd*dd)+(Is/2+s.Rational(49,50))**2/Ds
ck('highX exact printed pair margin',lower-upper==s.Rational(9014845822117732608757183663183337,97049415152918630536485588216000000))
ck('highX terminal strict margin',lower-upper>s.Rational(92,1000))

# Inspect every supplementary cell's actual entrance and positivity, without contours.
for filename,family in [('TMOMENT6251_PROPOSAL.json','M'),('T6251_HIGHX_PROPOSAL.json','X')]:
    data=json.loads((V/filename).read_text())
    for i,row in enumerate(data['cells']):
        p,r,l,u,H,de,Abar,U,J=[F(row[k]) for k in ['rho_lo','rho_hi','yl','yh','H','delta','Abar','U','J']]
        tag=f'{family}{i}'
        ck(tag+'/finite polar denominator positive',1-r*r*(1+de)**2/2>0)
        ck(tag+'/inverse Cayley actual entry',(r*(1+de)/2)**2<de)
        radial=(6250*p/2-J)/(r+u)
        if family=='M':
            ck(tag+'/independent field',J==4*u and (H<=(6250*p-F(2637,32))/(2*(1+r)) or H<=radial))
        else:
            ck(tag+'/independent field',J==F(23,5)*u and (H<=95 or H<=radial))
        br=(r/2+J/6250)/(1-r/2)
        ck(tag+'/lawful ratio upper',U>0 and (U>=u or U*U>=br or (l>0 and U*l>=br)))

result={'status':'PASS','checks':len(checks),'scope':'Fresh symbolic identities, degenerate finite Cayley example, normalized D25, independent Sturm highX entry, exact terminal margin, all 30 supplementary scalar entrances. No contour verifier imported or executed.','predicates':checks}
(OUT/'origin_extensions.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='predicates'},indent=2))
