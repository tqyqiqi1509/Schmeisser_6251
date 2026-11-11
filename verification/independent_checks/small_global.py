#!/usr/bin/env python3
"""Fresh targeted algebraic checks for the local/small/global final review."""
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as s
checks=[]
def ck(name, condition, value=None):
    ok=bool(condition)
    checks.append({'name':name,'pass':ok,'value':str(value) if value is not None else None})
    assert ok,name
x,a,R,v,D,k=s.symbols('x a R v D k', real=True)
lower=(x-a/2)**2+(1-R*a*a/4)**2
ck('normalized distance derivative',s.expand(s.diff(lower,a)-(-x+(s.Rational(1,2)-R)*a+R*R*a**3/4))==0)
ck('normalized distance endpoint',s.expand(lower.subs(a,x)-(1+(s.Rational(1,4)-R/2)*x*x+R*R*x**4/16))==0)
ck('normalized distance derivative negative on declared whole box', -1+F(1,2)-F(9,100)+F(9,100)**2/F(4)<0)
rho=F(1,100);cap=F(9,100);n=6251
b0=4*(1-cap)/5;f=1/(1-F(9,8)*cap**2);dp=1+(3*cap/4)**2
K=lambda r:(n*r*b0/(n*r+2*dp*b0)-r)/f
K0=min(K(rho),K(cap))
ck('normalized source exact printed K0',K0==F(239012578698741,379195048000000),K0)
ck('normalized first source field above 25',K0/(cap*F(7,25))>25)
ck('normalized second source squared field above 625',K0*6250/(2*F(31,10))>625)
Nk=x*(4*x-1)*(4*x+9);Dk=40*(x+1)
poly=s.expand(s.Rational(62,5)*Nk*Dk-(x*Dk+Nk)**2)
t=s.symbols('t');P=s.Poly(s.expand(poly.subs(x,s.Rational(7,25)+s.Rational(18,25)*t)),t);deg=P.degree()
coef=[sum(P.nth(j)*s.binomial(i,j)/s.binomial(deg,j) for j in range(i+1)) for i in range(deg+1)]
expected=[s.Rational(120508031,244140625),s.Rational(1411226329,1953125),s.Rational(185744659,78125),s.Rational(416205433,78125),s.Rational(30945547,3125),s.Rational(407837,25),s.Integer(24527)]
ck('all seven printed normalized Bernstein controls',coef==expected,coef)
ck('all normalized Bernstein controls strictly positive',all(c>0 for c in coef))
expr=3*(2+v)**2-(D-s.Rational(91,10))*(2+v)+2*D
ck('D70 inverse quadratic branch',s.simplify(expr.subs(D,s.Rational(151,5)/v+s.Rational(211,10)+3*v))==0)
expr=s.Rational(5,4)*(2+v)**2-(D-s.Rational(23,2))*(2+v)+2*D
ck('Medium139 inverse quadratic branch',s.simplify(expr.subs(D,28/v+s.Rational(33,2)+s.Rational(5,4)*v))==0)
k=F(9,4);C0=F(5,4)*k*k+F(23,2)*k
ck('Medium139 entire hole constant',2*C0+8*k==F(2637,32),2*C0+8*k)
poly=s.expand(6250**2*(s.Rational(3,50)**2+a*a/4)-(s.Rational(1899,20)+278*(1+a))**2)
ck('low-middle positive quadratic coefficients',poly==9688341*a*a-s.Rational(1036801,5)*a+s.Rational(613319,400),poly)
P=s.Poly(poly,a);A,B,C=P.all_coeffs();minimum=s.factor(C-B*B/(4*A))
ck('low-middle whole real-line minimum',minimum==s.Rational(65688734375,155013456),minimum)
ck('low-middle whole real-line positivity',minimum>0)
# Full rightmost half-disk expansion, without a sampling argument.
x,y,tau,u,w,f=s.symbols('x y tau u w f',real=True);b=x+s.I*y;gamma=b-tau;z=b+u+s.I*w
sq=lambda e:s.expand(e*s.conjugate(e))
pb=1-sq(b);Lg=1-sq(gamma);A0=pb+tau*s.conjugate(b);V=f*f*Lg*s.conjugate(b)+tau*(1-f*f)
formula=f*f*sq(A0)-tau*tau-(1-f*f*sq(gamma))*(u*u+w*w)-2*s.re(V*(u+s.I*w))
actual=f*f*sq(1-s.conjugate(z)*gamma)-sq(gamma-z)
ck('entire half-disk factor identity',s.simplify(actual-formula)==0)
# Positive Abel summation checks the telescoping identity for independent symbols.
W=s.symbols('w1:5');S=s.symbols('s0:5');S=[0,*S[1:]]
left=sum((W[j]-W[j+1])*S[j+1] for j in range(3))+W[-1]*S[4]
right=sum(W[j]*(S[j+1]-S[j]) for j in range(4))
ck('positive cumulative Abel telescoping',s.expand(left-right)==0)
# The sign conditions in the larger signed-curvature disk.
ck('phase-conditioned chord 4/3 gate',F(21,25)**2*(1+F(3011,4000))*F(4,3)**2>2)
rr=F(3,50);lips=F(4,3)*rr*(3+4*rr+29*rr**2+436*rr**3)
ck('phase-conditioned cancellation exact number',lips==F(214911,781250),lips)
ck('phase-conditioned cancellation strict gate',lips<F(7,22))
# Fine200 uses its full mu=1/200 rectangle; curvature constants independent of D270.
ck('Fine200 full-domain q-squared curvature',F(11,5)+F(120000,18)*F(7,625)**3<F(9,4))
ck('Fine200 full-domain qmu curvature',F(220,18)<13)
# Exact seam and endpoint comparisons in the high small-x branch.
r0=F(30501,500000)
ck('high small-x x lower seam',F(47,50)**2>F(433,10000)**2+(1-r0)**2)
ck('high small-x tau upper seam',F(47,50)**2-r0*r0>(1-F(31,500))**2)
ck('high small-x normalized X lower seam',F(433,10000)/F(31,500)>F(69,100))
ck('low-middle normalization overlap',F(47,50)**2<r0*r0+(1-r0)**2)
ck('low-middle moment extension upper seam',F(37,40)**2-r0*r0>(1-F(39,500))**2)
result={'status':'PASS','checks':len(checks),'scope':'Fresh targeted exact identities and domain seams; no previously generated PASS input.','results':checks}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'checks':len(checks)}))
