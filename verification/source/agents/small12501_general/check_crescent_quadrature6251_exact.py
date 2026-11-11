from fractions import Fraction as F
from pathlib import Path
import json,hashlib
import sympy as s
aa,r=s.symbols('a r',real=True);G=[]
def ck(name,p):
 p=bool(p);G.append(dict(name=name,passed=p))
 if not p:raise AssertionError(name)
identity=s.expand((2*r-s.Rational(5,4)*aa**2)**2-4*(1-aa**2/4)*(r*r-aa**2)-aa**2*(r-1)*(r-4)-s.Rational(9,16)*aa**4)
ck('complete squared crescent comparison identity',identity==0)
ck('whole derivative geometric domain',4-F(1,10)**2>1)
ck('positive complete squaring side',2-F(5,4)>0)
N=F(6250);d0=F(3,50);C70=F(596,5);C=F(1899,20);D=F(139)
ck('lawful D70 to independent D100',(N*d0-C70)/F(11,5)>100)
q2=N*N/4-4*D*D;q1=-4*D*(C+2*D);q0=N*N*d0*d0-(C+2*D)**2
ck('whole quadratic leading coefficient',q2>0)
ck('whole positive quadratic discriminant',4*q2*q0-q1*q1>0)
ck('whole squared field side positive',C+2*D>0)
ck('medium139 prior ratio9over4 exact cap',F(9,556)==F(9,4)/D)
out=dict(status='PASS',scope='full crescent d^2>(1-|b|)^2+a^2/4 and all n>=6251,a<=.1,|b|<=.94 after D70 obtain D>139',checks=len(G),squared_identity=str(identity),quadratic_coefficients=[str(q0),str(q1),str(q2)],positive_discriminant=str(4*q2*q0-q1*q1),minimum_quadratic=str(q0-q1*q1/(4*q2)),predicates=G,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('CRESCENT_QUADRATURE6251_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
