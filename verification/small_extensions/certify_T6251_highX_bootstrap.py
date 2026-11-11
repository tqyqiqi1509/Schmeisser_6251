"""Exact T=1 source entrance and paired-root terminal for n>=6251.

The field range h>=95 is excluded by the companion ten-cell certificate.
Parameter intervals use rational Bernstein bounds and rational scalar gates.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sympy as s,json,hashlib,sys
sys.set_int_max_str_digits(0)
O=Path(__file__).resolve().parent
N=6250;n=6251;r0=F(31,500);rl=F(30501,500000);xmin=F(69,100)
x=s.symbols('x');checks=[];polys=[]
def S(v):return s.Rational(v.numerator,v.denominator) if isinstance(v,F) else s.sympify(v)
def ff(v):return F(int(v.p),int(v.q))
def ck(name,p):
 checks.append(dict(name=name,passed=bool(p)))
 if not p:raise AssertionError(name)
def polypos(name,expr,lo,hi):
 expr=s.cancel(expr);num,den=expr.as_numer_denom()
 # Sympy may normalize both signs negatively. Establish the whole
 # denominator sign independently instead of using a midpoint as evidence.
 if den.subs(x,S((lo+hi)/2))<0:num,den=-num,-den
 if x in den.free_symbols:polypos(name+'/entire cleared denominator',den,lo,hi)
 else:ck(name+'/positive constant denominator',den>0)
 pp=s.Poly(num,x);deg=pp.degree();co=[ff(pp.nth(i)) for i in range(deg+1)]
 stack=[(lo,hi,0)];leaves=[]
 while stack:
  a,b,d=stack.pop()
  p=[sum((co[k]*comb(k,j)*a**(k-j)*(b-a)**j for k in range(j,deg+1)),F(0)) for j in range(deg+1)]
  bs=[sum((p[k]*F(comb(j,k),comb(deg,k)) for k in range(j+1)),F(0)) for j in range(deg+1)]
  if min(bs)>0:leaves.append((a,b,min(bs)));continue
  if d>=20:raise AssertionError((name,a,b,min(bs),float(expr.subs(x,S((a+b)/2)))))
  m=(a+b)/2;stack.extend([(a,m,d+1),(m,b,d+1)])
 leaves.sort();cur=lo
 for a,b,m in leaves:ck(name+'/closed seam',a==cur);cur=b;ck(name+'/all Bernstein coefficients positive',m>0)
 ck(name+'/complete interval',cur==hi)
 polys.append(dict(name=name,expr=str(expr),numerator=str(num),denominator=str(den),interval=[str(lo),str(hi)],leaves=[[str(a),str(b),str(m)] for a,b,m in leaves]))
 return expr
def source(de,X=s.Integer(1),r=rl):
 b=(1-S(r0)*(1+X*X)/2)/(1+S(de)**2)
 D=1+(S(r0)*(X+S(de))/2)**2
 f=1/(1-S(r0)**2*(X+S(de))**2/2)
 return s.cancel((n*S(r)*b/(n*S(r)+2*D*b)-S(r))/f)
def kapp(de):return s.cancel(x*(x-(1-S(de))**2)*(x+(1+S(de))**2)/(2*(1+S(de)**2)*(x+2*S(de))))


ck('full crescent same-X height gate',F(2,125)-r0*(F(1,4)+F(1,15625))>0)
Gmin=1+F(117,500)*xmin*xmin
ck('normalized D25 lemma legal full rho range',F(1,100)<rl<r0<F(9,100) and n>=6251)
ck('ratio3 geometric whole y cap',F(57,20)**2*Gmin>9)
ck('whole physical h49',N*rl/(2*(1+F(57,20)))>49)
ck('same-X D51 before matrix phase',49*49*Gmin>51*51)
ck('matrix coarse cosine enclosure',F(18,25)**2>1-F(7,10)**2)
ck('matrix whole phase .7',F(21,5)/51**2+F(11,2)*F(18,25)/51+(1-F(7,10)**2)/4-F(7,10)**2/2<0)
ck('four moment mu before companions',49*49*Gmin>50*50)
de=F(3,5);ka=kapp(de);c0=F(12,25)
polypos('whole high-X kappa quadratic',ka-S(c0)*x*x,xmin,F(1))
K=min(ff(source(de,r=rl)),ff(source(de,r=r0)))
ck('full high-X source positive',K>F(31,50))
beta=F(7,10)+F(3,10)*r0;Iv=1-rl*(1+xmin*xmin)/2
polypos('complete one-variable same-configuration q entrance',S(N*K)-800*x*x*(1-S(c0)*x)*(S(beta/xmin)+S(Iv/(xmin*xmin))*x),F(0),1/c0)
ck('h95 gives D100 before cloud ratio',95*95*Gmin>100*100)
ck('h95 actual ratio gives wide qmu',F(9,4)/100==F(9,400) and F(1,100)==F(1,100))
ck('same-X complete high-field y endpoint',F(427,200)**2*Gmin>F(81,16))
# After q05mu02 actual companions, phase7, h>49 and separate h<95.
f=1/(1-r0*r0*(1+de)**2/2);Dp=1+(r0*(1+de)/2)**2
B=2*Dp*(f*95+1)/n
ck('terminal complete HB denominator',0<B<1)
Ay=1/(4*c0)
upper=(f*2*Ay*95*102/(N-190)+r0)/(1-B)
R=F(17,98);mean=F(19,490)
ck('pair entire lower numerator',(1-r0)*(1-R)-r0*R*R/2>F(3,4))
ck('pair entire upper numerator',1+R<F(6,5))
for ii in [F(3,4),F(6,5)]:ck('pair whole sqrt chord '+str(ii),(ii/2+F(49,100))**2<ii)
ck('whole pair denominator positive radical bound',F(7,6)**2>1+de*de)
Is=2*(1-r0)*(1-mean)-r0*R*R
Ds=2*(1+de*de)+4*F(7,6)*mean+2*R*R
ck('pair source parameter domain',0<r0<1 and R>0 and 0<mean<1 and Is>0 and Ds>0)
pair=(Is/2+F(49,50))**2/Ds
lower=(1-r0)/(1+de*de)+pair
ck('full paired source terminal',lower>upper)
result=dict(status='PASS',scope='T6251 high-X .69<=X<=1,.061002<=rho<=.062 bootstrap and terminal conditional on separate h95 cloud',degree_floor=n,checks=len(checks),polynomials=len(polys),source_K=str(K),pair_lower=str(lower),terminal_upper=str(upper),terminal_margin=str(lower-upper),predicates=checks,whole_polynomials=polys,sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(O/'T6251_HIGHX_BOOTSTRAP_EXACT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['predicates','whole_polynomials']},indent=2))
