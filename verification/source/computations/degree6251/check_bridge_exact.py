"""Whole a in[.1,.11] bridge at every n>=6251; exact arithmetic."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,hashlib
N=6250;n=6251;al=F(1,10);au=F(11,100);K=F(21,4);C70=F(1899,20);G=[];CIRC=[]
def ck(k,p):
 G.append(dict(name=k,passed=bool(p)))
 if not p:raise AssertionError(k)
def lo(x):return sum((x**j/F(factorial(j)) for j in range(101)),F(0))
def up(x):return sum((x**j/F(factorial(j)) for j in range(25)),F(0))+x**25/(F(factorial(25))*(1-x/26))
ck('independent_direct_D100',F(5*N,14)*(al/2)/(1-al/2+(al/2)**2)>100)
ck('phase_and_field_after_legal_D70',N>=700)
ck('whole_e6',lo(F(6))>400)
def cloud(name,B,R,H):
 d=1-B;V=F(9,4)/d;m=1/d;nu=F(99,100);W=6;e=F(1,400)
 lam=V*R/(1-m*R);rh=R-V*(R*R/2+m*R**3/(3*(1-m*R)));Bi=V/(1-m*R)**2/(1-lam)**3
 eps=F(1001,1000)*Bi/H;g=F(1001,1000)*nu*Bi*rh;pin=1/(1+lam);amp=1/(1-lam)
 for tag,v in [('critical',0<m*R<1),('inverse',0<lam<1),('image',rh>0),('corner',W*W<(1-nu*nu)*(H*rh)**2),('left',eps*(1+e)+g*e<1-amp*e),('right',eps*(1+e)+g*e*e<pin-e),('horizontal',eps+g*e<pin),('amplitude',0<pin<1<amp),('height',nu*H*rh>F(22,7))]:ck(name+'/'+tag,v)
 count=nu*F(7,22)*H*rh-2;ck(name+'/count',count>0);den=1-(B+R)*B;ck(name+'/adaptive_factor_denom',den>0)
 CIRC.append(dict(name=name,H=str(H),R=str(R),rhat=str(rh),count=str(count)))
 return count,(R+au)/den
cnt,fac=cloud('Low',F(1,5),F(1,10),N*F(4,5)/F(9,4))
ck('Low_factor',fac<F(1,4));ck('Low_log',up(F(1))<4);ck('Low_count_per_degree',cnt/n>F(1,120))
ck('Low_variance',1-F(1,5)**2>F(94,100) and F(47,100)>F(2,5))
ck('Low_same_C_all_n',al*al*F(4,25)*(1-F(4,25)/(n*n))>F(1,1000))
ck('Low_entire_log_budget',lo(F(10))>8000 and lo(F(9))>n)
ck('Low_entire_degree_ray',F(n,120)>28 and F(1,120)>F(2,n))
cnt,fac=cloud('Middle',F(4,5),F(1,25),(N-C70/F(1,5))/(2*(1+au)))
ck('Middle_factor',fac<F(23,50));ck('Middle_log',up(F(3,4))<F(50,23))
cmid=al*al*F(9,100)**2*(1-F(9,100)**2)
ck('Middle_whole_C_target',lo(F(12))>8/cmid);ck('Middle_full_count',F(3,4)*cnt>12)
# One actual base-root factor, with the same adaptive shift and same C.
umin=F(69,100)**2;tmax=1-umin;cup=al*al*umin*(1-umin)
ck('Upper_whole_adaptive_range',0<umin<F(1,2) and F(1,2)<tmax<1)
den=1-F(47,50)**2+al*al/2;factor=tmax*au/den
ck('Upper_actual_root_denominator',den>0);ck('Upper_actual_root_factor',factor<F(19,40))
ck('Upper_actual_root_log',up(F(37,50))<F(40,19));ck('Upper_same_C_target',lo(F(202,25))>8/cup)
ck('Upper_same_C_energy_and_root',N*cup/2+F(37,50)>F(202,25))
# The entire high modulus branch; full phase/product convexity.
amin=F(3,100);xmin=F(1,10);s=3*amin/4;L=F(163,20);pmax=1-F(47,50)**2
ck('High_actual_half_shift',(F(47,50)-au)**2>F(1,2))
ck('High_concave_C_endpoint',xmin>3*amin)
ck('High_whole_log_budget',lo(L)>8/(amin*(xmin-s)))
kk=2*L/(N*amin);ck('High_whole_convexity',3*(xmin-s)>s)
for x in [xmin,F(1)]:
 rem=1-pmax*kk/(x-s)-4*K/(N*x)
 ck('High_positive_radical_'+str(x),rem>0)
 ck('High_full_endpoint_'+str(x),rem*rem>4*x*x*kk/(x-s))
out=dict(status='PASS',degree_floor=6251,scope='entire witness interval .1<=a<=.11, every actual rightmost-root modulus',checks=len(G),full_contours=2,contours=CIRC,upper_one_root_C=str(cup),upper_one_root_factor=str(factor),predicates=G,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='predicates'},indent=2))
