#!/usr/bin/env python3
"""New complete medium139 actual companions; all accepting gates rational.
The whole relative Gaussian path is rebuilt coefficient by coefficient.
The complete model exponential tail and full actual R3 tail are retained.
"""
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
import json,hashlib
G=[]
def ck(name,l,r,op='<'):
 ok={'<':l<r,'>':l>r,'==':l==r,'<=':l<=r,'>=':l>=r}[op]
 G.append(dict(name=name,left=str(l),right=str(r),relation=op,passed=ok))
 if not ok:print('FAIL',name,float(l),op,float(r))
def eu(x):return sum((x**k/F(factorial(k)) for k in range(25)),F(0))+x**25/(F(factorial(25))*(1-x/26))
def el(x):return sum((x**k/F(factorial(k)) for k in range(101)),F(0))
def mul(a,b,n):
 out=[F(0)]*(n+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<=n:out[i+j]+=x*y
 return out
def poly(a,q):return sum((x*q**i for i,x in enumerate(a)),F(0))
Z=[F(63,10),F(21),F(133),F(1030),F(9000)]
Mc=[F(0),F(63,10),F(45),F(400),F(3400),F(86000),F(320000),F(1700000),F(9100000),F(39000000)]
B=[F(82000),F(350000),F(1900000),F(13000000),F(120000000),F(3100000000),F(23000000000)]
power=[F(1)]+[F(0)]*11;expo=power.copy()
for k in range(1,12):
 power=mul(power,Mc,11)
 for j in range(12):expo[j]+=power[j]/factorial(k)
P11=mul(Z,expo,11)
q=F(9,556);mu=F(1,139)
model=poly(B,q)+(poly(Z,q)*eu(poly(Mc,q))-poly(P11,q))/q**5
# Full relative-path first coefficient, with all u in [0,1].
cubic=[F(1),F(1),F(-10),F(10)];control=[]
for lo,hi in [(F(0),F(1,2)),(F(1,2),F(3,4)),(F(3,4),F(1))]:
 pw=[sum((cubic[k]*comb(k,j)*lo**(k-j)*(hi-lo)**j for k in range(j,4)),F(0)) for j in range(4)]
 bb=[sum((pw[j]*F(comb(i,j),comb(3,j)) for j in range(i+1)),F(0)) for i in range(4)]
 for i,b in enumerate(bb):ck('entire relative leading Bernstein/'+str(lo)+'/'+str(i),b,F(0),'>')
 control.append(dict(lo=str(lo),hi=str(hi),coefficients=[str(x) for x in bb]))
ck('whole pi-square witness',F(44,7)**2,F(40))
ck('whole leading path norm bound',F(44,7),F(63,10))
sq=mul(Z,Z,8);Mr=[F(0),F(63,10)]
for k in range(2,10):Mr.append((Z[k] if k<len(Z) else 0)+sq[k-1]/2)
Z0=F(6681,1000);M0=F(29,250);MR0=F(183,1000);E0=F(57,50);modelcap=F(89000)
RA=F(120000);RB=F(225);R=RA*q**5+RB*q*mu;Z1=F(672,100);L=F(28,25);beta=F(3,200)
r3=q*mu*Z0**3/(3*(1-mu*Z0));pgeom=Z0**4/(12*(1-mu*Z0));actualtail=pgeom*eu(MR0+2*r3)
for name,l,r in [
 ('legal entire old Gaussian rectangle',q,F(1,20)),
 ('complete center polynomial',poly(Z,q),Z0),
 ('complete actual model exponent',poly(Mc,q),M0),
 ('complete relative Gaussian exponent',poly(Mr,q),MR0),
 ('full model residual including infinite tail',model,modelcap),
 ('critical-free entire center path',mu*Z0,F(1)),
 ('actual derivative full exponential',eu(M0+r3),E0),
 ('relative exponential tail denominator',MR0+2*r3,F(26)),
 ('whole model q5 Rouche budget',E0*modelcap,(1-beta)*RA),
 ('whole relative actual qmu Rouche budget',actualtail,(1-beta)*RB),
 ('whole actual new circle containment',Z0+R,Z1),
 ('whole circle critical-free containment',mu*Z1,F(1)),
 ('complete logarithmic derivative',1+q*Z1+q*mu*Z1**2/(1-mu*Z1),L),
 ('full nonlinear circle error',(eu(L*R)-1-L*R)/(L*R),beta),
 ('base original root excluded',21*q+133*q*q+1030*q**3+9000*q**4+R,F(6)),
 ('both actual circles disjoint',2*R,12*(1-q-17*q*q-137*q**3-313*q**4)),
 ('entire rightmost phase4',21+133*q+1030*q*q+9000*q**3+RA*q**4+RB*mu,F(157,25)*4),
 ('complete signed q2 coefficient',F(6,5)+F(1250,24649)*RA*q**3,F(5,4)),
 ('complete signed qmu coefficient',F(1250,24649)*RB,F(23,2)),
 ]:ck(name,l,r)
ck('full improved medium field',2*(F(5,4)*F(81,16)+F(23,2)*F(9,4))+18,F(2637,32),'==')
ck('whole medium139 smaller branch',F(5,4)*F(81,16)-(139-F(23,2))*F(9,4)+278,F(0))
ck('whole medium139 branch derivative',2*F(5,4)*F(9,4)-(139-F(23,2)),F(0))
# A rigorous failure receipt for the unchanged old two-absolute-envelope route.
za=poly(Z,q);ma=poly(Mc,q);rr=q*mu*za**3/(3*(1-mu*za));old_min=za**4/(12*(1-mu*za))*el(2*(ma+rr))
ck('old absolute-path certificate cannot reserve beta=.015 at B225',old_min,(1-beta)*RB,'>')
fails=[x['name'] for x in G if not x['passed']]
vals=dict(q=q,mu=mu,Z0=Z0,M0=M0,MR0=MR0,E0=E0,modelcap=modelcap,model_actual=model,relative_polynomial=poly(Mr,q),R=R,r3=r3,actual_normalized_qmu=actualtail,Z1=Z1,L=L,beta=beta,RA=RA,RB=RB,field=F(2637,32),old_absolute_route_lower=old_min,old_available_qmu=(1-beta)*RB)
out=dict(status='PASS' if not fails else 'FAIL',gate_count=len(G),failures=fails,values={k:str(v) for k,v in vals.items()},decimal={k:float(v) for k,v in vals.items()},relative_path_coefficients=[str(x) for x in Mr],relative_leading_Bernstein=control,rebuilt_full_model_P11=[str(x) for x in P11],predicates=G,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('MEDIUM139_RELATIVE_COMPLETE_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','gate_count','failures','decimal']},indent=2))
if fails:raise SystemExit(1)
