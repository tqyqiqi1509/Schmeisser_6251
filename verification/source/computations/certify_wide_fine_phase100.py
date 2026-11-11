"""Exact whole-rectangle quartic companion and phase bound at D>=100.

No floating arithmetic and no author verifier imports.  The unchanged exact
quartic center/residual coefficients are reconstructed in the archived
independent moment proof; this script rechecks their entire model-tail bound
and every scalar gate for the NEW parameter rectangle.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
rows=[]
def ck(label,lhs,rhs,kind='<'):
 ok=lhs<rhs if kind=='<' else lhs==rhs
 rows.append(dict(label=label,lhs=str(lhs),rhs=str(rhs),kind=kind,pass_=ok))
 if not ok: raise AssertionError(rows[-1])
def eup(x):
 assert 0<=x<26
 return sum((x**k/factorial(k) for k in range(25)),F(0))+x**25/(factorial(25)*(1-x/26))
def conv(a,b,K=11):
 c=[F(0)]*min(K+1,len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<len(c):c[i+j]+=x*y
 return c
def val(a,x):return sum((c*x**j for j,c in enumerate(a)),F(0))
Z=[F(63,10),F(21),F(133),F(1030),F(9000)]
M=[F(0),F(63,10),F(45),F(400),F(3400),F(86000),F(320000),F(1700000),F(9100000),F(39000000)]
# Rebuild the degree11 truncation by direct powers, not the author's recurrence.
power=[F(1)];ep=[F(0)]*12
for k in range(12):
 for j,c in enumerate(power):ep[j]+=c/factorial(k)
 power=conv(power,M)
p11=conv(Z,ep)
x=F(1,20)
res=sum((F(c)*x**j for j,c in enumerate([0,0,0,0,0,82000,350000,1900000,13000000,120000000,3100000000,23000000000])),F(0))
res+=val(Z,x)*eup(val(M,x))-val(p11,x)
ck('entire Gaussian residual at q<=1/20',res,132200*x**5)
q=F(9,400);mu=F(1,100);nu=q*mu
Z0=F(343,50);M0=F(43,250);E0=F(5,4);P0=F(250)
Z1=F(174,25);L=F(59,50);beta=F(3,50)
A=F(250000);B=F(400);R=A*q**5+B*nu
ck('rectangle inside model range',q,F(1,20))
ck('center bound',val(Z,q),Z0)
ck('whole path exponent bound',val(M,q),M0)
ck('positive full log disk',mu*Z0,F(1))
ck('full exponential path and inverse derivative',eup(M0+nu*Z0**3/(3*(1-mu*Z0))),E0)
ck('primitive perturbation coefficient',Z0**4*E0/(12*(1-mu*Z0)),P0)
ck('whole circle containing radius',Z0+R,Z1)
ck('whole containing disk is critical free',mu*Z1,F(1))
ck('whole disk logarithmic derivative',1+q*Z1+nu*Z1**2/(1-mu*Z1),L)
ck('entire nonlinear circle error',eup(L*R)-1-L*R,beta*L*R)
ck('model coefficient Rouche reserve',E0*132200+beta*A,A)
ck('finite product coefficient Rouche reserve',E0*P0+beta*B,B)
correction=21*q+133*q**2+1030*q**3+9000*q**4
ck('origin excluded by whole companion disk',correction+R,F(6))
ck('the two whole companion disks are disjoint',2*(correction+R),F(12))
phase=21+133*q+1030*q**2+9000*q**3+A*q**4+B*mu
ck('phase coefficient 23/5 via pi>157/50',phase,2*F(157,50)*F(23,5))
ck('D100 enters q rectangle',F(9,4)/100,q,'=')
ck('D100 enters mu rectangle',F(1,100),mu,'=')
out=dict(status='PASS',checks=len(rows),q=str(q),mu=str(mu),radius=str(R),phase_factor='23/5',model_tail_endpoint=float(res/x**5),diagnostic_phase_coefficient=float(phase),new_row_rows=rows)
(ROOT/'computations/WIDE_FINE_PHASE100_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='new_row_rows'},indent=2))
