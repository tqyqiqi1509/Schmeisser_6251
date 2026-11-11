#!/usr/bin/env python3
"""Independent phase-conditioned full signed bridge for the W8 circles.

The old q<=.05 right-semicircle Lipschitz constant is not extrapolated.
This file proves a new whole-q<=.06 chord gate using an already legal
actual phase8. All center terms are freshly evaluated; no author numeric
receipt is used. The root's independent signed scalar13/10 is an input.
"""
from fractions import Fraction as F
from pathlib import Path
from flint import arb,acb,ctx
import json,hashlib
ctx.prec=448
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def A(v):
 if isinstance(v,F):return arb(v.numerator)/v.denominator
 return arb(v)
G=[]
def ck(name,yes,margin=None):
 yes=bool(yes);G.append(dict(name=name,passed=yes,margin=None if margin is None else str(margin)))
 if not yes:raise AssertionError(name)
qwide=F(3,50);v=qwide+17*qwide**2+137*qwide**3+313*qwide**4;m=1-v
ck('full theta disk positive real B',m>F(21,25))
ck('full tau imaginary bound',v/m<F(1,5))
ck('whole tau real bound',F(97,100)**2+F(1,5)**2<1)
ck('already legal phase8 gives full sigma imaginary bound',8*qwide==F(12,25))
ck('whole sigma real bound',F(7,8)**2+F(12,25)**2<1)
dot=F(7,8)*F(97,100)-F(12,25)*F(1,5)
chord=F(4,3)
ck('whole phase-conditioned chord',F(21,25)**2*(1+dot)*chord**2>2)
lip=chord*qwide*(3+4*qwide+29*qwide**2+436*qwide**3)
ck('whole phase-conditioned Lipschitz cancellation',lip<F(7,22))
p=2*arb.pi();p2=p*p
ck('true pi upper for complete cancellation',arb.pi()<A(F(22,7)))
for label,value,bound in [('b2',3-p2/2,17),('b3',15-23*p2/6,137),('b4',105-88*p2/3+7*p2*p2/8,313),('D3',21-5*p2/4,29),('D4',144-44*p2/3,436)]:
 ck('whole exact coefficient '+label,abs(value)<bound)
signed=ROOT/'computations/degree6251/ROOT_SIGNED_THETA_13_10_Q06_EXACT.json'
sg=json.loads(signed.read_text())
ck('independent root scalar input scope',sg['status']=='PASS' and sg['scope']=='entire complex theta disk |theta|<=3/50' and sg['coefficient']=='13/10')
center=ROOT/'agents/tang_weighted_moment_center_W8.json';raw=center.read_bytes();cc=json.loads(raw)
q=F(273,5000);mu=F(91,5000);R=22*q*q+45000*q*mu**3
ck('whole second circle falls in new signed bridge',q<qwide)
ck('second actual phase implies phase8',F(33,5)<8)
ck('weighted D scale has exact ratio3 cap',q==3*mu)
ss=acb(0,p);shift=arb(0);terms=[];means=[arb(0),arb(0)]
for row in cc['coefficients']:
 aa,bb,dd=row['monomial'];w=aa+2*bb+3*dd
 coeff=[A(F(z)) for z in row['coefficients']]
 zz=sum((acb(z)*ss**j for j,z in enumerate(coeff)),acb(0))
 oddmean=sum((acb(z)*ss**j for j,z in enumerate(coeff) if j%2==0),acb(0))
 if bb+dd==0 and w<=4:continue
 ck('full nonquartic term weight at least2 '+str(row['monomial']),w>=2)
 term=abs(zz)*A(q**aa*(q*mu)**bb*(q*mu*mu)**dd);shift+=term
 terms.append(dict(monomial=row['monomial'],weight=w,norm=str(term)))
effective=shift+A(R)
curvature=A(F(13,10)*q*q)+effective/(2*arb.pi()**2)
C0=curvature/A(mu*mu);field=2*C0+6*A(F(33,5))
ck('full weighted D curvature85',C0<85,A(85)-C0)
ck('full same finite-polynomial field210',field<210,A(210)-field)
rho=F(30501,500000);N=F(6250)
# Independently reconstruct all scalar joins after the full normalized
# source/matrix D25 entrance, retaining the actual phase at each stage.
Doutside=N/(2/rho+96)
ck('first outside matrix phase floor48',Doutside>48)
matrixout=F(21,5)/48**2+F(11,2)*F(3,4)/48+F(241,1000)*F(5,9)-F(2,9)
ck('whole outside matrix excludes imaginary2over3',matrixout<0)
Doutfirst=N*rho/(2*rho+F(4,3)*(1-rho)+3*(2-rho))
ck('whole outside q06 actual-circle entrance',Doutfirst>53)
Dlow0=N*F(3,50)/8
matrixlow=F(21,5)/Dlow0**2+F(11,2)*F(143,200)/Dlow0+F(241,1000)*F(51,100)-F(49,200)
ck('whole low-middle matrix excludes imaginary.7',matrixlow<0)
Dlowfirst=N*F(3,50)/(2*rho+F(7,5)*F(47,50)+3*F(97,50))
ck('whole low-middle q06 actual-circle entrance',Dlowfirst>51)
Dout2=(N-48/rho)/98
Dlow2=(N-48/F(3,50))/(2*rho/F(3,50)+97)
ck('full outside Q06-to-Q0546 join',Dout2>1/mu)
ck('full low-middle Q06-to-Q0546 join',Dlow2>1/mu)
ck('full low-middle field210 enters D70',(N*F(3,50)-210)/F(11,5)>70)
ck('full outside field210 enters D70',(N*rho-210)/F(11,5)>70)
result=dict(status='PASS',scope='New whole q<=.06 phase8-conditioned signed chord; entire second W8 circle signed field210; all scalar n6251 bootstrap joins',precision_bits=ctx.prec,checks=len(G),theta_v=str(v),theta_m=str(m),dot_lower=str(dot),chord_constant=str(chord),Lipschitz_endpoint=str(lip),center_shift=str(shift),effective_quartic_radius=str(effective),curvature_at_cap=str(curvature),C0=str(C0),field_constant=str(field),outside_phase8_D=str(Dout2),low_middle_phase8_D=str(Dlow2),full_nonquartic_terms=terms,predicates=G,center_sha256=hashlib.sha256(raw).hexdigest(),signed_root_sha256=hashlib.sha256(signed.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'INDEPENDENT_W8_SIGNED_FIELD210_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['full_nonquartic_terms','predicates']},indent=2))
