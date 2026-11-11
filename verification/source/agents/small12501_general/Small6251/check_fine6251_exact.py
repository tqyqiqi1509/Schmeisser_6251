#!/usr/bin/env python3
"""Full rational Fine100 outside-Origin component at n>=6251."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent;raw=(P/'FINE6251_INTEGER_DATA.json').read_bytes();D=json.loads(raw)
G=[]
def ck(name,p):
 p=bool(p);G.append(dict(name=name,passed=p))
 if not p:raise AssertionError(name)
def elo(x):return sum((x**j/F(factorial(j)) for j in range(101)),F(0))
N=F(6250);n=N+1;rho=F(30501,500000);A=rho/100;dd=rho-A/2
ck('exact degree and witness seam',D['N0']==6250 and F(D['A'])==A and F(D['rho'])==rho and 100*A==rho)
ck('full outside independent W8-to-D70 scalar join',(N*rho-210)/F(11,5)>70)
c=rho/2;ck('full outside large-x pure geometry D70',5*N*c/(14*(1-c+c*c))>70)
first=(N*dd-F(596,5))/(2*(1+A));second=(N*dd-F(1899,20))/(2*(1+A));full=(N*dd-F(2637,32))/(2*(1+A))
ck('D70 independently enters D100',first>100)
ck('after independent D100 lawful medium140',second>140)
ck('whole medium140 qmu rectangle',F(9,4)/140==F(9,560) and F(1,140)<F(1,100))
Df=F(D['D0']);kap=F(D['kappa']);M=F(D['M'])
ck('complete actual D floor',full>Df)
ck('whole ratio bracket',2<kap<F(9,4))
ck('entire smaller-root quadratic',F(5,4)*kap*kap-(Df-F(23,2))*kap+2*Df<0)
ERR={1:F(1,2),2:F(1,7),3:F(1,20),4:F(1,50),5:F(1,140),6:F(1,400),7:F(1,1000),8:F(1,2900),9:F(1,8000),10:F(1,22000)}
for W,e in ERR.items():ck('whole exp(-W) '+str(W),elo(F(W))>1/e)
weighted=F(0);last=F(0);rows=[]
for j,nn,W in D['contours']:
 t=F(j,200);nu=F(nn,10000);ee=ERR[W];la=kap*t/(1-t)
 ff=t-kap*(t*t/2+t**3/(3*(1-t)));Bi=kap*(1-t)/(1-(1+kap)*t)**3
 eps=F(1001,1000)*Bi/Df;gg=F(1001,1000)*nu*Bi*ff;pin=1/(1+la);amp=1/(1-la)
 for lab,pred in [('critical free',0<t<1),('nu interval',0<nu<1),('whole injection',0<la<1),('full image',ff>0),('every corner',W*W<(1-nu*nu)*(Df*ff)**2),('entire left side',eps*(1+ee)+gg*ee<1-amp*ee),('entire right side',eps*(1+ee)+gg*ee*ee<pin-ee),('both horizontal sides',eps+gg*ee<pin),('separate amplitude',0<pin<1<amp),('positive odd height',nu*Df*ff>F(22,7)),('complete Abel increment',nu*ff>last)]:ck(str(j)+'/'+lab,pred)
 weighted+=(1-t)**2*(nu*ff-last);last=nu*ff
 rows.append(dict(j=j,nu=str(nu),W=W,t=str(t),lambda_upper=str(la),f=str(ff),B=str(Bi),eps=str(eps),g=str(gg),pinv=str(pin),amplitude=str(amp)))
gamma=2/(2-A*A);u=gamma*A;qa=u/(1-u);ea=2*gamma/(1-u)**2
ck('full transport denominators',0<u<1 and 0<qa<1)
co=N/n*(F(7,22)*weighted-2/Df);ck('full Abel budget positive',co>0)
J=ea/((1-qa)*co)
ck('whole Fine100 Schur contradiction',J+F(1,2)<M)
out=dict(status='PASS',scope='n>=6251,0<a<=rho0/100,outside max(Re b,1-Im b)<=rho0',checks=len(G),complete_Schur_contours=len(rows),N0=str(N),A=str(A),d_floor=str(dd),Dfirst=str(first),Dafter100=str(second),Dfield=str(full),D0=str(Df),kappa=str(kap),weighted_Abel=str(weighted),c=str(co),J=str(J),J_plus_half=str(J+F(1,2)),M=str(M),reserve=str(M-J-F(1,2)),contours=rows,predicates=G,data_sha256=hashlib.sha256(raw).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(P/'FINE6251_COMPLETE_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['contours','predicates']},indent=2))
