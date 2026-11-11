#!/usr/bin/env python3
"""Complete High6251, conditional only on the named Origin6251 components.
Rebuilds all rational witnesses, local same-c Schur entrances, whole
contours, phase/product caps and exact binary cover; no discovery imports.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent
DATA=json.loads((HERE/'HIGH6251_INTEGER_DATA.json').read_text())
GATES=[]
def ck(name,l,r,op='>'):
 ok={'<':l<r,'>':l>r,'<=':l<=r,'>=':l>=r,'==':l==r}[op]
 GATES.append(dict(name=name,left=str(l),right=str(r),relation=op,passed=ok))
 if not ok:print('FAIL',name,float(l),op,float(r))
def elo(t):return sum((t**j/F(factorial(j)) for j in range(101)),F(0))
def eup(t):return sum((t**j/F(factorial(j)) for j in range(25)),F(0))+t**25/(F(factorial(25))*(1-t/26))
def ser(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,list):return [ser(z) for z in x]
 if isinstance(x,dict):return {k:ser(v) for k,v in x.items()}
 return x
N=F(6250);n=N+1;rho=F(30501,500000);A=rho/100
C70=F(596,5);C70100=F(1899,20);CF=F(2637,32);ca=F(5,4);cb=F(23,2)
ck('declared degree',N,F(DATA['N0']),'==');ck('declared A',A,F(DATA['A']),'==');ck('declared rho',rho,F(DATA['rho']),'==')
# The entire small-x high-modulus branch maps into the independently proved T cap.
xsmall=F(433,10000);rhi=F(31,500)
ck('small-x exact x>.0433 seam',F(47,50)**2,xsmall*xsmall+(1-rho)**2)
ck('small-x exact tau<.062 seam',F(47,50)**2-rho*rho,(1-rhi)**2)
ck('small-x complete normalized X>.69',xsmall/rhi,F(69,100))
ck('small-x normalized rho lower belongs to cap',rho,F(1,100))
# The large-x route enters D70 by independent pure geometry.
c=rho/2
ck('large-x direct crescent D70',5*N*c/(14*(1-c+c*c)),F(70))
ck('whole high actual half shift',(F(47,50)-F(1,10))**2,F(1,2))
ck('fresh medium signed q2',F(6,5)+F(1250,24649)*120000*F(9,640)**3,ca,'<')
ck('fresh medium signed qmu',F(1250,24649)*225,cb,'<')
ck('complete medium field identity',2*(ca*F(81,16)+cb*F(9,4))+18,CF,'==')
ck('whole medium prior branch negative',ca*F(81,16)-(160-cb)*F(9,4)+320,F(0),'<')
ck('whole medium branch decreasing',2*ca*F(9,4)-(160-cb),F(0),'<')
ERR={1:F(1,2),2:F(1,7),3:F(1,20),4:F(1,50),5:F(1,140),6:F(1,400),7:F(1,1000),8:F(1,2900),9:F(1,8000),10:F(1,22000)}
for W,e in ERR.items():ck('entire exponential inverse/'+str(W),elo(F(W)),1/e)
CAPS=[]
def cap(name,Blo,xl,xu,a0,L,K,known,M=None):
 p=1-Blo*Blo;s=F(3,4)*a0;k=2*L/(N*a0);shift=F(0) if known is None else known/2;coef=4*K if known is None else 2*K
 ck(name+'/complete concave witness interval',xl,3*a0)
 ck(name+'/whole logarithm budget',elo(L),8/(a0*(xl-s)))
 ck(name+'/phase positive denominator',xl,shift)
 ck(name+'/whole convexity',3*(xl-s),s-shift)
 ck(name+'/complete linear interval',xu,xl,'>=')
 ends=[]
 for lab,x in [('left',xl),('right',xu)]:
  rem=1-p*k/(x-s)-coef/(N*(x-shift));rad=4*x*x*k/(x-s)
  ck(name+'/'+lab+' positive radical side',rem,F(0));ck(name+'/'+lab+' whole squared gate',rem*rem,rad)
  ends.append(dict(x=x,remaining=rem,radicand=rad,reserve=rem*rem-rad))
 if M is not None:
  ga=1/M-F(3,4)/(M*M);kq=2*L/(N*ga)
  ck(name+'/full same-C join',xu,M*a0,'==');ck(name+'/concave quadratic endpoint',ga,F(1,4),'<')
  ck(name+'/exact witness-Schur identity',ga*xu*xu,a0*(xu-s),'==')
  rem=1-p*kq/(xu*xu)-coef/(N*(xu-shift))
  ck(name+'/quadratic positive side',rem,F(0));ck(name+'/entire decreasing quadratic branch',rem*rem,4*kq)
 CAPS.append(dict(name=name,Blo=Blo,xmin=xl,xmax=xu,a0=a0,L=L,K=K,known=known,M=M,endpoints=ends))
initial=F(13,1000)
cap('initial full high a<.013',F(47,50),rho,F(1),initial,F(1879,200),F(21,4),None)
dd=rho-initial/2
ck('initial cap lawful independent D100',(N*dd-C70)/(2*(1+initial)),F(100))
ck('initial stronger D70 field positive',N-C70100/dd,F(0))
K0=F(21,4);entry=(N*dd*dd-2*K0*dd)/(2*rho*dd+1-F(473,500)**2)
ck('same-c entire derivative premise',N*dd,2*K0)
ck('whole modulus946 independent medium160 entry',entry,F(160))
SCHUR=[]
def schur(index,bl,aU,D0,kap,Mtarget,profile):
 name='local Schur '+str(index);cc=rho-aU/2;Dfull=(N*cc*cc-8*cc)/(2*rho*cc+1-bl*bl)
 ck(name+'/same-c whole derivative premise',N*cc,F(8))
 ck(name+'/whole actual D floor',Dfull,D0)
 ck(name+'/ratio above2',kap,F(2));ck(name+'/ratio below prior9/4',kap,F(9,4),'<')
 ck(name+'/entire smaller-root gate',ca*kap*kap-(D0-cb)*kap+2*D0,F(0),'<')
 last=F(0);weighted=F(0);rows=[]
 for j,nn,W in profile:
  t=F(j,200);nu=F(nn,10000);e=ERR[W];la=kap*t/(1-t);f=t-kap*(t*t/2+t**3/(3*(1-t)));B=kap*(1-t)/(1-(1+kap)*t)**3
  eps=F(1001,1000)*B/D0;g=F(1001,1000)*nu*B*f;pi=1/(1+la);amp=1/(1-la)
  vals=[('critical free',t,F(1),'<'),('lambda positive',la,F(0),'>'),('whole inverse injection',la,F(1),'<'),('positive image',f,F(0),'>'),('whole corner',W*W,(1-nu*nu)*(D0*f)**2,'<'),('complete left side',eps*(1+e)+g*e,1-amp*e,'<'),('complete right side',eps*(1+e)+g*e*e,pi-e,'<'),('both whole horizontal sides',eps+g*e,pi,'<'),('homotopy positive',pi,F(0),'>'),('homotopy upper',amp,F(1),'>='),('positive odd height',nu*D0*f,F(22,7),'>'),('positive Abel increment',nu*f,last,'>')]
  for lab,l,r,op in vals:ck(name+'/'+str(j)+'/'+lab,l,r,op)
  weighted+=(1-t)**2*(nu*f-last);last=nu*f
  rows.append(dict(j=j,t=t,nu=nu,W=W,lambda_=la,f=f,B=B,eps=eps,g=g,pinv=pi,amplitude=amp))
 G=2/(2-aU*aU);u=G*aU;qa=u/(1-u);ea=2*G/(1-u)**2
 ck(name+'/radial positive denominator',u,F(1),'<');ck(name+'/whole transport sign',qa,F(1),'<')
 co=N/n*(F(7,22)*weighted-2/D0);ck(name+'/positive Abel budget',co,F(0));J=ea/((1-qa)*co)
 ck(name+'/complete same-configuration Schur',J+F(1,2),Mtarget,'<')
 out=dict(name=name,Blo=bl,a_cap=aU,Dfull=Dfull,Dfloor=D0,kappa=kap,M=Mtarget,c=co,J=J,reserve=Mtarget-J-F(1,2),contours=rows)
 SCHUR.append(out);return out
known=initial;lastbl=F(473,500);MAP={F(47,50):(initial,None,F(21,4),C70100)}
ck('whole profile-cap matching lengths',F(len(DATA['profiles'])),F(len(DATA['caps'])),'==')
for i,(s0,c0) in enumerate(zip(DATA['profiles'],DATA['caps'])):
 bl=F(s0['Blo']);au=F(s0['a_cap']);D0=F(s0['D0']);kap=F(s0['kappa']);M=F(s0['M'])
 ck('sequential '+str(i)+'/prior witness cap',au,known,'==');ck('sequential '+str(i)+'/inherited modulus',bl,lastbl,'>=')
 schur(i+1,bl,au,D0,kap,M,s0['contours'])
 cblo=F(c0['Blo']);a0=F(c0['a0']);L=F(c0['L']);prior=F(c0['known'])
 ck('sequential '+str(i)+'/cap same modulus',cblo,bl,'==');ck('sequential '+str(i)+'/cap same Schur M',F(c0['M']),M,'==')
 ck('sequential '+str(i)+'/cap prior already established',prior,known,'==');ck('sequential '+str(i)+'/strict cap refinement',a0,prior,'<')
 cap('joint cap '+str(i+1),bl,rho,M*a0,a0,L,F(4),known,M)
 known=a0;lastbl=bl;MAP[bl]=(a0,M,F(4),CF)
# Full endpoint is purely quadratic, including p=0, with no Fine lower a needed.
E=DATA['endpoint'];eb=F(E['Blo']);em=F(E['M']);ek=F(E['known']);eL=F(E['L']);ga=1/em-F(3,4)/(em*em)
ck('endpoint inherits actual final witness cap',ek,known,'==');ck('endpoint inherits actual final Schur M',em,SCHUR[-1]['M'],'==')
ck('endpoint complete positive concave coefficient',ga,F(0));ck('endpoint correct concave endpoint',ga,F(1,4),'<')
ck('endpoint entire logarithm budget',elo(eL),8/(ga*rho*rho))
kq=2*eL/(N*ga);rem=1-(1-eb*eb)*kq/(rho*rho)-8/(N*(rho-ek/2))
ck('endpoint positive radical side',rem,F(0));ck('endpoint entire decreasing quadratic comparison',rem*rem,4*kq)
LEAVES=[];SPLITS=[]
def leaf(name,bl,bh,au,xl,xu,M,K,Cfield,r):
 H=F(r['H_floor']);budget=F(r['budget']);R=F(r['R']);fac=F(r['factor']);ell=F(r['ell']);nu=F(r['nu_numerator'],10000);W=r['W']
 d=xl-au/2;m=V=1/d;pmax=1-bl*bl;hh=(N-Cfield/d)/(2*(1+au));den=max(2*xl+pmax/(xl-au/2),2*xu+pmax/(xu-au/2));hp=(N-2*K/d)/den
 Hex=hh if r['H_mode']=='hole' else hp
 ck(name+'/legal field mode',F(int(r['H_mode'] in ['hole','phase'])),F(1),'==');ck(name+'/positive H floor',H,F(0));ck(name+'/whole actual field',Hex,H)
 ck(name+'/positive whole phase numerator',N-2*K/d,F(0));ck(name+'/entire field ray nonnegative intercept',Cfield/d if r['H_mode']=='hole' else 2*K/d,F(0))
 Cmin=A*(xl-F(3,4)*A)
 if M is not None:Cmin=max(Cmin,(1/M-F(3,4)/(M*M))*xl*xl)
 ck(name+'/actual Fine lower concave interval',xl,3*A);ck(name+'/positive full C lower',Cmin,F(0));ck(name+'/complete product logarithm',elo(budget),8/Cmin)
 la=V*R/(1-m*R);rim=R-V*(R*R/2+m*R**3/(3*(1-m*R)));B=V/((1-m*R)**2*(1-la)**3)
 eps=F(1001,1000)*B/H;g=F(1001,1000)*nu*B*rim;pi=1/(1+la);amp=1/(1-la);e=ERR[W]
 vals=[('critical free',m*R,F(1),'<'),('lambda positive',la,F(0),'>'),('whole inverse injection',la,F(1),'<'),('positive full image',rim,F(0),'>'),('nu positive',nu,F(0),'>'),('nu below1',nu,F(1),'<'),('full corner',W*W,(1-nu*nu)*(H*rim)**2,'<'),('complete left',eps*(1+e)+g*e,1-amp*e,'<'),('complete right',eps*(1+e)+g*e*e,pi-e,'<'),('both full horizontal sides',eps+g*e,pi,'<'),('homotopy positive',pi,F(0),'>'),('homotopy upper',amp,F(1),'>='),('positive odd height',nu*H*rim,F(22,7),'>')]
 for lab,l,rr,op in vals:ck(name+'/'+lab,l,rr,op)
 pm=1-bh*bh;T=au/(2*pm);U=R/pm;ll=1+2*T*xu;lb=fac*fac-(1-fac*fac*bl*bl)*(T*T+U*U)
 ck(name+'/entire half-disk radial sign',pm,R*R);ck(name+'/factor positive',fac,F(0));ck(name+'/factor below1',fac,F(1),'<')
 ck(name+'/entire rightmost half-disk support',lb,2*U*fac*fac*ll*bh);ck(name+'/complete factor logarithm',eup(ell),1/fac,'<')
 co=nu*F(7,22)*H*rim-2;bonus=Cmin*H*H/(2*N);prod=ell*co+bonus
 ck(name+'/positive original-root count',co,F(0));ck(name+'/entire same-C product plus field energy',prod,budget)
 LEAVES.append(dict(name=name,Blo=bl,Bhi=bh,a_cap=au,xmin=xl,xmax=xu,M=M,K=K,field_constant=Cfield,Hfloor=H,H_exact=Hex,H_mode=r['H_mode'],d=d,Cmin=Cmin,energy_bonus=bonus,R=R,factor=fac,ell=ell,nu=nu,W=W,count=co,product=prod,budget=budget,reserve=prod-budget))
def visit(name,bl,bh,au,xl,xu,M,K,Cfield,node):
 ck(name+'/positive full modulus width',bh,bl);ck(name+'/positive full x width',xu,xl)
 if 'axis' not in node:leaf(name,bl,bh,au,xl,xu,M,K,Cfield,node);return
 axis=node['axis'];ck(name+'/exact declared split dimension',F(int(axis in ['x','b'])),F(1),'==');ck(name+'/both full children present',F(len(node['children'])),F(2),'==')
 SPLITS.append(dict(name=name,Blo=bl,Bhi=bh,xmin=xl,xmax=xu,axis=axis))
 if axis=='x':
  mid=(xl+xu)/2;visit(name+'L',bl,bh,au,xl,mid,M,K,Cfield,node['children'][0]);visit(name+'R',bl,bh,au,mid,xu,M,K,Cfield,node['children'][1])
 else:
  mid=(bl+bh)/2;visit(name+'L',bl,mid,au,xl,xu,M,K,Cfield,node['children'][0]);visit(name+'R',mid,bh,au,xl,xu,M,K,Cfield,node['children'][1])
lastB=F(47,50)
for i,r in enumerate(DATA['mothers']):
 bl,bh,au,xl,xu=[F(r[k]) for k in ['Blo','Bhi','a_cap','xmin','xmax']];M=None if r['M'] is None else F(r['M']);K=F(r['K']);cf=F(r['CF']);expected=MAP[bl]
 ck('mother '+str(i)+'/entire closed modulus seam',bl,lastB,'==');ck('mother '+str(i)+'/lawful inherited cap',au,expected[0],'==')
 ck('mother '+str(i)+'/lawful same M',F(-1) if M is None else M,F(-1) if expected[1] is None else expected[1],'==')
 ck('mother '+str(i)+'/lawful phase',K,expected[2],'==');ck('mother '+str(i)+'/lawful field constant',cf,expected[3],'==')
 ck('mother '+str(i)+'/whole x lower',xl,rho,'==');ck('mother '+str(i)+'/whole x upper',xu,F(1) if M is None else M*au,'==')
 visit('B'+str(i+1),bl,bh,au,xl,xu,M,K,cf,r['tree']);lastB=bh
ck('full terminal endpoint seam',lastB,eb,'==');ck('complete frozen leaf count',F(len(LEAVES)),F(DATA['leaf_count']),'==');ck('complete frozen split count',F(len(SPLITS)),F(DATA['split_count']),'==')
fails=[r['name'] for r in GATES if not r['passed']]
out=dict(status='PASS' if not fails else 'FAIL',gate_count=len(GATES),full_Schur_profiles=len(SCHUR),full_Schur_contours=sum(len(r['contours']) for r in SCHUR),whole_mothers=len(DATA['mothers']),exact_split_nodes=len(SPLITS),complete_closed_leaves=len(LEAVES),full_physical_contours=len(LEAVES),phase_product_cap_rows=len(CAPS),endpoint=dict(Blo=eb,M=em,known=ek,L=eL,kq=kq,positive_side=rem,squared_reserve=rem*rem-4*kq),failures=fails,Schur_profiles=SCHUR,caps=CAPS,leaves=LEAVES,split_nodes=SPLITS,predicates=GATES,data_sha256=hashlib.sha256((HERE/'HIGH6251_INTEGER_DATA.json').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),small_x_input='Lemma sm:T-highX-extension; verification component small_extensions_highx')
(HERE/'HIGH6251_COMPLETE_EXACT.json').write_text(json.dumps(ser(out),indent=2)+'\n')
summary={k:out[k] for k in ['status','gate_count','full_Schur_profiles','full_Schur_contours','whole_mothers','exact_split_nodes','complete_closed_leaves','full_physical_contours','phase_product_cap_rows','failures']}
summary['minimum_product_reserve']=float(min(r['reserve'] for r in LEAVES));summary['endpoint_squared_reserve']=float(rem*rem-4*kq)
print(json.dumps(summary,indent=2))
if fails:raise SystemExit(1)
