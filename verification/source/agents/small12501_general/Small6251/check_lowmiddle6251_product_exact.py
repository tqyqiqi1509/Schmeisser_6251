#!/usr/bin/env python3
"""Complete rational Low/Middle6251 product complement of the T-source cells.
No discovery functions are imported. Rebuild every contour, whole factor,
full logarithm, witness interval and set-coverage seam from exact inputs.
"""
from fractions import Fraction as F
from math import factorial,isqrt
from pathlib import Path
import json,hashlib,sys
sys.set_int_max_str_digits(0)
O=Path(__file__).resolve().parent
raw=(O/'LOW_MIDDLE6251_PRODUCT_PROPOSAL.json').read_bytes();DATA=json.loads(raw)
N=F(6250);n=N+1;rho=F(30501,500000);A=rho/100;CF=F(2637,32)
G=[];CON=[];ROWS=[]
def ck(name,p):
 p=bool(p);G.append(dict(name=name,passed=p))
 if not p:raise AssertionError(name)
def elo(t):return sum((t**j/F(factorial(j)) for j in range(101)),F(0))
def eup(t):
 ck('whole upper exp domain',0<=t<26)
 return sum((t**j/F(factorial(j)) for j in range(25)),F(0))+t**25/(F(factorial(25))*(1-t/26))
def sqrtup(t,name,den=10**12):
 k=isqrt(t.numerator*den*den//t.denominator)
 if F(k,den)**2<t:k+=1
 k+=1;z=F(k,den)
 ck(name+'/whole strict outward radical',z*z>t and (z-F(2,den))**2<t)
 return z
ln2low=2*sum((F(1,3)**(2*j+1)/(2*j+1) for j in range(25)),F(0))
def logweight(f,name):
 ck(name+'/whole factor contracts',0<f<1)
 zz=f;k=0
 while zz<F(1,2):zz*=2;k+=1
 t=(1-zz)/(1+zz)
 ck(name+'/whole dyadic logarithm interval',F(1,2)<=zz<1 and 0<t<=F(1,3))
 low=k*ln2low+2*sum((t**(2*j+1)/(2*j+1) for j in range(16)),F(0))
 w=F(low.numerator*10**10//low.denominator,10**10)
 ck(name+'/whole logarithm strict lower',0<w<=low<w+F(1,10**10))
 return w
ERR={1:F(1,2),2:F(1,7),3:F(1,20),4:F(1,50),5:F(1,140),6:F(1,400),7:F(1,1000),8:F(1,2900),9:F(1,8000),10:F(1,22000)}
for W,e in ERR.items():ck('full exponential width '+str(W),elo(F(W))>1/e)
ck('exact declared global parameters',F(DATA['A'])==A and F(DATA['rho'])==rho and DATA['N0']==N and DATA['radius_denominator']==5000)
ck('independent full energy variance floors',n>=2048)
ck('whole lawfully entered D100',(N*F(3,50)-F(596,5))/F(11,5)>100)
qa=N*N/4-4*139**2;qb=-4*139*(F(1899,20)+278);qc=N*N*F(3,50)**2-(F(1899,20)+278)**2
ck('whole quadrature D139 join',qa>0 and 4*qa*qc>qb*qb)
ck('full medium139 rectangle matches prior ratio',F(9,4)/139==F(9,556))
ck('whole upper-strip actual half shift',(F(41,50)-F(1,10))**2>F(1,2))
def cloud(name,d,V,H,R,nu,W):
 m=1/d;la=V*R/(1-m*R);er=V*(R*R/2+m*R**3/(3*(1-m*R)));rim=R-er
 Bi=V/((1-m*R)**2*(1-la)**3);eps=F(1001,1000)*Bi/H;gg=F(1001,1000)*nu*Bi*rim;pin=1/(1+la);amp=1/(1-la);ee=ERR[W]
 for lab,pred in [('critical free',0<m*R<1),('full injection',0<la<1),('whole image positive',rim>0),('nu range',0<nu<1),('every corner',W*W<(1-nu*nu)*(H*rim)**2),('entire left side',eps*(1+ee)+gg*ee<1-amp*ee),('entire right side',eps*(1+ee)+gg*ee*ee<pin-ee),('both whole horizontals',eps+gg*ee<pin),('separate amplitude',0<pin<1<amp),('positive odd height',nu*H*rim>F(22,7))]:ck(name+'/'+lab,pred)
 sl=nu*F(7,22)*rim;count=sl*H-2
 ck(name+'/positive finite actual count',count>0)
 CON.append(dict(name=name,d=str(d),V=str(V),H=str(H),R=str(R),nu=str(nu),W=W,lambda_upper=str(la),whole_log_tail=str(er),image=str(rim),inverse_derivative=str(Bi),slope=str(sl),count=str(count)))
 return sl,count
# Low component: the degree-dependent adaptive shift is retained.
bl=F(11,100);dl=1-bl;Hl=N*dl/F(9,4);sl,L=cloud('Low',dl,F(9,4)/dl,Hl,F(1,10),F(99,100),8)
ck('Low whole adaptive denominator',1-(bl+F(1,10))*bl>0)
ck('Low whole adaptive factor',F(1,5)<F(1,4)*(1-(bl+F(1,10))*bl))
ck('Low complete factor logarithm',eup(F(1))<4)
ck('Low whole count/n degree ray',L/n>F(1,100))
ck('Low complete variance comparison',1-bl*bl>F(98,100) and F(49,100)>F(2,5))
ck('Low whole adaptive coefficient ray',A*A*F(4,25)*(1-F(4,25)/(n*n))>F(1,20000000))
ck('Low complete logarithm budget',elo(F(20))>8*20000000 and elo(F(10))>n)
ck('Low entire unbounded degree contradiction',n/100>40 and F(1,100)>2/n)
ROWS.append(dict(name='Low',mode='variance',count=str(L),reserve=str(L/n-F(1,100))))
# Middle: rightmostness improves the full adaptive numerator to R^2+a^2.
bm=F(41,50);dm=1-bm;Hm=(N-CF/dm)/F(11,5);Rm=F(11,250);facm=F(3,8);ellm=F(49,50)
sl,L=cloud('Middle',dm,F(9,4)/dm,Hm,Rm,F(99,100),6);denm=1-(bm+Rm)*bm
ck('Middle positive whole adaptive denominator',denm>0)
ck('Middle complete rightmost adaptive numerator',Rm*Rm+F(1,100)<facm*facm*denm*denm)
ck('Middle full factor logarithm',eup(ellm)<1/facm)
ck('Middle complete concave first endpoint',A*(F(11,100)-A)>F(1,15000))
ck('Middle complete concave other endpoint',F(1,10)*F(1,100)>F(1,15000))
ck('Middle whole radial factor',1-F(11,100)**2>F(9,10))
Cm=F(9,10)/15000**2
ck('Middle complete logarithm target',elo(F(43,2))>8/Cm)
ck('Middle complete product contradiction',ellm*L>F(43,2))
ROWS.append(dict(name='Middle',mode='adaptive rightmost numerator',count=str(L),factor=str(facm),ell=str(ellm),target='43/2',reserve=str(ellm*L-F(43,2))))
seams={'upper':A,'highx':A,'lowy':A};expected_counts={'upper':7,'highx':11,'lowy':5};counts={k:0 for k in seams}
for z in DATA['rows']:
 name=z['name'];mode=z['mode'];ck(name+'/declared complete mode',mode in seams);counts[mode]+=1
 al,au,bl,bh,xu,d,V,H,Cmin,target=[F(z[k]) for k in ['a_lo','a_hi','Blo','Bhi','xmax','d','V','H','Cmin','target']]
 ck(name+'/full closed witness seam',al==seams[mode] and al<au);seams[mode]=au
 ck(name+'/lawful whole witness range',A<=al<au<=F(1,10))
 if mode=='upper':
  ck(name+'/complete upper modulus cell',bl==F(41,50) and bh==F(37,40) and xu==bh)
  ck(name+'/complete radial d floor',d==1-bh)
  ck(name+'/whole prior ratio V',V>=F(9,4)/d)
  ck(name+'/whole source field',H<(N-CF/d)/(2*(1+au)))
  ck(name+'/same-C lower',Cmin==al*al/4)
 elif mode=='highx':
  ck(name+'/complete high-x modulus cell',bl==F(37,40) and bh==F(47,50) and xu==bh)
  ck(name+'/whole crescent d floor',d>0 and (d<=1-bh or d<=rho-au/2 or d*d<=(1-bh)**2+al*al/4))
  ck(name+'/same-x positive horizontal gap',rho>au/2)
  ck(name+'/whole actual ratio and direct V',V>=F(9,4)/d or V>=1/(rho-au/2))
  ck(name+'/whole source field',H<(N-CF/d)/(2*(1+au)))
  if au<=rho:cm=min(al*rho-F(3,4)*al*al,au*rho-F(3,4)*au*au)
  elif al>=rho:cm=al*al/4
  else:cm=min(al*rho-F(3,4)*al*al,rho*rho/4)
  ck(name+'/full piecewise-concave same C',Cmin==cm)
 else:
  yy=F(5,4)
  ck(name+'/complete low-y cell',bl==F(37,40) and bh==F(47,50) and xu==rho and au<=rho)
  ck(name+'/whole low-y radius and reciprocal bound',d==F(3,50) and V>=yy/d)
  ck(name+'/whole same-source low-y phase field',H<(N-8*yy/d)/(2*rho+(2-d)*yy))
  ck(name+'/same-C lower',Cmin==al*al/4)
 ck(name+'/positive complete field',H>0)
 ck(name+'/full positive same-C and logarithm',Cmin>0 and elo(target)>8/Cmin)
 pm=1-bh*bh;T=au/(2*pm);LL=1+2*T*xu
 ws=[];cs=[];ss=[];prevR=F(0);prevcount=F(0);factors=[]
 for j,nn,W in z['contours']:
  R=F(j,5000);nu=F(nn,10000);label=name+'/'+str(j);slope,count=cloud(label,d,V,H,R,nu,W)
  ck(label+'/whole positive half-disk radial term',R*R<pm)
  U=R/pm;den=1+bl*bl*(T*T+U*U)-2*U*LL*bh
  ck(label+'/positive full factor denominator',den>0)
  ff=sqrtup((T*T+U*U)/den,label)
  lb=ff*ff-(1-ff*ff*bl*bl)*(T*T+U*U)
  ck(label+'/entire rightmost half disk',lb>2*U*ff*ff*LL*bh)
  w=logweight(ff,label)
  ck(label+'/nested whole disks and increasing counts',R>prevR and count>prevcount)
  if ws:ck(label+'/nonincreasing complete weights',ws[-1]>=w)
  prevR,prevcount=R,count;ws.append(w);cs.append(count);ss.append(slope);factors.append(str(ff))
 ck(name+'/nonempty entire profile',bool(ws))
 budget=sum((w*(c-p) for w,c,p in zip(ws,cs,[F(0)]+cs[:-1])),F(0))
 abel=ws[-1]*cs[-1]+sum(((ws[j]-ws[j+1])*cs[j] for j in range(len(ws)-1)),F(0))
 ck(name+'/complete finite Abel identity',budget==abel)
 Wstar=sum((w*(s-p) for w,s,p in zip(ws,ss,[F(0)]+ss[:-1])),F(0))
 ck(name+'/single complete initial count loss',budget==H*Wstar-2*ws[0] and Wstar>0)
 bonus=N*Cmin/2
 ck(name+'/complete product plus same-C energy contradiction',budget+bonus>target)
 ROWS.append(dict(name=name,mode=mode,a_lo=str(al),a_hi=str(au),Blo=str(bl),Bhi=str(bh),xmax=str(xu),d=str(d),V=str(V),H=str(H),Cmin=str(Cmin),contours=len(ws),budget=str(budget),bonus=str(bonus),target=str(target),reserve=str(budget+bonus-target),weights=list(map(str,ws)),counts=list(map(str,cs)),factors=factors))
for mode,last in [('upper',F(1,10)),('highx',F(1,10)),('lowy',rho)]:
 ck(mode+'/entire final witness endpoint',seams[mode]==last and counts[mode]==expected_counts[mode])
# Complete normalized-source remainder mapping; the source contradiction is
# the explicitly separate new full certificate, not inferred by this file.
ck('full small-x normalized tau lower',1-F(47,50)==F(3,50))
ck('full small-x normalized tau upper',F(37,40)**2-rho*rho>(1-F(39,500))**2)
ck('full small-x normalized X at tau<=rho',F(47,50)**2<rho*rho+(1-rho)**2)
ck('full normalized rho interval monotonicity',rho<F(1,2) and F(3,50)<rho<F(39,500))
ck('whole low-y/high-y seam',F(5,4)<F(9,4))
out=dict(status='PASS',scope='All n>=6251,A=rho0/100<=a<=.1,|b|<=.94 except precisely the complete normalized T source rectangle rho[.06,.078],y[1.25,2.25],X<=1; the latter is an explicit separate contradiction input',checks=len(G),complete_physical_contours=len(CON),whole_product_rows=len(DATA['rows']),N0=str(N),A=str(A),rows=ROWS,contours=CON,predicates=G,proposal_sha256=hashlib.sha256(raw).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(O/'LOW_MIDDLE6251_PRODUCT_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['rows','contours','predicates']},indent=2))
print('reserves',[(r['name'],float(F(r['reserve']))) for r in ROWS])
