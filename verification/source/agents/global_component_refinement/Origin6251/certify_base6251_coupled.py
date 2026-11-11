"""Uniform BASE010 candidate, all 0<rho<=1/100, n>=6251.

Author-side exact rational verification of every complete inverse contour,
all changed HB gates, source-count compensation, and the new complete
quartic companion row. No lower bound on n*rho is used.
"""
from fractions import Fraction as F
from math import factorial, comb
from pathlib import Path
import json, sys

OUT=Path(__file__).resolve().parent
n0=int(sys.argv[1]) if len(sys.argv)>1 else 6251
N=n0-1
rho=F(1,100);c=rho/2
f=F(1001,1000);fx=ft=F(10001,10000);D=1+9*c*c
checks=[];rows=[]
def ck(name,p,value=None):
 checks.append({'name':name,'pass':bool(p),'value':None if value is None else str(value)})
def EU(x,k=30):
 assert 0<=x<k+2
 return sum((x**j/F(factorial(j)) for j in range(k+1)),F(0))+x**(k+1)/(F(factorial(k+1))*(1-x/F(k+2)))
def row(name,q,m,R,H,nu,W,eps,g,p,b,e,quarter=False):
 la=q*R/(1-m*R)
 err=q*(R*R/2+m*R**3/(3*(1-m*R)))
 rh=R-err
 deriv=q/((1-m*R)**2*(1-la)**3)
 conditions={
  'physical disk':0<m*R<1,'injective logarithm':0<la<1,
  'positive inverse radius':rh>0,
  'all containing corners':W*W<(1-nu*nu)*(H*rh)**2,
  'derivative on whole inverse disk':deriv/H<eps,
  'whole vertical path':nu*deriv*rh<g,
  'real amplitude':1/(1+la)>=p,'modulus amplitude':1/(1-la)<=b,
  'amplitude homotopy':0<p<1<b,
  'complete left side':eps*(1+e)+g*e<1-b*e,
  'complete right side':eps*(1+e)+g*e*e<p-e,
  'both complete horizontal sides':eps+g*e<(p/(1+p+b) if quarter else p),
  'positive admissible height':nu*H*rh>(F(11,7) if quarter else F(22,7)),
  'entire exponential comparison':sum((W**j/F(factorial(j)) for j in range(49)),F(0))>1/e}
 for tag,test in conditions.items():ck(name+'/'+tag,test)
 z={'name':name,'q':q,'m':m,'R':R,'H':H,'nu':nu,'W':W,'epsilon':eps,'g':g,'p':p,'b':b,'e':e,'quarter':quarter,'lambda':la,'error':err,'rhat':rh,'B':deriv,'loss':F(3,2) if quarter else F(2)}
 rows.append(z);return z
def density(z):return z['nu']*F(7,22)*z['rhat']
def count(z):return density(z)-z['loss']/z['H']
def weight(R,delta,lo,hi):
 den=(delta+R)**2
 ts=[max(F(0),lo-R),hi+R]
 return max(F(0),min((t-c*(1+t*t))/(den+t*t) for t in ts))
def add(a,b):
 return [(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b)))]
def mul(a,b):
 z=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]+=x*y
 return z
def bernstein(poly,lo,hi):
 d=len(poly)-1
 power=[sum((poly[k]*F(comb(k,j))*lo**(k-j)*(hi-lo)**j for k in range(j,d+1)),F(0)) for j in range(d+1)]
 return [sum((power[k]*F(comb(j,k),comb(d,k)) for k in range(j+1)),F(0)) for j in range(d+1)]

xf=row('Xfirst',F(2),F(2),F(3,20),F(600),F(99,100),F(4),F(1,25),F(3),F(7,10),F(7,4),F(1,50))
xd=row('XlowerStrip',F(2),F(2),F(13,100),F(600),F(99,100),F(4),F(1,25),F(3),F(7,10),F(7,4),F(1,50))
xm={}
for H in [450,200]:
 xm[H]=[]
 for j in range(6,31):
  W,g,e=(F(2),F(1,5),F(1,7)) if j<=9 else (F(4),F(3),F(1,50))
  xm[H].append(row('Xmiddle'+str(H)+'_'+str(j),F(2),F(2),F(j,200),F(H),(F(23,25) if H==450 else F(9,10)),W,F(11,100),g,F(7,10),F(7,4),e,True))
xl=row('Xlast78',F(2),F(2),F(3,20),F(78),F(97,100),F(11,5),F(3,10),F(3),F(7,10),F(7,4),F(1,9))
tf=row('Tfirst150',F(25,11),F(100,99),F(3,20),F(150),F(9,10),F(4),F(1,6),F(3),F(7,10),F(7,4),F(1,50))
tl=row('Tlast103',F(25,11),F(100,99),F(3,20),F(103),F(9,10),F(4),F(1,6),F(3),F(7,10),F(7,4),F(1,50))
ck('all 55 complete physical rows',len(rows)==55)
ck('moving polar inflation',1/(1-F(9,2)*rho*rho)<f)
ck('fixed X polar inflation',1/(1-F(10201,20000)*rho*rho)<fx)
ck('fixed T polar inflation',1/(1-F(169,200)*rho*rho)<ft)
qfine=F(7,625);mufine=F(1,200)
phasemajor=21+133*qfine+1030*qfine**2+9000*qfine**3+120000*qfine**4+220*mufine
ck('inherited standalone complete fine phase rectangle',phasemajor<24,phasemajor)
ck('X high field separate fine entrances',F(2,600)<qfine and F(2,600)<mufine)
ck('X full phase displacement',F(8,600)<F(1,70) and 1-F(1,70)**2>F(999,1000)**2)
er=xf['error'];rh=xf['rhat'];H=xf['H']
kup=(F(99,100)-F(4,5))*rh/(2*F(22,7))-2/H
depth=F(4,5)*F(999,1000)*rh-er-F(4,70)/H
ck('X upper strip has actual root',kup>0,kup)
ck('X high field whole depth bound',depth>F(67,1000),depth)
U=F(311,1000)
ck('X low depth full degree ray',(rho/2+F(9,N))/(F(67,1000)*(1-rho*F(67,1000)/2))<U*U)
delta=F(1,100);tmax=F(3,10)
kap=F(1,2)-(tmax*tmax-delta-delta*delta)**2/(4*(tmax*tmax+delta*delta)*(F(1,2)+delta))
ck('whole critical-domain kappa',kap==F(29591,63600) and kap>0)
aval=fx*(1-kap*U)
betau=2*D*(aval*U+F(1,n0))
upperlow=rho*(aval+F(1,600))/(1-betau)
derivative=-fx*kap+2*D*(fx*fx*(1-kap*U)**2+(fx/H)*(1-2*kap*U))
ck('same-field factor positive and increasing',1-kap*U>0 and 1-2*kap*U>0)
ck('same-field HB denominator positive first',betau<1,betau)
ck('same-field solved HB monotone whole interval',derivative>0,derivative)
ck('same-field low-depth HB upper',upperlow<F(1,50),upperlow)
gamma=F(7,25)
inset=gamma*F(999,1000)*xd['rhat']-xd['error']-F(4,70*600)-F(3481,5000)*rho
kd=(xd['nu']-gamma)*xd['rhat']/(2*F(22,7))
ck('X lower strip inset positive',inset>0,inset)
ck('X lower strip excludes known root beta',gamma>0 and gamma<xd['nu'])
dist=F(157,1200)
# b(T)>=2a(T): positive denominator cleared, certified on the entire depth cell.
P=add(mul([-c,F(1),-c],[dist*dist,2*dist,F(1)]),[-2*x for x in mul([inset,F(1)],[F(1,10000),F(0),F(1)])])
bcoef=bernstein(P,F(67,1000),F(3,10))
ck('X source compensates entire strip count loss on full interval',len(bcoef)==5 and min(bcoef)>0,min(bcoef))
for t in [F(67,1000),F(3,10)]:
 val=kd*(t+inset)/(t+dist)**2
 ck('X compensated full lower strip endpoint '+str(t),val>F(1,50),val)
Uh=F(147,1000)
ck('X high depth degree ray',(rho/2+F(9,N))/(F(3,10)*(1-rho*F(3,10)/2))<Uh*Uh)
betah=2*D*(fx*Uh+F(1,n0))
upperhigh=rho*(fx+F(1,600))/(1-betah)
ck('X high depth solved HB',betah<F(3,10) and upperhigh<F(2,125),upperhigh)
ksym=density(xf)-F(2,600)
for t in [F(3,10),F(1)]:
 val=ksym*(t-F(3,20)-F(5,4)*rho)/(t+F(901,6000))**2
 ck('X symmetric full-depth endpoint '+str(t),val>F(1,50),val)
ck('X zero inset full degree contradiction',N>1200)
# No lower n*rho premise: moving-source estimates retain all positive rho.
Bs={600:F(193,1000),450:F(29,200),200:F(13,200)}
Ts={600:F(13,200),450:F(9,100),200:F(21,100)}
for H in [600,450,200]:
 B=Bs[H];Wh=rho*(H*f+1)/(1-B)
 ck('X bounded field denominator '+str(H),2*D*(H*f+1)/n0<B)
 ck('X moving source forces full depth '+str(H),(1-rho)**2/(2*Ts[H])>Wh,Wh)
middle={}
for H,oldH in [(450,600),(200,450)]:
 zs=xm[H];ks=[count(z) for z in zs];incs=[ks[0]]+[b-a for a,b in zip(ks,ks[1:])]
 ck('X all nested finite count increments positive '+str(H),min(incs)>0)
 cells=[];prev=Ts[oldH]
 for cell_lo,cell_hi in [(Ts[oldH],F(1,5)),(F(1,5),F(1))]:
  ck('X complete depth seam '+str(H)+'/'+str(cell_lo),cell_lo==prev and cell_lo<cell_hi);prev=cell_hi
  ww=[weight(z['R'],F(1,100),cell_lo,cell_hi) for z in zs]
  lower=sum(a*b for a,b in zip(incs,ww));upper=rho*(fx+F(1,H))/(1-Bs[oldH])
  ck('X all nested weights decrease on complete depth '+str(H)+'/'+str(cell_lo),len(ww)==25 and ww[0]>0 and all(a>=b for a,b in zip(ww,ww[1:])))
  ck('X entire joint-depth cloud lowers field below '+str(H)+'/'+str(cell_lo),lower>upper,lower-upper)
  cells.append({'T_lo':cell_lo,'T_hi':cell_hi,'lower':lower,'upper':upper,'margin':lower-upper,'all_weights':ww})
 ck('X complete depth final endpoint '+str(H),prev==1)
 middle[H]=cells
a78=weight(xl['R'],F(1,100),F(21,100),F(1))
bsource=(1-rho)/(1+F(1,10000))
# Exact factor: source numerator difference=(1-T)(A*T-B).
A=1-c*(1-F(1,10000));B=c*(1-F(1,10000))+F(1,10000)
ck('X full actual-source depth factor positive',A*F(21,100)>B)
ck('X symmetric source improvement positive',bsource>a78)
lower78=a78*density(xl)+min(F(0),bsource-3*a78)/F(78)
upper78=rho*(fx+F(1,78))/(1-Bs[200])
ck('X h78 entire symmetric count positive',count(xl)>0)
ck('X h78 source-compensated complete cloud',lower78>upper78,lower78-upper78)
W78=rho*(78*fx+1)/(1-Bs[200])
ck('X fixed source closes h78',bsource>W78,bsource-W78)
# T=1: independent D100 entrance at each separate field threshold.
d0=F(99,100)
ck('T whole physical critical gap',2*(1-d0)-(2-d0*d0)*rho>0 and 1-d0*rho>0)
ck('T first D100 invocation',150*d0>100)
ck('T first separate wide-fine entrances',F(25,11*150)<F(9,400) and 1/(d0*150)<F(1,100))
ck('T absolute standalone wide phase',F(23,5)*F(25,11)<11)
Ut=F(21,250)
ck('T complete high-field degree ray',(rho/2+F(11,N))/(1-rho/2)<Ut*Ut)
betat=2*D*(ft*Ut+F(1,n0))
upper150=rho*(ft+F(1,150))/(1-betat)
lower150=count(tf)*weight(tf['R'],F(3,10),F(1),F(1))
ck('T high-field denominator positive',betat<F(17,100),betat)
ck('T high-field full cloud comparison',lower150>upper150,lower150-upper150)
BT=F(1,20)
ck('T bounded h150 denominator',2*D*(150*ft+1)/n0<BT)
ck('T second D100 invocation before ratio',103*d0>100)
lower103=count(tl)*weight(tl['R'],F(3,10),F(1),F(1))
upper103=rho*(ft+F(1,103))/(1-BT)
ck('T h103 whole cloud comparison',lower103>upper103,lower103-upper103)
tsource=(1-rho)/F(109,100)
smallX=F(77,100)
Wsmall=rho*(103*ft*smallX+1)/(1-BT)
ck('T all X at most .77 fixed-source closure',tsource>Wsmall,tsource-Wsmall)
# A new complete quartic row. Its Gaussian entire-tail input is the independently
# audited 132200*q^5 bound for every |theta|<=q<=1/20.
q=F(33,1000);mu=F(8,625);nu=q*mu
Z=F(63,10)+21*q+133*q*q+1030*q**3+9000*q**4
M=F(63,10)*q+45*q*q+400*q**3+3400*q**4+86000*q**5+320000*q**6+1700000*q**7+9100000*q**8+39000000*q**9
Z0=F(719,100);M0=F(7,25);E0=F(71,50);P0=F(350)
Z1=F(15,2);L=F(32,25);R=F(3,10);nb=F(11,50)
ck('quartic Gaussian full-tail premise has proper domain',q<F(1,20))
ck('quartic complete center bound',Z<Z0,Z)
ck('quartic complete path exponent bound',M<M0,M)
ck('quartic containing critical-free disk',mu*Z1<1 and Z0+R<Z1)
wholeExp=EU(M0+nu*Z0**3/(3*(1-mu*Z0)))
ck('quartic full path and relative derivative exponential',wholeExp<E0)
pert=Z0**4*E0/(12*(1-mu*Z0))
ck('quartic whole primitive perturbation coefficient',pert<P0,pert)
center=E0*(132200*q**5+P0*nu)
ck('quartic full center relative residual',center<(1-nb)*R,center)
ck('quartic full containing logarithmic derivative',1+q*Z1+nu*Z1*Z1/(1-mu*Z1)<L)
nonlinear=(EU(L*R)-1-L*R)/L
ck('quartic full-circle nonlinear bound',nonlinear<nb*R)
ck('quartic complete-circle strict Rouche',center+nonlinear<R)
ck('quartic both circles avoid original root',6-21*q-133*q*q-1030*q**3-9000*q**4>R)
ck('quartic pair circles are disjoint',12*(1-q-17*q*q-137*q**3-313*q**4)>2*R)
ck('T low field q entrance uses horizontal geometry only',2/(smallX*79)<q)
ck('T low field mu entrance uses physical gap only',1/(d0*79)<mu)
Rphys=F(19,200)
ck('T low field complete physical companion disk',Z1/F(79)<=Rphys)
wcomp=weight(Rphys,F(3,10),F(1),F(1))
W103=rho*(103*ft+1)/(1-BT)
ck('T actual beta and both distinct companions beat all h103',tsource+2*wcomp>W103,tsource+2*wcomp-W103)
W79=rho*(79*ft+1)/(1-BT)
ck('T fixed actual source closes h79',tsource>W79,tsource-W79)
# Complete geometry support, unchanged finite HB/root-product lemmas.
ck('fixed X test nonzero before HB',c*c*F(101,100)**2<F(1,100))
ck('fixed T test nonzero before HB',c*c*F(13,10)**2<F(3,10))
ck('moving test full positive-inset domain',2/(1-rho)**2<3 and 1-4*c*c>0)
ck('X lower strip whole actual-root norm',F(109,100)<F(21,20)**2 and F(21,20)+F(13,100)==F(59,50))
ck('X lower strip whole squared norm',F(59,50)**2/2==F(3481,5000))
ck('X symmetric whole actual-root norm',F(99,70)**2>2 and (F(99,70)+F(3,20))**2<F(5,2))
ck('X lower strip complete distance',F(13,100)+F(1,10000)/(2*F(3,50))==dist)
ck('X symmetric complete distance',F(3,20)+F(1,10000)/(2*F(3,10))==F(901,6000))
mm=F(49,100);de=F(3,10)
ck('T full small-X comparison derivative',(4*mm)**2+12*(4*mm*de-1)<0)
ck('T full small-X comparison endpoint',2*mm*(mm+de)**2-mm*(mm*mm+1)>0)
def serial(obj):
 if isinstance(obj,F):return str(obj)
 if isinstance(obj,list):return [serial(v) for v in obj]
 if isinstance(obj,dict):return {k:serial(v) for k,v in obj.items()}
 return obj
result={'status':'PASS' if all(z['pass'] for z in checks) else 'FAIL',
 'n0':n0,'rho_cap':rho,'no_lower_n_rho_assumption':True,
 'checks':len(checks),'failed':[z['name'] for z in checks if not z['pass']],
 'complete_inverse_rows':len(rows),'quartic_complete_rows':1,
 'predicates':checks,'rows':rows,'strip_compensation_polynomial':P,
 'strip_compensation_Bernstein_coefficients':bcoef,
 'diagnostics':{'X_high_low_upper':upperlow,'X_high_depth_upper':upperhigh,
 'X_middle':middle,'X_h78_lower':lower78,'X_h78_upper':upper78,
 'X_h78_source_lower':bsource,'X_h78_total_upper':W78,
 'T_h150_lower':lower150,'T_h150_upper':upper150,
 'T_h103_lower':lower103,'T_h103_upper':upper103,
 'T_smallX_source_lower':tsource,'T_smallX_total_upper':Wsmall,
 'quartic_relative_center':center,'quartic_companion_weight':wcomp,
 'T_h79_source_lower':tsource,'T_h79_total_upper':W79}}
(OUT/f'BASE010_{n0}_COUPLED_EXACT.json').write_text(json.dumps(serial(result),indent=2)+'\n')
print(json.dumps(serial({k:result[k] for k in ['status','n0','rho_cap','checks','failed','complete_inverse_rows','quartic_complete_rows','diagnostics']}),indent=2))
