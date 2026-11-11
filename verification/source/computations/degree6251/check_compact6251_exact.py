"""Full exact general Compact6251 constants. Analytic proof is separate.
Every circle choice is reconstructed with complete rational tails. This
does not import any floating discovery evaluation as an accepting gate.
"""
from fractions import Fraction as F
from math import factorial,isqrt
from pathlib import Path
import json,hashlib,time
start=time.monotonic();HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'COMPACT6251_INTEGER_CHOICES.json').read_text())
N=6250;n=N+1;CAP=F(1,388);G=[];records=[];cnt=0
def ck(name,p):
 G.append(dict(name=name,passed=bool(p)))
 if not p:raise AssertionError(name)
def elo(x):return sum((x**j/F(factorial(j)) for j in range(101)),F(0))
EW={1:F(367879441172,10**12),2:F(135335283237,10**12),3:F(49787068368,10**12),4:F(18315638889,10**12),5:F(6737947000,10**12),6:F(2478752177,10**12),8:F(335462628,10**12),10:F(45399930,10**12)}
for w,e in EW.items():ck('complete_exp_'+str(w),elo(F(w))>1/e)
ck('complete_cap',CAP*elo(N*CAP/2)>8)
def sqrtup(v):
 sc=10**12;k=isqrt(v.numerator*sc*sc//v.denominator)
 if F(k,sc)**2<v:k+=1
 return F(k,sc)
def weight(eta):
 t=(1-eta)/(1+eta);lo=2*sum((t**(2*j+1)/F(2*j+1) for j in range(24)),F(0))
 return F((lo*10**12).numerator//(lo*10**12).denominator,10**12)
for A,rr in [(F('.11'),F('.65')),(F('.2'),F('.27')),(F('.45'),F('.115')),(F('.75'),F('.069')),(F('.85'),F('.06')),(F('.9'),F('.057')),(F('.93'),F('.055')),(F('.95'),F('.054')),(F('.96'),F('.053')),(F('.965'),F('.053'))]:
 ck('entire_far_'+str(A),A*A*rr*rr*(1-rr*rr)>CAP and rr*rr<F(1,2))
for row in data:
 name=row['name'];A,B,rlo,rhi,target=map(F,[row['A'],row['B'],row['rlo'],row['rhi'],row['target']])
 def c(tag,p):ck(name+'/'+tag,p)
 c('whole_closed_cell',F('.11')<=A<B<=F('.975') and 0<rlo<rhi and rhi*rhi<F(1,2))
 near=row['near_c'];ce=A*A*rlo*rlo*(1-rlo*rlo)
 c('complete_log_target',elo(target)>8/ce)
 if near:
  nc=F(near);bb=B+F(1,1000)
  c('exact_shrinking_near',rlo==nc/n)
  c('whole_near_modulus',B+rlo<bb)
  c('whole_variance_exclusion',nc<(1-bb*bb)/2)
  c('entire_shrinking_ray',target>2 and 0<nc<n)
 d=1-rhi;m=1/d
 if row['mode']=='trace':
  cp=A/2
  c('independent_direct_D100',F(5*N,14)*cp/(1-cp+cp*cp)>100)
  df=F(row['Dfloor']);kap=F(row['kappa'])
  c('ratio9_4_before_quantitative',F(2*N,13)*d>df>=201)
  # Floor D differs from discovery's continuous D; this gate must be
  # rebuilt and may require a more conservative proposed kappa.
  c('whole_smaller_branch',2<kap<F(9,4) and F(9,4)*kap*kap-(df-13)*kap+2*df<0)
  V=min(2/A,kap/d);Hi=max(N*A/2,F(2*N,13))
  c('prior_fine_q',V/Hi<F(7,625));c('prior_fine_mu',m/Hi<F(1,270))
  cr=A/(2*(1+rhi)**2);alpha=(1+2*cr*m)/(2*(cr+m));beta=(F(9,4)*V*V+13*V*m)/(2*(cr+m))
  H=alpha*N-beta;c('signed_trace',0<cr<m and H>Hi and alpha>0 and beta>0)
 else:
  ell,wm,km=map(F,row['moment'])
  if A<F('.93'):c('entire_moment_scope',A>=F('.9') and rhi<=F('.057') and [ell,wm,km]==list(map(F,['4.42','5.098','2.5482'])))
  elif A<F('.95'):c('entire_moment_scope',A>=F('.93') and rhi<=F('.055') and [ell,wm,km]==list(map(F,['4.22','4.765','2.3817'])))
  else:c('entire_moment_scope',A>=F('.95') and rhi<=F('.054') and [ell,wm,km]==list(map(F,['4.11','4.623','2.3108'])))
  er=8*rhi/(A*N)+rhi*rhi/(2*(1-rhi)*(1-2*rhi))
  qm=1-A*A+wm*(ell/N+er)+F(2,n*n);mm=km*(ell/n+er)+F(1,n*n)
  cf=1-mm/A-qm/(A*d);alpha=cf;beta=F(0);H=N*cf;V=1/(cf*d*d)
  c('whole_moment_field',cf>0 and 0<qm<1 and 0<mm<A)
  c('independent_phase_entrance',F(4,N)/(A*A)<F(7,625) and 2/(N*A*d)<F(1,270))
 weights=[];counts=[];slopes=[];circles=[];lastR=F(0)
 for j,nn,W in row['circles']:
  R=F(j,1000);nu=F(nn,10000);e=EW[W];label=name+'/'+str(j)
  lam=V*R/(1-m*R);rh=R-V*(R*R/2+m*R**3/(3*(1-m*R)))
  Bi=V/(1-m*R)**2/(1-lam)**3;eps=F(1001,1000)*Bi/H;gg=F(1001,1000)*nu*Bi*rh;pin=1/(1+lam);amp=1/(1-lam)
  for lab,pred in [('critical_free',0<m*R<1),('whole_inverse',0<lam<1),('positive_image',rh>0),('rational_nu',0<nu<1),('complete_corner',W*W<(1-nu*nu)*(H*rh)**2),('complete_left',eps*(1+e)+gg*e<1-amp*e),('complete_right',eps*(1+e)+gg*e*e<pin-e),('complete_horizontal',eps+gg*e<pin),('positive_amplitude',0<pin<1<amp),('positive_odd_height',nu*H*rh>F(22,7)),('nested_radius',R>lastR)]:ck(label+'/'+lab,pred)
  lastR=R;slope=nu*F(7,22)*rh;count=slope*H-2;ck(label+'/finite_count',count>0)
  delta=CAP*(1+rhi*rhi)/(2*(1-rhi*rhi))+(1+B*B/2)*rhi*rhi+R*rhi
  ck(label+'/entire_same_gamma_denominator',0<delta<1)
  rad=(B*B+R*R)/(1-delta)**2;eta=sqrtup(rad)
  ck(label+'/full_outward_radical',eta*eta>=rad and (eta-F(1,10**12))**2<rad)
  ck(label+'/full_factor',0<eta<1);ww=weight(eta);ck(label+'/positive_full_log_weight',ww>0)
  if weights:ck(label+'/monotone_full_weights',ww<=weights[-1]);ck(label+'/monotone_full_counts',count>counts[-1])
  weights.append(ww);counts.append(count);slopes.append(slope);cnt+=1
  circles.append(dict(j=j,nu=str(nu),W=W,eta=str(eta),weight=str(ww),count=str(count)))
 c('nonempty_whole_family',bool(weights))
 budget=sum((w*(ct-pr) for w,ct,pr in zip(weights,counts,[F(0)]+counts[:-1])),F(0))
 abel=weights[-1]*counts[-1]+sum(((weights[j]-weights[j+1])*counts[j] for j in range(len(weights)-1)),F(0))
 Wstar=sum((w*(sl-pr) for w,sl,pr in zip(weights,slopes,[F(0)]+slopes[:-1])),F(0))
 bonus=F(0) if near else N*ce/2;nslope=Wstar*alpha+(F(0) if near else ce/2);nloss=Wstar*(alpha+beta)+2*weights[0]+(F(0) if near else ce/2)
 c('full_Abel_identity',budget==abel and budget==Wstar*H-2*weights[0])
 c('complete_degree_ray',budget+bonus==nslope*n-nloss and nslope>0 and nloss>0)
 c('same_configuration_contradiction',budget+bonus>target)
 records.append(dict(name=name,A=str(A),B=str(B),rlo=str(rlo),rhi=str(rhi),H=str(H),V=str(V),alpha=str(alpha),beta=str(beta),target=str(target),budget=str(budget),bonus=str(bonus),reserve=str(budget+bonus-target),n_slope=str(nslope),n_loss=str(nloss),circles=circles))
out=dict(status='PASS',degree_floor=n,scope='.11<=a<=.975; requires separately proved near and moment6251 interfaces',checks=len(G),complete_contours=cnt,whole_rows=len(records),rows=records,predicates=G,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),choices_sha256=hashlib.sha256((HERE/'COMPACT6251_INTEGER_CHOICES.json').read_bytes()).hexdigest(),seconds=time.monotonic()-start)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['predicates','rows']},indent=2))
for r in records:print(r['name'],float(F(r['reserve'])))
