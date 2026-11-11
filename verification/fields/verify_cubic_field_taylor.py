"""Root's independent full-cube local-field verifier.

No imports from either discovery or either preceding verifier.  Every input
float is interpreted as its exact binary rational.  Gaussian primitives and
positive perturbation integrals use independently implemented Taylor recurrences
with full complex Cauchy tails on |t|=2, rather than Acb analytic quadrature.
Only outward-enclosed strict inequalities accept a complete closed leaf.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import argparse,json,time,hashlib
from flint import arb,acb,ctx
TARGET=F(3);ENTRY=30;PREC=224;ORDER=256;TAIL_DEGREE=8

def ar(x):
 if isinstance(x,F):return arb(x.numerator)/x.denominator
 if isinstance(x,float):
  a,b=x.as_integer_ratio();return arb(a)/b
 return arb(x)

def closed_box(rec):
 u=[[F(0),F(1)] for _ in range(4)]
 for c in rec['path']:
  axis,side=divmod(int(c),2);mid=(u[axis][0]+u[axis][1])/2
  u[axis][1-side]=mid
 endpoints=[TARGET/270,TARGET/ENTRY,F(13*(ENTRY+1),ENTRY**2)]
 lo,hi=endpoints[rec['branch']:rec['branch']+2]
 q=[ar(lo+(hi-lo)*v) for v in u[0]]
 ph=[ar(v)*arb.pi()/2 for v in u[1]]
 aa=[ar(2*v/TARGET) for v in u[2]]
 bb=[ar(2*(2*v-1)/TARGET) for v in u[3]]
 return u,q,ph,aa,bb

def coverage(leaves):
 out=[]
 for b in [0,1]:
  raw=[r['path'] for r in leaves if r['branch']==b]
  assert raw and len(set(raw))==len(raw)
  current=set(raw);depth=max(map(len,raw))
  assert all(all(c in '01234567' for c in p) for p in raw)
  for k in range(depth,0,-1):
   groups={}
   for p in current:
    if len(p)==k:groups.setdefault(p[:-1],[]).append(int(p[-1]))
   for parent,digits in groups.items():
    digits.sort();assert len(digits)==2 and digits[0]%2==0 and digits[1]==digits[0]+1
    assert parent not in current
    for d in digits:current.remove(parent+str(d))
    current.add(parent)
  assert current=={''}
  vol=sum((F(1,2**len(p)) for p in raw),F(0));assert vol==1
  out.append(dict(branch=b,leaves=len(raw),maxdepth=depth,exact_volume=str(vol)))
 assert sum(r['leaves'] for r in out)==len(leaves)
 return out

def exp_polynomial_integral(coeff,complex_output=False):
 # Entire exp(sum_j c_j t^j) integrated over[0,1]. Different from Acb quadrature.
 typ=acb if complex_output else arb
 coeff=[typ(v) for v in coeff]
 seq=[typ(1)];total=typ(1)
 for k in range(1,ORDER+1):
  z=typ(0)
  for j in range(1,min(k,len(coeff))+1):z+=j*coeff[j-1]*seq[k-j]
  z/=k;seq.append(z);total+=z/(k+1)
 # Two complete Cauchy circles; the better proved enclosure is selected.
 tails=[]
 for rr in (F(3,2),F(2)):
  rad=ar(rr);M=sum((abs(z).upper()*rad**(j+1) for j,z in enumerate(coeff)),arb(0)).exp()
  tails.append((M/(rad**(ORDER+1)*(ORDER+2)*(1-1/rad))).upper())
 tail=min(tails,key=lambda x:float(x))
 if complex_output:return total+acb(arb(0,tail),arb(0,tail)),tail
 return total+arb(0,tail),tail

def third_moment_ball(qlo,qhi,plo,phi,alo,ahi,blo,bhi,mucap,kappa,strip=None):
 # Independent weighted-covariance enclosure of sum(v_j^3), before phase.
 cmin=phi.cos().lower();cmax=plo.cos().upper()
 if not cmin>0:raise ValueError('whole covariance denominator not positive')
 amin=(qlo*alo/2).lower();amax=(qhi*ahi/2).upper()
 vals=[-q*b/2 for q in (qlo,qhi) for b in (blo,bhi)]
 cxylo=min(v.lower() for v in vals);cxyhi=max(v.upper() for v in vals)
 if cxylo<=0<=cxyhi:cxy2=arb(0)
 else:cxy2=min((cxylo*cxylo).lower(),(cxyhi*cxyhi).lower())
 real_lo=(4*amin*amin/cmax-3*mucap*mucap*cmax).lower()
 real_hi=min((mucap*amax).upper(),(mucap*mucap*cmax-4*cxy2/cmax).upper())
 if not real_hi>=real_lo:raise ValueError('empty covariance rectangle requires separate prune')
 centers=[3*A*C/c for A in (amin,amax) for C in (cxylo,cxyhi) for c in (cmin,cmax)]
 variance=(mucap*mucap*cmax-(amin*amin+cxy2)/cmax).upper()
 variance=max(arb(0),variance)
 imagerr=(3*variance/2+mucap*(qhi-amin)).upper()
 imag_lo=(min(z.lower() for z in centers)-imagerr).lower()
 imag_hi=(max(z.upper() for z in centers)+imagerr).upper()
 if strip is not None:
  i,count=strip['index'],strip['count']
  assert type(i) is int and type(count) is int and count in (2,4,8,16) and 0<=i<count
  left,span=imag_lo,imag_hi-imag_lo
  imag_lo=(left+span*ar(F(i,count))).lower()
  imag_hi=(left+span*ar(F(i+1,count))).upper()
 zmid=acb((real_lo+real_hi)/2,(imag_lo+imag_hi)/2)
 radius=(((real_hi-real_lo)/2)**2+((imag_hi-imag_lo)/2)**2).sqrt().upper()
 pm=(plo+phi)/2;cmid=zmid*acb(0,-3*pm).exp()
 angular=2*abs(zmid).upper()*(3*(phi-plo)/4).sin()
 return (radius+angular+abs(cmid-kappa).upper()).upper()

def verify_model(rec):
 try:
  unit,(ql,qh),(pl,ph),(al,ah),(bl,bh)=closed_box(rec)
  mu=(qh/ar(TARGET) if rec['branch']==0 else arb(1)/ENTRY).upper()
  bm2=arb(0) if unit[3][0]<=F(1,2)<=unit[3][1] else (bh*bh if unit[3][1]<F(1,2) else bl*bl)
  radius=2*pl.cos()/ar(TARGET) if rec['branch']==0 else 2*pl.cos()/(ENTRY*ql)
  # a ranges below 1, so min (1-a)^2 occurs at a=ah.
  assert ah<1
  tests=[('moment_lens',al*al+bm2-radius*radius),('unit_moment_disk',(1-ah)**2+bm2-1)]
  pi=arb.pi();de=3*pi*pi-8
  envelope=9*pi*pi*mu*mu/de+(2+24*pi/de)*mu*pl.cos()+(40-3*pi*pi)/(2*de)*pl.cos()**2-pl.sin()**2/2
  tests.append(('matrix',ql-envelope))
  for label,margin in tests:
   if margin>0:return True,dict(method=label,margin_lower=float(margin.lower()))
  if 'point' not in rec:return False,dict(reason='no independently certified exclusion')
  w=acb(ar(rec['point'][0]),ar(rec['point'][1]));th=acb(ar(rec['theta'][0]),ar(rec['theta'][1]))
  qm=(ql+qh)/2;pm=(pl+ph)/2;am=(al+ah)/2;bm=(bl+bh)/2
  core=acb(1-am,bm);dc=(((ah-al)/2)**2+((bh-bl)/2)**2).sqrt()
  center=qm*core*acb(0,-2*pm).exp()
  # Expand q*C*exp(-2i phi) about qm,Cm,pm; every cross-term retained.
  hr=((qh-ql)/2*(abs(core).upper()+dc)+qm*dc+2*qm*abs(core).upper()*((ph-pl)/2).sin()+abs(center-th).upper()).upper()
  W=abs(w).upper();disk=(1-mu*W).lower()
  if not disk>0:return False,dict(reason='full critical logarithm disk')
  hb=(hr*W*W/2).upper()
  cubic='cubic' in rec
  kap=acb(ar(rec['cubic'][0]),ar(rec['cubic'][1])) if cubic else acb(0)
  h3=third_moment_ball(ql,qh,pl,ph,al,ah,bl,bh,mu,kap,rec.get('cubic_slice')) if cubic else arb(0)
  hc=(h3*W**3/3).upper()
  B=th*w*w/2;C=kap*w*w*w/3
  integ,tail=exp_polynomial_integral([w,B,C],True)
  J0=w*integ;g0=(w+B+C).exp();JA=abs(J0);GA=abs(g0)
  cp=[w.real.upper(),B.real.upper()+hb,C.real.upper()+hc]+[arb(0)]*(TAIL_DEGREE-2)
  tailone=arb(0)
  for k in range(4 if cubic else 3,TAIL_DEGREE+1):
   coeff=(qh*mu**(k-2)*W**k/k).upper();cp[k-1]+=coeff;tailone+=coeff
  k=TAIL_DEGREE+1
  coeff=(qh*mu**(k-2)*W**k/(k*disk)).upper();cp[k-1]+=coeff;tailone+=coeff
  Iplus,tplus=exp_polynomial_integral(cp)
  Ibase,tbase=exp_polynomial_integral([w.real.upper(),B.real.upper(),C.real.upper()])
  EJ=(W*(Iplus-Ibase)).upper()
  if not EJ>=0:return False,dict(reason='primitive error enclosure sign')
  derivative_error=(GA.upper()*((hb+hc+tailone).exp()-1)).upper()
  if not JA.lower()>EJ:return False,dict(reason='actual primitive disk may contain zero',primitive_error=float(EJ),modelabs=float(JA.lower()))
  angle=((ph-pl)/4).sin();sig=acb(0,pm).exp()
  P=((w*sig.conjugate()).real-2*W*angle).lower()
  if not P>0:return False,dict(reason='physical right half-plane')
  numerator=(sig*g0*J0.conjugate()).real+2*GA.upper()*JA.upper()*angle+GA.upper()*EJ
  upper=numerator/(JA*JA-EJ*EJ)+derivative_error/(JA-EJ)-P/(W*W)
  return bool(upper<0),dict(method='cubic_Taylor_Cauchy' if cubic else 'Gaussian_Taylor_Cauchy',upper=float(upper.upper()),projection=float(P),domain=float(disk),model_tail=float(tail),positive_integral_tail=float((tplus+tbase).upper()),primitive_error=float(EJ))
 except Exception as e:return False,dict(reason=repr(e))

def verify(rec):
 if 'cubic_slice' in rec:return False,dict(reason='standalone conditional slice is forbidden')
 if 'family' not in rec:return verify_model(rec)
 try:
  fam=rec['family'];count=len(fam)
  assert type(fam) is list and count in (2,4,8,16)
  details=[]
  for i,member in enumerate(fam):
   assert set(member)=={'theta','cubic','point','cubic_slice'}
   strip=member['cubic_slice'];assert set(strip)=={'index','count'}
   assert type(strip['index']) is int and type(strip['count']) is int
   assert strip['index']==i and strip['count']==count
   r={'branch':rec['branch'],'path':rec['path'],'status':'point',**member}
   ok,detail=verify_model(r)
   if not ok:return False,dict(reason='complete family member failed',index=i,count=count,detail=detail)
   details.append(detail)
  points=[d for d in details if 'upper' in d]
  if not points:return True,details[0]
  assert len(points)==count
  return True,dict(method='complete_cubic_family_Taylor',family_members=count,
   upper=max(d['upper'] for d in points),domain=min(d['domain'] for d in points),
   projection=min(d['projection'] for d in points),
   positive_integral_tail=max(d['positive_integral_tail'] for d in points))
 except Exception as e:return False,dict(reason=repr(e))

def init(target,entry,prec,order,tail_degree):
 global TARGET,ENTRY,PREC,ORDER,TAIL_DEGREE
 TARGET=F(target);ENTRY=entry;PREC=prec;ORDER=order;TAIL_DEGREE=tail_degree;ctx.prec=prec

def main():
 ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--prec',type=int,default=224);ap.add_argument('--order',type=int,default=256);a=ap.parse_args()
 raw=Path(a.source).read_bytes();data=json.loads(raw);assert not data['unresolved']
 target=F(str(data['target']));entry=int(data['entrance']);tail_degree=int(data.get('tail_degree',8));assert tail_degree>=4;assert F(9,4)<=target<=3 and 0<entry<=270
 part=coverage(data['leaves']);bad=[];counts=Counter();maxupper=-float('inf');mint=float('inf');minp=float('inf');maxtail=0.;started=time.time()
 with ProcessPoolExecutor(max_workers=a.workers,initializer=init,initargs=(str(target),entry,a.prec,a.order,tail_degree)) as pool:
  for j,(ok,detail) in enumerate(pool.map(verify,data['leaves'],chunksize=24)):
   counts[detail.get('method','error')]+=1
   if not ok:bad.append(dict(index=j,record=data['leaves'][j],detail=detail))
   if 'upper' in detail:maxupper=max(maxupper,detail['upper']);mint=min(mint,detail['domain']);minp=min(minp,detail['projection']);maxtail=max(maxtail,detail['positive_integral_tail'])
   if (j+1)%1000==0:print(json.dumps(dict(checked=j+1,failed=len(bad),seconds=time.time()-started)),flush=True)
 result=dict(status='PASS' if not bad else 'FAIL',source_sha256=hashlib.sha256(raw).hexdigest(),target=str(target),entrance=entry,precision_bits=a.prec,Taylor_order=a.order,method='independent full tree/boxes, full P3 strip families, degree9 Taylor integrals with complete Cauchy tails; no author verifier imports or quadrature',partition=part,counts=dict(counts),failed_count=len(bad),worst_quotient_upper=maxupper,min_log_disk=mint,min_projection=minp,max_integral_tail=maxtail,seconds=time.time()-started,failures=bad)
 Path(a.output).write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='failures'},indent=2))
if __name__=='__main__':main()
