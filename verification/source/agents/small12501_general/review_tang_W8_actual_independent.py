#!/usr/bin/env python3
"""Independent W8 actual companions: exact Laurent cancellation and full tails.
No author function or numeric receipt is imported. QQ[s] input coefficients
are revalidated by a different exact relative-primitive construction.
Numeric parameter exp uses weighted Euler; integration uses Ik recurrence;
full paths use Taylor disks on 128 closed intervals, not Bernstein controls.
"""
from fractions import Fraction as F
from math import factorial,comb,ceil
from pathlib import Path
from flint import fmpq_poly as QP,fmpq as QQ,acb,arb,acb_poly,ctx
import json,time,hashlib,argparse
ctx.prec=448
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
ZERO=(0,0,0)
def wt(k):return k[0]+2*k[1]+3*k[2]
def plus(k,l):return tuple(a+b for a,b in zip(k,l))
def minus(k,l):return tuple(a-b for a,b in zip(k,l))
def mons(W):return sorted(((a,b,c) for a in range(W+1) for b in range(W//2+1) for c in range(W//3+1) if a+2*b+3*c<=W),key=lambda k:(wt(k),k))
def A(x):
 if isinstance(x,F):return arb(x.numerator)/x.denominator
 return arb(x)
def C(x):return acb(A(F(x)))
PRED=[]
def ck(name,p):
 ok=bool(p);PRED.append(dict(name=name,passed=ok))
 if not ok:raise AssertionError(name)
def pmul(p,q,W=None):
 out={}
 for k,v in p.items():
  for l,z in q.items():
   m=plus(k,l)
   if W is None or wt(m)<=W:out[m]=out.get(m,v*0)+v*z
 return out
def tadd(p,q):
 r=p.copy()+[QP()]*(max(len(p),len(q))-len(p))
 for j,v in enumerate(q):r[j]+=v
 while len(r)>1 and not r[-1]:r.pop()
 return r
def tmul(p,q):
 r=[QP() for _ in range(len(p)+len(q)-1)]
 for i,v in enumerate(p):
  for j,z in enumerate(q):r[i+j]+=v*z
 while len(r)>1 and not r[-1]:r.pop()
 return r
def exact_check(Z):
 W=8;zp=[{ZERO:QP([1])},Z]
 for j in range(2,5):zp.append(pmul(zp[-1],Z,W))
 AA={}
 for k,v in Z.items():
  if wt(k)>0:AA[k]=[QP(),v]
 for j,shift,den in [(2,(1,0,0),2),(3,(0,1,0),3),(4,(0,0,1),4)]:
  for k,v in zp[j].items():
   kk=plus(k,shift)
   if wt(kk)<=W:AA[kk]=tadd(AA.get(kk,[QP()]),[QP()]*j+[v*QQ(1,den)])
 for k,v in AA.items():v[0]-=sum(v,QP())
 E={ZERO:[QP([1])]}
 for n in mons(W)[1:]:
  value=[QP()]
  for k,v in AA.items():
   if all(k[j]<=n[j] for j in range(3)):
    term=tmul(v,E[minus(n,k)]);term=[x*QQ(wt(k),wt(n)) for x in term];value=tadd(value,term)
  E[n]=value
 JC={}
 for n in mons(W):
  value=[QP()]
  for k,z in Z.items():
   if all(k[j]<=n[j] for j in range(3)):value=tadd(value,[x*z for x in E[minus(n,k)]])
  JC[n]=value
 degree=max(len(p)-1 for p in JC.values());IP=[QP()]
 for j in range(1,degree+1):
  poly=[QQ(0)]*(degree+1)
  for l in range(j):poly[degree-l-1]=QQ((-1)**l*factorial(j),factorial(j-l))
  IP.append(QP(poly))
 for n,p in JC.items():
  integrated=sum((p[j]*IP[j] for j in range(len(p))),QP())
  ck('QQ[s] full relative coefficient vanishes '+str(n),not integrated)
 return dict(weight=W,monomials=len(JC),common_Laurent_shift=degree,method='weighted Euler in QQ[s][t], multiply full antiderivative by a common s power')
def upper(x,name,scale=10**14):
 f=F(ceil(float(x.upper())*scale)+3,scale);ck(name,A(f)>x);return f
def path_bound(p,name,parts=128):
 co=[p[j] for j in range(5)];h=F(1,2*parts);hh=[h**j for j in range(5)];vals=[]
 for i in range(parts):
  mid=F(2*i+1,2*parts)
  jets=[sum((co[k]*C(comb(k,j)*mid**(k-j)) for k in range(j,5)),acb(0)) for j in range(5)]
  bound=sum((abs(jets[j])*A(hh[j]) for j in range(5)),arb(0));vals.append(bound)
 b=F(ceil(max(float(x.upper()) for x in vals)*10**14)+3,10**14)
 ck(name+' entire 128 Taylor disks',all(A(b)>x for x in vals))
 return b
def evaluate(label,q,mu,Ar,Br,beta,phase,Z1,L,exactZ,exact_info,source_hash):
 global PRED
 start=time.monotonic();before=len(PRED);s=acb(0,2*arb.pi());nu=q*mu;lam=q*mu*mu
 Z={}
 for k,p in exactZ.items():
  v=acb(0)
  for a in reversed(p.coeffs()):v=v*s+C(F(str(a)))
  Z[k]=v*C(q**k[0]*nu**k[1]*lam**k[2])
 zp=[{ZERO:acb(1)},Z]
 for j in range(2,5):zp.append(pmul(zp[-1],Z))
 AA={}
 for k,z in Z.items():
  if wt(k)>0:AA[k]=acb_poly([0,z])
 for j,shift,fac in [(2,(1,0,0),q/2),(3,(0,1,0),nu/3),(4,(0,0,1),lam/4)]:
  for k,z in zp[j].items():
   kk=plus(k,shift);AA[kk]=AA.get(kk,acb_poly())+acb_poly([0]*j+[z*C(fac)])
 for k,p in AA.items():AA[k]=p-acb_poly([sum(p.coeffs(),acb(0))])
 ck(label+'/full exponent weight35',max(map(wt,AA))==35)
 ck(label+'/whole exponent degree4',max(len(p)-1 for p in AA.values())==4)
 ck(label+'/no constant parameter exponent',ZERO not in AA)
 W=28;MM=mons(W);E={ZERO:acb_poly([1])}
 entries=list(AA.items())
 for n in MM[1:]:
  p=acb_poly()
  for k,v in entries:
   if all(k[j]<=n[j] for j in range(3)):p+=v*E[minus(n,k)]*wt(k)
  E[n]=p*C(F(1,wt(n)))
 JC={}
 for n in MM:
  p=acb_poly()
  for k,z in Z.items():
   if all(k[j]<=n[j] for j in range(3)):p+=E[minus(n,k)]*z
  JC[n]=p
 degree=max(len(p)-1 for p in JC.values());IP=[acb(0)]
 for j in range(1,degree+1):IP.append((1-j*IP[-1])/s)
 finite=arb(0);coefficients=[]
 for n,p in JC.items():
  value=sum((p[j]*IP[j] for j in range(len(p))),acb(0))
  if wt(n)<=8:ck(label+'/numeric consistency with independent exact identity '+str(n),value.contains(0))
  else:
   ck(label+'/finite coefficient total q degree >=2 '+str(n),sum(n)>=2);finite+=abs(value)
  coefficients.append(dict(monomial=n,weight=wt(n),value=str(value)))
 print(json.dumps(dict(event='independent finite residual',label=label,finite=str(finite),seconds=time.monotonic()-start)),flush=True)
 am=[F(0)]*36;paths=[]
 for k,p in sorted(AA.items()):
  bound=path_bound(p,label+'/'+str(k));am[wt(k)]+=bound;paths.append(dict(monomial=k,bound=str(bound)))
 zm=[upper(sum((abs(z) for k,z in Z.items() if wt(k)==j),arb(0)),label+'/center weight '+str(j)) for j in range(9)]
 # A second distinct construction: scalar weighted differential recurrence.
 ep=[F(1)]
 for j in range(1,W+1):ep.append(sum((F(k)*am[k]*ep[j-k] for k in range(1,j+1)),F(0))/j)
 major=[sum((zm[j]*ep[k-j] for j in range(min(8,k)+1)),F(0)) for k in range(W+1)]
 Zb=sum(zm);Mrel=sum(am);tail=A(Zb)*A(Mrel).exp()-A(sum(major))
 ck(label+'/full positive infinite model tail',tail>0)
 ck(label+'/entire infinite tail q degree >=10',(W+1+2)//3>=10)
 ck(label+'/whole center critical-free',mu*Zb<1)
 r5=q*mu**3*Zb**5/(5*(1-mu*Zb))
 model=(finite+tail)*A(r5).exp()
 actual=A(q*mu**3*Zb**6/(30*(1-mu*Zb)))*A(Mrel+2*r5).exp()
 R=Ar*q*q+Br*q*mu**3;nonlinear=(A(L*R).exp()-1-A(L*R))/A(L)
 ck(label+'/full actual containing disk',Zb+R<Z1)
 ck(label+'/full critical-free containing disk',mu*Z1<1)
 ck(label+'/complete logarithmic derivative',1+q*Z1+q*mu*Z1**2/(1-mu*Z1)<L)
 ck(label+'/whole variable model q2',model/A(q*q)<A((1-beta)*Ar))
 ck(label+'/whole actual variable product qmu3',actual/A(q*mu**3)<A((1-beta)*Br))
 ck(label+'/full-circle nonlinear fraction',nonlinear/A(R)<A(beta))
 ck(label+'/complete actual Rouche domination',A(R)-model-actual-nonlinear>0)
 displacement=sum(zm[1:]);phase_margin=2*arb.pi()*A(phase*q)-A(displacement+R)
 ck(label+'/full homogeneous phase',phase_margin>0)
 ck(label+'/both circles exclude source',2*arb.pi()>A(displacement+R))
 ck(label+'/two full circles disjoint',4*arb.pi()>A(2*(displacement+R)))
 result=dict(status='PASS',label=label,q_max=str(q),mu_max=str(mu),center_weight=8,tail_weight=28,precision_bits=ctx.prec,independent_exact_center=exact_info,center_sha256=source_hash,parameter_recursion='weighted Euler',integral_recursion='I0=0; Ik=(1-k I(k-1))/s',path_enclosure='128 complete Taylor disks per full quartic path polynomial',path_polynomials=len(AA),path_intervals=len(AA)*128,center_upper=str(Zb),Mrelative=str(Mrel),R5=str(r5),finite_residual=str(finite),infinite_model_tail=str(tail),model_relative=str(model),actual_product_relative=str(actual),radius_q2=str(Ar),radius_qmu3=str(Br),radius=str(R),nonlinear_fraction=str(beta),containing=str(Z1),logder=str(L),whole_nonlinear=str(nonlinear),whole_circle_margin=str(A(R)-model-actual-nonlinear),phase=str(phase),phase_margin=str(phase_margin),complete_center_weights=list(map(str,zm)),weighted_path_majorants=list(map(str,am)),path_coefficients=paths,integrated_coefficients=coefficients,predicates=PRED[before:],checks=len(PRED)-before,seconds=time.monotonic()-start)
 (HERE/('INDEPENDENT_'+label+'_ACTUAL_RECEIPT.json')).write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['complete_center_weights','weighted_path_majorants','path_coefficients','integrated_coefficients','predicates']},indent=2),flush=True)
 return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--row',choices=['both','q06','q0546'],default='both');args=ap.parse_args()
 source=ROOT/'agents/tang_weighted_moment_center_W8.json';raw=source.read_bytes();d=json.loads(raw)
 Z={tuple(r['monomial']):QP([QQ(v) for v in r['coefficients']]) for r in d['coefficients'] if wt(tuple(r['monomial']))<=8}
 ck('input contains exactly41 full center coefficients',len(Z)==41)
 info=exact_check(Z);print(json.dumps(dict(event='independent exact relative cancellation',**info)),flush=True)
 exact_result=dict(status='PASS',exact_center=info,predicates=PRED.copy(),center_sha256=hashlib.sha256(raw).hexdigest())
 (HERE/'INDEPENDENT_W8_EXACT_RELATIVE_CANCELLATIONS.json').write_text(json.dumps(exact_result,indent=2)+'\n')
 if args.row in ['both','q06']:evaluate('W8_Q06',F(3,50),F(1,50),F(64),F(95000),F(1,3),F(8),F(19,2),F(7,4),Z,info,hashlib.sha256(raw).hexdigest())
 if args.row in ['both','q0546']:evaluate('W8_Q0546',F(273,5000),F(91,5000),F(22),F(45000),F(1,5),F(33,5),F(9),F(17,10),Z,info,hashlib.sha256(raw).hexdigest())
if __name__=='__main__':main()
