"""Independent entire-polydisk four-moment companion verification.

No root-author function imports.  Parameter exponentials use a coordinate
partial-derivative recurrence on Acb polynomials in the complete path
variable.  Entire path bounds use complex Bernstein convex hulls on 128
closed subintervals, not sampled paths or the author's interval-Horner grid.
Exponential integrals use complete finite antiderivative sums.
"""
from fractions import Fraction as F
from pathlib import Path
from flint import acb,arb,acb_poly,ctx
from math import comb,factorial,ceil
import json,time,hashlib

ctx.prec=352
PRED=[]
def A(x):
 x=F(x);return arb(x.numerator)/x.denominator
def C(x):return acb(A(x))
def ck(name,p):
 ok=bool(p);PRED.append({'name':name,'pass':ok})
 if not ok:raise AssertionError(name)
def wt(m):return m[0]+2*m[1]+3*m[2]
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def minus(a,b):return tuple(x-y for x,y in zip(a,b))
def monomials(wmax):
 return sorted(((a,b,c) for a in range(wmax+1) for b in range(wmax//2+1)
  for c in range(wmax//3+1) if a+2*b+3*c<=wmax),key=lambda x:(wt(x),x))
def pmul(a,b):
 out={}
 for k,v in a.items():
  for l,w in b.items():
   m=plus(k,l);out[m]=out.get(m,acb(0))+v*w
 return out
def cup(x,name,scale=10**12):
 f=F(ceil(float(x.upper())*scale)+2,scale)
 ck(name,A(f)>x);return f

def full_path_bernstein(p,end,name,parts=128):
 """Complex polynomial lies in convex hull of its Bernstein coefficients."""
 co=[p[j] for j in range(5)]
 den=F(1,parts);allabs=[];relabs=[]
 for cell in range(parts):
  lo=F(cell,parts)
  shifted=[sum((co[j]*C(comb(j,k)*lo**(j-k)*den**k) for j in range(k,5)),acb(0))
           for k in range(5)]
  bern=[sum((shifted[j]*C(F(comb(k,j),comb(4,j))) for j in range(k+1)),acb(0))
        for k in range(5)]
  allabs.extend(abs(x) for x in bern)
  relabs.extend(abs(x-end) for x in bern)
 pa=max(float(x.upper()) for x in allabs)
 pr=max(float(x.upper()) for x in relabs)
 ba=F(ceil(pa*10**12)+2,10**12);br=F(ceil(pr*10**12)+2,10**12)
 ck(name+'/whole_Bernstein_control_polygon',all(A(ba)>x for x in allabs))
 ck(name+'/relative_Bernstein_control_polygon',all(A(br)>x for x in relabs))
 return ba,br

def main():
 start=time.time();s=acb(0,2*arb.pi());q=F(1,20);nu=F(1,1000);lam=F(1,50000)
 # Independently transcribed exact coefficient list, already symbolically
 # checked against the full formal primitive by this agent's separate file.
 c=s-s**2/2;d=3*s-2*s**2+s**3/2
 e=15*s-C(F(21,2))*s**2+C(F(23,6))*s**3-C(F(5,8))*s**4
 f=105*s-72*s**2+C(F(88,3))*s**3-C(F(22,3))*s**4+C(F(7,8))*s**5
 u=-s**3/3+s**2-2*s
 v=C(F(5,6))*s**4-C(F(13,3))*s**3+12*s**2-20*s
 w=-C(F(7,4))*s**5+C(F(157,12))*s**4-51*s**3+131*s**2-210*s
 r=s**5/3-C(F(13,6))*s**4+C(F(26,3))*s**3-22*s**2+40*s
 l=-s**4/4+s**3-3*s**2+6*s
 m=C(F(3,4))*s**5-C(F(19,4))*s**4+18*s**3-51*s**2+90*s
 Z={(0,0,0):s,(1,0,0):c*C(q),(2,0,0):d*C(q*q),
    (3,0,0):e*C(q**3),(4,0,0):f*C(q**4),(0,1,0):u*C(nu),
    (1,1,0):v*C(q*nu),(2,1,0):w*C(q*q*nu),(0,2,0):r*C(nu*nu),
    (0,0,1):l*C(lam),(1,0,1):m*C(q*lam)}
 Zpow=[{(0,0,0):acb(1)},Z]
 for j in range(2,5):Zpow.append(pmul(Zpow[-1],Z))
 AA={}
 for key,z in Z.items():
  if wt(key)>0:AA[key]=acb_poly([0,z])
 for j,shift,fac in [(2,(1,0,0),q/2),(3,(0,1,0),nu/3),(4,(0,0,1),lam/4)]:
  for key,z in Zpow[j].items():
   k=plus(key,shift);AA[k]=AA.get(k,acb_poly())+acb_poly([0]*j+[z*C(fac)])
 ck('full_A_max_weight19',max(map(wt,AA))==19)
 ck('full_A_has_no_constant_parameter',not (0,0,0) in AA)
 ck('whole_path_degree4',max(len(p)-1 for p in AA.values())==4)
 # For each parameter monomial n choose a nonzero coordinate j.  The
 # coefficient of d exp(A)/d x_j = (d A/d x_j) exp(A) gives this recursion.
 degree=14;E={(0,0,0):acb_poly([1])}
 for n in monomials(degree)[1:]:
  coord=next(i for i in range(3) if n[i])
  accum=acb_poly()
  for k,a in AA.items():
   if k[coord] and all(k[i]<=n[i] for i in range(3)):
    accum+=a*E[minus(n,k)]*k[coord]
  E[n]=accum*C(F(1,n[coord]))
 Jcoef={}
 for n in monomials(degree):
  p=acb_poly()
  for k,z in Z.items():
   if all(k[i]<=n[i] for i in range(3)):p+=E[minus(n,k)]*z
  Jcoef[n]=p
 maxdeg=max(len(p)-1 for p in Jcoef.values())
 # exp(s)=1.  The endpoint degree-m antiderivative cancels exactly with
 # the lower endpoint, leaving the full m-term expression below.
 Ip=[acb(0)]
 for k in range(1,maxdeg+1):
  val=acb(0)
  for j in range(k):
   val+=(-1)**j*factorial(k)//factorial(k-j)*s**(-j-1)
  Ip.append(val)
 finite=arb(0);integrated=[];weighted=[arb(0) for _ in range(degree+1)]
 for n,p in Jcoef.items():
  val=sum((p[j]*Ip[j] for j in range(len(p))),acb(0))
  integrated.append({'monomial':list(n),'weight':wt(n),'integral':str(val)})
  if wt(n)<=4:ck('exact_identity_consistency_'+str(n),val.contains(0))
  else:finite+=abs(val);weighted[wt(n)]+=abs(val)
 print(json.dumps({'event':'independent_finite_residual','value':str(finite),'seconds':time.time()-start}),flush=True)
 am=[F(0)]*20;mc=F(0);mr=F(0);paths=[]
 for key,p in sorted(AA.items()):
  end=sum(p.coeffs(),acb(0))
  b,br=full_path_bernstein(p,end,str(key))
  ce=cup(abs(end),'center_'+str(key))
  am[wt(key)]+=b;mc+=ce;mr+=br
  paths.append({'monomial':list(key),'whole_path':str(b),'relative_path':str(br),'center':str(ce)})
 zm=[cup(sum((abs(z) for k,z in Z.items() if wt(k)==j),arb(0)),'center_weight_'+str(j)) for j in range(5)]
 # Reconstruct the scalar majorant exponential by direct powers/factorials,
 # independently from the author's weighted differential recurrence.
 def scalar_mul(a,b):
  return [sum((a[j]*b[k-j] for j in range(k+1)),F(0)) for k in range(degree+1)]
 ap=am[:degree+1];power=[F(1)]+[F(0)]*degree;ec=power.copy()
 for j in range(1,degree+1):
  power=scalar_mul(power,ap)
  for k in range(degree+1):ec[k]+=power[k]/factorial(j)
 major=[sum((zm[j]*ec[k-j] for j in range(min(4,k)+1)),F(0)) for k in range(degree+1)]
 tail=A(sum(zm))*A(sum(am)).exp()-A(sum(major))
 ck('all_positive_majorant_tail',tail>0)
 Zbound=sum(zm);mub=F(1,50);k5=q*mub**3
 r5=k5*Zbound**5/(5*(1-mub*Zbound))
 model=(finite+tail)*A(mc+r5).exp()
 actual_tail=A(k5*Zbound**6/(30*(1-mub*Zbound)))*A(mr+2*r5).exp()
 R=F(9,20);Z1=F(17,2);L=F(31,20)
 ck('center_plus_circle_contained',Zbound+R<Z1)
 ck('entire_containing_disk_critical_free',mub*Z1<1)
 logder=1+q*Z1+q*mub*Z1**2/(1-mub*Z1)
 ck('whole_log_derivative_bound',logder<L)
 nonlinear=(A(L*R).exp()-1-A(L*R))/A(L)
 margin=A(R)-model-actual_tail-nonlinear
 ck('whole_circle_Rouche',margin>A(F(1,100)))
 displacement=sum(zm[1:])
 ck('each_circle_excludes_source',2*arb.pi()-A(displacement)>A(R))
 ck('two_complete_disks_disjoint',4*arb.pi()-A(2*displacement)>A(2*R))
 out={'status':'PASS','scope':'entire q<=1/20,mu<=1/50 and both actual companion disks',
      'precision_bits':ctx.prec,'parameter_expansion_method':'coordinate_partial_derivative_recurrence',
      'integration_method':'finite_full_antiderivative_endpoint_sum',
      'path_method':'complex_Bernstein_convex_hull_128_closed_intervals',
      'path_monomials':len(AA),'whole_path_control_points':len(AA)*128*5*2,
      'center_upper':str(Zbound),'finite_residual':str(finite),'infinite_model_tail':str(tail),
      'Mcenter':str(mc),'Mrelative':str(mr),'Mpath':str(sum(am)),'R5':str(r5),
      'model_relative':str(model),'finite_product_relative':str(actual_tail),
      'whole_circle_nonlinear':str(nonlinear),'strict_margin':str(margin),
      'weighted_finite_bounds':list(map(str,weighted)),
      'path_coefficients':paths,'integrated_coefficients':integrated,'predicates':PRED,
      'checks':len(PRED),'seconds':time.time()-start,
      'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path(__file__).with_name('tang_four_moment_independent_receipt.json').write_text(json.dumps(out,indent=2))
 print(json.dumps({k:v for k,v in out.items() if k not in ['weighted_finite_bounds','path_coefficients','integrated_coefficients','predicates']},indent=2))

if __name__=='__main__':main()
