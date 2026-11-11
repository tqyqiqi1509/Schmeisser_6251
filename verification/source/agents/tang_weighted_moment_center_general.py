"""Exact higher-weight implicit center for the full four-moment primitive.

The coefficient ring is QQ[s]. The only transcendental identity used is
exp(s)=1 for s=+-2*pi*i. Every center coefficient is obtained by cancelling
the full primitive's next weighted coefficient; no analytic remainder is
discarded as a theorem. Actual roots still need a separate full-tail circle.
"""
from flint import fmpq_poly as P,fmpq as Q
from math import factorial
from pathlib import Path
import json,hashlib,time,argparse

ZERO=(0,0,0);S=P([0,1])
def wt(k):return k[0]+2*k[1]+3*k[2]
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def add(a,b):
    z=dict(a)
    for k,v in b.items():
        z[k]=z.get(k,P())+v
        if not z[k]:del z[k]
    return z
def scale(a,v):return {k:p*v for k,p in a.items() if p*v}
def mul(a,b,W):
    z={}
    for k,p in a.items():
        for l,q in b.items():
            m=plus(k,l)
            if wt(m)<=W:z[m]=z.get(m,P())+p*q
    return {k:p for k,p in z.items() if p}
def monomials(W):
    return sorted(((a,b,c) for a in range(W+1) for b in range(W//2+1)
      for c in range(W//3+1) if a+2*b+3*c<=W),key=lambda k:(wt(k),k))
def residual(Z,W):
    delta=dict(Z);delta.pop(ZERO,None)
    power={ZERO:P([1])};expdelta=dict(power)
    for j in range(1,W+1):
        power=mul(power,delta,W);expdelta=add(expdelta,scale(power,Q(1,factorial(j))))
    zp=[{ZERO:P([1])}]
    for j in range(1,2*W+1):zp.append(mul(zp[-1],Z,W))
    out={}
    for k in monomials(W):
        a,b,c=k;weight=wt(k);order=2*a+3*b+4*c
        antider={}
        for j in range(order+1):
            antider=add(antider,scale(zp[j],Q((-1)**(order-j)*factorial(order),factorial(j))))
        term=mul(expdelta,antider,W)
        term=add(term,{ZERO:P([-(-1)**order*factorial(order)])})
        coefficient=Q(1,2**a*3**b*4**c*factorial(a)*factorial(b)*factorial(c))
        for kk,v in term.items():
            key=plus(k,kk)
            if wt(key)<=W:out[key]=out.get(key,P())+v*coefficient
    return {k:p for k,p in out.items() if p}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--weight',type=int,default=6)
    ap.add_argument('--output',default='agents/tang_weighted_moment_center_W6.json');a=ap.parse_args()
    started=time.monotonic();Z={ZERO:S};steps=[]
    for w in range(1,a.weight+1):
        r=residual(Z,w);assert all(wt(k)==w for k in r),r
        new={k:-v for k,v in r.items()};Z=add(Z,new)
        steps.append(dict(weight=w,added=len(new),seconds=time.monotonic()-started))
        print(json.dumps(steps[-1]),flush=True)
    final=residual(Z,a.weight);assert not final,final
    coefficients=[dict(monomial=list(k),weight=wt(k),polynomial=str(p),coefficients=list(map(str,p.coeffs())))
      for k,p in sorted(Z.items(),key=lambda z:(wt(z[0]),z[0]))]
    out=dict(status='PASS_EXACT_FORMAL_CANCELLATIONS',warning='No actual-root assertion before complete remainder and whole-circle verification',
        center_weight=a.weight,primitive='integral exp(t+theta*t^2/2+kappa*t^3/3+lambda*t^4/4)dt',
        weights=[1,2,3],assumption='exp(s)=1',coefficients=coefficients,steps=steps,
        exact_remaining_coefficients=0,seconds=time.monotonic()-started,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('coefficients','steps')},indent=2))
if __name__=='__main__':main()
