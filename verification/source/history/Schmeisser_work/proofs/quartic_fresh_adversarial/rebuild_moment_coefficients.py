"""Independent straight-path moment reconstruction: no author imports or checker calls.
I_0=0, I_m=(1-m I_{m-1})/s at s=+/-2 pi i. Center corrections are
recovered successively from the vanishing conditions, then a_0,...,a_11
are computed from the exponential of the path polynomial.
"""
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
import sympy as sp
from fractions import Fraction
from pathlib import Path
from math import factorial
import json,time,hashlib
HERE=Path(__file__).resolve().parent
start=time.monotonic()
RR,s,t=ring('s,t',QQ)
zero,one=RR.zero,RR.one

def integrate_moments(poly):
    # Each monomial s^a t^m is integrated with the exact Laurent polynomial
    # I_m(s)=sum_{k=0}^{m-1} (-1)^k m!/(m-k)! s^(-k-1); m=0 gives zero.
    laurent={}
    for (ds,dt),v in poly.items():
        for k in range(dt):
            deg=ds-k-1
            c=v*((-1)**k)*QQ(factorial(dt),factorial(dt-k))
            laurent[deg]=laurent.get(deg,QQ.zero)+c
    laurent={k:v for k,v in laurent.items() if v}
    assert min(laurent,default=0)>=0, ('uncancelled Laurent denominator',laurent)
    return RR.from_dict({(k,0):v for k,v in laurent.items()})

def series(center,N):
    A=[zero]*(2*len(center))
    for j in range(1,len(A)):
        if j<len(center):A[j]+=t*center[j]
        A[j]+=(t*t/2)*sum((center[k]*center[j-1-k] for k in range(len(center)) if 0<=j-1-k<len(center)),zero)
    B=[one]+[zero]*N
    for n in range(1,N+1):
        B[n]=sum((j*A[j]*B[n-j] for j in range(1,min(n,len(A)-1)+1)),zero)/n
    coefficients=[]
    for n in range(N+1):
        C=sum((center[k]*B[n-k] for k in range(min(n,len(center)-1)+1)),zero)
        coefficients.append(integrate_moments(C))
    return A[1:],coefficients

center=[s]
implicit_steps=[]
for n in range(1,5):
    _,coeff=series(center+[zero],n)
    correction=-coeff[n]
    implicit_steps.append({'order':n,'uncorrected_residual':str(coeff[n].as_expr()),'derived_correction':str(correction.as_expr())})
    center.append(correction)
A,coeff=series(center,11)
assert all(x==zero for x in coeff[:5])
ss,tt,pp=sp.symbols('s t p',real=True)
# JSON comparison happens after the independent construction, and is evidence
# of matching the submitted candidate, not a source of the construction.
submission=json.loads((HERE.parent/'root_quartic_companions'/'QUARTIC_COEFFICIENTS.json').read_text())
compare={}
for key,computed in [('center',center),('coefficients',coeff),('path_exponent_coefficients',A)]:
    expected=submission[key]
    compare[key]=len(expected)==len(computed) and all(sp.expand(x.as_expr()-sp.sympify(y,locals={'s':sp.Symbol('s'),'t':sp.Symbol('t')}))==0 for x,y in zip(computed,expected))
assert all(compare.values()),compare
out={'method':'straight-path exponential coefficient recursion and I_m=(1-m I_(m-1))/s, I_0=0; center solved implicitly in orders 1..4','implicit_center_steps':implicit_steps,'center':[str(x.as_expr()) for x in center],'coefficients':[str(sp.factor(x.as_expr())) for x in coeff],'path_exponent_coefficients':[str(sp.factor(x.as_expr())) for x in A],'submission_comparisons':compare,'all_match':all(compare.values()),'elapsed_seconds':time.monotonic()-start}
(HERE/'INDEPENDENT_MOMENT_COEFFICIENTS.json').write_text(json.dumps(out,indent=2))
# Emit a self-contained genuine-WL program from our derived expressions.
from sympy.printing.mathematica import mathematica_code
W=lambda x:mathematica_code(x)
checks=[]
for j,x in enumerate(coeff[:5]):checks.append(f'"moment_zero_{j}" -> ({W(x.as_expr())}===0)')
for i,(name,seq,bounds) in enumerate([
 ('path',A,[sp.Rational(63,10),45,400,3400,86000,320000,1700000,9100000,39000000]),
 ('residual',coeff[5:],[82000,350000,1900000,13000000,120000000,3100000000,23000000000]),
 ('center',center[1:],[21,133,1030,9000])]):
    for j,(poly,bd) in enumerate(zip(seq,bounds)):
        ex=poly.as_expr().subs({sp.Symbol('s'):sp.I*pp,sp.Symbol('t'):tt})
        mod=sp.expand(ex*sp.conjugate(ex))
        assert not mod.has(sp.I)
        dom='6<=p<=44/7 && 0<=t<=1' if name=='path' else '6<=p<=44/7'
        var='{p,t}' if name=='path' else 'p'
        label=f'{name}_{j+5 if name=="residual" else j+1}'
        checks.append(f'"{label}" -> (Reduce[{dom} && ({W(mod)}) >= ({W(sp.sympify(bd))})^2,{var},Reals]===False)')
program='Clear["Global`*"];\nchecks=<|\n'+',\n'.join(checks)+'\n|>;\n<|"method"->"independent moment reconstruction; full real domains", "count"->Length[checks],"all_true"->And@@Values[checks],"checks"->checks|>\n'
(HERE/'INDEPENDENT_COEFFICIENT_DOMAIN.wl').write_text(program)
print(json.dumps({'implicit_center_steps':implicit_steps,'submission_comparisons':compare,'all_match':out['all_match'],'Wolfram_checks':len(checks),'elapsed_seconds':out['elapsed_seconds']},indent=2))
