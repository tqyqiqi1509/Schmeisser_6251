"""Independent exact rational three-variable Bernstein proof of signed13/10 on radius3/50.

No author/reviewer imports. SymPy builds exact real polynomials in r,c,p^2.
Machin alternating sums enclose p^2=(2*pi)^2 in an explicit rational box;
every tensor Bernstein coefficient on all three whole boxes is positive.
The radius interval is never subdivided. No floating-point accepting gate.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sympy as s
import json,hashlib,time
start=time.time();r,c,z=s.symbols('r c z',real=True)
def rational(q):return s.Rational(q.numerator,q.denominator)
def frac(q):return F(int(s.numer(q)),int(s.denom(q)))
def atan_sum(t,n):return sum(((-1)**k*t**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
lo=16*atan_sum(F(1,5),20)-4*atan_sum(F(1,239),21)
hi=16*atan_sum(F(1,5),21)-4*atan_sum(F(1,239),20)
zl=F(3947841760435743,10**14);zh=F(3947841760435744,10**14)
assert 0<lo<hi and zl<4*lo*lo<4*hi*hi<zh
b=[s.Integer(1),s.Integer(1),3-z/2,15-23*z/6,105-88*z/3+7*z*z/8]
d=[s.Integer(0),s.Integer(1),s.Integer(4),21-5*z/4,144-44*z/3]
for j,bd in [(2,17),(3,137),(4,313)]:
    # b2,b3 decrease; b4 increases throughout this narrow whole interval.
    deriv=s.Poly(s.diff(b[j],z),z)
    vals=[frac(b[j].subs(z,rational(t))) for t in (zl,zh)]
    assert max(map(abs,vals))<bd
    if j==4:assert frac(deriv.as_expr().subs(z,rational(zl)))>0
T=s.expand(sum(b[i]*b[j]*r**(i+j)*s.chebyshevt(abs(i-j),c) for i in range(5) for j in range(5)))
U=s.expand(sum(b[i]*b[j]*r**(i+j+1)*s.chebyshevt(i+j+1,c) for i in range(5) for j in range(5)))
V=s.expand(sum(d[i]*b[j]*r**(i+j)*s.chebyshevt(abs(i-j),c) for i in range(5) for j in range(5)))
assert s.Poly(U-V,r).nth(0)==0 and s.Poly(U-V,r).nth(1)==0
assert V.subs(r,0)==0 and T.subs(r,0)==1
P=s.cancel((U-V)/r**2);W=s.cancel(V*(T-1)/r**2)
assert s.Poly(P,r,c,z) and s.Poly(W,r,c,z)
assert s.expand(P.subs(r,0)-W.subs(r,0)/2)==1-5*c*c
cap=F(3,50);m=1-cap-17*cap**2-137*cap**3-313*cap**4
assert m>0
Ap=s.Poly(s.Rational(13,10)*T-P,r,c,z)
Qp=s.Poly((1+rational(m))*Ap.as_expr()+W,r,c,z)

def convert(poly,bounds):
    """Successive exact affine transformations then tensor Bernstein maps."""
    deg=poly.degree_list()
    coefficients={mon:frac(co) for mon,co in poly.terms()}
    for axis,(left,right) in enumerate(bounds):
        new={};width=right-left
        for powers,value in coefficients.items():
            j=powers[axis]
            for k in range(j+1):
                vv=value*comb(j,k)*left**(j-k)*width**k
                if not vv:continue
                pp=list(powers);pp[axis]=k;pp=tuple(pp)
                new[pp]=new.get(pp,F(0))+vv
        coefficients=new
    # A coordinatewise power basis conversion. Explicit zeros are added by
    # the sums; all Bernstein control coefficients, including edges, remain.
    for axis,n in enumerate(deg):
        new={}
        for powers,value in coefficients.items():
            k=powers[axis]
            for j in range(k,n+1):
                pp=list(powers);pp[axis]=j;pp=tuple(pp)
                vv=value*F(comb(j,k),comb(n,k))
                new[pp]=new.get(pp,F(0))+vv
        coefficients=new
    count=1
    for n in deg:count*=n+1
    assert len(coefficients)==count
    minimum=min(coefficients.values())
    assert minimum>0
    return {'degree':list(map(int,deg)),'controls':count,'minimum_exact':str(minimum),
            'minimum_decimal':float(minimum)}
cells=[]
for l,h in [(F(-1),F(0)),(F(0),F(1,2)),(F(1,2),F(1))]:
    bounds=[(F(0),cap),(l,h),(zl,zh)]
    row={'r':['0',str(cap)],'cos':[str(l),str(h)],'p_squared':[str(zl),str(zh)],
         'A':convert(Ap,bounds),'combination':convert(Qp,bounds)}
    cells.append(row)
out={'status':'PASS','coefficient':'13/10','scope':'entire complex theta disk |theta|<=3/50',
 'arithmetic':'exact rational; no floating point acceptance',
 'method':'three-variable tensor Bernstein, including full rational pi-squared interval',
 'pi_enclosure':'Machin identity with alternating atan sums,20 and21 terms',
 'global_B_lower':str(m),'closed_boxes':len(cells),
 'total_controls':sum(x['A']['controls']+x['combination']['controls'] for x in cells),
 'cells':cells,'seconds':time.time()-start,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('ROOT_SIGNED_THETA_13_10_Q06_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='cells'}|{'minima':[
 [x['A']['minimum_decimal'],x['combination']['minimum_decimal']] for x in cells]},indent=2))
