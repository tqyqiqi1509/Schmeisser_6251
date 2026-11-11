"""Exact weighted companion center for a full cubic/quartic log model.

theta has weight1, kappa=+P3 weight2, lambda=-P4 weight3.
This verifies all formal cancellations and whole imaginary-s parameter
coefficient bounds.  Actual roots still require a full-tail Rouche proof.
"""
from pathlib import Path
import sympy as S
import json,math,hashlib,time

s,X,Y,V,p,t=S.symbols('s X Y V p t',real=True)
N=4
c=-s**2/2+s
d=s**3/2-2*s**2+3*s
e=-5*s**4/8+S.Rational(23,6)*s**3-S.Rational(21,2)*s**2+15*s
f=S.Rational(7,8)*s**5-S.Rational(22,3)*s**4+S.Rational(88,3)*s**3-72*s**2+105*s
u=-s**3/3+s**2-2*s
v=S.Rational(5,6)*s**4-S.Rational(13,3)*s**3+12*s**2-20*s
w=-S.Rational(7,4)*s**5+S.Rational(157,12)*s**4-51*s**3+131*s**2-210*s
r=s**5/3-S.Rational(13,6)*s**4+S.Rational(26,3)*s**3-22*s**2+40*s
l=-s**4/4+s**3-3*s**2+6*s
m=S.Rational(3,4)*s**5-S.Rational(19,4)*s**4+18*s**3-51*s**2+90*s
Z=[s,c*X,d*X**2+u*Y,e*X**3+v*X*Y+l*V,
   f*X**4+w*X**2*Y+r*Y**2+m*X*V]

def mul(a,b):
    return [S.expand(sum(a[j]*b[k-j] for j in range(k+1))) for k in range(N+1)]

def power(a,k):
    out=[S.Integer(1)]+[S.Integer(0)]*N
    for _ in range(k):out=mul(out,a)
    return out

def bound_ip(poly):
    val=S.expand(poly.subs(s,S.I*p));ns=S.expand(val*S.conjugate(val))
    q=S.Poly(S.expand(ns.subs(p,6+S.Rational(2,7)*t)),t);degree=q.degree()
    b=[sum(q.nth(j)*S.binomial(k,j)/S.binomial(degree,j) for j in range(k+1))
       for k in range(degree+1)]
    hi=max(b);num,den=hi.as_numer_denom();integer=math.isqrt(int(num)//int(den))+1
    assert S.Integer(integer)**2>hi
    return {'polynomial':str(poly),'bound':integer,'s_domain':'s=+-ip, 6<=p<=44/7',
            'Bernstein_squared_upper':str(hi)}

def main():
    started=time.time();delta=Z.copy();delta[0]=0
    expd=[S.Integer(1)]
    for j in range(1,N+1):
        expd.append(S.expand(sum(k*delta[k]*expd[j-k] for k in range(1,j+1))/j))
    powers=[power(Z,j) for j in range(9)]
    J=[S.Integer(0)]*(N+1)
    for aa in range(5):
      for bb in range(3):
       for cc in range(2):
        weight=aa+2*bb+3*cc
        if weight>N:continue
        order=2*aa+3*bb+4*cc
        P=[S.expand(sum((-1)**(order-j)*S.factorial(order)/S.factorial(j)*powers[j][k]
                        for j in range(order+1))) for k in range(N+1)]
        term=mul(expd,P);term[0]-=(-1)**order*S.factorial(order)
        coeff=X**aa*Y**bb*V**cc/(2**aa*3**bb*4**cc*S.factorial(aa)*S.factorial(bb)*S.factorial(cc))
        for k in range(weight,N+1):J[k]+=coeff*term[k-weight]
    J=[S.expand(x) for x in J];assert all(x==0 for x in J),J
    bounds={name:bound_ip(poly) for name,poly in [('c',c),('d',d),('e',e),('f',f),
              ('u_kappa',u),('v_theta_kappa',v),('w_theta2_kappa',w),('r_kappa2',r),
              ('l_lambda',l),('m_theta_lambda',m)]}
    result={'status':'PASS_FORMAL_CANCELLATIONS_AND_FULL_COEFFICIENT_BOUNDS',
            'warning':'Does not by itself prove actual companion roots; whole log tail and Rouche circle still required',
            'model':'J(z)=integral_0^z exp(t+theta*t^2/2+kappa*t^3/3+lambda*t^4/4)dt',
            'weights':{'theta':1,'kappa':2,'lambda':3},'assumption':'exp(s)=1, s=+-2*pi*i',
            'vanishing_weighted_coefficients':[str(x) for x in J],
            'coefficient_bounds':bounds,'seconds':time.time()-started,
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('tang_moment_companion_center_receipt.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':result['status'],'seconds':result['seconds'],
                      'bounds':{k:v['bound'] for k,v in bounds.items()}},indent=2))

if __name__=='__main__':main()
