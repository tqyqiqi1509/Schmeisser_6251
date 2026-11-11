#!/usr/bin/env python3
"""Exact full-domain polynomial identities; no numerical acceptance."""
import json
from pathlib import Path
import sympy as s
x,y,tau,f,u,v,R,Blo,Bhi,T,U,L,N,N0,c,F=s.symbols('x y tau f u v R Blo Bhi T U L N N0 c F',real=True)
b=x+s.I*y;g=b-tau;z=b+u+s.I*v;p=1-x*x-y*y;Lg=1-(x-tau)**2-y*y
A0=1-s.conjugate(b)*g
V=f*f*A0*s.conjugate(g)+tau
abs2=lambda z:s.expand(z*s.conjugate(z))
checks={}
checks['exact_cross_coefficient']=s.expand(V-(f*f*Lg*s.conjugate(b)+tau*(1-f*f)))==0
exact=f*f*abs2(1-s.conjugate(z)*g)-abs2(g-z)
rhs=f*f*abs2(A0)-tau*tau-(1-f*f*abs2(g))*(u*u+v*v)-2*s.re(V*(u+s.I*v)).expand(complex=True)
checks['entire_real_factor_identity']=s.expand(exact-rhs)==0
Cb=f*f*p*p-(1-f*f*(x*x+y*y))*(tau*tau+R*R)+2*f*f*tau*x*(p-R*R)+f*f*tau*tau*R*R
checks['entire_radial_constant']=s.expand(Cb-(f*f*abs2(A0)-tau*tau-(1-f*f*abs2(g))*R*R))==0
lb=f*f-(1-f*f*Blo*Blo)*(T*T+U*U)
checks['whole_half_disk_linear_gate']=s.expand(lb-2*U*f*f*L*Bhi-(f*f*(1+Blo*Blo*(T*T+U*U)-2*U*L*Bhi)-(T*T+U*U)))==0
checks['full_degree_ray_ratio']=s.factor((N-c)/(F*N)-(N0-c)/(F*N0)-c*(N-N0)/(F*N*N0))==0
checks['positive_support_real_part']=s.expand(s.re(V).expand(complex=True)-(f*f*Lg*x+tau*(1-f*f)))==0
checks['complete_support_imaginary_part']=s.expand(s.im(V).expand(complex=True)+f*f*Lg*y)==0
out={'status':'PASS' if all(checks.values()) else 'FAIL','identities':checks,'method':'exact SymPy polynomial and rational reduction'}
path=Path(__file__).with_name('RIGHTMOST_HALF_DISK_MOBIUS_SYMBOLIC_EXACT.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
assert all(checks.values())
