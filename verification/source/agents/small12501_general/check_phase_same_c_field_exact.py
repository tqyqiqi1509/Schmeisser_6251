#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import sympy as s,json
N,K,x,c,p,ss=s.symbols('N K x c p ss',real=True)
num=N*c*c-2*K*c;den=2*x*c+p
ids={
 'fixed_x_derivative':s.expand(s.diff(num,c)*den-num*s.diff(den,c)-2*(N*x*c*c+p*(N*c-K)))==0,
 'same_x_derivative':s.expand((2*N*c-2*K)*(2*c*c+2*ss*c+p)-(N*c*c-2*K*c)*(4*c+2*ss)-2*((N*ss+2*K)*c*c+p*(N*c-K)))==0,
 'uncapped_derivative':s.expand((N*x/2-K)*(x*x+p)-(N*x*x/4-K*x)*2*x-(K*x*x+p*(N*x/2-K)))==0,
}
N0=F(6250);rho=F(30501,500000);au=F(13,1000);K0=F(21,4)
checks=[]
def ck(name,l,r):
 checks.append(dict(name=name,left=str(l),right=str(r),passed=l>r));assert l>r,name
ck('uncapped positive entire numerator',N0*rho,4*K0)
uncap=(N0*rho*rho/4-K0*rho)/(rho*rho+1-F(197,200)**2)
ck('uncapped entire modulus985 enters160',uncap,F(160))
cc=rho-au/2
ck('capped positive entire numerator',N0*cc,2*K0)
dentry=(N0*cc*cc-2*K0*cc)/(2*rho*cc+1-F(473,500)**2)
ck('capped entire modulus946 enters160',dentry,F(160))
rows=[]
for bl in [F(473,500),F(191,200),F(97,100),F(98,100),F(99,100),F(199,200),F(997,1000),F(998,1000)]:
 D=(N0*cc*cc-8*cc)/(2*rho*cc+1-bl*bl)
 ck('medium same-c field/'+str(bl),D,F(160))
 rows.append(dict(Blo=str(bl),a_cap=str(au),K='4',D=str(D),D_decimal=float(D)))
out=dict(status='PASS' if all(ids.values()) else 'FAIL',symbolic_identities=ids,predicates=checks,uncapped985=str(uncap),capped946_entry=str(dentry),medium_rows=rows)
Path(__file__).with_name('PHASE_SAME_C_FIELD_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));assert all(ids.values())
