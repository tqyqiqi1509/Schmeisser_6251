"""Exact complete scalar gates for the enlarged fine quartic row mu<=1/200.

The analytic center and Gaussian residual polynomial package is unchanged
from the previously audited q<=7/625 fine row. This verifies every scalar
condition that changes when mu is increased, and the phase corollary.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
checks=[]
def ck(name,condition):
    checks.append((name,bool(condition)))
def eu(x):
    assert 0<=x<26
    return sum((x**k/F(factorial(k)) for k in range(25)),F(0))+x**25/(F(factorial(25))*(1-x/26))
q=F(7,625);mu=F(1,200)
Z0=F(657,100);M0=F(2,25);E0=F(109,100);Z1=F(659,100);L=F(11,10);beta=F(1,100)
Z=F(63,10)+21*q+133*q*q+1030*q**3+9000*q**4
M=F(63,10)*q+45*q*q+400*q**3+3400*q**4+86000*q**5+320000*q**6+1700000*q**7+9100000*q**8+39000000*q**9
nu=q*mu
ck('q lies in inherited fine Gaussian-model interval',q==F(7,625))
ck('complete quartic center bound',Z<Z0)
ck('complete quartic path exponent bound',M<M0)
ck('center logarithm domain',1-mu*Z0>0)
eps=nu*Z0**3/(3*(1-mu*Z0))
ck('entire exponential upper bound',eu(M0+eps)<E0)
P0=Z0**4*E0/(12*(1-mu*Z0))
ck('complete primitive perturbation coefficient',P0<175)
model_relative=86200*E0
pert_relative=175*E0
ck('complete relative model coefficient',model_relative==93958)
ck('complete relative perturbation coefficient',pert_relative==F(763,4))
R0=120000*q**5+220*nu
ck('containing physical circle',Z0+R0<Z1)
ck('critical-free containing disk',mu*Z1<1)
logL=1+q*Z1+nu*Z1*Z1/(1-mu*Z1)
ck('whole containing-disk logarithmic derivative',logL<L)
nonlinear=(eu(L*R0)-1-L*R0)/(L*R0)
ck('entire nonlinear Taylor error',nonlinear<beta)
ck('Rouche q5 coefficient',120000-model_relative-beta*120000>0)
ck('Rouche qmu coefficient',220-pert_relative-beta*220>0)
correction=21*q+133*q*q+1030*q**3+9000*q**4
ck('both disks avoid zero',6-correction-R0>0)
ck('both disks disjoint',12-2*correction-2*R0>0)
phase=21+133*q+1030*q*q+9000*q**3+120000*q**4+220*mu
ck('entire phase coefficient',phase<24)
Dexact=F(5625,28)
ck('exact combined entrance exceeds200',Dexact>200)
ck('ratio nine-fourths q entrance',F(9,4)/Dexact==q)
ck('convenient D201 entrance',F(9,4)/201<q and F(1,201)<mu)
fail=[name for name,ok in checks if not ok]
out={'status':'PASS' if not fail else 'FAIL','check_count':len(checks),'failures':fail,
     'q_upper':str(q),'mu_upper':str(mu),'P0_exact':str(P0),'P0_decimal':float(P0),
     'R0_exact':str(R0),'R0_decimal':float(R0),'nonlinear_upper_exact':str(nonlinear),
     'nonlinear_upper_decimal':float(nonlinear),'phase_coefficient_exact':str(phase),
     'phase_coefficient_decimal':float(phase),'D_phase_exact':str(Dexact),
     'changed_scalar_checks':checks,
     'inherited_input':'Previously audited quartic coefficient package and Gaussian primitive residual86200q^5 for q<=7/625'}
Path(__file__).with_name('fine200_phase201_exact_receipt.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['changed_scalar_checks','nonlinear_upper_exact']},indent=2))

