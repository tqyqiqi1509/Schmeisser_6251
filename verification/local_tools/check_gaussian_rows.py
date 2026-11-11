#!/usr/bin/env python3
"""All five complete Gaussian companion rows in the article, over Q."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
gates = []

def check(name, lower, upper):
    margin = upper - lower
    assert margin > 0, (name, lower, upper)
    gates.append({"name": name, "positive_margin": str(margin)})

def exp_upper(x):
    assert 0 <= x < 26
    return sum((x**j / factorial(j) for j in range(25)), F(0)) + \
           x**25 / (factorial(25) * (1-x/26))

def multiply(a, b):
    ans = [F(0)] * 12
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i+j < 12:
                ans[i+j] += x*y
    return ans

def value(a, q):
    return sum((x*q**j for j, x in enumerate(a)), F(0))

Z = list(map(F, [F(63,10), 21, 133, 1030, 9000]))
M = list(map(F, [0, F(63,10), 45, 400, 3400, 86000, 320000,
                1700000, 9100000, 39000000]))
residual = list(map(F, [82000, 350000, 1900000, 13000000,
                       120000000, 3100000000, 23000000000]))
power, exponential = [F(1)], [F(0)] * 12
for k in range(12):
    for j, x in enumerate(power):
        exponential[j] += x / factorial(k)
    power = multiply(power, M)
P11 = multiply(Z, exponential)

def full_model(q):
    return sum((x*q**(j+5) for j, x in enumerate(residual)), F(0)) + \
           value(Z,q)*exp_upper(value(M,q))-value(P11,q)

for q, bound in [(F(1,20),132200), (F(97,2000),126000),
                 (F(7,625),86200), (F(9,556),89000)]:
    check("full positive model tail at " + str(q), full_model(q), bound*q**5)

# name, q, mu, Z0, M0, E0, model, perturb, A, B, fixed R, Z1, L, beta, phase.
rows = [
 ("coarse","97/2000","1/270","39/5","101/200","43/25",126000,550,
  300000,1300,None,"41/5","143/100","27/100",None),
 ("fine","7/625","1/200","657/100","2/25","109/100",86200,175,
  120000,220,None,"659/100","11/10","1/100","4"),
 ("wide","9/400","1/100","343/50","43/250","5/4",132200,250,
  250000,400,None,"174/25","59/50","3/50","23/5"),
 ("origin","33/1000","8/625","719/100","7/25","71/50",132200,350,
  None,None,"3/10","15/2","32/25","11/50",None),
 ("medium","9/556","1/139","6681/1000","29/250","57/50",89000,215,
  120000,225,None,"168/25","28/25","3/200","4"),
]
for row in rows:
    name,q,mu,Z0,M0,E0,model,perturb,ca,cb,fixed,Z1,L,beta,phase = row
    q,mu,Z0,M0,E0,Z1,L,beta = map(F,(q,mu,Z0,M0,E0,Z1,L,beta))
    R = F(fixed) if fixed else ca*q**5+cb*q*mu
    r3 = q*mu*Z0**3/(3*(1-mu*Z0))
    ck = lambda label,left,right:check(name+"/"+label,left,right)
    ck("center",value(Z,q),Z0)
    ck("path",value(M,q),M0)
    ck("critical-free log",mu*Z0,F(1))
    ck("full center exponential",exp_upper(M0+r3),E0)
    ck("full residual model",full_model(q),model*q**5)
    if name == "medium":
        Mr = list(map(F,[0,F(63,10),F(2653,10),F(10442,5),18282,
                       F(174349,2),325990,1727450,9270000,40500000]))
        Mr0=F(183,1000)
        ck("relative path",value(Mr,q),Mr0)
        coefficient=Z0**4*exp_upper(Mr0+2*r3)/(12*(1-mu*Z0))
        ck("both actual logarithmic losses",coefficient,F(perturb))
        ck("model radius reserve",E0*model,(1-beta)*ca)
        ck("actual radius reserve",F(perturb),(1-beta)*cb)
    else:
        ck("primitive perturbation coefficient",
           Z0**4*E0/(12*(1-mu*Z0)),F(perturb))
        if fixed:
            ck("total radius reserve",E0*(model*q**5+perturb*q*mu),(1-beta)*R)
        else:
            ck("model radius reserve",E0*model,(1-beta)*ca)
            ck("actual radius reserve",E0*perturb,(1-beta)*cb)
    ck("containing circle",Z0+R,Z1)
    ck("critical-free circle",mu*Z1,F(1))
    ck("whole logarithmic derivative",1+q*Z1+q*mu*Z1**2/(1-mu*Z1),L)
    ck("complete exponential circle remainder",exp_upper(L*R)-1-L*R,beta*L*R)
    correction=value([0,21,133,1030,9000],q)
    ck("source exclusion",correction+R,F(6))
    ck("pair separation",2*R,12*(1-q-17*q*q-137*q**3-313*q**4))
    if phase:
        ck("homogeneous phase",correction/q+R/q,F(157,25)*F(phase))

result = {"status":"PASS","rows":len(rows),"gates":len(gates),
          "arithmetic":"exact rational", "predicates":gates}
(HERE/"GAUSSIAN_ROWS_COMPLETE.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:v for k,v in result.items() if k!="predicates"},indent=2))
