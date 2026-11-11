#!/usr/bin/env python3
"""Complete coefficient sums for both same-pair statements in the article.

Existence and the variable radius are supplied by the W8 Q06 complete-circle
verification. This verifies only its explicitly stated pair consequences.
The center coefficients are independently checked over Q[s] by that verifier.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from flint import arb, acb, ctx

ctx.prec = 448
HERE = Path(__file__).resolve().parent
source = HERE.parent / "source/agents/tang_weighted_moment_center_W8.json"
raw = source.read_bytes()
data = json.loads(raw)

def A(x):
    x = F(x)
    return arb(x.numerator) / x.denominator

s = acb(0, 2 * arb.pi())
results = []
for weight,q,mu,radius,root_cap,mean_cap,phase_cap in [
    (4,F(1,20),F(1,50),F(9,20),F(17,2),F(19,10),F(7)),
    (8,F(57,1000),F(1,51),
     64*F(57,1000)**2+95000*F(57,1000)*F(1,51)**3,
     F(8917,1000),F(2109,1000),F(37,5)),
]:
    displacement, mean, terms = arb(0), arb(0), []
    for row in data["coefficients"]:
        a,b,c = row["monomial"]
        if a+b+c == 0 or a+2*b+3*c > weight:
            continue
        scale=q**a*(q*mu)**b*(q*mu*mu)**c
        coefficients=[F(x) for x in row["coefficients"]]
        value=sum((acb(A(x))*s**j for j,x in enumerate(coefficients)),acb(0))
        even=sum((acb(A(x))*s**j for j,x in enumerate(coefficients)
                  if j % 2 == 0),acb(0))
        displacement += abs(value)*A(scale)
        mean += abs(even)*A(scale)
        terms.append({"monomial":row["monomial"],"scale":str(scale),
                      "modulus":str(abs(value)),"even_modulus":str(abs(even))})
    root=2*arb.pi()+displacement+A(radius)
    same_pair_mean=mean+A(radius)
    phase=(displacement+A(radius))/(2*arb.pi()*A(q))
    assert root < A(root_cap)
    assert same_pair_mean < A(mean_cap)
    assert phase < A(phase_cap)
    results.append({"weight":weight,"q":str(q),"mu":str(mu),"radius":str(radius),
                    "absolute_root_upper":str(root),
                    "same_pair_mean_upper":str(same_pair_mean),
                    "phase_upper":str(phase),"nonconstant_terms":len(terms),
                    "terms":terms})
result={"status":"PASS","rows":results,"coefficient_sha256":hashlib.sha256(raw).hexdigest()}
(HERE/"PAIR_CONSEQUENCES_COMPLETE.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({**result,"rows":[{k:v for k,v in row.items() if k!="terms"}
                                  for row in results]},indent=2))
