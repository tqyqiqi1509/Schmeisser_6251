"""Exact Bernstein verification of the complete Gaussian coefficient domains.

Input is the independently reconstructed exact polynomials, not a numerical
receipt.  All accepting predicates use fractions.  The fixed boxes are the
larger original domain 6 <= p <= 44/7, 0 <= t <= 1, not samples of 2*pi.
Bernstein subdivision uses exact de Casteljau convex combinations.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'source/history/Schmeisser_work/proofs/root_quartic_companions/QUARTIC_COEFFICIENTS.json'
data = json.loads(SOURCE.read_text())
s, p, t = sp.symbols('s p t', real=True)


def bernstein(poly, bounds):
    degrees = poly.degree_list()
    coeff = {powers: F(int(v.p), int(v.q)) for powers, v in poly.terms()}
    for axis, (left, right) in enumerate(bounds):
        new = {}
        for powers, value in coeff.items():
            j = powers[axis]
            for k in range(j + 1):
                v = value * comb(j, k) * left ** (j-k) * (right-left) ** k
                if v:
                    powers1 = list(powers); powers1[axis] = k; powers1 = tuple(powers1)
                    new[powers1] = new.get(powers1, F(0)) + v
        coeff = new
    for axis, degree in enumerate(degrees):
        new = {}
        for powers, value in coeff.items():
            k = powers[axis]
            for j in range(k, degree + 1):
                powers1 = list(powers); powers1[axis] = j; powers1 = tuple(powers1)
                new[powers1] = new.get(powers1, F(0)) + value * F(comb(j, k), comb(degree, k))
        coeff = new
    return degrees, {powers: coeff.get(powers, F(0)) for powers in product(*(range(d+1) for d in degrees))}


def subdivide(degrees, coeff, axis):
    left, right = {}, {}
    other = [i for i in range(len(degrees)) if i != axis]
    for values in product(*(range(degrees[i]+1) for i in other)):
        key = [0] * len(degrees)
        for i, v in zip(other, values): key[i] = v
        row = []
        for j in range(degrees[axis]+1):
            key[axis] = j; row.append(coeff[tuple(key)])
        a = [row[0]]; b = [row[-1]]
        while len(row) > 1:
            row = [(x+y)/2 for x,y in zip(row,row[1:])]
            a.append(row[0]); b.append(row[-1])
        b.reverse()
        for j, (x,y) in enumerate(zip(a,b)):
            key[axis] = j; left[tuple(key)] = x; right[tuple(key)] = y
    return left, right


def certify(name, expression, bound, path=False):
    expr = sp.sympify(expression, locals={'s':s, 't':t}).subs(s, sp.I*p)
    squared = sp.expand(expr*sp.conjugate(expr))
    assert not squared.has(sp.I)
    variables = (p,t) if path else (p,)
    poly = sp.Poly(sp.Rational(bound.numerator,bound.denominator)**2-squared, *variables, domain=sp.QQ)
    boxes = [(F(6),F(44,7))]+([(F(0),F(1))] if path else [])
    degrees, controls = bernstein(poly,boxes)
    queue = [(controls,0)]; leaves=0; smallest=None; depthmax=0
    while queue:
        c,depth=queue.pop(); minimum=min(c.values())
        if minimum>0:
            leaves+=1; depthmax=max(depthmax,depth)
            smallest=minimum if smallest is None else min(smallest,minimum)
            continue
        assert depth<20, (name,depth,str(minimum))
        # Path dependency is predominantly in t. Alternate axes only once
        # t has received its first split, preserving complete closed boxes.
        axis=(1 if depth%3!=2 else 0) if path else 0
        left,right=subdivide(degrees,c,axis)
        queue.extend([(left,depth+1),(right,depth+1)])
    return {'name':name,'bound':str(bound),'degrees':list(map(int,degrees)),
            'positive_closed_boxes':leaves,'maximum_depth':depthmax,
            'minimum_bernstein_coefficient':str(smallest)}


rows=[]
for label,seq,bounds,path in [
    ('path',data['path_exponent_coefficients'],[F(63,10),F(45),F(400),F(3400),F(86000),F(320000),F(1700000),F(9100000),F(39000000)],True),
    ('residual',data['coefficients'][5:],[F(82000),F(350000),F(1900000),F(13000000),F(120000000),F(3100000000),F(23000000000)],False),
    ('center',data['center'][1:],[F(21),F(133),F(1030),F(9000)],False)]:
    for j,(expr,bound) in enumerate(zip(seq,bounds)):
        rows.append(certify(f'{label}_{j+(5 if label=="residual" else 1)}',expr,bound,path))
assert all(sp.sympify(x)==0 for x in data['coefficients'][:5])
out={'status':'PASS','arithmetic':'exact rational tensor Bernstein with exact dyadic de Casteljau subdivision',
     'domain':'6 <= p <= 44/7 and, for paths, 0 <= t <= 1',
     'bounds':len(rows),'closed_boxes':sum(x['positive_closed_boxes'] for x in rows),
     'coefficient_zero_identities':5,'rows':rows}
(HERE/'GAUSSIAN_DOMAINS_NEW_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
