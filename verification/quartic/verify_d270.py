"""New author general C13/coupled-curvature D270 certificate.
Transparently adapted from accepted D299 and author D278 parameter proofs.
Author computation only: a fresh independent reviewer is required.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
P=Path(__file__).resolve().parent
check={};evidence={}
def ck(name,value,margin=None):
    check[name]=bool(value)
    if margin is not None:evidence[name]={'exact':str(margin),'diagnostic':float(margin)}
    assert value,name

def EU(x):
    assert 0<=x<26
    return sum((x**k/F(factorial(k)) for k in range(25)),F(0))+x**25/(F(factorial(25))*(1-x/26))
def pol(a,x):
    value=F(0)
    for c in a[::-1]:value=value*x+c
    return value

def mul(a,b):
    c=[F(0)]*12
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<12:c[i+j]+=x*y
    return c
Z=list(map(F,[F(63,10),21,133,1030,9000]));M=list(map(F,[0,F(63,10),45,400,3400,86000,320000,1700000,9100000,39000000]));B=list(map(F,[82000,350000,1900000,13000000,120000000,3100000000,23000000000]))
ex=[F(0)]*12;power=[F(1)]+[F(0)]*11
for k in range(12):
    for i in range(12):ex[i]+=power[i]/factorial(k)
    power=mul(power,M)
P11=mul(Z,ex)
rows=[dict(name='coarse',q=F(97,2000),mu=F(1,270),Z0=F(39,5),M0=F(101,200),E0=F(43,25),model=F(126000),perturb=F(550),A=F(300000),B=F(1300),Z1=F(41,5),L=F(143,100),beta=F(27,100)),dict(name='fine',q=F(7,625),mu=F(1,270),Z0=F(657,100),M0=F(2,25),E0=F(109,100),model=F(86200),perturb=F(175),A=F(120000),B=F(220),Z1=F(659,100),L=F(11,10),beta=F(1,100))]
receipts=[]
for row in rows:
    nm=row['name'];q=row['q'];mu=row['mu'];nu=q*mu;z=row['Z0'];m=row['M0'];e=row['E0'];rad=row['A']*q**5+row['B']*nu;z1=row['Z1'];L=row['L'];beta=row['beta']
    fullmodel=sum((B[k]*q**(k+5) for k in range(7)),F(0))+pol(Z,q)*EU(pol(M,q))-pol(P11,q)
    r3=nu*z**3/(3*(1-mu*z));pert=z**4*e/(12*(1-mu*z));betaactual=(EU(L*rad)-1-L*rad)/(L*rad)
    margins={'q_positive':q,'q_within_shared_coefficient_domain':F(1,20)-q,'Z':z-pol(Z,q),'M':m-pol(M,q),'model':row['model']-fullmodel/q**5,'log_domain':1-mu*z,'exponential':e-EU(m+r3),'perturbation':row['perturb']-pert,'containing_disk':z1-z-rad,'critical_free':1-mu*z1,'log_derivative':L-(1+q*z1+nu*z1*z1/(1-mu*z1)),'nonlinear':beta-betaactual,'rouche_q5':row['A']-row['model']*e-beta*row['A'],'rouche_nu':row['B']-row['perturb']*e-beta*row['B'],'avoids_origin':6-pol([0]+Z[1:],q)-rad,'disjoint':12*(1-q-17*q*q-137*q**3-313*q**4)-2*rad}
    for label,value in margins.items():ck(nm+'_'+label,value>0,value)
    receipts.append({k:str(v) if isinstance(v,F) else v for k,v in row.items()}|{'nu_max':str(nu),'radius_max':str(rad),'complete_model_ratio_exact':str(fullmodel/q**5),'complete_model_ratio_diagnostic':float(fullmodel/q**5),'nonlinear_fraction_exact':str(betaactual),'nonlinear_fraction_diagnostic':float(betaactual),'relative_model_coefficient':str(row['model']*e),'relative_nu_coefficient':str(row['perturb']*e)})
def phaseC(r):
    v=r+17*r*r+137*r**3+313*r**4;m=1-v;a=(1+17*r+137*r*r+313*r**3)/m;b=(18+154*r+450*r*r+313*r**3)/m
    return 1+(2+F(100,199))*r*a*a+3*r*b+F(600,199)*r*r*a*a+4*r*a+29*r+436*r*r
cq=F(97,2000);fq=F(7,625)
for name,gap in {'coarse_phase':F(9)-phaseC(cq),'fine_phase':F(11,5)-phaseC(fq),'coarse_curvature_q2':F(11)-F(9)-F(300000,18)*cq**3,'coarse_curvature_nu':73-F(1300,18),'fine_curvature_q2':F(9,4)-F(11,5)-F(120000,18)*fq**3,'fine_curvature_nu':13-F(220,18)}.items():ck(name,gap>0,gap)
D=F(270);upper=13*(1+1/D)
ck('D270_mu_entry',1/D<=F(1,270),F(1,270)-1/D)
ck('D270_q_entry',13*(D+1)/D**2<cq,cq-13*(D+1)/D**2)
ck('fine_entry_after_y3',3/D<fq,fq-3/D)
def coarseP(y,D):return y*y*phaseC(y/D)+F(300000,18)*y**5/D**3+F(1300,18)*y-D*(y-2)
for nm,y in [('coarse_left',F(3)),('coarse_right',upper)]:
    v=coarseP(y,D);ck(nm,v<0,v)
for nm,y in [('fine_left',F(9,4)),('fine_right',F(3))]:
    v=F(9,4)*y*y+13*y-D*(y-2);ck(nm,v<0,v)
ck('convex_power_series_denominator_positive',1-pol([0,1,17,137,313],cq)>0)
ck('coarse_lower_gt2',F(3)>2);ck('fine_lower_gt2',F(9,4)>2)
result={'status':'AUTHOR_DISCOVERY_PASS_NOT_INDEPENDENT_ACCEPTANCE','count':len(check),'all_true':all(check.values()),'checks':check,'evidence':evidence,'rows':receipts,'phase_C_coarse_exact':str(phaseC(cq)),'phase_C_fine_exact':str(phaseC(fq)),'new_mechanism':'retain C(y/D), use nonnegative power series convexity in y and monotonicity in D, C13 entry','C13_entry':'13(D+1)/D^2 <= 13*271/270^2 <97/2000 for D>=270'}
(P/'D270_EXACT_CERTIFICATE.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'count':len(check),'all_true':all(check.values()),'rows':receipts,'endpoint_diagnostics':{k:evidence[k] for k in ['coarse_left','coarse_right','fine_left','fine_right']}},indent=2))
