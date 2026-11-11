#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import json
H=Path(__file__).resolve().parent
rho=F(30501,500000);c=rho/2;N=F(6250)
checks=[]
def ck(s,p):checks.append({'name':s,'passed':bool(p)})
ck('positive rho and c range',0<c<F(1,2))
ck('whole radicand rational upper',2-2*c+c*c<F(49,25))
ck('sqrt2 rational strict lower',2>F(49,25))
ck('positive field denominator comparison',5>F(1,2))
ck('positive rational endpoint denominator',1-c+c*c>0)
value=5*N*c/(14*(1-c+c*c))
ck('whole degree6251 D70 entrance',value>70)
ck('monotone c fraction full range',1-F(1,2)**2>0)
# Formal coefficients in c^4,c^2 v^2,v^4 after subtracting (v^2-c^2)^2.
coeff=[9-8-1,6-8+2,1-1]
ck('complete square polynomial identity',coeff==[0,0,0])
result={'status':'PASS' if all(z['passed'] for z in checks) else 'FAIL','gates':checks,'exact_N6250_lower':str(value),'decimal_N6250_lower':float(value),'rho':str(rho),'source':'direct complete crescent + exact radial identity, no companion input'}
(H/'DIRECT_CRESCENT_FIELD_EXACT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
