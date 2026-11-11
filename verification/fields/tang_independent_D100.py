"""Verify the D100 closed cover by the direct point-obstruction criterion.

Fraction arithmetic reconstructs the two boxes and checks their dyadic covers.
Arb/Acb analytic quadrature encloses the Gaussian primitive. Every logarithmic
tail and phase error is included, and acceptance requires a strict ball sign.
Floating-point witnesses are interpreted as exact binary rationals.
"""
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import json, hashlib, time, argparse
from flint import arb, acb, ctx
ctx.prec=144
TOL=arb(2)**-96
TARGET=F(9,4)
ENTRANCE=100
TAIL_DEGREE=0

def ar(x):
    if isinstance(x,F):return arb(x.numerator)/x.denominator
    if isinstance(x,float):
        a,b=x.as_integer_ratio();return arb(a)/b
    return arb(x)

def repartition(records):
    counts=[]
    for branch in [0,1]:
        paths=[r['path'] for r in records if r['branch']==branch]
        assert paths and len(set(paths))==len(paths)
        pending=set(paths)
        assert all(all('0'<=s<='7' for s in path) for path in paths)
        deepest=max(map(len,pending))
        for level in range(deepest,0,-1):
            groups={}
            for path in tuple(pending):
                if len(path)==level:groups.setdefault(path[:-1],[]).append(int(path[-1]))
            for parent,digits in groups.items():
                digits.sort()
                assert len(digits)==2 and digits[0]%2==0 and digits[1]==digits[0]+1, (branch,parent,digits)
                assert parent not in pending
                for digit in digits:pending.remove(parent+str(digit))
                pending.add(parent)
        assert pending=={''},(branch,pending)
        counts.append({'branch':branch,'leaves':len(paths),'maximum_depth':deepest,
                       'exact_binary_volume':str(sum((F(1,2**len(p)) for p in paths),F(0)))})
    assert len(records)==sum(x['leaves'] for x in counts)
    return counts

def box(rec):
    unit=[[F(0),F(1)] for _ in range(4)]
    for s in rec['path']:
        j,side=divmod(int(s),2);mid=sum(unit[j])/2
        unit[j][1-side]=mid
    end=[TARGET/270,TARGET/ENTRANCE,F(13*(ENTRANCE+1),ENTRANCE**2)]
    lo,hi=end[rec['branch']:rec['branch']+2]
    qs=[ar(lo+(hi-lo)*x) for x in unit[0]]
    ps=[arb.pi()*ar(x)/2 for x in unit[1]]
    aas=[ar(2*x/TARGET) for x in unit[2]]
    bbs=[ar(2*(2*x-1)/TARGET) for x in unit[3]]
    return unit,qs,ps,aas,bbs

def integrate(f):
    out=acb.integral(f,0,1,rel_tol=TOL,abs_tol=TOL,deg_limit=256,eval_limit=10000)
    assert out.is_finite()
    return out

def analyze(rec):
    try:
        unit,(ql,qh),(pl,ph),(al,ah),(bl,bh)=box(rec)
        mu=(qh/ar(TARGET) if rec['branch']==0 else arb(1)/ENTRANCE).upper()
        if unit[3][0]<=F(1,2)<=unit[3][1]:min_beta2=arb(0)
        elif unit[3][1]<F(1,2):min_beta2=bh*bh
        else:min_beta2=bl*bl
        radius=(2*pl.cos()/ar(TARGET) if rec['branch']==0 else 2*pl.cos()/(ENTRANCE*ql))
        mp=al*al+min_beta2-radius*radius
        qp=(1-ah)**2+min_beta2-1
        pi=arb.pi();D=3*pi*pi-8
        env=9*pi*pi*mu*mu/D+(2+24*pi/D)*mu*pl.cos()+(40-3*pi*pi)/(2*D)*pl.cos()**2-pl.sin()**2/2
        mm=ql-env
        for name,margin in [('moment',mp),('qdisk',qp),('matrix',mm)]:
            if margin>0:return True,{'method':name,'margin':float(margin.lower())}
        if 'point' not in rec:return False,{'reason':'no independently valid pruning predicate','author_status':rec['status']}
        w=acb(ar(rec['point'][0]),ar(rec['point'][1]))
        theta0=acb(ar(rec['theta'][0]),ar(rec['theta'][1]))
        qm=(ql+qh)/2;pm=(pl+ph)/2;am=(al+ah)/2;bm=(bl+bh)/2
        dq=(qh-ql)/2
        core=acb(1-am,bm)
        dc=(((ah-al)/2)**2+((bh-bl)/2)**2).sqrt()
        mismatch=abs(qm*core*acb(0,-2*pm).exp()-theta0).upper()
        # Different decomposition: (q-qm) C E + qm(C-Cm) E + qm Cm(E-Em).
        h=(dq*(abs(core).upper()+dc)+qm*dc+
           2*qm*abs(core).upper()*((ph-pl)/2).sin()+mismatch).upper()
        W=abs(w).upper();W2=W*W;domain=(1-mu*W).lower()
        if not domain>0:return False,{'reason':'full logarithm disk'}
        eps=(qh*mu*W**3/(3*domain)).upper();hb=(h*W2/2).upper()
        Bc=theta0*w*w/2
        A=w.real;B=Bc.real
        phase_radius=((ph-pl)/4).sin()
        sigmid=acb(0,pm).exp()
        Plo=((w*sigmid.conjugate()).real-2*W*phase_radius).lower()
        if not Plo>0:return False,{'reason':'whole-leaf physical halfplane'}
        # Validated analytic integration, independent of Taylor coefficients.
        U0=integrate(lambda t,_:((t-1)*w+(t*t-1)*Bc).exp())
        U0=w*U0
        g0=(w+Bc).exp()
        J0=U0*g0
        JA=abs(J0);GA=abs(g0)
        Au=A.upper();Bu=B.upper()
        if TAIL_DEGREE:
            # Every term of the infinite logarithmic tail is retained.
            # Orders3..K are kept individually; k>=K+1 uses the geometric
            # majorant 1/k<=1/(K+1), t^k<=t^(K+1) on[0,1].
            tail_coeff=[(k,(qh*mu**(k-2)*W**k/k).upper())
                        for k in range(3,TAIL_DEGREE+1)]
            kk=TAIL_DEGREE+1
            tail_coeff.append((kk,(qh*mu**(kk-2)*W**kk/(kk*domain)).upper()))
        else:
            tail_coeff=[(3,eps)]
        def tail_bound(t):
            return sum((c*t**k for k,c in tail_coeff),acb(0))
        tail_one=sum((c for k,c in tail_coeff),arb(0)).upper()
        EJ=(W*integrate(lambda t,_:(Au*t+Bu*t*t).exp()*
                       ((hb*t*t+tail_bound(t)).exp()-1)).real).upper()
        dJ=(GA.upper()*((hb+tail_one).exp()-1)).upper()
        direct_ok=False;direct_upper=None
        if JA.lower()>EJ:
            phase_correction=2*GA.upper()*JA.upper()*phase_radius
            re_num=(sigmid*g0*J0.conjugate()).real+phase_correction+GA.upper()*EJ
            den=JA*JA-EJ*EJ
            up=re_num/den+dJ/(JA-EJ)-Plo/W2
            direct_ok=bool(up<0);direct_upper=float(up.upper())
        return direct_ok, {
            'method': 'direct' if direct_ok else 'unresolved',
            'direct_upper': direct_upper,
            'projection_lower': float(Plo),
            'theta_radius': float(h),
            'primitive_error': float(EJ),
            'domain_lower': float(domain),
        }
    except Exception as e:return False,{'reason':repr(e)}

def initialize_worker(target,entrance,tail_degree=0):
    global TARGET,ENTRANCE,TOL,TAIL_DEGREE
    TARGET=F(target);ENTRANCE=int(entrance);TAIL_DEGREE=int(tail_degree)
    ctx.prec=144;TOL=arb(2)**-96

def main():
    global TARGET,ENTRANCE,TAIL_DEGREE
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',default=str(Path(__file__).with_name('point_D100_y225.json')))
    ap.add_argument('--output',default=str(Path(__file__).with_name('d100_receipt.json')))
    ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args()
    p=Path(args.source);raw=p.read_bytes();data=json.loads(raw)
    TARGET=F(str(data['target']));ENTRANCE=int(data['entrance'])
    TAIL_DEGREE=int(data.get('tail_degree',0))
    assert TARGET==F(9,4) and ENTRANCE==100 and TAIL_DEGREE==0
    assert not data['unresolved']
    leaves=data['leaves'];part=repartition(leaves)
    started=time.time();fail=[];counts=Counter();max_direct=-float('inf')
    min_proj=float('inf')
    with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize_worker,
                             initargs=(str(TARGET),ENTRANCE,TAIL_DEGREE)) as pool:
        for i,(ok,detail) in enumerate(pool.map(analyze,leaves,chunksize=24)):
            counts[detail.get('method','error')]+=1
            if not ok:fail.append({'index':i,'record':leaves[i],'detail':detail})
            if detail.get('direct_upper') is not None:max_direct=max(max_direct,detail['direct_upper'])
            if 'projection_lower' in detail:min_proj=min(min_proj,detail['projection_lower'])
            if (i+1)%500==0:print(json.dumps({'processed':i+1,'unresolved':len(fail),'seconds':time.time()-started}),flush=True)
    result={'status':'PASS' if not fail else 'FAIL','target_ratio':str(TARGET),'entrance':ENTRANCE,
            'independent_methods':'Fraction partition, Acb analytic quadrature, direct Jprime/J obstruction',
            'source_sha256':hashlib.sha256(raw).hexdigest(),'precision_bits':ctx.prec,'quadrature_tolerance_bits':96,'tail_degree':TAIL_DEGREE,
            'partition':part,'counts':dict(counts),'failed_count':len(fail),'failures':fail,
            'diagnostics':{'max_direct_upper':max_direct,'min_projection':min_proj},
            'seconds':time.time()-started}
    Path(args.output).write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='failures'},indent=2))
if __name__=='__main__':main()

