"""Complete D25 covariance-family checker on the explicit C11A20/D30 domain.
Fresh exact cube/phase reconstruction and Acb whole-tail integrals.
Only the direct Jprime/J disk obstruction is an acceptance criterion.
No source/discovery proof module is imported.
"""
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import json, hashlib, time, argparse
from flint import arb, acb, ctx
ctx.prec=144
TOL=arb(2)**-96
TARGET=F(3)
ENTRANCE=100
TAIL_DEGREE=0

def ar(x):
    if isinstance(x,F):return arb(x.numerator)/x.denominator
    if isinstance(x,float):
        a,b=x.as_integer_ratio();return arb(a)/b
    return arb(x)

def interval(lo,hi):
    return arb((lo+hi)/2,(hi-lo)/2)

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
    end=[TARGET/30,TARGET/ENTRANCE,F(20+11*ENTRANCE,ENTRANCE**2)]
    lo,hi=end[rec['branch']:rec['branch']+2]
    qs=[ar(lo+(hi-lo)*x) for x in unit[0]]
    ps=[arb.pi()*ar(x)/4 for x in unit[1]]
    aas=[ar(2*x/TARGET) for x in unit[2]]
    bbs=[ar(2*(2*x-1)/TARGET) for x in unit[3]]
    return unit,qs,ps,aas,bbs

def cubic_coefficient_ball(qs,ps,aas,bbs,mu,kappa0,slice_data=None):
    ql,qh=qs;pl,ph=ps;al,ah=aas;bl,bh=bbs
    clo=ph.cos().lower();chi=pl.cos().upper()
    if not clo>0:raise ValueError('cubic covariance requires cos(phi)>0')
    Alo=(ql*al/2).lower();Ahi=(qh*ah/2).upper()
    corners=[-q*b/2 for q in (ql,qh) for b in (bl,bh)]
    Clo=min(c.lower() for c in corners);Chi=max(c.upper() for c in corners)
    Cmin2=arb(0) if Clo<=0<=Chi else (Clo*Clo).lower() if Clo>0 else (Chi*Chi).lower()
    rl=(4*Alo*Alo/chi-3*mu*mu*chi).lower()
    rh=min((mu*Ahi).upper(),(mu*mu*chi-4*Cmin2/chi).upper())
    if not rh>=rl:raise ValueError('empty cubic real interval')
    centers=[3*a*c/d for a in (Alo,Ahi) for c in (Clo,Chi) for d in (clo,chi)]
    variance=(mu*mu*chi-(Alo*Alo+Cmin2)/chi).upper()
    if variance<0:variance=arb(0)
    rv=(3*variance/2+mu*(qh-Alo)).upper()
    il=(min(c.lower() for c in centers)-rv).lower()
    ih=(max(c.upper() for c in centers)+rv).upper()
    if slice_data is not None:
        i=slice_data['index'];m=slice_data['count']
        assert type(i) is int and type(m) is int and m in (2,4,8,16) and 0<=i<m
        left=il;span=ih-il
        il=(left+span*ar(F(i,m))).lower()
        ih=(left+span*ar(F(i+1,m))).upper()
    mid=acb((rl+rh)/2,(il+ih)/2)
    radius=(((rh-rl)/2)**2+((ih-il)/2)**2).sqrt().upper()
    pm=(pl+ph)/2
    mismatch=abs(mid*acb(0,-3*pm).exp()-kappa0).upper()
    out=(radius+2*abs(mid).upper()*(3*(ph-pl)/4).sin()+mismatch).upper()
    return out

def integrate(f):
    out=acb.integral(f,0,1,rel_tol=TOL,abs_tol=TOL,deg_limit=256,eval_limit=10000)
    assert out.is_finite()
    return out

def analyze_model(rec):
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
        is_cubic='cubic' in rec
        kappa0=acb(ar(rec['cubic'][0]),ar(rec['cubic'][1])) if is_cubic else acb(0)
        h3=cubic_coefficient_ball((ql,qh),(pl,ph),(al,ah),(bl,bh),mu,kappa0,rec.get('cubic_slice')) if is_cubic else arb(0)
        W=abs(w).upper();W2=W*W;domain=(1-mu*W).lower()
        if not domain>0:return False,{'reason':'full logarithm disk'}
        eps=(qh*mu*W**3/(3*domain)).upper();hb=(h*W2/2).upper()
        Bc=theta0*w*w/2
        Cc=kappa0*w*w*w/3
        hc=(h3*W**3/3).upper()
        A=w.real;B=Bc.real;C=Cc.real
        phase_radius=((ph-pl)/4).sin()
        sigmid=acb(0,pm).exp()
        Plo=((w*sigmid.conjugate()).real-2*W*phase_radius).lower()
        if not Plo>0:return False,{'reason':'whole-leaf physical halfplane'}
        # Validated analytic integration, independent of Taylor coefficients.
        U0=integrate(lambda t,_:((t-1)*w+(t*t-1)*Bc+(t*t*t-1)*Cc).exp())
        U0=w*U0
        g0=(w+Bc+Cc).exp()
        J0=U0*g0
        JA=abs(J0);GA=abs(g0)
        Au=A.upper();Bu=B.upper();Cu=C.upper()
        if TAIL_DEGREE:
            # Every term of the infinite logarithmic tail is retained.
            # Orders3..K are kept individually; k>=K+1 uses the geometric
            # majorant 1/k<=1/(K+1), t^k<=t^(K+1) on[0,1].
            tail_coeff=[(k,(qh*mu**(k-2)*W**k/k).upper())
                        for k in range(4 if is_cubic else 3,TAIL_DEGREE+1)]
            kk=TAIL_DEGREE+1
            tail_coeff.append((kk,(qh*mu**(kk-2)*W**kk/(kk*domain)).upper()))
        else:
            assert not is_cubic
            tail_coeff=[(3,eps)]
        def tail_bound(t):
            return sum((c*t**k for k,c in tail_coeff),acb(0))
        tail_one=sum((c for k,c in tail_coeff),arb(0)).upper()
        EJ=(W*integrate(lambda t,_:(Au*t+Bu*t*t+Cu*t*t*t).exp()*
                       ((hb*t*t+hc*t*t*t+tail_bound(t)).exp()-1)).real).upper()
        dJ=(GA.upper()*((hb+hc+tail_one).exp()-1)).upper()
        direct_ok=False;direct_upper=None
        if JA.lower()>EJ:
            phase_correction=2*GA.upper()*JA.upper()*phase_radius
            re_num=(sigmid*g0*J0.conjugate()).real+phase_correction+GA.upper()*EJ
            den=JA*JA-EJ*EJ
            up=re_num/den+dJ/(JA-EJ)-Plo/W2
            direct_ok=bool(up<0);direct_upper=float(up.upper())
        return direct_ok,{'method':('cubic_direct' if is_cubic else 'direct') if direct_ok else 'unresolved',
                          'direct_upper':direct_upper,'projection_lower':float(Plo),
                          'theta_radius':float(h),'cubic_radius':float(h3),
                          'primitive_error':float(EJ),'domain_lower':float(domain)}
    except Exception as e:return False,{'reason':repr(e)}

def analyze(rec):
    # A single conditional P3 strip is never an unconditioned outer leaf.
    if 'cubic_slice' in rec:return False,{'reason':'unaccompanied cubic slice forbidden'}
    if 'family' not in rec:return analyze_model(rec)
    try:
        family=rec['family'];m=len(family)
        assert type(family) is list and m in (2,4,8,16)
        details=[]
        for i,child in enumerate(family):
            assert set(child)=={'theta','cubic','point','cubic_slice'}
            sl=child['cubic_slice']
            assert type(sl['index']) is int and type(sl['count']) is int
            assert sl['index']==i and sl['count']==m
            item={'branch':rec['branch'],'path':rec['path'],'status':'accept_point',**child}
            ok,detail=analyze_model(item)
            if not ok:return False,{'reason':'one cubic strip failed','strip_index':i,'strip_count':m,'detail':detail}
            details.append(detail)
        point_details=[d for d in details if d.get('method')=='cubic_direct']
        if not point_details:
            # Every member independently pruned the same entire outer box.
            return True,details[0]
        assert len(point_details)==m
        return True,{'method':'cubic_family','family_members':m,
                     'direct_upper':max(d['direct_upper'] for d in point_details),
                     'projection_lower':min(d['projection_lower'] for d in point_details),
                     'theta_radius':max(d['theta_radius'] for d in point_details),
                     'cubic_radius':max(d['cubic_radius'] for d in point_details),
                     'primitive_error':max(d['primitive_error'] for d in point_details),
                     'domain_lower':min(d['domain_lower'] for d in point_details)}
    except Exception as e:return False,{'reason':'family structure: '+repr(e)}

def initialize_worker(target,entrance,tail_degree=0):
    global TARGET,ENTRANCE,TOL,TAIL_DEGREE
    TARGET=F(target);ENTRANCE=int(entrance);TAIL_DEGREE=int(tail_degree)
    ctx.prec=144;TOL=arb(2)**-96

def main():
    global TARGET,ENTRANCE,TAIL_DEGREE
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',default='agents/tang_D25_cubic_family_cover.json')
    ap.add_argument('--output',default='agents/tang_D25_cubic_family_144bit_receipt.json')
    ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args()
    p=Path(args.source);raw=p.read_bytes();data=json.loads(raw)
    assert hashlib.sha256(raw).hexdigest()=='873eed333c984f84e9a5c4ab02f7a89d280c9830c4f7733440f926506059a7a9', 'frozen D25 source mismatch'
    TARGET=F(str(data['target']));ENTRANCE=int(data['entrance'])
    assert TARGET==3 and ENTRANCE==25
    domain={'closure_entrance':30,'closure_source_sha256':'4fc2368104d287b2d5391ebba7d16ecc6aa3f053fed221d9afef47b2a6b7eb3f','affine_intercept':20,'affine_coefficient':11,'phase_denominator':4}
    assert all(data.get(k)==v for k,v in domain.items())
    pp=F(157,50);dd=3*pp*pp-8
    assert 9*pp*pp/dd<F(21,5) and 2+24*pp/dd<F(11,2) and (40-3*pp*pp)/(2*dd)<F(1,4)
    assert ar(pp)<arb.pi() and F(71,100)**2>F(1,2)
    assert F(21,5)/25**2+F(181,200)/25-F(1,8)==-F(513,6250)
    assert F(20+11*25,25**2)==F(59,125)
    TAIL_DEGREE=int(data.get('tail_degree',0))
    assert TAIL_DEGREE==8
    assert not any('cubic' in r or 'family' in r for r in data['leaves']) or TAIL_DEGREE>=4
    assert F(9,4)<=TARGET<=3 and 0<ENTRANCE<=270 and not data['unresolved']
    leaves=data['leaves'];part=repartition(leaves)
    started=time.time();fail=[];counts=Counter();max_direct=-float('inf');max_corr=-float('inf')
    min_proj=float('inf');point_certificates=0
    with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize_worker,
                             initargs=(str(TARGET),ENTRANCE,TAIL_DEGREE)) as pool:
        for i,(ok,detail) in enumerate(pool.map(analyze,leaves,chunksize=24)):
            counts[detail.get('method','error')]+=1
            point_certificates+=detail.get('family_members',1 if detail.get('method') in ('direct','direct_and_correlated','cubic_direct') else 0)
            if not ok:fail.append({'index':i,'record':leaves[i],'detail':detail})
            if detail.get('direct_upper') is not None:max_direct=max(max_direct,detail['direct_upper'])
            if 'correlated_upper' in detail:max_corr=max(max_corr,detail['correlated_upper'])
            if 'projection_lower' in detail:min_proj=min(min_proj,detail['projection_lower'])
            if (i+1)%500==0:print(json.dumps({'processed':i+1,'unresolved':len(fail),'seconds':time.time()-started}),flush=True)
    result={'status':'PASS' if not fail else 'FAIL','target_ratio':str(TARGET),'entrance':ENTRANCE,
            'independent_methods':'Fraction partition, Acb analytic quadrature, direct Jprime/J with complete checked P3 strip families, independent cubic covariance rectangle, and full tail; no author code imported',
            'source_sha256':hashlib.sha256(raw).hexdigest(),'domain_metadata':domain,'precision_bits':ctx.prec,'quadrature_tolerance_bits':96,'tail_degree':TAIL_DEGREE,
            'partition':part,'counts':dict(counts),'point_certificates':point_certificates,'failed_count':len(fail),'failures':fail,
            'diagnostics':{'max_direct_upper':max_direct,'min_projection':min_proj},
            'seconds':time.time()-started}
    Path(args.output).write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='failures'},indent=2))
if __name__=='__main__':main()

