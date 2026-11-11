"""Fresh standard-library replay of affine C11,A20 witnesses; adapted from the earlier independent C13 verifier. Imports only Python standard-library modules.
Reads only the author's six-column rational witnesses and four-column skips.
No author checker, discoveries, numerical output, or precomputed checks enter.
"""
from fractions import Fraction as Q
from math import isqrt, lcm, factorial
from pathlib import Path
from collections import Counter
import csv,json,hashlib,time
HERE=Path(__file__).resolve().parent
INPUT=HERE
C=Q(11); INTERCEPT=Q(20); PI0=Q(6283,2000)
DEN=3*PI0**2-8
KU=9*PI0**2/DEN; AU=2+24*PI0/DEN; GU=(40-3*PI0**2)/(2*DEN)
EDGES=tuple(map(Q,[0,1,2,3,4,5,6,8,10,12]))
counts=Counter(); failures=[]; records=[]
def need(label,truth,where='global'):
    counts[label]+=1
    if not truth: failures.append({'where':where,'predicate':label})
    assert truth,(where,label)
def lower(x,scale=10**4):
    return Q((x*scale).__floor__(),scale)
def upper(x,scale=10**4): return -lower(-x,scale)
def sqrt_low(x,scale=10**6):
    assert x>=0
    return Q(isqrt(x.numerator*scale*scale//x.denominator),scale)
def sqrt_high(x,scale=10**6):
    q=sqrt_low(x,scale)
    return q if q*q==x else q+Q(1,scale)
def integer_coeffs(*xs):
    den=lcm(*(v.denominator for v in xs))
    return tuple(int(den*v) for v in xs)
def root_below_mesh(a,b,c,scale):
    # Positive upper root; integer radical gives its exact floor on the mesh.
    A,B,D=integer_coeffs(a,b*scale,c*scale*scale)
    disc=B*B-4*A*D
    n=(-B+isqrt(disc))//(2*A)
    if A*n*n+B*n+D==0:n-=1
    return Q(n,scale)
def first_inside_lower_root(a,b,c,scale):
    # Quadratic a*z^2+b*z+c<0; select first mesh point above its lower root.
    A,B,D=integer_coeffs(a,b*scale,c*scale*scale)
    disc=B*B-4*A*D
    sd=isqrt(disc); sud=sd if sd*sd==disc else sd+1
    n=(-B-sud)//(2*A)+1
    return Q(n,scale)
def exp_poly(x,n):
    # Horner evaluation differs from the author's term-accumulation code.
    v=Q(1)
    for k in range(n,0,-1):v=1+x*v/k
    return v
def reciprocal_exp_interval(x,where):
    need('exponential whole-tail domain',0<=x<82,where)
    e_upper=exp_poly(x,80)+x**81/Q(factorial(81))/(1-x/82)
    e_lower=exp_poly(x,40)
    need('exponential rational ordering',0<e_lower<=e_upper,where)
    return 1/e_upper,1/e_lower
def log_coefficients(r,where):
    need('complete log radius',0<r<1,where)
    # Integrated finite geometric sums, then outward rounding.
    hl=sum(((-r)**i/Q(i+2) for i in range(16)),Q(0))
    ju=sum((r**i/Q(i+2) for i in range(16)),Q(0))+r**16/(18*(1-r))
    h,j=lower(hl,10**6),upper(ju,10**6)
    need('rational log coefficient signs',0<h and 0<j,where)
    need('rational log coefficient enclosure',h<=hl and j>=ju,where)
    return h,j
def tails(alpha,a,where):
    need('positive Mills domain',alpha>0 and a>0,where)
    emin,emax=reciprocal_exp_interval(alpha*a*a,where)
    low=emin*a/(1+2*alpha*a*a)
    high=emax/(2*alpha*a)
    ans={0:(low,high)}
    for n in (2,4,6):
        low=(a**(n-1)*emin+(n-1)*low)/(2*alpha)
        high=(a**(n-1)*emax+(n-1)*high)/(2*alpha)
        need('positive recursive tail interval',0<=low<=high,where)
        ans[n]=(low,high)
    return ans
def one_case(row,phase,idx):
    tl,th,yl,yhi,rho,k=row;where=f'{phase}:{idx}'
    need('input rectangle domain',0<=tl<th<=Q(7,10) and 0<=yl<yhi<=Q(1,4),where)
    need('positive witness',rho>0 and k>0,where)
    dl=sqrt_low(1+tl*tl);dh=sqrt_low(1+th*th)
    need('slope roots directed',dl*dl<=1+tl*tl<(dl+Q(1,10**6))**2 and dh*dh<=1+th*th<(dh+Q(1,10**6))**2,where)
    eps=2/(C*dl);need('centered disk fiber',yl<=eps,where)
    r0=1-sqrt_high(eps*eps-yl*yl);g0=GU-tl*tl/2
    yh=min(yhi,eps,sqrt_high(1-r0*r0))
    need('retained disk and component domain',r0>yh>=yl>=0 and r0*r0+yl*yl<1 and g0>0,where)
    h0=root_below_mesh(g0,AU-C*dl,KU-INTERCEPT,100)
    P=lambda v:g0*v*v+(AU-C*dl)*v+KU-INTERCEPT
    need('entrance witness exact',h0>0 and P(0)<0 and P(h0)<0 and P(h0+Q(1,100))>=0,where)
    q0=INTERCEPT+C*dl*h0;rq=sqrt_low(q0,10**4)
    bm=sqrt_low(r0*r0+yl*yl)
    need('positive curvature and mass roots',0<bm and bm*bm<=r0*r0+yl*yl and 0<rq and rq*rq<=q0,where)
    a=first_inside_lower_root(bm+1,-(C*bm+1),C,10**6)
    need('saddle witness exact',1<=a<=2 and (a*bm-1)*(C-a)>a*a,where)
    need('saddle first mesh witness',( (a-Q(1,10**6))*bm-1)*(C-a+Q(1,10**6))<=(a-Q(1,10**6))**2,where)
    s=a/C;U=1/(1-s);p=U*U-1;e=a*a/(bm*(C-a))
    qu=1/dl;cu=th/dh;qlow=1/sqrt_high(1+th*th);clow=tl/sqrt_high(1+tl*tl)
    rf=r0/(r0*r0+yl*yl);yf=yh/(r0*r0+yh*yh);rl=sqrt_low(1-yh*yh)
    need('phase scalar bounds ordered',0<qlow<=qu and 0<=clow<=cu and 0<rl<=rf and 0<=yl<=yf,where)
    if phase==1:
        dc=qu*rf+cu*yf+e
        ic=max(Q(0),cu*rf-qlow*yl,qu*yf-clow*rl)+e
    else:
        dc=max(Q(0),qu*rf-clow*yl)+e
        ic=cu*rf+qu*yf+e
    X=min(s,dc/C);Y=min(s,ic/C)
    hcrit=1/(C*dl);LL=(1-s)**2
    xdual=(sqrt_high(LL*LL+4*hcrit*hcrit*(1+Y)**2)-LL)/(2*hcrit)
    need('pointwise saddle dual upper root',hcrit*xdual*xdual+LL*xdual-hcrit*(1+Y)**2>=0,where)
    X=min(X,xdual)
    old=X*((1+X)/(1-s*s)+X*U*U)
    ratio=X/(1-Y)
    new=X*(1-Y+X)/((1-Y)**2+X*X)
    beta=min(old,new)
    need('convex variance whole-cell domain',0<=ratio<1 and new==ratio*(1+ratio)/(1+ratio*ratio),where)
    qden=1-2*X*X-2*Y
    need('saddle moment positive denominators',0<s<1 and 0<=X<=s and 0<=Y<=s and qden>0 and 0<=beta<Q(1,2),where)
    qm=min(U*U,1/qden)
    R=lower(max(r0-p,(1-2*beta)/(1+2*Y+s*s)))
    I=upper(min(yh+p,2*qm*sqrt_high(beta*(1-beta))))
    M=upper(qm);G=sqrt_low(2*PI0/M,10**4)
    need('curvature and Gaussian bounds',R>0 and I>=0 and M>=qm>0 and G>0 and G*G<=2*PI0/M,where)
    lam=lower(rho/(k/rq+dc/C));T0=lower(rho*rq)
    delta=sqrt_high(rho*rho*(1+1/(lam*lam)),10**5)
    need('full contour radius and aspect',lam>1 and T0>0 and T0<=rho*rq and lam*(k/rq+dc/C)<=rho and delta*delta>=rho*rho*(1+1/(lam*lam)) and U*delta<1,where)
    hv,jv=log_coefficients(U*rho,where);hc,jc=log_coefficients(U*delta,where)
    ql=lower(1/(1+2*Y+s*s),10**6)
    av=lower(hv-(hv+jv)*beta,10**6);ac=lower(hc-(hc+jc)*beta,10**6);bc=upper(jc-(hc+jc)*beta,10**6)
    need('connector coefficient signs',ql>0 and av>0 and bc>=0,where)
    alpha=lower(min(R/2,ql*av),10**5)
    damp=lower(ql*(ac-bc/(lam*lam))-(hc+jc)*I/lam,10**5)
    # Reconstruct the new ORIGINAL-coordinate connector bound from its pointwise proof.
    invlam=1/lam;extent=1+ic/(lam*dc);shift=extent-1
    original_radius=sqrt_high(rho*rho*(invlam*invlam+extent*extent),10**5)
    direct_used=False
    if original_radius<1 and shift<1:
        ho,jo=log_coefficients(original_radius,where+'/original connector')
        vertical_coefficient=ho-(ho+jo)/(C*dl)
        if vertical_coefficient>0:
            saddle_series=upper(1/(1-s)-sum((s**i/Q(i+2) for i in range(16)),Q(0)),10**6)
            coupling=(ho+jo)*yh*invlam+cu/(lam*dc)
            left=1-shift;vertex=coupling/(2*vertical_coefficient)
            candidates=[-vertical_coefficient*left*left+coupling*left,-vertical_coefficient*extent*extent+coupling*extent]
            if left<=vertex<=extent:candidates.append(coupling*coupling/(4*vertical_coefficient))
            maximum=max(candidates)
            direct_damping=lower(-maximum-jo*invlam*invlam-qu*k/(lam*dc*T0)-saddle_series*a*a/(lam*lam*dc*dc),10**5)
            if direct_damping>damp:damp=direct_damping;direct_used=True
    need('complete vertical and horizontal damping',alpha>0 and damp>0,where)
    total=Q(0);shellinfo=[]
    for j,lo in enumerate(EDGES):
        hi=EDGES[j+1] if j+1<len(EDGES) else None
        loc=f'{where}/shell{j}'
        rad=rho if hi is None else min(rho,hi/rq)
        radU=U*rad;hh,jj=log_coefficients(radU,loc)
        al=lower(min(R/2,ql*(hh-(hh+jj)*beta)),10**5)
        need('shell real damping',0<al<Q(1,2) and al<=R/2 and al<=ql*(hh-(hh+jj)*beta),loc)
        odd=Q(1,3)+radU*radU/(5*(1-radU*radU))
        even=Q(1,4)+radU*radU/(6*(1-radU*radU))
        need('shell entire series signs',odd>0 and even>0 and 0<rad<=rho and (hi is None or rad<=hi/rq),loc)
        if lo==0:
            ar=sqrt_low(al);pr=sqrt_high(Q(22,7))
            need('origin moment square roots',ar>0 and ar*ar<=al and pr*pr>=Q(22,7),loc)
            top4=3*pr/(8*al*al*ar);top6=15*pr/(16*al**3*ar)
        else:
            tail=tails(al,lo,loc+'/left');top4=tail[4][1];top6=tail[6][1]
        if hi is None:bottom4=bottom6=Q(0)
        else:
            tail=tails(al,hi,loc+'/right');bottom4=tail[4][0];bottom6=tail[6][0]
        d4=top4-bottom4;d6=top6-bottom6
        need('shell entire positive integrals',d4>=0 and d6>=0,loc)
        err=(2*qm*U*U*even*d4+qm*qm*U*U*odd*odd*d6)/(rq*rq)
        need('shell complete paired error',err>=0,loc)
        total+=err;shellinfo.append((lo,hi,rad,al,odd,even,d4,d6,err))
    error=upper(total);F0=lower((G-error-Q(1,1000)-Q(1,100))/2)
    J0=upper(1+ic/(lam*dc));x=damp*T0*T0
    need('complete primitive positive lower bound',0<F0<=(G-error-Q(1,1000)-Q(1,100))/2 and error>=total,where)
    need('endpoint ordinate enclosure',dc>0 and J0>=1+ic/(lam*dc),where)
    need('entire Gaussian tail small',exp_poly(R*T0*T0/2,2)>2000/(R*T0),where)
    need('whole connector ray monotonic',x>Q(1,2),where)
    need('both complete connectors small',exp_poly(x,24)>200*T0/lam,where)
    need('whole quotient ray monotonic',damp*(k*k+J0*J0*T0*T0)>J0*J0,where)
    ratio_margin=exp_poly(x,24)*k*F0/(k*k+J0*J0*T0*T0)
    need('actual original-root quotient contradiction',ratio_margin>Q(11,10),where)
    # Small exact rounded diagnostic margins; no floating diagnostic is a premise.
    rec={'phase':phase,'index':idx,'input':[str(v) for v in row],
         'Q0_lower':str(lower(q0)),'beta_upper':str(upper(beta,10**6)),
         'gaussian_lower':str(G),'paired_error_upper':str(error),
         'primitive_lower':str(F0),'damping_lower':str(damp),
         'quotient_factor_lower':str(lower(ratio_margin,10**6)),
         'shells':len(shellinfo)}
    records.append(rec)
    return rec

def cover(phase,rows,excluded):
    rects=[r[:4] for r in rows]+excluded
    where=f'cover {phase}'
    need('cover nonempty',len(rects)>0,where)
    for t0,t1,y0,y1 in rects:
        need('cover rectangle domain',0<=t0<t1<=Q(7,10) and 0<=y0<y1<=Q(1,4),where)
    for t0,t1,y0,y1 in excluded:
        eps=2/(C*sqrt_low(1+t0*t0))
        need('closed exclusion',y0>eps or (y0<=eps and (1-sqrt_high(eps*eps-y0*y0))**2+y0*y0>=1) or GU-t0*t0/2<=0,where)
    xs=sorted({t for r in rects for t in r[:2]})
    need('cover extreme fibers',xs[0]==0 and xs[-1]==Q(7,10),where)
    def fiber(x,closed):
        spans=sorted((r[2],r[3]) for r in rects if ((r[0]<=x<=r[1]) if closed else (r[0]<x<r[1])))
        stop=Q(0)
        for lo,hi in spans:
            need('fiber gap absent',lo<=stop,where)
            stop=max(stop,hi)
        need('fiber full upper edge',stop==Q(1,4),where)
    for x in xs:fiber(x,True)
    for a,b in zip(xs,xs[1:]):fiber((a+b)/2,False)
    area=sum(((b-a)*(d-c) for a,b,c,d in rects),Q(0))
    need('cover exact area',area==Q(7,40),where)
    # Rectangles also have disjoint interiors, independently of the fiber proof.
    for i,a in enumerate(rects):
        for b in rects[i+1:]:
            need('partition interior disjoint',min(a[1],b[1])<=max(a[0],b[0]) or min(a[3],b[3])<=max(a[2],b[2]),where)
    return {'phase':phase,'retained':len(rows),'excluded':len(excluded),'open_slabs':len(xs)-1,'closed_fibers':len(xs),'area':str(area)}

def main():
    began=time.monotonic()
    need('compact entrance constants',KU<INTERCEPT and AU<C and 0<2*GU<Q(49,100))
    need('compact curvature domain',Q(2,11)<Q(1,4))
    hashes={}; covers=[]
    for phase in (1,-1):
        src=INPUT/f'C11_A20_PHASE_{phase}_DIRECT_QUADRATIC_METRIC_CELLS.csv'
        skip=INPUT/f'C11_A20_PHASE_{phase}_DIRECT_QUADRATIC_METRIC_SKIPS.csv'
        hashes[src.name]=hashlib.sha256(src.read_bytes()).hexdigest();hashes[skip.name]=hashlib.sha256(skip.read_bytes()).hexdigest()
        with src.open() as f:rows=[tuple(Q(x) for x in row) for row in list(csv.reader(f))[1:]]
        with skip.open() as f:excluded=[tuple(Q(x) for x in row[:4]) for row in list(csv.reader(f))[1:]]
        for i,row in enumerate(rows,1):
            one_case(row,phase,i)
            if i%25==0:print(json.dumps({'phase':phase,'verified':i,'seconds':round(time.monotonic()-began,2)}),flush=True)
        covers.append(cover(phase,rows,excluded))
    result={'status':'PASS' if not failures else 'FAIL','method':'fresh same-agent exact Fraction replay with integer algebraic root choices and Horner exponential sums','checks':sum(counts.values()),'counts':dict(sorted(counts.items())),'failures':failures,'covers':covers,'source_sha256':hashes,'elapsed_seconds':round(time.monotonic()-began,3),'retained_cells':len(records),'symmetric_shells':sum(r['shells'] for r in records),'least_quotient_factor_lower':str(min(Q(r['quotient_factor_lower']) for r in records))}
    (HERE/'FRESH_AFFINE_REPLAY_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'FRESH_AFFINE_CELL_DIAGNOSTICS.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
