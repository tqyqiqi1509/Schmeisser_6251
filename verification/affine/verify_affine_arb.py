"""Independent whole-cell audit of affine Q <= 20+11H.

Only the fixed six-rational-column CSV witnesses and rectangle exclusions are
read. No author Python is imported. All parameter envelopes are reconstructed
with Fraction; complete logarithms, Gaussian moments and exponential tests are
cross-checked by 192-bit Arb, rather than the author's Taylor/Mills scheme.
This establishes a separate certificate for the same closed parameter cover.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt, factorial
from collections import Counter
import csv, json, hashlib, time
from flint import arb, ctx

ctx.prec = 192
OUT = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve().parent
C = F(11)
A = F(20)
PILO = F(6283,2000)
DEN = 3*PILO**2-8
K = 9*PILO**2/DEN
LIN = 2+24*PILO/DEN
QUAD = (40-3*PILO**2)/(2*DEN)
COUNTS = Counter()
DETAILS = []

def check(tag, condition):
    COUNTS[tag] += 1
    if not condition:
        raise AssertionError(tag)

def ball(x):
    if isinstance(x,F):return arb(x.numerator)/x.denominator
    return arb(x)

def floor_mesh(x, den=10**6):
    if isinstance(x,F):return F((x*den).__floor__(),den)
    k=(x*den).floor().unique_fmpz()
    if k is None:raise ArithmeticError('Arb floor not unique')
    return F(int(k),den)

def ceil_mesh(x, den=10**6):return -floor_mesh(-x,den)

def sqrtlo(x, den=10**6):
    ans=F(isqrt(x.numerator*den*den//x.denominator),den)
    check('directed lower square root', ans*ans<=x<(ans+F(1,den))**2)
    return ans

def sqrthi(x, den=10**6):
    ans=sqrtlo(x,den)
    ans=ans if ans*ans==x else ans+F(1,den)
    check('directed upper square root',ans*ans>=x)
    return ans

def log_coeffs(radius):
    check('full logarithm domain',0<radius<1)
    r=ball(radius)
    trueh=(r-(1+r).log())/(r*r)
    truej=(-r-(1-r).log())/(r*r)
    h=floor_mesh(trueh,10**12)
    j=ceil_mesh(truej,10**12)
    check('Arb full logarithm directed coefficients',ball(h)<=trueh and ball(j)>=truej and h>0)
    return h,j

def upper_gaussian_shell(alpha, lo, hi, degree):
    """Exact positive integral represented by a rigorous Arb enclosure."""
    aa=ball(alpha)
    def tail(x):
        xx=ball(x)
        m=arb.pi().sqrt()*((aa.sqrt()*xx).erfc())/(2*aa.sqrt())
        e=(-aa*xx*xx).exp()
        for d in range(2,degree+1,2):
            m=(xx**(d-1)*e+(d-1)*m)/(2*aa)
        return m
    result=tail(lo) if hi is None else tail(lo)-tail(hi)
    check('positive complete Gaussian shell',result>0)
    return result

def reconstruct(row, phase, index):
    tlo,thi,ylo,yhi,rho,k=row
    check('input rational rectangle',0<=tlo<thi<=F(7,10) and 0<=ylo<yhi<=F(1,4))
    check('positive contour data',rho>0 and k>0)
    dl=sqrtlo(1+tlo*tlo)
    dh=sqrtlo(1+thi*thi)
    eps=2/(C*dl)
    check('retained centered disk fiber',ylo<=eps)
    xr=1-sqrthi(eps*eps-ylo*ylo)
    yu=min(yhi,eps,sqrthi(1-xr*xr))
    check('positive real curvature fiber',xr>yu>=ylo>=0 and xr*xr+ylo*ylo<1)
    ga=QUAD-tlo*tlo/2
    pb=LIN-C*dl
    pc=K-A
    check('entrance quadratic signs',ga>0 and pb<0 and pc<0)
    positive_root=(-ball(pb)+(ball(pb*pb-4*ga*pc)).sqrt())/(2*ball(ga))
    eta0=floor_mesh(positive_root,100)
    poly=lambda x:ga*x*x+pb*x+pc
    if poly(eta0)==0:eta0-=F(1,100)
    check('whole entrance lower witness',eta0>0 and poly(eta0)<0 and poly(eta0+F(1,100))>=0)
    qmin=A+C*dl*eta0
    rootq=sqrtlo(qmin,10**4)
    curvature=sqrtlo(xr*xr+ylo*ylo)
    qa=curvature+1
    qb=-(C*curvature+1)
    qc=C
    lowerroot=(-ball(qb)-(ball(qb*qb-4*qa*qc)).sqrt())/(2*ball(qa))
    a=ceil_mesh(lowerroot,10**6)
    rouche=lambda z:(z*curvature-1)*(C-z)-z*z
    if rouche(a)==0:a+=F(1,10**6)
    check('complete finite-product Rouche circle',1<=a<=2 and a<C and rouche(a)>0)
    s=a/C
    umax=1/(1-s)
    edelta=a*a/(curvature*(C-a))
    # Independent phase ranges. The lower t/sqrt(1+t^2) endpoint uses tlo.
    qhi=1/dl
    chi=thi/dh
    qlow=1/sqrthi(1+thi*thi)
    clow=tlo/sqrthi(1+tlo*tlo)
    realhi=xr/(xr*xr+ylo*ylo)
    imaghi=yu/(xr*xr+yu*yu)
    reallow=sqrtlo(1-yu*yu)
    check('whole phase comparison domain',0<qlow<=qhi and 0<=clow<=chi and xr>yu)
    if phase==1:
        dc=qhi*realhi+chi*imaghi+edelta
        ic=max(F(0),chi*realhi-qlow*ylo,qhi*imaghi-clow*reallow)+edelta
    else:
        dc=max(F(0),qhi*realhi-clow*ylo)+edelta
        ic=chi*realhi+qhi*imaghi+edelta
    xx=min(s,dc/C)
    yy=min(s,ic/C)
    hmax=1/(C*dl)
    l=(1-s)**2
    dual=(sqrthi(l*l+4*hmax*hmax*(1+yy)**2)-l)/(2*hmax)
    check('saddle displacement dual root enclosure',hmax*dual*dual+l*dual-hmax*(1+yy)**2>=0)
    xx=min(xx,dual)
    vold=xx*((1+xx)/(1-s*s)+xx*umax*umax)
    vnew=xx*(1-yy+xx)/((1-yy)**2+xx*xx)
    var=min(vold,vnew)
    check('variance full unbounded-domain entrance',0<=xx<1-yy and 0<=var<F(1,2))
    den=1-2*xx*xx-2*yy
    check('saddle moment mass denominator',den>0)
    masshi=min(umax*umax,1/den)
    masslo=1/(1+2*yy+s*s)
    change=umax*umax-1
    realcurv=max(xr-change,(1-2*var)*masslo)
    imagcurv=min(yu+change,2*masshi*sqrthi(var*(1-var)))
    check('complete saddle moments',realcurv>0 and imagcurv>=0 and 0<masslo<=masshi)
    aspect=floor_mesh(rho/(k/rootq+dc/C),10**4)
    tmin=floor_mesh(rho*rootq,10**4)
    check('contour aspect and ray entrance',aspect>1 and tmin>0 and aspect*(k/rootq+dc/C)<=rho)
    radius=sqrthi(rho*rho*(1+1/(aspect*aspect)),10**5)
    check('entire translated contour domain',radius*umax<1)
    hv,jv=log_coeffs(umax*rho)
    hc,jc=log_coeffs(umax*radius)
    vertical=hv-(hv+jv)*var
    oldneg=hc-(hc+jc)*var
    oldpos=jc-(hc+jc)*var
    check('vertical damping sign',vertical>0)
    kappas=[]
    if oldpos>=0 and oldneg-oldpos/(aspect*aspect)>0:
        old=masslo*(oldneg-oldpos/(aspect*aspect))-(hc+jc)*imagcurv/aspect
        if old>0:kappas.append(('translated',old))
    # Direct complete-log estimate in the original critical coordinates.
    inv=1/aspect
    shift=ic/(aspect*dc)
    J=1+shift
    orig_radius=sqrthi(rho*rho*(inv*inv+J*J),10**5)
    direct=None
    if shift<1 and orig_radius<1:
        ho,jo=log_coeffs(orig_radius)
        neg=ho-(ho+jo)*hmax
        if neg>0:
            linear=(ho+jo)*yu*inv+chi/(aspect*dc)
            vertex=linear/(2*neg)
            left=1-shift
            polyv=lambda x:-neg*x*x+linear*x
            values=[polyv(left),polyv(J)]
            if left<=vertex<=J:values.append(polyv(vertex))
            maximum=max(values)
            true_f=((ball(s)/(1-ball(s)))+(1-ball(s)).log())/(ball(s)**2)
            fupper=ceil_mesh(true_f,10**12)
            check('complete saddle normalization series',ball(fupper)>=true_f and fupper>0)
            direct=-maximum-jo*inv*inv-qhi*k/(aspect*dc*tmin)-fupper*a*a/(aspect*aspect*dc*dc)
            if direct>0:kappas.append(('original',direct))
    check('at least one valid whole-connector bound',bool(kappas))
    method,kappa=max(kappas,key=lambda z:z[1])
    # The entire paired logarithmic error, including the final infinite shell.
    edges=[F(z) for z in [0,1,2,3,4,5,6,8,10,12]]
    error=arb(0)
    for j,start in enumerate(edges):
        stop=edges[j+1] if j+1<len(edges) else None
        local=rho if stop is None else min(rho,stop/rootq)
        rr=umax*local
        hh,jj=log_coeffs(rr)
        alpha=min(realcurv/2,masslo*(hh-(hh+jj)*var))
        check('whole shell damping sign',alpha>0)
        odd=F(1,3)+rr*rr/(5*(1-rr*rr))
        even=F(1,4)+rr*rr/(6*(1-rr*rr))
        m4=upper_gaussian_shell(alpha,start,stop,4)
        m6=upper_gaussian_shell(alpha,start,stop,6)
        error+=ball(2*masshi*umax*umax*even/(rootq*rootq))*m4
        error+=ball(masshi*masshi*umax*umax*odd*odd/(rootq*rootq))*m6
    gaussian=(2*arb.pi()/ball(masshi)).sqrt()
    primitive=(gaussian-error-ball(F(1,1000))-ball(F(1,100)))/2
    check('full finite-polynomial primitive lower bound',primitive>0)
    # Direct Arb exponential comparisons, independent of E2 and E24 witnesses.
    tail=ball(2/(realcurv*tmin))*(-ball(realcurv*tmin*tmin/2)).exp()
    check('both complete Gaussian tails',tail<ball(F(1,1000)))
    x=kappa*tmin*tmin
    check('whole connector ray monotone',x>F(1,2))
    connectors=ball(2*tmin/aspect)*(-ball(x)).exp()
    check('both complete horizontal connectors',connectors<ball(F(1,100)))
    check('whole original-root quotient ray monotone',kappa*(k*k+J*J*tmin*tmin)>J*J)
    factor=ball(x).exp()*ball(k)*primitive/ball(k*k+J*J*tmin*tmin)
    check('strict actual original-root quotient contradiction',factor>ball(F(11,10)))
    rec={'phase':phase,'index':index,'row':[str(z) for z in row],
         'connector':method,'quotient_factor_lower':str(F(10**6) if factor>10**6 else floor_mesh(factor,10**8)),
         'paired_error_upper':str(ceil_mesh(error,10**8)),
         'primitive_lower':str(floor_mesh(primitive,10**8)),
         'damping_lower':str(floor_mesh(kappa,10**8)),
         'Q_lower':str(qmin),'saddle_radius_upper':str(s),
         'used_true_Gaussian_tails':True}
    DETAILS.append(rec)

def cover(rows, excluded, phase):
    rects=[r[:4] for r in rows]+excluded
    for a,b,c,d in rects:
        check('cover rectangle domain',0<=a<b<=F(7,10) and 0<=c<d<=F(1,4))
    for a,b,c,d in excluded:
        eps=2/(C*sqrtlo(1+a*a))
        impossible=c>eps
        if not impossible:
            r=1-sqrthi(eps*eps-c*c)
            impossible=(r*r+c*c>=1) or (QUAD-a*a/2<=0)
        check('analytic closed exclusion',impossible)
    xs=sorted({x for r in rects for x in r[:2]})
    ys=sorted({y for r in rects for y in r[2:]})
    check('cover rectangle extreme boundaries',xs[0]==ys[0]==0 and xs[-1]==F(7,10) and ys[-1]==F(1,4))
    # A separate atomic-grid verification: all open cells, all vertical and
    # horizontal open edges, and all vertices are covered by closed rectangles.
    xpoints=sorted(set(xs+[(a+b)/2 for a,b in zip(xs,xs[1:])]))
    ypoints=sorted(set(ys+[(a+b)/2 for a,b in zip(ys,ys[1:])]))
    for x in xpoints:
        for y in ypoints:
            covering=sum(1 for a,b,c,d in rects if a<=x<=b and c<=y<=d)
            check('complete closed atomic-grid cover',covering>=1)
    check('exact rectangle area',sum((b-a)*(d-c) for a,b,c,d in rects)==F(7,40))
    return {'phase':phase,'retained':len(rows),'exclusions':len(excluded),
            'atomic_x_points':len(xpoints),'atomic_y_points':len(ypoints)}

def main():
    start=time.monotonic()
    check('matrix entrance rational constants',K<A and LIN<C and 0<2*QUAD<F(49,100))
    check('complete phase rectangle',F(2,11)<F(1,4))
    check('certified pi enclosure',ball(PILO)<arb.pi()<ball(F(22,7)))
    summaries=[]
    hashes={}
    for phase in [1,-1]:
        csvpath=SOURCE/f'C11_A20_PHASE_{phase}_DIRECT_QUADRATIC_METRIC_CELLS.csv'
        skipath=SOURCE/f'C11_A20_PHASE_{phase}_DIRECT_QUADRATIC_METRIC_SKIPS.csv'
        with csvpath.open() as f:rows=[tuple(F(v) for v in row) for row in list(csv.reader(f))[1:]]
        with skipath.open() as f:skips=[tuple(F(v) for v in row[:4]) for row in list(csv.reader(f))[1:]]
        hashes[csvpath.name]=hashlib.sha256(csvpath.read_bytes()).hexdigest()
        hashes[skipath.name]=hashlib.sha256(skipath.read_bytes()).hexdigest()
        for j,row in enumerate(rows):
            reconstruct(row,phase,j)
            if j%25==0:print(json.dumps({'phase':phase,'done':j,'elapsed':time.monotonic()-start}),flush=True)
        summaries.append(cover(rows,skips,phase))
    result={'status':'PASS','arithmetic':'Fraction parameter enclosures and 192-bit Arb complete log / erfc / exp',
            'checks':sum(COUNTS.values()),'counts':dict(COUNTS),'closed_covers':summaries,
            'retained_cells':len(DETAILS),'paired_shells':10*len(DETAILS),'source_hashes':hashes,
            'least_quotient_factor_lower':str(min(F(r['quotient_factor_lower']) for r in DETAILS)),
            'connector_choices':dict(Counter(r['connector'] for r in DETAILS)),
            'elapsed_seconds':time.monotonic()-start}
    (OUT/'INDEPENDENT_ARB_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'INDEPENDENT_ARB_CELLS.json').write_text(json.dumps(DETAILS,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
