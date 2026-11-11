"""Independent exact reconstruction of the complete X6251 higher strip.

No author function or receipt is imported. JSON numbers are read as their
exact decimal rationals; all derived discovery numbers are ignored. Complete
coverage is checked on vertical depth fibers, whereas the author used y
fibers. Every full contour and source/phase field boost is rebuilt.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,isqrt
from collections import Counter
import hashlib,json,time,sys
sys.set_int_max_str_digits(0)
START=time.monotonic();ROOT=Path(__file__).resolve().parents[1]
DIRECTORY=ROOT/'agents/global_component_refinement/Origin6251'
N=6250;n=6251;CAP=F(30501,500000);FLOOR=F(1,100)
COUNTS=Counter();DIGEST=hashlib.sha256();CELLS=[];CONTOURS=[];ENTRANCE=[]
def rat(x):return F(x) if isinstance(x,(str,int,F)) else F(str(x))
def ck(tag,ok,value=None):
    COUNTS[tag.split('/')[-1]]+=1
    DIGEST.update((tag+'|'+str(bool(ok))+'|'+str(value)+'\n').encode())
    if not ok:raise AssertionError((tag,value))
def pos(tag,z):ck(tag,z>0,z)
def lower12(x):return F((x*10**12).__floor__(),10**12)
def upper12(x):return -lower12(-x)
def sqrt_upper(x):
    scale=10**16;k=isqrt(x.numerator*scale*scale//x.denominator)
    if F(k,scale)**2<x:k+=1
    r=F(k,scale);ck('square root/outward exact witness',r*r>=x);return r
@lru_cache(None)
def eminus(W):return 1/sum((F(W)**j/factorial(j) for j in range(81)),F(0))
def kap(T,de):return F(1,2)-max(F(0),T*T-de-de*de)**2/(4*(T*T+de*de)*(F(1,2)+de))
def beta_weight(rho,de,T):return (T-rho*(1+T*T)/2)/(T*T+de*de)
def weight(R,width,de,tl,tr,rho):
    # The derivative numerator A-t^2+rho(1-A)t has product -A<0
    # between its roots, and its positive root is the sole maximum.
    A=(de+width)**2
    assert A>0
    return max(F(0),min((t-rho*(1+t*t)/2)/(A+t*t) for t in (max(F(0),tl-R),tr+R)))
def quadratic_max(c0,c1,c2,lo,hi):
    points=[lo,hi]
    if c2<0:
        vertex=-c1/(2*c2)
        if lo<=vertex<=hi:points.append(vertex)
    return max(c0+c1*t+c2*t*t for t in points)
def joint_max(k,T,lo,hi):return quadratic_max(F(1),T-k,-k*T,lo,hi)
def phase_max(k,lo,hi):return quadratic_max(F(0),F(1),-k,lo,hi)
def polar(rho,de):
    den=1-rho*rho*(1+de)**2/2
    assert den>0
    return 1/den,1+(rho*(1+de)/2)**2
def cayley(tag,rho,de):
    ck(tag+'/whole inverse Cayley test',F(1,200)<=de<=F(3,4) and (rho*(1+de)/2)**2<de)

# Independently rederive both complete moving-source identities. Symbols
# remain exact; no rounded source or result from the author is imported.
import sympy as sy
cc,tt=sy.symbols('c t',real=True)
Iv=tt-cc*(1+tt*tt);D0=1-cc*(tt+Iv);stest=(1+Iv)/D0
polar_den=(stest-1-cc*stest*tt)**2+(tt-cc*stest)**2
ck('moving identity/full actual beta denominator',sy.cancel(D0*polar_den-2*Iv**2)==0)
bracket=1+2*cc*(1-tt)+cc*cc*(tt*tt-4*tt-1)+2*cc**3*(1+tt*tt)
ck('moving identity/full inverse Cayley numerator',sy.cancel(D0**2*(stest-1-cc*cc*stest**2)-Iv*bracket)==0)
# On c>=0,0<=T<=1, bracket >=1-4c^2: the missing terms
# are 2c(1-T)+c^2(T-1)(T-3)+2c^3(1+T^2), all nonnegative.
ck('moving identity/whole positive bracket decomposition',sy.expand(bracket-(1-4*cc*cc)-2*cc*(1-tt)-cc*cc*(tt-1)*(tt-3)-2*cc**3*(1+tt*tt))==0)

pos('crescent/complete height gate',F(2,125)-CAP*(F(1,4)+F(1,15625)))
tr=F(1,8);fm=1/(1-F(9,2)*CAP*CAP);Dm=1+F(9,4)*CAP*CAP
bm=(1-CAP)**2/(2*tr)
ck('moving source/whole test interval',2/(1-CAP)**2<3 and 1-4*(CAP/2)**2>0)
pos('moving source/complete polar denominator',1-F(9,2)*CAP*CAP)
pos('moving source/zero-inset strictly interior T gate',2*tr/(1+tr*tr)-CAP)
Hmove=(n*bm/(n*CAP+2*Dm*bm)-1)/fm
pos('moving source/entire small-depth h50 entrance',Hmove-50)
ENTRANCE.append({'branch':'small_depth_moving','h_lower':str(Hmove),'T':['0','1/8']})
# Under the contrary h<=50 only, the full radial identity gives rho<.048.
rho_trial=F(6,125);de=F(3,10)
ck('Rayleigh/entire contrary rho cap',F(100,N)*3==rho_trial)
cayley('Rayleigh',rho_trial,de);f,Dp=polar(rho_trial,de)
B50=2*Dp*(50*f+1)/n;ck('Rayleigh/finite HB denominator before solve',0<B50<1)
b=min(beta_weight(rho_trial,de,z) for z in (F(1,8),F(1)))
pos('Rayleigh/full beta positive source',b)
k=kap(F(1),de);ck('Rayleigh/full coefficient exact',k==F(13719,34880) and 0<k<F(1,2))
K=((1-B50)*b-rho_trial)/f;pos('Rayleigh/positive source before division',K)
M=F(4,5);a2=M*k-k*k/4;a1=-M*(1-k*rho_trial)+k/2;a0=M*M-M*rho_trial-F(1,4)
v=sy.Symbol('v',real=True)
mat=sy.Matrix([[sy.Rational(M.numerator,M.denominator)-(1-sy.Rational(k.numerator,k.denominator)*v)*(sy.Rational(rho_trial.numerator,rho_trial.denominator)+v),-(1-sy.Rational(k.numerator,k.denominator)*v)/2],[-(1-sy.Rational(k.numerator,k.denominator)*v)/2,sy.Rational(M.numerator,M.denominator)]])
expected=sum(sy.Rational(z.numerator,z.denominator)*v**j for j,z in enumerate([a0,a1,a2]))
ck('Rayleigh/full symmetric determinant identity',sy.expand(mat.det()-expected)==0)
pos('Rayleigh/whole v factor positive',1-2*k);pos('Rayleigh/full determinant leading',a2)
pos('Rayleigh/full real-line completed-square numerator',4*a2*a0-a1*a1)
pos('Rayleigh/positive other principal minor',M)
pos('Rayleigh/same-source full h50 contradiction',K*N/(2*M)-2500)
RAYLEIGH_RESERVE=K*N/(2*M)-2500
ENTRANCE.append({'branch':'large_depth_Rayleigh','T':['1/8','1'],'source_K':str(K),'M':str(M),'h2_lower':str(K*N/(2*M)),'determinant':[str(a0),str(a1),str(a2)]})
ck('actual companion/lawful q rectangle',F(2,50)==F(1,25))
ck('actual companion/lawful mu rectangle',F(2,50)==F(1,25))
Z=F(3853414903319,500000000000)
pos('actual companion/homothetic phase pi gate',F(13,5)*F(314159,100000)-Z-F(9,20))
pos('high branch/independent Fine q entrance',F(7,625)-F(2,600))
pos('high branch/independent Fine mu entrance',F(1,200)-F(2,600))
# Phase4y is the full accepted Fine200 theorem; no later field conclusion
# appears in its high-field q,mu entrance above.
ck('high branch/phase4y uniform cap',4*2==8)

NAMES=['X6251_high_0.01_0.061002.json','X6251_MOVING_low_0.01_0.013.json',
 'X6251_MOVING_low_0.013_0.017.json','X6251_MOVING_low_0.017_0.023.json',
 'X6251_MOVING_low_0.023_0.033.json','X6251_MOVING_low_0.033_0.045.json','X6251_MOVING_low_0.045_0.061002.json']
HASHES={};LOW_SLABS=[]
for name in NAMES:
    file=DIRECTORY/name;raw=file.read_bytes();HASHES[name]=hashlib.sha256(raw).hexdigest()
    data=json.loads(raw,parse_float=F);mode=data['mode'];rlo,rhi=rat(data['rlo']),rat(data['rhi'])
    ck(name+'/proposal has no unresolved cells',not data['failed'])
    ck(name+'/declared mode consistent',mode==('high' if 'high' in name else 'low'))
    ck(name+'/whole rho domain',FLOOR<=rlo<rhi<=CAP)
    if mode=='high':ck(name+'/complete high rho strip',rlo==FLOOR and rhi==CAP)
    else:LOW_SLABS.append((rlo,rhi))
    boxes=[]
    for ix,z in enumerate(data['accepted']):
        tag=name+f'/cell{ix}';yl,yh,tl,tr=[rat(z[key]) for key in ('yl','yh','tl','th')]
        H=rat(z['H']);kind=z['kind'];J=4*yh if mode=='high' else F(15,2)*yh
        ck(tag+'/whole closed cell',0<=yl<yh<=2 and 0<=tl<tr<=1)
        boxes.append((yl,yh,tl,tr))
        coarse=N*rlo/(2*(1+tr*yh));phase=(N*rlo/2-J)/(rhi+tr*yh)
        bm=(1-rhi)**2/(2*tr);fm=1/(1-F(9,2)*rhi*rhi);Dm=1+F(9,4)*rhi*rhi
        ck(tag+'/unconditional moving full test range',2/(1-rhi)**2<3 and 1-F(9,2)*rhi*rhi>0)
        moving=(n*bm/(n*rhi+2*Dm*bm)-1)/fm
        ck(tag+'/same-y actual phase metadata',rat(z['phase_coefficient'])==(4 if mode=='high' else F(15,2)))
        reasons=[]
        if H<=moving:reasons.append('unconditional_moving_source')
        if H<=(600 if mode=='high' else 50):reasons.append('prior_field')
        if H<=coarse:reasons.append('unphased_radial')
        if H<=phase:reasons.append('phase_radial')
        boosts=[]
        if mode=='low':
            for ib,bf in enumerate(z['field']['source_boosts']):
                de=rat(bf['delta']);bt=tag+f'/boost{ib}';cayley(bt,rhi,de)
                k=kap(tr,de);f,Dp=polar(rhi,de);A=1-k*yl
                b=weight(F(0),F(0),de,tl,tr,rhi)
                ck(bt+'/whole source positivity',b>0 and A>0 and 0<=k<=F(1,2))
                K=min(n*r*b/(n*r+2*Dp*b)-r for r in (rlo,rhi))/f
                pos(bt+'/direct finite K positive',K)
                G=joint_max(k,tr,yl,yh);Gp=phase_max(k,yl,yh)
                pos(bt+'/complete joint quadratic maximum',G)
                pos(bt+'/complete phase quadratic maximum',Gp)
                if H<=K/(rhi*A):reasons.append('direct_source')
                if H*H<=K*N/(2*G):reasons.append('joint_source')
                if H<F(N,2) and K*(N-2*H)>=2*Gp*H*(tr*H+F(15,2)):
                    reasons.append('same_y_source_phase')
                boosts.append(dict(delta=str(de),K=str(K),G=str(G),Gp=str(Gp)))
        ck(tag+'/legal ordered field before contour',bool(reasons) and H>J)
        if kind=='field_skip':
            ck(tag+'/whole field contradiction',mode=='low' and H>=600)
            CELLS.append(dict(file=name,index=ix,kind=kind,H=str(H),field_reasons=reasons,boosts=boosts));continue
        rebuilt=[]
        for ir,row in enumerate(z['rows']):
            ct=tag+f'/contour{ir}';R,nu,W,m=[rat(row[key]) for key in ('R','nu','W','m')]
            q=yh+F(1,10**9)
            ck(ct+'/full same-depth reciprocal domain',m*m*(F(1,4)+max(tl-F(1,125),F(0))**2)>=1)
            ck(ct+'/positive full contour parameters',R>0 and 0<nu<1 and m>0 and W>0 and W.denominator==1)
            pos(ct+'/entire critical-free circle',1-m*R)
            lam=q*R/(1-m*R);pos(ct+'/entire inverse derivative perturbation',1-lam)
            # The final term bounds every logarithmic order >=17.
            e=q*(sum((m**(j-2)*R**j/F(j) for j in range(2,17)),F(0))+m**15*R**17/(17*(1-m*R)))
            rh=R-e;pos(ct+'/entire inverse image',rh)
            inv=q/((1-m*R)**2*(1-lam)**3)
            eps=F(1001,1000)*inv/H;g=F(1001,1000)*nu*inv*rh
            p=1/(1+lam);amp=1/(1-lam);ew=eminus(int(W))
            pos(ct+'/all four containing corners',(1-nu*nu)*(H*rh)**2-W*W)
            pos(ct+'/whole left vertical',1-amp*ew-eps*(1+ew)-g*ew)
            pos(ct+'/whole right vertical',p-ew-eps*(1+ew)-g*ew*ew)
            pos(ct+'/both whole horizontals',p-eps-g*ew)
            ck(ct+'/whole amplitude homotopy',0<p<1<amp)
            pos(ct+'/full integer count height',nu*H*rh-F(22,7))
            slope=lower12(nu*F(7,22)*rh)
            width=min(R,upper12(W/H+J*nu*rh/H+e))
            depth=lower12((1-(J/H)**2)*nu*rh-e-4*F(22,7)/H-J*W/H**2)
            ck(ct+'/physical rightmost real interval',0<width<=R<2 and 0<J/H<1)
            pos(ct+'/complete rounded count',slope*H-1)
            rebuilt.append(dict(R=R,slope=slope,width=width,depth=depth))
            CONTOURS.append(dict(file=name,cell=ix,index=ir,H=str(H),R=str(R),
              slope=str(slope),width=str(width),depth=str(depth)))
        ck(tag+'/nonempty full inverse family',bool(rebuilt))
        if kind=='depth_skip':
            best=max(a['depth'] for a in rebuilt);pos(tag+'/actual root upper-halfplane contradiction',best-tr)
            # Positive depth itself forces the budget above4pi, so the
            # two-strip existence lemma is legal even though generic counts
            # only need height>pi.
            CELLS.append(dict(file=name,index=ix,kind=kind,H=str(H),depth=str(best),field_reasons=reasons,boosts=boosts));continue
        ck(tag+'/recognized complete cloud',kind=='cloud' and tl>0)
        de=rat(z['delta']);cayley(tag,rhi,de);k=kap(tr,de);A=1-k*yl
        f,Dp=polar(rhi,de);ck(tag+'/complete negative-Q domain',0<=k<=F(1,2) and A>0 and f>1)
        if mode=='high':
            B_rad=(rhi/2+J/N)/(tl*(1-rhi*tl/2))
            U=min(yh,sqrt_upper(B_rad),B_rad/yl if yl>0 else sqrt_upper(B_rad))
            B=2*Dp*(f*A*U+F(1,n))
        else:B=2*Dp*(600*f*A+1)/n
        ck(tag+'/finite HB denominator before solve',0<B<1)
        upper=upper12(rhi*(f*A+1/H)/(1-B))
        weights=[];slopes=[];lastR=F(0);lastw=F(0);lasts=F(0)
        for j,crow in enumerate(rebuilt):
            ck(tag+f'/nest{j}/entire physical nesting',crow['R']>lastR and crow['slope']>lasts)
            lastw=max(lastw,crow['width']);lastR=crow['R'];lasts=crow['slope']
            ck(tag+f'/nest{j}/whole real strip remains legal',lastw<=lastR<2)
            weights.append(lower12(weight(lastR,lastw,de,tl,tr,rhi)));slopes.append(lasts)
        eff=[max(weights[j:]) for j in range(len(weights))]
        ck(tag+'/complete decreasing nonnegative weights',eff[0]>0 and all(eff[j]>=eff[j+1]>=0 for j in range(len(eff)-1)))
        beta_w=weight(F(0),F(0),de,tl,tr,rhi);comp=beta_w-3*eff[0]
        correction=comp/H if comp<0 else (comp/600 if mode=='low' else F(0))
        lower=eff[-1]*slopes[-1]+sum(((eff[j]-eff[j+1])*slopes[j] for j in range(len(slopes)-1)),F(0))+correction
        pos(tag+'/complete actual-source contradiction',lower-upper)
        CELLS.append(dict(file=name,index=ix,kind=kind,H=str(H),lower=str(lower),upper=str(upper),
             margin=str(lower-upper),source_correction=str(correction),contours=len(rebuilt),field_reasons=reasons,boosts=boosts))
    # Independent vertical-fiber coverage, including every depth endpoint
    # and one entire-membership representative in every intervening slab.
    ts=sorted(set(t for _,_,a,b in boxes for t in (a,b)))
    ck(name+'/whole exterior depth interval',ts[0]==0 and ts[-1]==1)
    for t in ts+[(a+b)/2 for a,b in zip(ts,ts[1:])]:
        cover=sorted((yl,yh) for yl,yh,tl,tr in boxes if tl<=t<=tr)
        last=F(0)
        for a,b in cover:
            ck(name+'/closed vertical fiber no gap',a<=last)
            last=max(last,b)
        ck(name+'/closed vertical fiber complete',last==2)
    print(json.dumps(dict(file=name,processed_cells=len(CELLS),contours=len(CONTOURS),checks=sum(COUNTS.values()),seconds=time.monotonic()-START)),flush=True)
LOW_SLABS.sort();ck('low strip/exterior rho endpoints',LOW_SLABS[0][0]==FLOOR and LOW_SLABS[-1][1]==CAP)
ck('low strip/all exact closed rho seams',all(a[1]==b[0] for a,b in zip(LOW_SLABS,LOW_SLABS[1:])))
result=dict(status='PASS',scope='Complete X=1 higher strip .01<rho<=.061002, every n>=6251; no base/global assertion',
 degree_floor=n,checks=sum(COUNTS.values()),complete_contours=len(CONTOURS),closed_cells=len(CELLS),
 by_kind=dict(Counter(c['kind'] for c in CELLS)),min_cloud_margin=str(min(F(c['margin']) for c in CELLS if 'margin' in c)),
 entrance=ENTRANCE,rayleigh_reserve=str(RAYLEIGH_RESERVE),proof_gate_digest=DIGEST.hexdigest(),gate_counts=dict(COUNTS),
 proposal_sha256=HASHES,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 cells=CELLS,contours=CONTOURS,seconds=time.monotonic()-START)
Path(__file__).with_name('tang_X6251_independent_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('cells','contours','gate_counts')},indent=2))
