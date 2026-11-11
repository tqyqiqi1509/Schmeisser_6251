"""Independent T6251 full-strip validation.

No author function, bound, receipt, or polynomial certificate is imported.
The proposal contributes only exact rational choices of cells and contours.
Whole rational functions are certified by exact Sturm root exclusion and an
endpoint sign, rather than the author's Bernstein subdivision. Scalar and
all complete contour gates use exact Fractions; exp(-W) is bounded by the
reciprocal of an independently generated positive 80-term exponential sum.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from functools import lru_cache
import sympy as sy
import json,hashlib,time,sys
sys.set_int_max_str_digits(0)
T0=time.monotonic();ROOT=Path(__file__).resolve().parents[1]
DATA_DIR=ROOT/'agents/global_component_refinement/Origin6251'
N=6250;n=N+1;rho=F(30501,500000);rho_min=F(1,100)
X=sy.Symbol('X',real=True);checks=[];polycert=[];contours=[];cells=[]
def S(v):
    if isinstance(v,F):return sy.Rational(v.numerator,v.denominator)
    return sy.sympify(v)
def FQ(v):
    v=sy.cancel(v);assert v.is_Rational
    return F(int(v.p),int(v.q))
def ck(name,test,value=None):
    ok=bool(test);checks.append({'name':name,'pass':ok,'value':None if value is None else str(value)})
    if not ok:raise AssertionError((name,value))
def pos(name,value):ck(name,value>0,value)
def below(name,value,upper):pos(name,upper-value)
def rational_positive(name,expr,a,b):
    """Both polynomial numerator and denominator have no real zero on [a,b]."""
    expr=sy.cancel(expr);num,den=expr.as_numer_denom()
    assert F(a)<F(b)
    vals=[]
    for kind,p in [('numerator',num),('denominator',den)]:
        pp=sy.Poly(p,X,domain=sy.QQ);va=pp.eval(S(a));vb=pp.eval(S(b))
        ck(name+'/'+kind+' endpoint nonzero',va!=0 and vb!=0)
        roots=pp.count_roots(S(a),S(b)) if pp.degree()>0 else 0
        ck(name+'/'+kind+' entire Sturm zero exclusion',roots==0,roots)
        vals.append((va,vb))
    ck(name+'/same endpoint sign',vals[0][0]*vals[1][0]>0)
    polycert.append(dict(name=name,interval=[str(a),str(b)],expression=str(expr),
                         numerator=str(num),denominator=str(den),
                         numerator_root_count=0,denominator_root_count=0,
                         endpoint_sign=1 if vals[0][0]>0 else -1))
    return expr
def active_kappa(delta):
    d=S(delta)
    # Derived first from the full unbounded-label extremum, independently
    # of the author's factored form.
    return sy.cancel(X*X/2-X*(1-X*d-d*d)**2/(4*(1+d*d)*(X/2+d)))
def source(delta,x,r):
    delta=S(delta);x=S(x);r=S(r)
    b=(1-S(rho)*(1+x*x)/2)/(1+delta*delta)
    polar=1+(S(rho)*(x+delta)/2)**2
    invf=1-S(rho)**2*(x+delta)**2/2
    return sy.cancel((n*r*b/(n*r+2*polar*b)-r)*invf)

# Independent algebra of the general critical comparison and direct HB
# source, including the inactive positive-part case.
B,t,d,kap,x=sy.symbols('B t d kap x',real=True)
Ac=B-d;b=x/2+d;L=1+d*d
delta_clear=(x*Ac-kap)*(B*B+t*t)-x*(Ac-x/2)*(Ac*Ac+(1-t)**2)
rhs=(x*x/2-kap)*B*B-x*(B-b)*(L-2*d*B)+(x*x/2-kap)*t*t+2*x*(Ac-x/2)*t
ck('complete unbounded-label cleared identity',sy.expand(delta_clear-rhs)==0)
w=sy.symbols('w',real=True)
f=(1-b*w)*(L*w-2*d)
vertex=(L+2*d*b)/(2*b*L)
ck('full active extremum identity',sy.cancel(f.subs(w,vertex)-(L-2*d*b)**2/(4*b*L))==0)
ck('inactive endpoint zero',sy.expand(f.subs(w,1/b))==0)
kk=x*x/2-x*(1-x*d-d*d)**2/(4*(1+d*d)*(x/2+d))
fact=x*(x-(1-d)**2)*(x+(1+d)**2)/(2*(1+d*d)*(x+2*d))
ck('complete active factor identity',sy.cancel(kk-fact)==0)
rr,bb,DD,nn=sy.symbols('r b D n',positive=True)
src=nn*rr*bb/(nn*rr+2*DD*bb)-rr
ck('direct source full concavity identity',sy.cancel(sy.diff(src,rr,2)+4*nn**2*DD*bb**2/(nn*rr+2*DD*bb)**3)==0)
ck('direct source monotone degree',sy.cancel(sy.diff(src,nn)-2*DD*rr*bb**2/(nn*rr+2*DD*bb)**2)==0)

# Complete geometric entrance, with X retained.
pos('whole crescent height comparison',F(2,125)-rho*(F(1,4)+F(1,15625)))
pos('positive squared height comparator',1-rho/F(125))
G=1+sy.Rational(117,500)*X**2
Y=3*(1-sy.Rational(11,100)*X**2)
rational_positive('whole same-X y envelope after exact zero factor',
  sy.cancel(((Y/3)**2*G-1)/X**2),F(0),F(49,100))
pos('positive envelope before squaring',1-F(11,100)*F(49,100)**2)

half=F(1,2);kh=active_kappa(half)
half_endpoints=[FQ(source(half,1,r)) for r in (rho_min,rho)]
Kh=min(half_endpoints);pos('direct half source positive',Kh)
rational_positive('half kappa divided ratio derivative',sy.diff(sy.cancel(kh/X),X),F(0),F(1))
below('all direct-source inverse Cayley tests',(rho*F(171,100)/2)**2,F(1,200))
pos('first D25 small-X entrance',Kh/(rho*F(2,5))-25)
rational_positive('first D25 full large-X quadratic maximum',4*S(F(227,200))*kh-(X+kh)**2,F(2,5),F(1))
pos('first D25 independent large-X entrance',Kh*N/(2*F(227,200))-25**2)

below('small-X A(1+y) divided gap',F(1,4)*4,F(21,20))
half_ratio=sy.cancel(kh/X)
below('full y derivative sign on half small-X',7*FQ(half_ratio.subs(X,S(F(41,100)))),1)
rational_positive('whole same-gap small-X upper',S(F(21,20))*G-(X-kh*Y)*(1+Y),F(1,4),F(41,100))
rational_positive('whole same-gap large-X upper',4*S(F(21,20))*G*kh-(X+kh)**2,F(41,100),F(1))
pos('ordered first D45 entrance',Kh*N/(2*F(21,20))-45**2)

# Universal matrix coefficients and the first lawful full-angle exclusion.
pl=F(157,50);pl2=pl*pl
below('matrix K0 rational pi lower bound',9*pl2/(3*pl2-8),F(21,5))
below('matrix A0 rational pi lower bound',2+24*pl/(3*pl2-8),F(11,2))
below('matrix gamma rational pi lower bound',(40-3*pl2)/(2*(3*pl2-8)),F(1,4))
pos('matrix45 complete cosine upper',F(18,25)**2-(1-F(7,10)**2))
below('matrix45 entire forbidden-angle sector',F(21,5)/45**2+F(11,2)*F(18,25)/45+(1-F(7,10)**2)/4-F(7,10)**2/2,0)

# Reconstruct the two new whole bivariate source bounds with an independent
# one-dimensional method. For each X interval, the complete quadratic in y
# is either monotone, concave, or has positive completed-square minimum.
# All sign claims are exact Sturm exclusions, not author tensor controls.
def has_positive(expr,a,b):
    expr=sy.cancel(expr);num,den=expr.as_numer_denom()
    if expr==0:return False
    for p in (num,den):
        pp=sy.Poly(p,X,domain=sy.QQ)
        if pp.eval(S(a))==0 or pp.eval(S(b))==0:return False
        if pp.degree()>0 and pp.count_roots(S(a),S(b))!=0:return False
    return num.subs(X,S(a))*den.subs(X,S(a))>0

def quadratic_rectangle(name,P,yvar):
    pol=sy.Poly(P,yvar);a=sy.cancel(pol.coeff_monomial(yvar**2));b=sy.cancel(pol.coeff_monomial(yvar));c=sy.cancel(pol.coeff_monomial(1))
    ck(name+'/genuine complete quadratic',pol.degree()==2)
    endpoint=c+3*b+9*a
    regions=[];pending=[(F(0),F(1),0)]
    while pending:
        lo,hi,depth=pending.pop()
        opts=[('increasing',sy.cancel(b/X)),('decreasing',sy.cancel(-(6*a+b)/X)),('concave',sy.cancel(-a/X)),('positive_square',sy.cancel((4*a*c-b*b)/X))]
        chosen=None
        for kind,expr in opts:
            if has_positive(expr,lo,hi) and (kind!='positive_square' or has_positive(sy.cancel(a/X),lo,hi)):
                chosen=(kind,expr);break
        if chosen is None:
            ck(name+'/bounded independent subdivision',depth<14)
            mid=(lo+hi)/2;pending.extend([(lo,mid,depth+1),(mid,hi,depth+1)]);continue
        kind,expr=chosen
        rational_positive(name+'/'+kind+' sign',expr,lo,hi)
        rational_positive(name+'/y0',c,lo,hi)
        rational_positive(name+'/y3',endpoint,lo,hi)
        if kind=='positive_square':rational_positive(name+'/positive quadratic coefficient',sy.cancel(a/X),lo,hi)
        regions.append((lo,hi,kind))
    regions.sort();cur=F(0)
    for lo,hi,kind in regions:ck(name+'/exact X seam',lo==cur);cur=hi
    ck(name+'/complete X endpoint',cur==1)
    # At X=0 the quadratic coefficients vanish; positive c(0) handles
    # that closed boundary without dividing by a vanishing quantity.
    ck(name+'/Xzero full y edge',a.subs(X,0)==0 and b.subs(X,0)==0 and c.subs(X,0)>0)
    return [{'X':[str(lo),str(hi)],'y':['0','3'],'method':kind} for lo,hi,kind in regions]

Yv=sy.Symbol('y',real=True);de=F(3,5);k06=active_kappa(de)
K06=source(de,X,rho);K06low=source(de,X,rho_min)
rational_positive('delta06 source endpoint order',K06low-K06,F(0),F(1))
rational_positive('delta06 source positive',K06,F(0),F(1))
beta06=F(7,10)+F(3,10)*rho;I06=1-rho_min/2
A06=X-k06*Yv;R06=S(beta06)+S(I06)*Yv
TENSOR=[]
for label,hgoal,gap in [('same-X D51',51,G),('physical h45',45,sy.Integer(1))]:
    expr=sy.cancel(N*K06*gap-2*hgoal**2*A06*R06)
    TENSOR.append({'name':label,'expression':str(expr),'regions':quadratic_rectangle(label,expr,Yv)})
print(json.dumps({'event':'independent_entry_feedback','quadratic_rectangles':len(TENSOR),'X_regions':sum(len(t['regions']) for t in TENSOR),'seconds':time.monotonic()-T0}),flush=True)
pos('matrix51 complete cosine upper',F(3,4)**2-F(5,9))
below('matrix51 entire forbidden-angle sector',F(21,5)/51**2+F(11,2)*F(3,4)/51+F(5,36)-F(2,9),0)
below('lawful actual companion mu rectangle',F(1,51),F(1,50))

delta=F(71,100);kd=active_kappa(delta);xc=(1-delta**2)/delta
rational_positive('whole .71 kappa ratio derivative',sy.diff(sy.cancel(kd/X),X),F(0),xc)
K=source(delta,X,rho);Klow=source(delta,X,rho_min)
rational_positive('entire source endpoint choice',Klow-K,F(0),F(1))
rational_positive('same-X entire source positivity',K,F(0),F(1))
Kstar=min(FQ(source(delta,1,r)) for r in (rho_min,rho))
pos('global .71 direct source positive',Kstar)
rational_positive('entire same-X lower source with endpoint factor',sy.cancel((K-S(Kstar))/(1-X)),F(0),F(1))
split=F(49,100);gam=F(2,3);qcap=F(57,1000)
beta=S(gam*(1-rho_min))+S(rho_min)*X;In=1-S(rho_min)*(1+X*X)/2
ratio=sy.cancel(kd/X);rmax=FQ(ratio.subs(X,S(split)))
Imin=1-rho_min*(1+split**2)/2;Imax=1-rho_min/2
bmin=gam*(1-rho_min);bmax=bmin+rho_min*split
pos('low-X radial coefficient decreases in rho',gam-split)
pos('common concave derivative lower at y0',2*bmin)
pos('common concave derivative lower at y3',2*bmin+9*(Imin-rmax*bmax)-36*rmax*Imax)
pos('concavity of common lower polynomial',rmax)
rational_positive('full same-X q057 entrance',N*K-S(2/qcap**2)*Y**2*(X-kd*Y)*(beta+In*Y),F(0),split)
rational_positive('whole active large-X kappa lower quadratic',kd-S(F(12,25))*X**2,split,xc)
pos('inactive large-X quadratic lower',F(1,2)-F(12,25))
pos('same-gap high-X y bound',F(73,25)**2*(1+F(117,500)*split**2)-9)
high_q2=2*F(25,48)*F(73,25)*(gam+(1-gam)*rho+Imax*F(73,25))/(N*Kstar)
below('large-X full q057 rectangle',high_q2,qcap*qcap)
below('q057 belongs to independently checked actual Q06',qcap,F(3,50))

# Complete general terminal Ay bound, independent of degree and ceilings.
m=(1-delta)**2
pos('uniform small-X m',m)
ck('uniform large-A discriminant zero',(1-delta**2)**2-m*(4*delta+m)==0)
below('uniform small-A derivative discriminant',16*m*m+12*(4*m*delta-1),0)
ck('uniform small-A endpoint identity',2*m*(m+delta)**2-m*(m*m+1)==m*delta**4)
pos('uniform small-A endpoint positive',m*delta**4)
below('small-X Ay terminal',3*m,F(11,20))
below('intermediate Ay full y derivative',6*FQ(ratio.subs(X,S(F(35,100)))),1)
rational_positive('intermediate entire Ay bound',S(F(11,20))-3*(X-3*kd),m,F(35,100))
rational_positive('active whole Ay quadratic',kd-S(F(5,11))*X*X,F(35,100),xc)
pos('inactive whole Ay quadratic',F(1,2)-F(5,11))

# Recompute the new same-actual-pair consequence from all 40 exact center
# terms. The two actual circles, including full tails, are a stated input;
# no signed/curvature theorem is used anywhere in this T proof.
from flint import arb,acb,ctx
ctx.prec=384
center_file=ROOT/'agents/tang_weighted_moment_center_W8.json'
center=json.loads(center_file.read_text());ck('frozen exact W8 center hash',hashlib.sha256(center_file.read_bytes()).hexdigest()=='2a0f8a1d121649c20e1a12e0fb045a7b0914b8794d8a0a7ce93e51e697005fe2')
def ar(v):
    v=F(v);return arb(v.numerator)/v.denominator
ss=acb(0,2*arb.pi());Zsum=arb(0);Even=arb(0);Zdelta=arb(0)
for row in center['coefficients']:
    a,b,c=row['monomial'];factor=qcap**(a+b+c)*F(1,51)**(b+2*c)
    co=list(map(F,row['coefficients']));v=acb(0);ev=acb(0)
    for j,cc in enumerate(co):
        v+=acb(ar(cc))*ss**j
        if j%2==0:ev+=acb(ar(cc))*ss**j
    bound=abs(v)*ar(factor);Zsum+=bound;Even+=abs(ev)*ar(factor)
    if sum(row['monomial']):Zdelta+=bound
Rnew=64*qcap*qcap+95000*qcap*F(1,51)**3
ck('new actual pair entire displacement',Zsum+ar(Rnew)<ar(F(8917,1000)))
ck('new same-parameter pair even mean',Even+ar(Rnew)<ar(F(2109,1000)))
ck('new whole smaller-rectangle homothetic phase',Zdelta+ar(Rnew)<2*arb.pi()*ar(qcap*F(37,5)))
PAIR_CONSEQUENCE={'q':str(qcap),'mu':'1/51','center_bound':str(Zsum),'even_center_bound':str(Even),'actual_radius':str(Rnew),'w_upper':'8917/1000','mean_upper':'2109/1000','phase':'37/5'}

# Complete ordered ceilings. First h>=105 enters D100 independently.
# Second 90<=h<105 enters D70 independently, then the smaller ratio branch.
ck('h105 independent D100 entrance',105>100)
ck('h90 independent D70 entrance',90>70>25)
P70=lambda z,D:3*z*z-(D-F(91,10))*z+2*D
below('D90 excludes prior large branch',P70(F(8,3),90),0)
below('D90 strict ratio49/20',P70(F(49,20),90),0)
below('D90 ratio49/20 decreases over entire D ray',2-F(49,20),0)
below('D90 prior8/3 decreases over entire D ray',2-F(8,3),0)
below('h90 old q05 companion entrance',F(49,20)/90,F(1,20))
below('h90 old mu02 companion entrance',F(1,90),F(1,50))
@lru_cache(None)
def exp_upper(W):
    assert W>0
    return 1/sum((F(W)**j/factorial(j) for j in range(81)),F(0))
def floor12(v):return F((v*10**12).__floor__(),10**12)
def ceil12(v):return -floor12(-v)
from math import isqrt
def sqrtup(v):
    k=isqrt(v.numerator*10**32//v.denominator)+1;u=F(k,10**16)
    ck('complete squared radical enclosure',u*u>v and (u-F(1,10**16))**2<=v)
    return u
pos('all inverse Cayley domains',F(1,200)-rho*rho)
PROPOSAL_HASHES={};STAGES=[]
for stage,Hfloor,ratio_cap,phase in [('H105',105,F(9,4),F(23,5)),('H90',90,F(49,20),F(21,4))]:
    DATA=DATA_DIR/('T6251_'+stage+'_PROPOSAL.json');proposal=json.loads(DATA.read_text())
    PROPOSAL_HASHES[stage]=hashlib.sha256(DATA.read_bytes()).hexdigest()
    ck(stage+'/proposal degree',proposal['degree_floor']==n)
    ck(stage+'/proposal no failed cells',not proposal['failed'])
    boxes=[];count0=len(contours);cell0=len(cells)
    for iz,z in enumerate(proposal['cells']):
        tag=stage+f'/cell{iz}'
        lo,hi,yl,yh,H,de,A,U,Xcap,J=[F(z[key]) for key in ('rho_lo','rho_hi','yl','yh','H','delta','Abar','U','Xcap','J')]
        ck(tag+'/complete closed domain',rho_min<=lo<hi<=rho and 0<=yl<yh<=ratio_cap)
        ck(tag+'/same-y lawful phase',J==phase*yh)
        ck(tag+'/field lower without new ceiling',H>=Hfloor and (H==Hfloor or H<=(N*lo/2-J)/(hi+yh)))
        ck(tag+'/positive actual phase direction',0<J<H)
        ck(tag+'/same-X crescent bound',0<Xcap<=1 and (Xcap==1 or yl**2*(1+F(117,500)*Xcap**2)>=ratio_cap**2))
        ck(tag+'/strict inverse Cayley test',F(1,200)<=de<=1 and (hi*(1+de)/2)**2<de)
        transition=(1-de*de)/de;ka=active_kappa(de)
        for yend in (yl,yh):
            if transition>0:rational_positive(tag+'/entire active X y='+str(yend),S(A)-X+ka*S(yend),F(0),min(Xcap,transition))
            if transition<Xcap:rational_positive(tag+'/entire inactive X y='+str(yend),S(A)-X+X*X*S(yend)/2,max(F(0),transition),Xcap)
        ff=1/(1-hi*hi*(1+de)**2/2);dp=1+(hi*(1+de)/2)**2
        if z['kind']=='pair':
            ck(tag+'/prior h105 ceiling only in second stage',stage=='H90')
            hb=2*dp*(ff*A*105+1)/n;ck(tag+'/pair HB denominator',0<hb<1 and A>0)
            upper=hi*(ff*A*105+1)/(1-hb)
            R=F(17,180);mean=F(19,900)
            pos(tag+'/full pair numerator lower',(1-hi)*(1-R)-hi*R*R/2-F(3,4))
            below(tag+'/full pair numerator upper',1+R,F(6,5))
            for vv in (F(3,4),F(6,5)):pos(tag+'/whole square-root chord '+str(vv),vv-(vv/2+F(49,100))**2)
            Is=2*(1-hi)*(1-mean)-hi*R*R;Ds=2*(1+de*de)+4*sqrtup(1+de*de)*mean+2*R*R
            ck(tag+'/pair source parameter domain',0<hi<1 and R>0 and 0<mean<1 and Is>0 and Ds>0)
            lower=(1-hi)/(1+de*de)+(Is/2+F(49,50))**2/Ds
            pos(tag+'/full actual pair contradiction',lower-upper)
            cells.append(dict(stage=stage,index=iz,kind='pair',rho=[str(lo),str(hi)],y=[str(yl),str(yh)],contours=0,lower=str(lower),upper=str(upper),margin=str(lower-upper)))
            boxes.append((lo,hi,yl,yh));continue
        ck(tag+'/explicit cloud kind',z['kind']=='cloud')
        br=(hi/2+J/N)/(1-hi/2)
        ck(tag+'/complete radial upper h/n',U>0 and (U>=yh or U*U>=br or (yl>0 and U*yl>=br)))
        hb=2*dp*(ff*A*U+F(1,n));ck(tag+'/finite HB denominator before solve',0<hb<1 and A>0)
        upper=hi*(ff*A+1/H)/(1-hb)
        slopes=[];weights=[];last_R=F(0);last_width=F(0);last_slope=F(0)
        for ir,choice in enumerate(z['contours']):
            R,nu,W=F(choice[0]),F(choice[1]),int(choice[2]);label=tag+f'/contour{ir}'
            ck(label+'/physical disk and trim',0<R<1 and 0<nu<1 and W>0)
            invden=1-R;lam=yh*R/invden
            logtail=yh*(sum((R**j/F(j) for j in range(2,17)),F(0))+R**17/(17*invden));image=R-logtail
            ck(label+'/whole analytic inverse domain',0<lam<1 and image>0)
            C=yh/(invden**2*(1-lam)**3)
            eps=F(1001,1000)*C/H;gg=F(1001,1000)*nu*C*image
            p=1/(1+lam);amp=1/(1-lam);eW=exp_upper(W)
            pos(label+'/complete containing corners',(1-nu**2)*(H*image)**2-W**2)
            pos(label+'/entire left vertical',1-amp*eW-eps*(1+eW)-gg*eW)
            pos(label+'/entire right vertical',p-eW-eps*(1+eW)-gg*eW**2)
            pos(label+'/both entire horizontals',p-eps-gg*eW)
            ck(label+'/whole amplitude homotopy',0<p<1<amp)
            pos(label+'/complete finite-height count',nu*H*image-F(22,7))
            slope=floor12(nu*F(7,22)*image)
            width=max(last_width,min(R,ceil12(W/H+J*nu*image/H+logtail)))
            ck(label+'/entire nested disks and real strips',R>last_R and slope>last_slope and last_width<=width<=R)
            aa=(de+width)**2;pos(label+'/whole weight derivative root product',aa)
            endpoints=[(tt-hi*(1+tt*tt)/2)/(aa+tt*tt) for tt in (1-R,1+R)]
            weight=floor12(max(F(0),min(endpoints)));slopes.append(slope);weights.append(weight)
            contours.append(dict(stage=stage,cell=iz,index=ir,radius=str(R),image_radius=str(image),H=str(H),lambda_=str(lam),width=str(width),slope=str(slope),weight=str(weight)))
            last_R,last_slope,last_width=R,slope,width
        ck(tag+'/nonempty complete cloud',bool(slopes))
        effective=[max(weights[j:]) for j in range(len(weights))]
        ck(tag+'/nonnegative decreasing full weights',effective[0]>0 and all(effective[j]>=effective[j+1]>=0 for j in range(len(effective)-1)))
        bsource=(1-hi)/(1+de*de);compensation=min(F(0),bsource-3*effective[0])/H
        lower=effective[-1]*slopes[-1]+sum(((effective[j]-effective[j+1])*slopes[j] for j in range(len(slopes)-1)),F(0))+compensation
        pos(tag+'/complete finite cloud contradiction',lower-upper)
        cells.append(dict(stage=stage,index=iz,kind='cloud',rho=[str(lo),str(hi)],y=[str(yl),str(yh)],H=str(H),contours=len(slopes),lower=str(lower),upper=str(upper),margin=str(lower-upper)))
        boxes.append((lo,hi,yl,yh))
    # Vertical fibers, independent of the author's slab traversal. Every
    # y seam is checked, then its active rho intervals tile the full strip.
    cuts=sorted(set([F(0),ratio_cap]+[c for z in boxes for c in z[2:]]))
    for y0,y1 in zip(cuts,cuts[1:]):
        probe=(y0+y1)/2;active=sorted((a,b) for a,b,l,u in boxes if l<probe<u);end=rho_min
        for a,b in active:ck(stage+'/independent rho fiber seam',a==end);end=b
        ck(stage+'/independent rho fiber complete endpoint',end==rho)
    ck(stage+'/exact closed y endpoint',cuts[0]==0 and cuts[-1]==ratio_cap)
    STAGES.append({'stage':stage,'cells':len(cells)-cell0,'contours':len(contours)-count0})
    print(json.dumps({'event':'independent_complete_ceiling','stage':stage,'cells':len(cells)-cell0,'contours':len(contours)-count0,'seconds':time.monotonic()-T0}),flush=True)

# Two continuous actual-field terminal intervals, with their own lawful
# actual pair/phase theorem and exact common seam h=70.
TERMINALS=[];delta=F(71,100)
for label,hlo,hhi,Z,M,phase in [('low',45,70,F(8917,1000),F(2109,1000),F(37,5)),('high',70,90,F(17,2),F(19,10),F(7))]:
    if label=='high':
        below('h70 old companion q entry',F(3,70),F(1,20));below('h70 old companion mu entry',F(1,70),F(1,50))
    polar=1+(rho*(1+delta)/2)**2;ff=1/(1-rho*rho*(1+delta)**2/2)
    hb=2*polar*(ff*hhi+1)/n;ck(label+'/terminal finite HB denominator',0<hb<1)
    pos(label+'/terminal whole degree denominator',N-2*hhi)
    upper=(ff*2*F(11,20)*hhi*(hhi+phase)/(N-2*hhi)+rho)/(1-hb)
    R=Z/hlo;mean=M/hlo
    pos(label+'/entire pair numerator lower',(1-rho)*(1-R)-rho*R*R/2-F(3,4))
    below(label+'/entire pair numerator upper',1+R,F(6,5))
    for vv in (F(3,4),F(6,5)):pos(label+'/whole square-root chord '+str(vv),vv-(vv/2+F(49,100))**2)
    pos(label+'/paired denominator norm bound',F(123,100)**2-(1+delta*delta))
    Is=2*(1-rho)*(1-mean)-rho*R*R;Ds=2*(1+delta*delta)+4*F(123,100)*mean+2*R*R
    ck(label+'/pair source parameter domain',0<rho<1 and R>0 and 0<mean<1 and Is>0 and Ds>0)
    lower=(1-rho)/(1+delta*delta)+(Is/2+F(49,50))**2/Ds
    pos(label+'/complete paired-source terminal',lower-upper)
    TERMINALS.append({'name':label,'h':[hlo,hhi],'pair_lower':str(lower),'HB_upper':str(upper),'margin':str(lower-upper)})
ck('two full field intervals have exact seam',TERMINALS[0]['h'][1]==TERMINALS[1]['h'][0])
result=dict(status='PASS',scope='Complete T=1 higher strip .01<rho<=.061002, all n>=6251; no complementary X/base claim; independently accepted Q06 actual-circle input only',checks=len(checks),whole_sturm_certificates=len(polycert),whole_quadratic_rectangles=len(TENSOR),quadratic_X_regions=sum(len(t['regions']) for t in TENSOR),closed_cells=len(cells),complete_contours=len(contours),pair_cells=sum(z['kind']=='pair' for z in cells),stages=STAGES,min_cell_margin=str(min(F(row['margin']) for row in cells)),terminals=TERMINALS,half_source_K=str(Kh),delta071_source_K=str(Kstar),high_X_q2_upper=str(high_q2),pair_consequence=PAIR_CONSEQUENCE,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),proposal_hashes=PROPOSAL_HASHES,seconds=time.monotonic()-T0,predicates=checks,polynomials=polycert,quadratic_certificates=TENSOR,cells=cells,contours=contours)
Path(__file__).with_name('tang_T6251_independent_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['predicates','polynomials','quadratic_certificates','cells','contours']},indent=2))
