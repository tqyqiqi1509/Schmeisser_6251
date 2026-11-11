"""Independent arithmetic review of BASE010 at n>=6251, retaining the same depth in every cloud.

No author proof/checker module is imported. Rational witnesses are rebuilt
from the manuscripts. Transcendentals use Arb at256bits, not an author's
exponential Taylor routine. Every inverse contour and companion rectangle
is retained, together with complete depth intervals and entrance order.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib,json,time
from flint import arb,acb,ctx
import sympy as sp
ctx.prec=384
ROOT=Path(__file__).resolve().parents[2]
checks=[];contours=[];notes={};tic=time.time()
def A(x):
 if isinstance(x,F):return arb(x.numerator)/x.denominator
 return arb(x)
def ck(name,test,margin=None):
 ok=bool(test);checks.append({'name':name,'pass':ok,'margin':None if margin is None else str(margin)})
 if not ok:raise AssertionError(name+': '+str(margin))
def gt(name,a,b):ck(name,a>b,a-b)
def exp_lt(name,x,y):gt(name,A(y),A(x).exp())
def exprow(name,x,y):gt(name,A(x).exp(),A(y))
def profile(name,q,m,R,H,nu,W,eps,g,p,b,e,quarter=False):
 q,m,R,H,nu,W,eps,g,p,b,e=map(F,(q,m,R,H,nu,W,eps,g,p,b,e))
 z=m*R; lam=q*R/(1-z)
 err=q*R**2*(F(1,2)+z/(3*(1-z))); rr=R-err
 B=q/((1-z)**2*(1-lam)**3)
 tests={'distance':1-z,'inverse':1-lam,'radius':rr,
 'corners':(H*rr)**2*(1-nu**2)-W**2,'derivative':eps-B/H,
 'vertical':g-nu*B*rr,'left':1-b*e-eps*(1+e)-g*e,
 'right':p-e-eps*(1+e)-g*e**2,
 'horizontal':(p/(1+p+b) if quarter else p)-eps-g*e,
 'height':nu*H*rr-(F(11,7) if quarter else F(22,7))}
 for key,val in tests.items():gt(name+'/'+key,val,0)
 ck(name+'/amplitude',0<p<1<b and p<=1/(1+lam) and b>=1/(1-lam))
 exprow(name+'/full_exponential',W,1/e)
 out=dict(name=name,q=q,m=m,R=R,H=H,nu=nu,W=W,epsilon=eps,g=g,p=p,b=b,e=e,quarter=quarter,rr=rr,err=err,B=B,density=nu*F(7,22)*rr,loss=F(3,2) if quarter else F(2))
 contours.append(out);return out

def weight(rho,R,delta,lo,hi):
 ts=(max(F(0),lo-R),hi+R)
 return max(F(0),min((t-rho*(1+t*t)/2)/(t*t+(R+delta)**2) for t in ts))
def sourcecloud(rho,row,delta,lo,hi):
 a=weight(rho,row['R'],delta,lo,hi);b=weight(rho,F(0),delta,lo,hi)
 return a*row['density']+min(F(0),b-3*a)/row['H']
def qcenter(q):return F(63,10)+21*q+133*q**2+1030*q**3+9000*q**4
def pathmajor(q):return sum(c*q**j for j,c in enumerate([0,F(63,10),45,400,3400,86000,320000,1700000,9100000,39000000]))
def companion(name,q,mu,Z0,M0,E0,P0,Z1,L,R,share,coeff=None):
 q,mu,Z0,M0,E0,P0,Z1,L,R,share=map(F,(q,mu,Z0,M0,E0,P0,Z1,L,R,share));nu=q*mu
 gt(name+'/model_domain',F(1,20),q)
 gt(name+'/center',Z0,qcenter(q));gt(name+'/path',M0,pathmajor(q))
 gt(name+'/critical_free',1,mu*Z1);gt(name+'/containment',Z1,Z0+R)
 expo=M0+nu*Z0**3/(3*(1-mu*Z0))
 exp_lt(name+'/full_path_exp',expo,E0)
 gt(name+'/primitive_tail',P0,Z0**4*E0/(12*(1-mu*Z0)))
 gt(name+'/log_derivative',L,1+q*Z1+nu*Z1**2/(1-mu*Z1))
 nonlinear=(A(L*R).exp()-1-A(L*R))/A(L)
 gt(name+'/whole_circle_nonlinear',A(share*R),nonlinear)
 center=E0*(132200*q**5+P0*nu)
 if coeff:
  aq,bmu=map(F,coeff)
  gt(name+'/coefficient_q5_reserve',(1-share)*aq,132200*E0)
  gt(name+'/coefficient_qmu_reserve',(1-share)*bmu,P0*E0)
  ck(name+'/variable_radius_endpoint',R==aq*q**5+bmu*nu)
 else:gt(name+'/center_residual',(1-share)*R,center)
 gt(name+'/strict_Rouche',A(R),A(center)+nonlinear)
 correction=qcenter(q)-F(63,10)
 gt(name+'/source_avoidance',F(6),correction+R)
 gt(name+'/pair_separation',F(12),2*(correction+R))
 notes[name]={'q':q,'mu':mu,'center_upper':center,'radius':R,'nonlinear_upper':str(nonlinear)}

# Independent reconstruction of the inherited full path majorants at the actual pi.
s=acb(0,2*arb.pi())
z=[s,s-s*s/2,3*s-2*s*s+s**3/2,15*s-21*s*s/2+23*s**3/6-5*s**4/8,
   105*s-72*s*s+88*s**3/3-22*s**4/3+7*s**5/8]
for j,cap in enumerate([21,133,1030,9000],1):gt('center_coefficient_'+str(j),A(cap),abs(z[j]))
path_leaves=[]
for j,cap in enumerate([F(63,10),45,400,3400,86000,320000,1700000,9100000,39000000],1):
 zz=sum((z[a]*z[j-1-a] for a in range(5) if 0<=j-1-a<5),acb(0))
 pending=[(F(0),F(1),0)];leaves=0
 while pending:
  lo,hi,dep=pending.pop()
  t=arb(A((lo+hi)/2),A((hi-lo)/2))
  val=t*(z[j] if j<=4 else acb(0))+t*t*zz/2
  if A(cap)>abs(val):leaves+=1;continue
  ck('path_subdivision_finite_'+str(j),dep<25)
  mid=(lo+hi)/2;pending.extend([(lo,mid,dep+1),(mid,hi,dep+1)])
 path_leaves.append(leaves)
notes['whole_path_coefficient_cover_leaves']=path_leaves

# Independently sum the entire model majorant at q=1/20. The seven
# residual coefficient bounds are the explicitly retained historical input.
mc=[F(0),F(63,10),F(45),F(400),F(3400),F(86000),F(320000),F(1700000),F(9100000),F(39000000)]
zc=list(map(F,[F(63,10),21,133,1030,9000]))
ecoef=[F(1)]
for j in range(1,12):
 ecoef.append(sum((k*mc[k]*ecoef[j-k] for k in range(1,min(j,9)+1)),F(0))/j)
pc=[sum((zc[k]*ecoef[j-k] for k in range(min(j,4)+1)),F(0)) for j in range(12)]
qtest=F(1,20)
P11=sum((v*qtest**j for j,v in enumerate(pc)),F(0))
finite=sum((F(v)*qtest**j for j,v in enumerate([82000,350000,1900000,13000000,120000000,3100000000,23000000000],5)),F(0))
model=(A(finite)+A(qcenter(qtest))*A(pathmajor(qtest)).exp()-A(P11))/A(qtest**5)
gt('entire_Gaussian_model_tail_at_full_q_domain',A(132200),model)
notes['Gaussian_entire_tail_ratio_upper']=str(model)

# BASE011: all55 contours rebuilt from the manuscript, with no removed rows.
n=6251;N=6250;r=F(1,100);c=r/2;f=F(1001,1000);fx=F(10001,10000);D=1+9*c*c
x1=profile('B11_Xfirst',2,2,F(3,20),600,F(99,100),4,F(1,25),3,F(7,10),F(7,4),F(1,50))
xd=profile('B11_Xstrip',2,2,F(13,100),600,F(99,100),4,F(1,25),3,F(7,10),F(7,4),F(1,50))
xmid={}
for H in [450,200]:
 xmid[H]=[]
 for j in range(6,31):
  W,g,e=(F(2),F(1,5),F(1,7)) if j<10 else (F(4),F(3),F(1,50))
  xmid[H].append(profile('B11_Xmid_'+str(H)+'_'+str(j),2,2,F(j,200),H,(F(23,25) if H==450 else F(9,10)),W,F(11,100),g,F(7,10),F(7,4),e,True))
x78=profile('B11_X78',2,2,F(3,20),78,F(97,100),F(11,5),F(3,10),3,F(7,10),F(7,4),F(1,9))
t150=profile('B11_T150',F(25,11),F(100,99),F(3,20),150,F(9,10),4,F(1,6),3,F(7,10),F(7,4),F(1,50))
t103=profile('B11_T103',F(25,11),F(100,99),F(3,20),103,F(9,10),4,F(1,6),3,F(7,10),F(7,4),F(1,50))
ck('B11_exact_55_rows',len(contours)==55)
for tag,ratio,fac in [('moving',F(9,2),f),('X',F(10201,20000),fx),('T',F(169,200),fx)]:gt('B11_polar_'+tag,fac,1/(1-ratio*r*r))
ck('B11_direct_Fine200_entry',F(2,600)<F(7,625) and F(2,600)<F(1,200))
ck('B11_phase_angle',F(8,600)<F(1,70) and 1-F(1,70)**2>F(999,1000)**2)
gt('B11_upper_strip_actual_count',(F(99,100)-F(4,5))*x1['rr']*F(7,44),F(2,600))
gt('B11_upper_strip_depth',F(4,5)*F(999,1000)*x1['rr']-x1['err']-F(4,42000),F(67,1000))

def radial_cap(tag,r,J,lo,U):gt(tag,U**2,(r/2+F(J,N))/(lo*(1-r*lo/2)))
def samefield(tag,r,fx,D,kappa,U,H):
 a=fx*(1-kappa*U);be=2*D*(a*U+F(1,n))
 gt(tag+'/positive_factor',1,2*kappa*U);gt(tag+'/denominator',1,be)
 gt(tag+'/monotonicity',-fx*kappa+2*D*(a*a+fx*F(1,H)*(1-2*kappa*U)),0)
 return r*(a+F(1,H))/(1-be)
kap=F(29591,63600);U=F(311,1000)
radial_cap('B11_low_depth_radial',r,9,F(67,1000),U)
up=samefield('B11_low_depth',r,fx,D,kap,U,600);gt('B11_low_depth_HB',F(1,50),up)
Delta=F(7,25)*F(999,1000)*xd['rr']-xd['err']-F(4,42000)-F(3481,5000)*r
gt('B11_strip_inset',Delta,0)
t=sp.Symbol('t');x=sp.Symbol('x')
S=lambda a:sp.Rational(a.numerator,a.denominator) if isinstance(a,F) else sp.Rational(a)
dist=F(157,1200)
poly=(t-S(c)*(1+t*t))*(t+S(dist))**2-2*(t+S(Delta))*(t*t+sp.Rational(1,10000))
expanded=sp.Poly(poly.subs(t,S(F(67,1000))+(S(F(3,10))-S(F(67,1000)))*x),x)
deg=expanded.degree();co=[expanded.nth(j) for j in range(deg+1)]
bern=[sum(co[j]*sp.Rational(comb(k,j),comb(deg,j)) for j in range(k+1)) for k in range(deg+1)]
ck('B11_whole_strip_Bernstein',all(b>0 for b in bern));notes['B11_Bernstein']=[str(v) for v in bern]
for tt in [F(67,1000),F(3,10)]:gt('B11_strip_endpoint_'+str(tt),(F(99,100)-F(7,25))*xd['rr']*F(7,44)*(tt+Delta)/(tt+dist)**2,F(1,50))
radial_cap('B11_high_depth_radial',r,9,F(3,10),F(147,1000))
beta=2*D*(fx*F(147,1000)+F(1,n));gt('B11_high_depth_HB',F(2,125),r*(fx+F(1,600))/(1-beta))
for tt in [F(3,10),F(1)]:gt('B11_high_depth_cloud_'+str(tt),(x1['density']-F(2,600))*(tt-F(3,20)-F(5,4)*r)/(tt+F(901,6000))**2,F(1,50))
ck('B11_zero_inset',N>1200)
Bs={600:F(193,1000),450:F(29,200),200:F(13,200)};Ts={600:F(13,200),450:F(9,100),200:F(21,100)}
for H in [600,450,200]:
 gt('B11_bounded_den_'+str(H),Bs[H],2*D*(H*f+1)/n)
 gt('B11_moving_depth_'+str(H),(1-r)**2/(2*Ts[H]),r*(H*f+1)/(1-Bs[H]))
for H,old in [(450,600),(200,450)]:
 for lo,hi in [(Ts[old],F(1,5)),(F(1,5),F(1))]:
  total=F(0);prev=F(0);lastweight=None;weights=[];ks=[]
  for z in xmid[H]:
   k=z['density']-z['loss']/H;w=weight(r,z['R'],F(1,100),lo,hi)
   gt('coupled_count_increment_'+z['name']+'_'+str(lo),k,prev)
   if lastweight is not None:ck('coupled_nested_weight_'+z['name']+'_'+str(lo),lastweight>=w)
   total+=w*(k-prev);prev=k;lastweight=w;weights.append(w);ks.append(k)
  abel=weights[-1]*ks[-1]+sum(((weights[j]-weights[j+1])*ks[j] for j in range(len(weights)-1)),F(0))
  ck('full_coupled_Abel_'+str(H)+'_'+str(lo),abel==total)
  gt('whole_coupled_depth_cloud_'+str(H)+'_'+str(lo),total,r*(fx+F(1,H))/(1-Bs[old]))
  notes['coupled_'+str(H)+'_'+str(lo)]={'depth':[lo,hi],'lower':total,'upper':r*(fx+F(1,H))/(1-Bs[old])}
 ck('closed_coupled_depth_cover_'+str(H),Ts[old]<F(1,5)<1)
a=weight(r,F(3,20),F(1,100),F(21,100),F(1));b=(1-r)/(1+F(1,10000))
gt('B11_beta_source_full_depth',(1-c*(1-F(1,10000)))*F(21,100),c*(1-F(1,10000))+F(1,10000))
gt('B11_X78_source_cloud',a*x78['density']+min(F(0),b-3*a)/78,r*(fx+F(1,78))/(1-Bs[200]))
gt('B11_X78_final',b,r*(78*fx+1)/(1-Bs[200]))
d0=F(99,100);BT=F(1,20)
gt('B11_T_gap',2*(1-d0),(2-d0*d0)*r)
for H in [150,103]:gt('B11_T_D100_'+str(H),d0*H,100)
ck('B11_Twide_independent_entries',F(25,1650)<F(9,400) and 1/(150*d0)<F(1,100))
gt('B11_Twide_phase',11,F(23,5)*F(25,11))
radial_cap('B11_T_radial',r,11,F(1),F(21,250))
beta=2*D*(fx*F(21,250)+F(1,n));gt('B11_T150_denominator',1,beta)
for z,den in [(t150,1-beta),(t103,1-BT)]:gt('B11_T_cloud_'+z['name'],(z['density']-2/z['H'])*weight(r,z['R'],F(3,10),F(1),F(1)),r*(fx+1/z['H'])/den)
gt('B11_T_bounded150',BT,2*D*(150*fx+1)/n)
bT=(1-r)/F(109,100);gsmall=F(77,100)
gt('B11_T_smallX',bT,r*(103*fx*gsmall+1)/(1-BT))
ck('B11_T_low_geometry_entry',2/(gsmall*79)<F(33,1000) and 1/(d0*79)<F(8,625))
companion('B11_companion',F(33,1000),F(8,625),F(719,100),F(7,25),F(71,50),350,F(15,2),F(32,25),F(3,10),F(11,50))
gt('B11_companion_physical',F(19,200),F(15,2)/79)
gt('B11_T_source_pair',bT+2*weight(r,F(19,200),F(3,10),F(1),F(1)),r*(103*fx+1)/(1-BT))
gt('B11_T79_final',bT,r*(79*fx+1)/(1-BT))


# New quadratic polar identity, independently expanded as real polynomials.
rhos,ss,rr,tt,al=sp.symbols('rho s r t alpha',real=True)
cc=rhos/2
lhs=(ss-rr-cc*ss*tt)**2+(tt-cc*ss*rr)**2-(1-rhos**2*ss**2/2)*((ss-rr)**2+tt**2)
rhs=rhos*ss**2*(cc*((rr-al)**2+tt**2)-tt)+rhos**2*ss**2*(ss-al)*(ss+al-2*rr)/2+rhos**2*ss**2*(rr**2+tt**2)/4
ck('new_whole_quadratic_polar_identity',sp.expand(lhs-rhs)==0)
ck('whole_quadratic_loss_positive',1-F(9,2)*r*r>0)
ck('moving_inverse_Cayley_positive_bracket',1-4*c*c>0 and 2/(1-r)**2<3)
qf=F(7,625);mf=F(1,200)
ck('full_fine_phase_before_use',21+133*qf+1030*qf*qf+9000*qf**3+120000*qf**4+220*mf<24)
ck('complete55_contours',len(contours)==55)
paths=['agents/global_component_refinement/Origin6251/certify_base6251_coupled.py']
def serial(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
 if isinstance(x,list):return list(map(serial,x))
 return x
receipt={'status':'PASS','n0':n,'precision_bits':ctx.prec,'author_function_imports':False,'method':'independent manuscript constants, Fractions, Arb exponentials, adaptive full path intervals, SymPy Bernstein coefficient construction','checks':len(checks),'full_inverse_contours':len(contours),'complete_companion_rectangles':1,'whole_path_majorants':9,'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},'notes':notes,'predicates':checks,'seconds':time.time()-tic}
out=ROOT/'computations/origin6251/ROOT_BASE010_6251_INDEPENDENT.json';out.write_text(json.dumps(serial(receipt),indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','n0','precision_bits','checks','full_inverse_contours','complete_companion_rectangles','whole_path_majorants','seconds']},indent=2))
