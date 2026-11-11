"""Whole finite T source cloud at n>=6251.
No discovery functions imported; all continuous source extrema and complete
critical-product logarithm tails are reconstructed with exact rationals.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sympy as s,json,hashlib,sys
sys.set_int_max_str_digits(0)
O=Path(__file__).resolve().parent
N=6250;n=6251;r0=F(39,500);rl=F(3,50)
x=s.symbols('x');checks=[];polys=[]
def S(v):return s.Rational(v.numerator,v.denominator) if isinstance(v,F) else s.sympify(v)
def ff(v):return F(int(v.p),int(v.q))
def ck(name,p):
 checks.append(dict(name=name,passed=bool(p)))
 if not p:raise AssertionError(name)
def polypos(name,expr,lo,hi):
 expr=s.cancel(expr);num,den=expr.as_numer_denom()
 # Sympy may normalize both signs negatively. Establish the whole
 # denominator sign independently instead of using a midpoint as evidence.
 if den.subs(x,S((lo+hi)/2))<0:num,den=-num,-den
 if x in den.free_symbols:polypos(name+'/entire cleared denominator',den,lo,hi)
 else:ck(name+'/positive constant denominator',den>0)
 pp=s.Poly(num,x);deg=pp.degree();co=[ff(pp.nth(i)) for i in range(deg+1)]
 stack=[(lo,hi,0)];leaves=[]
 while stack:
  a,b,d=stack.pop()
  p=[sum((co[k]*comb(k,j)*a**(k-j)*(b-a)**j for k in range(j,deg+1)),F(0)) for j in range(deg+1)]
  bs=[sum((p[k]*F(comb(j,k),comb(deg,k)) for k in range(j+1)),F(0)) for j in range(deg+1)]
  if min(bs)>0:leaves.append((a,b,min(bs)));continue
  if d>=20:raise AssertionError((name,a,b,min(bs),float(expr.subs(x,S((a+b)/2)))))
  m=(a+b)/2;stack.extend([(a,m,d+1),(m,b,d+1)])
 leaves.sort();cur=lo
 for a,b,m in leaves:ck(name+'/closed seam',a==cur);cur=b;ck(name+'/all Bernstein coefficients positive',m>0)
 ck(name+'/complete interval',cur==hi)
 polys.append(dict(name=name,expr=str(expr),numerator=str(num),denominator=str(den),interval=[str(lo),str(hi)],leaves=[[str(a),str(b),str(m)] for a,b,m in leaves]))
 return expr
def kapp(de):return s.cancel(x*(x-(1-S(de))**2)*(x+(1+S(de))**2)/(2*(1+S(de)**2)*(x+2*S(de))))


from math import factorial,isqrt
rows=[];cells=[];terminals=[]
def fl(v,den=10**12):return F(v.numerator*den//v.denominator,den)
def ce(v,den=10**12):return -fl(-v,den)
def sqrtup(v,den=10**12):
 k=isqrt(v.numerator*den*den//v.denominator)
 if F(k,den)**2<v:k+=1
 z=F(k,den);ck('whole positive radical bound',z*z>=v and (z-F(1,den))**2<v);return z
DEN={1:2,2:7,3:20,4:50,5:148,6:403,8:2900,10:22000}
for W,dd in DEN.items():ck('full exponential '+str(W),sum((F(W)**j/factorial(j) for j in range(70)),F(0))>dd)

data=json.loads((O/'TMOMENT6251_PROPOSAL.json').read_text())
ck('complete rational proposal',data['status']=='PROPOSAL_ONLY' and data['degree_floor']==n and not data['failed'])
ck('whole crescent gap positive at rho upper',F(1,4)-r0/2>0)
ck('all actual quarter loss conventions',F(22,7)>F(314159,100000)>3)
bounds=[]
for iz,z in enumerate(data['cells']):
 tag='Tmoment/'+str(iz)
 pl,ph,yl,yh,H,de,A,U,Xcap,J=[F(z[k]) for k in ['rho_lo','rho_hi','yl','yh','H','delta','Abar','U','Xcap','J']]
 ck(tag+'/entire closed cell domain',rl<=pl<ph<=r0 and F(5,4)<=yl<yh<=F(9,4))
 ck(tag+'/legal independent medium139 phase',J==4*yh)
 cf=F(2637,32)
 field=(N*pl-cf)/(2*(1+ph));radfield=(N*pl/2-J)/(ph+yh)
 ck(tag+'/same-rho full medium139 field',H>0 and (H<=field or H<=radfield))
 ck(tag+'/full phase real direction',0<J/H<1)
 ck(tag+'/no unproved same-X gap',Xcap==1)
 ck(tag+'/full inverse Cayley source test',F(1,200)<=de<=1 and (ph*(1+de)/2)**2<de)
 ka=kapp(de);xc=(1-de*de)/de
 for yy in [yl,yh]:
  if xc>0:polypos(tag+'/whole active X source '+str(yy),S(A)-x+ka*S(yy),F(0),min(Xcap,xc))
  if xc<Xcap:polypos(tag+'/whole inactive X source '+str(yy),S(A)-x+x*x*S(yy)/2,max(F(0),xc),Xcap)
 f=1/(1-ph*ph*(1+de)**2/2);Dp=1+(ph*(1+de)/2)**2
 Br=(ph/2+J/N)/(1-ph/2)
 ck(tag+'/whole radial field ratio',U>0 and (U>=yh or U*U>=Br or (yl>0 and U*yl>=Br)))
 B=2*Dp*(f*A*U+F(1,n))
 ck(tag+'/HB denominator before solve',0<B<1 and A>0)
 upper=ph*(f*A+1/H)/(1-B)
 ws=[];ss=[];prevR=F(0);prevsl=F(0);prevell=F(0)
 for ir,rr in enumerate(z['contours']):
  label=tag+'/contour'+str(ir);R,nu,W=F(rr[0]),F(rr[1]),int(rr[2]);qq=yh
  la=qq*R/(1-R)
  er=qq*(sum((R**j/j for j in range(2,17)),F(0))+R**17/(17*(1-R)))
  rh=R-er;Bi=qq/((1-R)**2*(1-la)**3)
  eps=F(1001,1000)*Bi/H;g=F(1001,1000)*nu*Bi*rh;p=1/(1+la);amp=1/(1-la);ee=F(1,DEN[W])
  for name,pred in [('full critical-free disk',0<R<1),('whole log injection',0<la<1 and rh>0),('all containing corners',W*W<(1-nu*nu)*(H*rh)**2),('entire left side',eps*(1+ee)+g*ee<1-amp*ee),('entire right side',eps*(1+ee)+g*ee*ee<p-ee),('both whole horizontals',eps+g*ee<p),('whole amplitude homotopy',0<p<1<amp),('finite odd height',nu*H*rh>F(22,7))]:ck(label+'/'+name,pred)
  slope=fl(nu*F(7,22)*rh);ell=min(R,ce(W/H+J*nu*rh/H+er));ell=max(prevell,ell)
  ck(label+'/nested whole root regions',R>prevR and slope>prevsl and prevell<=ell<=R)
  w=max(F(0),min((t-ph*(1+t*t)/2)/((de+ell)**2+t*t) for t in [1-R,1+R]));w=fl(w)
  ws.append(w);ss.append(slope)
  rows.append(dict(label=label,R=str(R),H=str(H),q=str(qq),nu=str(nu),W=W,lambda_upper=str(la),full_log_error=str(er),image_radius=str(rh),inverse_derivative=str(Bi),slope=str(slope),real_width=str(ell),raw_weight=str(w)))
  prevR,prevsl,prevell=R,slope,ell
 ck(tag+'/full nonempty family',bool(ws))
 for j in range(len(ws)-2,-1,-1):ws[j]=max(ws[j],ws[j+1])
 ck(tag+'/whole nested weights',ws[0]>0 and all(a>=b>=0 for a,b in zip(ws,ws[1:])))
 b=(1-ph)/(1+de*de);comp=min(F(0),b-3*ws[0])/H
 lower=sum((w*(a-bb) for w,a,bb in zip(ws,ss,[F(0)]+ss[:-1])),F(0))+comp
 abel=ws[-1]*ss[-1]+sum(((ws[j]-ws[j+1])*ss[j] for j in range(len(ws)-1)),F(0))+comp
 ck(tag+'/finite full Abel identity',abel==lower)
 ck(tag+'/whole cloud contradiction',lower>upper)
 cells.append(dict(label=tag,kind='cloud',rho=[str(pl),str(ph)],y=[str(yl),str(yh)],H=str(H),delta=str(de),Abar=str(A),U=str(U),Xcap=str(Xcap),J=str(J),upper=str(upper),lower=str(lower),margin=str(lower-upper),contours=len(ws),source_compensation=str(comp)))

 bounds.append((pl,ph,yl,yh))
slabs=sorted(set((p,r) for p,r,l,u in bounds));prev=rl
for pl,ph in slabs:
 ck('full closed rho seam',pl==prev);prev=ph;cur=F(5,4)
 for aa,bb in sorted((l,u) for p,r,l,u in bounds if p==pl and r==ph):
  ck('full closed y seam',aa==cur and bb>aa);cur=bb
 ck('full y endpoint',cur==F(9,4))
ck('whole rho endpoint',prev==r0)
out=dict(status='PASS',scope='Whole normalized T=1 source cloud, rho in [.06,.078], y in [5/4,9/4], n>=6251, under independently legal D139/phase4/CF2637/32',degree_floor=n,checks=len(checks),whole_polynomials=len(polys),closed_cells=len(cells),full_contours=len(rows),minimum_margin=str(min(F(c['margin']) for c in cells)),cells=cells,contours=rows,predicates=checks,polynomials=polys,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),proposal_sha256=hashlib.sha256((O/'TMOMENT6251_PROPOSAL.json').read_bytes()).hexdigest())
(O/'TMOMENT6251_COMPLETE_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['cells','contours','predicates','polynomials']},indent=2))
