"""Full rational n>=6251 T=1 strip certificate, including complete contours.

Imports only author's exact analytic-bootstrap implementation, never any
floating discovery function. Proposal data are rational choices only.
"""
from certify_T6251_highX_bootstrap import *
from math import factorial
rows=[];cells=[]
def fl(v,den=10**12):return F(v.numerator*den//v.denominator,den)
def ce(v,den=10**12):return -fl(-v,den)
DEN={1:2,2:7,3:20,4:50,5:148,6:403,8:2900,10:22000}
expbound={W:sum((F(W)**j/factorial(j) for j in range(70)),F(0)) for W in DEN}
for W,dd in DEN.items():ck('entire exponential W'+str(W),expbound[W]>dd)
ck('h95 lawful D100 entry',95*95*Gmin>10000)
ck('D100 wide q entry',F(9,4)/100==F(9,400))
ck('D100 wide mu entry',F(1,100)==F(1,100))
ck('all test inverse Cayley entry',r0*r0<F(1,200))
data=json.loads((O/'T6251_HIGHX_PROPOSAL.json').read_text())
ck('whole proposal complete',data['status']=='PROPOSAL_ONLY' and not data['failed'] and data['degree_floor']==n)
bounds=[]
for iz,z in enumerate(data['cells']):
 tag='cloud'+str(iz)
 pl,ph,yl,yh,H,de,A,U,Xcap,J=[F(z[k]) for k in ['rho_lo','rho_hi','yl','yh','H','delta','Abar','U','Xcap','J']]
 ck(tag+'/entire cell domain',rl<=pl<ph<=r0 and 0<=yl<yh<=F(427,200))
 ck(tag+'/entire same-y phase',J==F(23,5)*yh)
 ck(tag+'/ordered field',H>=95 and (H==95 or H<=(N*pl/2-J)/(ph+yh)))
 ck(tag+'/positive phase real direction',0<J/H<1)
 ck(tag+'/same-X full-gap cap',xmin<=Xcap<=1 and (Xcap==1 or yl*yl*(1+F(117,500)*Xcap*Xcap)>=F(81,16)))
 ck(tag+'/full finite source test',F(1,200)<=de<=1 and (ph*(1+de)/2)**2<de)
 ka=kapp(de);xc=(1-de*de)/de
 for yy in [yl,yh]:
  if xc>xmin:polypos(tag+'/whole active X source y'+str(yy),S(A)-x+ka*S(yy),xmin,min(Xcap,xc))
  if xc<Xcap:polypos(tag+'/whole inactive X source y'+str(yy),S(A)-x+x*x*S(yy)/2,max(xmin,xc),Xcap)
 f=1/(1-ph*ph*(1+de)**2/2);Dp=1+(ph*(1+de)/2)**2
 Br=(ph/2+J/N)/(1-ph/2)
 ck(tag+'/whole radial field ratio',U>0 and (U>=yh or U*U>=Br or (yl>0 and U*yl>=Br)))
 B=2*Dp*(f*A*U+F(1,n));ck(tag+'/HB denominator before solve',0<B<1 and A>0)
 upper=ph*(f*A+1/H)/(1-B)
 ws=[];ss=[];prevR=F(0);prevsl=F(0);prevell=F(0)
 for ir,rr in enumerate(z['contours']):
  label=tag+'/contour'+str(ir);R,nu,W=F(rr[0]),F(rr[1]),int(rr[2]);q=yh
  la=q*R/(1-R)
  er=q*(sum((R**j/j for j in range(2,17)),F(0))+R**17/(17*(1-R)))
  rh=R-er;Bi=q/((1-R)**2*(1-la)**3)
  eps=F(1001,1000)*Bi/H;g=F(1001,1000)*nu*Bi*rh
  p=1/(1+la);amp=1/(1-la);ee=F(1,DEN[W])
  ck(label+'/entire critical free disk',0<R<1)
  ck(label+'/entire logarithm injection',0<la<1 and rh>0)
  ck(label+'/all four enclosing corners',W*W<(1-nu*nu)*(H*rh)**2)
  ck(label+'/entire left vertical',eps*(1+ee)+g*ee<1-amp*ee)
  ck(label+'/entire right vertical',eps*(1+ee)+g*ee*ee<p-ee)
  ck(label+'/both whole horizontals',eps+g*ee<p)
  ck(label+'/whole amplitudes',0<p<1<amp)
  ck(label+'/finite height to count all strips',nu*H*rh>F(22,7))
  slope=fl(nu*F(7,22)*rh)
  ell=min(R,ce(W/H+J*nu*rh/H+er));ell=max(prevell,ell)
  ck(label+'/whole nested source region',R>prevR and slope>prevsl and prevell<=ell<=R)
  w=max(F(0),min((t-ph*(1+t*t)/2)/((de+ell)**2+t*t) for t in [1-R,1+R]))
  w=fl(w);ws.append(w);ss.append(slope)
  rows.append(dict(label=label,R=str(R),H=str(H),q=str(q),nu=str(nu),W=W,lambda_upper=str(la),log_error=str(er),image_radius=str(rh),inverse_derivative=str(Bi),slope=str(slope),real_width=str(ell),raw_weight=str(w)))
  prevR,prevsl,prevell=R,slope,ell
 ck(tag+'/nonempty full family',bool(ws))
 for j in range(len(ws)-2,-1,-1):ws[j]=max(ws[j],ws[j+1])
 ck(tag+'/nested nonnegative all-root weights',ws[0]>0 and all(a>=b>=0 for a,b in zip(ws,ws[1:])))
 b=(1-ph)/(1+de*de);comp=min(F(0),b-3*ws[0])/H
 lower=sum((w*(a-bb) for w,a,bb in zip(ws,ss,[F(0)]+ss[:-1])),F(0))+comp
 abel=ws[-1]*ss[-1]+sum(((ws[j]-ws[j+1])*ss[j] for j in range(len(ws)-1)),F(0))+comp
 ck(tag+'/exact finite Abel identity',abel==lower)
 ck(tag+'/strict whole cloud contradiction',lower>upper)
 cells.append(dict(label=tag,rho=[str(pl),str(ph)],y=[str(yl),str(yh)],H=str(H),delta=str(de),Abar=str(A),U=str(U),Xcap=str(Xcap),J=str(J),upper=str(upper),lower=str(lower),margin=str(lower-upper),contours=len(ws),source_compensation=str(comp)))
 bounds.append((pl,ph,yl,yh))
slabs=sorted(set((p,r) for p,r,l,u in bounds));prev=rl
for pl,ph in slabs:
 ck('full closed rho seam',pl==prev);prev=ph;cur=F(0)
 for a,b in sorted((l,u) for p,r,l,u in bounds if p==pl and r==ph):
  ck('full closed y seam',a==cur and b>a);cur=b
 ck('complete y ray endpoint',cur==F(427,200))
ck('complete rho final endpoint',prev==r0)
out=dict(status='PASS',scope='Complete T=1 higher strip .061002<=rho<=.062,.69<=X<=1 for all n>=6251',degree_floor=n,checks=len(checks),whole_polynomials=len(polys),closed_cells=len(cells),complete_contours=len(rows),min_cloud_margin=str(min(F(z['margin']) for z in cells)),paired_terminal_margin=result['terminal_margin'],cells=cells,contours=rows,predicates=checks,polynomial_certificates=polys,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),bootstrap_sha256=hashlib.sha256((O/'certify_T6251_highX_bootstrap.py').read_bytes()).hexdigest(),proposal_sha256=hashlib.sha256((O/'T6251_HIGHX_PROPOSAL.json').read_bytes()).hexdigest())
(O/'T6251_HIGHX_COMPLETE_EXACT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['cells','contours','predicates','polynomial_certificates']},indent=2))
