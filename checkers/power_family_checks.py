#!/usr/bin/env python3
"""Exact-rational continuous-parameter certificates; no sampled-law proof."""
from fractions import Fraction as F
from pathlib import Path
import json,time
import sympy as S
from certificate_checks import mandatory, decimal_bound
class I:
 def __init__(self,a,b=None): self.a=F(a);self.b=F(a if b is None else b);mandatory(self.a<=self.b, 'power_family_checks.py:8 mandatory self.a<=self.b')
 def __add__(self,o):
  o=cv(o);return I(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return I(-self.b,-self.a)
 def __sub__(self,o):return self+-cv(o)
 def __rsub__(self,o):return cv(o)+-self
 def __mul__(self,o):
  o=cv(o);q=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b];return I(min(q),max(q))
 __rmul__=__mul__
 def __truediv__(self,o):
  o=cv(o);mandatory(not o.a<=0<=o.b, (o.a,o.b));return self*I(1/o.b,1/o.a)
 def __rtruediv__(self,o):return cv(o)/self
 def __pow__(self,n):
  mandatory(n>=0, 'power_family_checks.py:22 mandatory n>=0')
  if n==0:return I(1)
  if self.a>=0:return I(self.a**n,self.b**n)
  if self.b<=0:return I(self.b**n,self.a**n) if n%2==0 else I(self.a**n,self.b**n)
  return I(0,max(abs(self.a),abs(self.b))**n) if n%2==0 else I(self.a**n,self.b**n)
 def absmax(self):return max(abs(self.a),abs(self.b))
 def union(self,o):return I(min(self.a,o.a),max(self.b,o.b))
 def js(self):return {'lower':str(self.a),'upper':str(self.b),'decimal_lower':decimal_bound(self.a),'decimal_upper':decimal_bound(self.b,True)}
def cv(o):return o if isinstance(o,I) else I(o)
def frac(x):return F(int(S.numer(x)),int(S.denom(x)))
def poly(expr,vars):return [(powers,frac(coef)) for powers,coef in S.Poly(S.expand(expr),*vars).terms()]
def ev(p,vals):
 out=I(0)
 for powers,c in p:
  term=I(c)
  for val,k in zip(vals,powers):term=term*cv(val)**k
  out=out+term
 return out
A=F(1,20); B=F(3,20); K=F(1,50); Z0=F(13,20);Z1=F(3,4);L0=F(-1,4);L1=F(1,2)
t,z,la,c=S.symbols('t z la c');f=10+la*(200*t-20);W=2*f+(t+S.Rational(1,50))*S.diff(f,t)
pmin=F(2401,25000) # .09604

def make(n):
 R=(t+z)**n;Q=(t+S.Rational(1,50))*R-(t+z)**(n+1)/(n+1)
 mean=S.integrate(Q*f,(t,S.Rational(1,20),S.Rational(3,20)))
 P=Q-mean
 J=S.integrate(((t+S.Rational(1,50))*S.diff(R,t)-R)*f,(t,c,S.Rational(3,20)))
 J2=S.integrate(((t+S.Rational(1,50))*S.diff(R,t,2)-S.diff(R,t))*f,(t,c,S.Rational(3,20)))
 mass=S.integrate(f,(t,c,S.Rational(3,20)))
 LL=(S.Rational(17,100))*f.subs(t,S.Rational(3,20))*R.subs(t,S.Rational(3,20))-S.integrate(W*R,(t,c,S.Rational(3,20)))
 N1=S.expand((1+S.Rational(1,50)-z)*mass*(c+z)**n-(c+z-1)*J)
 N2=S.expand((c+z)*J2-n*J)
 def pack(ex):return {'p':poly(ex,[c,z,la]),'derivs':[poly(S.diff(ex,v),[c,z,la]) for v in [c,z,la]]}
 return {'N1':pack(N1),'N2':pack(N2),'LL':pack(LL),'P':poly(P,[t,z,la]),'Pz':poly(S.diff(P,z),[t,z,la]),'Pla':poly(S.diff(P,la),[t,z,la]),'J':poly(J,[c,z,la]),'J2':poly(J2,[c,z,la]),'mass':poly(mass,[c,z,la]),'Lambda':poly(LL,[c,z,la]),'mean':poly(mean,[z,la])}

def centered(pack,vals):
 mids=[(v.a+v.b)/2 for v in vals];base=ev(pack['p'],mids);err=sum(ev(dp,vals).absmax()*(v.b-v.a)/2 for dp,v in zip(pack['derivs'],vals));return I(base.a-err,base.b+err)

def root_at(n,obj,zz,ll):
 lo=A;hi=B
 mandatory(ev(obj['P'],[lo,zz,ll]).b<0<ev(obj['P'],[hi,zz,ll]).a,'midpoint root initially bracketed')
 for _ in range(44):
  mid=(lo+hi)/2;v=ev(obj['P'],[mid,zz,ll]);mandatory(v.a==v.b, 'power_family_checks.py:63 mandatory v.a==v.b')
  if v.a<0:lo=mid
  else:hi=mid
 return I(lo,hi)

def cell(n,obj,zz,ll):
 zm=(zz.a+zz.b)/2;lm=(ll.a+ll.b)/2;center=root_at(n,obj,zm,lm)
 dz=ev(obj['Pz'],[I(A,B),zz,ll]).absmax();dl=ev(obj['Pla'],[I(A,B),zz,ll]).absmax()
 error=(dz*(zz.b-zz.a)+dl*(ll.b-ll.a))/(2*pmin)
 rt=I(max(A,center.a-error),min(B,center.b+error));yy=rt+zz
 ws=(yy-1)/(yy**n);rs=n/yy
 N1=centered(obj['N1'],[rt,zz,ll]);N2=centered(obj['N2'],[rt,zz,ll]);LL=centered(obj['LL'],[rt,zz,ll])
 C1=N1/yy**n;C2=N2/yy
 rho=2*(n-(n-1)*(B+zz))/(n*(B+zz)**(n-1))
 if C2.b>=0:return None, {'c2':C2.js(),'root':rt.js()}
 star=-C1/C2;gap=rho-star
 return {'root':rt,'C1':C1,'C2':C2,'Astar':star,'rho':rho,'implementation_gap':gap,'Lambda':LL}, {'midpoint_root':center.js(),'midpoint_root_endpoint_values':[ev(obj['P'],[center.a,zm,lm]).js(),ev(obj['P'],[center.b,zm,lm]).js()],'bisections':44,'Pz_abs_bound':str(dz),'Plambda_abs_bound':str(dl),'root_enlargement':str(error),'numerator_N1':N1.js(),'numerator_N2':N2.js(),'Lambda_polynomial':LL.js()}

def run(n,nz,nl):
 obj=make(n);agg={};bad=0;records=[]
 for iz in range(nz):
  zz=I(Z0+(Z1-Z0)*iz/nz,Z0+(Z1-Z0)*(iz+1)/nz)
  for il in range(nl):
   ll=I(L0+(L1-L0)*il/nl,L0+(L1-L0)*(il+1)/nl)
   out,details=cell(n,obj,zz,ll)
   if out is None:return None,details
   for k,v in out.items():agg[k]=v if k not in agg else agg[k].union(v)
   if not(out['C1'].a>0 and out['C2'].b<0 and out['Lambda'].a>0 and out['implementation_gap'].a>0):bad+=1
   records.append({'z_index':iz,'lambda_index':il,'z_box':zz.js(),'lambda_box':ll.js(),'root_certificate':details,'bounds':{k:v.js() for k,v in out.items()},'status':'pass' if out['C1'].a>0 and out['C2'].b<0 and out['Lambda'].a>0 and out['implementation_gap'].a>0 else 'fail'})
 return {'nu':n,'z_cells':nz,'lambda_cells':nl,'all_cells_verified':bad==0,'bad_cells':bad,'bounds':{k:v.js() for k,v in agg.items()},'cells':records},None

PRINTED={
 2:{'C1':('.0675','.1142'),'C2':('-.1454','-.1071'),'Lambda':('.4926','.8666'),'Astar':('.5655','.8822'),'gap':'.4919'},
 3:{'C1':('.0698','.1450'),'C2':('-.1999','-.1262'),'Lambda':('.4002','.8348'),'Astar':('.4216','.9702'),'gap':'.2926'},
 4:{'C1':('.0661','.1819'),'C2':('-.2461','-.1313'),'Lambda':('.3253','.8008'),'Astar':('.3244','1.1828'),'gap':'.1235'}}
def verify_printed(result):
 checks=[];targets=PRINTED[result['nu']]
 for name in ('C1','C2','Lambda','Astar'):
  box=result['bounds'][name];lo,hi=map(F,targets[name])
  mandatory(lo<=F(box['lower']) and F(box['upper'])<=hi,'Table 1 '+str(result['nu'])+' '+name)
  checks.append('Table 1 '+name+' interval contained')
 mandatory(F(result['bounds']['implementation_gap']['lower'])>=F(targets['gap']),'Table 1 joint implementation gap')
 for name,condition in [('C1 > .066',F(result['bounds']['C1']['lower'])>F('.066')),('C2 < -.107',F(result['bounds']['C2']['upper'])<F('-.107')),('Lambda > .325',F(result['bounds']['Lambda']['lower'])>F('.325')),('joint gap > .123',F(result['bounds']['implementation_gap']['lower'])>F('.123'))]:
  mandatory(condition,name);checks.append(name)
 mandatory(len(result['cells'])==256 and result['z_cells']==16 and result['lambda_cells']==16,'exact original 256 full parameter cells')
 result['printed_table_checks']=checks+['joint gap supports printed lower bound','all 256 original cells retained']
if __name__=='__main__':
 start=time.time();results=[]
 for n in [2,3,4]:
  for grid in [16]:
   result,error=run(n,grid,grid)
   print(n,grid,'error' if error else result['bad_cells'],flush=True)
   if result and result['all_cells_verified']:verify_printed(result);results.append(result);break
  else:raise AssertionError((n,error,result))
 mandatory(F('.32')*F('.12')*F(16,25)>F(1,50),'uniform normalized leading slope exceeds 1/50')
 out={'parameter_box':{'D':'1','a':str(A),'b':str(B),'kappa':str(K),'z':[str(Z0),str(Z1)],'lambda':[str(L0),str(L1)],'H1':'3/2','H2':'7/4','e':'1/100'},'method':'Exact Fraction intervals; centered root plus mean-value parameter bound; every rectangular cell covered.','results':results,'uniform_leading_slope_lower':str(F('.32')*F('.12')*F(16,25)),'status':'pass'}
 Path(__file__).with_name('power_family_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'pass','verified_cells':sum(len(v['cells']) for v in results),'bad_cells':0,'elapsed_seconds':time.time()-start},indent=2))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
