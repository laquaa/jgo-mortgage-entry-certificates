#!/usr/bin/env python3
"""Analytic-bound arithmetic and finite localization, exact rationals only."""
from fractions import Fraction as F
from pathlib import Path
import json
from certificate_checks import mandatory, decimal_bound
S=10**12
def floor(x):return F(x.numerator*S//x.denominator,S)
def ceil(x):return -floor(-x)
def js(x):return {'rational':str(x),'decimal_lower':decimal_bound(x),'decimal_upper':decimal_bound(x,True)}
L=F(7,10);Y=F(9,10);U=F(11,10);kap=F(1,50);a=F(1,20);b=F(3,20);eta=F(1,100)
rows=[]
printed={2:('.058209','.967431','.151517'),3:('.148838','.982628','.225967'),4:('.390947','.866621','.361303')}
for n in [2,3,4]:
 low=L**n/(F(n+1)*Y**2);interior=L**(n-2)/F(n+1);high=min(interior,F(n)/(Y*(U+1-L)))
 alpha=(a+kap)*min(low,interior,high)
 if n==2:cint=(1+2/L)/3
 elif n==3:cint=(Y+2+3/L)/4
 else:cint=(Y**2+2*Y+3+4/L)/5
 defect=cint+F(n)/(4*L*(1-Y))
 mandatory(alpha>F(1,250), 'power_family_bounds.py:19 mandatory alpha>F(1,250)');mandatory(defect<17, 'power_family_bounds.py:19 mandatory defect<17')
 r=F(n)/L;k1=r-n+1;fy=L**(n-2)/F(n+1)
 H=ceil(eta*k1/(fy-eta*r/2));mandatory(0<H<L, 'power_family_bounds.py:21 mandatory 0<H<L')
 iterates=[]
 for j in range(6):
  q=1-H/L;m=floor(L**(n-2)*sum(F(k+1,n+1)*q**k for k in range(n))-eta*r/2)
  mandatory(m>0, 'power_family_bounds.py:25 mandatory m>0')
  cmax=ceil((b+kap)*k1**2/(4*m))
  iterates.append({'n':j,'H_upper':js(H),'m_lower':js(m),'C_upper':js(cmax)})
  H=min(H,ceil(eta*k1/m))
 mandatory(cmax<F(181,500), 'power_family_bounds.py:29 mandatory cmax<F(181,500)')
 h0=F(iterates[0]['H_upper']['rational']);m5=F(iterates[5]['m_lower']['rational'])
 ph,pm,pc=map(F,printed[n])
 mandatory(h0<=ph and m5>=pm and cmax<=pc,'Section 5.2 printed H0/m5/C5 bounds, exponent '+str(n))
 rows.append({'nu':n,'alpha_lower':js(alpha),'defect_upper':js(defect),'bootstrap':iterates})
pmin=min(F(n)*(a+kap)*L**(n-1) for n in [2,3,4]);mandatory(pmin==F(2401,25000), 'power_family_bounds.py:31 mandatory pmin==F(2401,25000)')
cm=F(181,500);boundary=15*cm**2/(2*pmin);mandatory(boundary<11, 'power_family_bounds.py:32 mandatory boundary<11')
outside=L**5/5-F(1,100);mandatory(outside==F(11807,500000), 'power_family_bounds.py:33 mandatory outside==F(11807,500000)')
# Primitive positive lower funding: barV >=13/8-4/5; Phi<=.9^3/3.
funding=F(13,8)-F(4,5)-kap-Y**3/3;mandatory(funding==F(281,500), 'power_family_bounds.py:35 mandatory funding==F(281,500)')
# Independently evaluate the displayed common primitive constants.
mandatory(min(L**n for n in (2,3,4))==F(2401,10000),'family Rmin')
mandatory(min(F(n)/Y**2 for n in (2,3,4))==F(200,81),'family curvature minimum')
beta_max=max(F(n)*(b+kap)/(b+F(13,20))-1 for n in (2,3,4))
gamma_min=min(F(n)*(a+kap)/(a+F(3,4)) for n in (2,3,4))
mandatory(beta_max==-F(3,20) and gamma_min==F(7,40),'family beta and gamma margins')
mandatory(min(20+ll*(600*tt-36) for ll in (F(-1,4),F(1,2)) for tt in (a,b))==F(13,2),'family W minimum')
mandatory(min(10+ll*(200*tt-20) for ll in (F(-1,4),F(1,2)) for tt in (a,b))==5,'family density minimum')
mandatory(max(10+ll*(200*tt-20) for ll in (F(-1,4),F(1,2)) for tt in (a,b))==15,'family density maximum')
mandatory(min(F(n)*F(1,4)**(n-1) for n in (2,3,4))==F(1,16),'family buffered direction density floor')
mandatory(min((n-(n-1)*Y)/Y**(n+1) for n in (2,3,4))>1,'family wprime lower bound')
ybmin=b+F(13,20)
mandatory(max(2*(n-(n-1)*ybmin)/(n*ybmin**(n-1)) for n in (2,3,4))==F(25,16),'family endpoint ratio maximum')
mandatory(F('.32')*F('.12')*F(16,25)>F(1,50),'family normalized leading slope arithmetic')
# Two elementary failed members: positive primitives persist.
W_at_b=20+54*F(-1,2);beta_n6=6*(b+kap)/(b+F(7,10))-1
mandatory(W_at_b==-7 and beta_n6==F(1,5), 'power_family_bounds.py:38 mandatory W_at_b==-7 and beta_n6==F(1,5)')
out={'rows':rows,'uniform':{'pmin':js(pmin),'fmin':js(F(5)),'Wmin':js(F(13,2)),'beta0':js(F(3,20)),'gamma0':js(F(7,40)),'alpha':js(F(1,250)),'defect':js(F(17)),'funding_lower':js(funding),'participation_slack':js(outside),'direct_radius':js(eta),'direct_Cmax':js(cm),'boundary_constant':js(F(11))},'failed_cases':{'lambda_minus_half_Wb':js(W_at_b),'nu6_z_point_beta_b':js(beta_n6)},'assertions':'All computations exact Fraction; decimals are display only.'}
Path(__file__).with_name('power_family_bounds.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'direct_Cmax':'181/500','boundary_constant':11,'global_alpha':'1/250','global_defect':17,'bootstrap_rows':len(rows),'all_assertions_passed':True}))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
