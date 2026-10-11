#!/usr/bin/env python3
"""Verify the paper's compact rational root boxes and rounded inequalities.
Run from this directory with Python 3 and SymPy. This imports only the supplied
standalone exact arithmetic definitions, not a previous output certificate.
"""
from pathlib import Path
import importlib.util,sys,json
from fractions import Fraction as Q
from certificate_checks import mandatory, decimal_bound
p=Path(__file__).resolve().parent/'numeric'/'precision_certificate.py'
spec=importlib.util.spec_from_file_location('pc',p);m=importlib.util.module_from_spec(spec);sys.modules['pc']=m;spec.loader.exec_module(m)
zb=m.Box(Q('0.393103890949911708598915'),Q('0.393103890949911708598916'))
rb=m.Box(Q('0.0763775839493957928969'),Q('0.0763775839493957928970'))
kb=m.Box(Q('0.0067476090616167526307'),Q('0.0067476090616167526308'))
boxes={m.z:zb,m.r:rb,m.s:m.Box.point(Q(655999774,10**11)),m.k:kb,m.eta:m.Box.point(Q(1,20))}
checks=[]
root_checks=[]
def require(label,condition):
 mandatory(condition, label);checks.append(label)
def opposing(expr,var,box):
 v0=m.evaluate(expr,{**boxes,var:m.Box.point(box.lo)});v1=m.evaluate(expr,{**boxes,var:m.Box.point(box.hi)})
 require('root opposing '+str(var),(v0.hi<0 and v1.lo>0) or(v0.lo>0 and v1.hi<0))
 root_checks.append({'polynomial':str(expr),'variable':str(var),'root_box':box.output(),'left_value':v0.output(),'right_value':v1.output(),'status':'pass'})
opposing(m.pricing,m.z,zb);opposing(m.P.subs(m.t,m.r),m.r,rb);opposing(m.g.subs(m.t,m.r),m.k,kb)
require('unique pricing root',m.evaluate(m.S.diff(m.pricing,m.z),{m.z:m.Box(Q('.393'),Q('.394'))}).lo>0)
require('unique margin root',m.bernstein(m.S.diff(m.P,m.t),m.t,m.a,m.b,boxes).lo>0)
require('unique compensating root',m.bernstein(m.S.diff(m.g.subs(m.t,m.r),m.k),m.k,Q('0.0067'),Q('0.0068'),boxes).hi<0)
require('gain derivative above .00018',m.bernstein(m.S.diff(m.g,m.t),m.t,m.a,m.b,boxes).lo>Q('.00018'))
require('gain derivative below .00023',m.bernstein(m.S.diff(m.g,m.t),m.t,m.a,m.b,boxes).hi<Q('.00023'))
require('margin derivative positive',m.bernstein(m.S.diff(m.p,m.t),m.t,m.a,m.b,boxes).lo>0)
actual=m.evaluate(m.primitive,{**boxes,m.t:m.Box.point(m.b)})-m.evaluate(m.primitive,{**boxes,m.t:rb})
I=m.evaluate(m.baseline_primitive,{**boxes,m.t:m.Box.point(m.b)})-m.evaluate(m.baseline_primitive,{**boxes,m.t:rb})
require('printed profit lower',actual.lo>Q('.00247240900689302'))
require('printed profit upper',actual.hi<Q('.00247240900689303'))
gain_b=m.evaluate(m.g,{**boxes,m.t:m.Box.point(m.b)})
require('printed endpoint gain interval',Q('.00000453776882578')<gain_b.lo and gain_b.hi<Q('.00000453776882579'))
bernstein_rows={}
for name,expr,targets in [
 ('gain_derivative',m.S.diff(m.g,m.t),[('.00022333850256','.00022333850257'),('.00020344024233','.00020344024234'),('.00018213489718','.00018213489720')]),
 ('margin_derivative',m.S.diff(m.p,m.t),[('.26851116350379','.26851116350381'),('.34564453410350','.34564453410352'),('.43255800555610','.43255800555613'),('.53000157786159','.53000157786163')])]:
 cc=[m.evaluate(c,boxes) for c in m.coefficients(expr,m.t,m.a,m.b)]
 require(name+' coefficient count',len(cc)==len(targets))
 for index,(val,(lower,upper)) in enumerate(zip(cc,targets)):
  require(name+' printed coefficient '+str(index),Q(lower)<=val.lo and val.hi<=Q(upper))
 bernstein_rows[name]=[v.output() for v in cc]
require('shifted thresholds including maximum tested errors stay buffered',(zb+m.a-kb-Q('0.00001')).lo>Q('.4') and (zb+m.b+boxes[m.s]+Q('0.00001')).hi<Q('.5'))
require('implemented promises remain positive',(m.Box.point(Q(43,32))-zb-boxes[m.s]-Q('0.00001')).lo>0 and (m.Box.point(Q(247,160))-zb+kb-Q('0.00001')).lo>0)
L=Q(773,500)
Cf=20*(Q(1,20)+Q(1,50))+20*(Q(1,10)+Q(1,50))+40*Q(1,20)
MP=m.bernstein(m.S.diff(m.P,m.t),m.t,m.a,m.b,boxes)
Fmax=Q(9,2)+Q(1,20)
tail_mass=20*(m.Box.point(m.b)-rb)
require('C_f equals 5.8',Cf==Q('5.8'))
require('M_P bound .54',0<MP.lo and MP.hi<Q('.54'))
require('type density fmax equals 20',1/(Q(1,10)-Q(1,20))==20)
require('severity density F bound 4.55',Fmax==Q('4.55'))
require('exact-client grant mass printed prefix',Q('.47244832101208414')<tail_mass.lo and tail_mass.hi<Q('.47244832101208415'))
require('receipt Lipschitz bound exact',L==1+Q(3,25)*Q(91,20))
coarse=actual-Q('5.8')*Q('0.000001')-Q('.54')*20*(Q('0.000001')/Q('.00018'))**2
require('Eq 46 printed analytic lower',coarse.lo>Q('.0021332756'))
require('Eq 46 fixed cost .0021 permitted',coarse.lo>Q('.0021'))
# Compact rational cutoff and margin-root boxes appearing in the paper.
rows=[('0','0.0000001',
       '.07588291727878614','.07588291727878615',
       '.07687325301479577','.07687325301479578',
       '.07632915816785827','.07632915816785828',
       '.00247037026058'),
      ('0','0.000001',
       '.07147493902527849','.07147493902527850',
       '.08138053055345920','.08138053055345921',
       '.07633264873629321','.07633264873629322',
       '.00227577718525'),
      ('.000001','.0000001',
       '.08087552439376768','.08087552439376769',
       '.08188661072343721','.08188661072343722',
       '.07632915816785827','.07632915816785828',
       '.00234620611201')]
out=[]
contribution_targets=[('.00247116245533','-.00000079219476'),('.00236784441240','-.00009206722716'),('.00234620611201','0')]
for tau,rho,lo0,lo1,hi0,hi1,q00,q01,bound in rows:
 tau,rho=Q(tau),Q(rho)
 lower=m.Box(Q(lo0),Q(lo1));upper=m.Box(Q(hi0),Q(hi1));q0=m.Box(Q(q00),Q(q01))
 opposing(m.g-m.S.Rational(tau-rho),m.t,lower)
 opposing(m.g-m.S.Rational(tau+rho),m.t,upper)
 adjusted=m.p-m.S.Rational(L*rho)
 opposing(adjusted,m.t,q0)
 require('gain threshold roots ordered',lower.hi<upper.lo)
 if tau==0:require('negative margin zero inside ambiguous band',lower.hi<q0.lo and q0.hi<upper.lo)
 else:require('negative margin zero below ambiguous band',q0.hi<lower.lo)
 F=m.primitive-20*m.S.Rational(L*rho)*m.t
 def value(q):return m.evaluate(F,{**boxes,m.t:q})
 guaranteed=value(m.Box.point(m.b))-value(upper)
 if lower.lo>=q0.hi:negative=m.Box.point(0)
 elif upper.hi<=q0.lo:negative=value(upper)-value(lower)
 else:negative=value(q0)-value(lower)
 B=guaranteed+negative
 require('paper robust lower '+str(tau)+' '+str(rho),B.lo>Q(bound))
 tail_target,loss_target=map(Q,contribution_targets[len(out)])
 require('printed guaranteed-tail lower '+str(tau)+' '+str(rho),guaranteed.lo>=tail_target)
 require('printed ambiguous-loss lower '+str(tau)+' '+str(rho),negative.lo>=loss_target)
 if tau==0 and rho==Q('0.0000001'):
  require('K I*+6e-6',B.lo>I.hi+Q('0.000006'))
  require('cost 1e-6 grant compensated K I*+5e-6',B.lo-Q('0.000001')>I.hi+Q('0.000005'))
 out.append({'tau':str(tau),'rho':str(rho),'r_minus':lower.output(),'r_plus':upper.output(),'adjusted_margin_zero':q0.output(),'guaranteed':guaranteed.output(),'ambiguous_negative':negative.output(),'lower':B.output(),'printed_total_lower':bound,'grant_compensated_lower':(B-Q('0.000001')).output() if tau==0 and rho==Q('0.0000001') else None,'interpretation':'fixed feasible law, every admitted payment-error pair; tau row is uncompensated; grant-compensated field offsets the separate known cost 1/1000000'})
# Both payments increased by the specified admissible error.
rho=Q('0.00001')
advg=m.S.expand((m.Phi(m.t+m.z+m.s-m.S.Rational(rho))+m.Phi(m.t+m.z-m.k-m.S.Rational(rho))-2*m.Phi(m.t+m.z)+m.eta*(m.Psi(m.t+m.z+m.s-m.S.Rational(rho))-m.Psi(m.t+m.z-m.k-m.S.Rational(rho))))/2)
advbox=m.bernstein(advg,m.t,m.a,m.b,boxes)
require('paper adverse rounding upper',advbox.hi<Q('-0.00000515635360'))
advcoeff=[m.evaluate(c,boxes) for c in m.coefficients(advg,m.t,m.a,m.b)]
advtargets=[('-.000013243054951','-.000013243054950'),('-.000010160049439','-.000010160049438'),('-.000007455831852','-.000007455831850'),('-.000005156353605','-.000005156353603')]
require('adverse gain cubic has four coefficients',len(advcoeff)==4)
for index,(value,(lo,hi)) in enumerate(zip(advcoeff,advtargets)):
 require('adverse printed Bernstein interval '+str(index),Q(lo)<=value.lo and value.hi<=Q(hi) and value.hi<0)
result={'status':'pass','checks':checks,'nuisance_boxes':{'z':zb.output(),'t_star':rb.output(),'compensation_k':kb.output()},'root_checks':root_checks,'Bernstein_coefficients':bernstein_rows,'actual_profit':actual.output(),'endpoint_gain':gain_b.output(),'I_star':I.output(),'precision_constants':{'C_f':str(Cf),'M_P_actual_derivative':MP.output(),'M_P_bound':'.54','fmax':'20','severity_density_bound':str(Fmax),'receipt_Lipschitz':str(L)},'exact_client_mass':tail_mass.output(),'Eq46_analytic_lower':coarse.output(),'rows':out,'adverse_gain':advbox.output(),'adverse_Bernstein_coefficients':[v.output() for v in advcoeff],'scope':'Compact paper boxes certify all displayed Appendix D signs, coefficient intervals, receipts and three robust lower bounds. Analytic propositions establish uniform error-pair coverage under one fixed feasible law; this is not law-wise guaranteed profit.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print('all',len(checks),'compact-box assertions passed')
for row in out:print(row['tau'],row['rho'],row['lower']['decimal_lower'])

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
