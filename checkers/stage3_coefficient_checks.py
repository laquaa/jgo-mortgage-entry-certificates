#!/usr/bin/env python3
"""Exact coefficient and one feasible finite-radius witness.

This checks the finite-dimensional tangent reduction and the stated witness.
It does not certify the full optimized prelimit upper remainder. Interval
arithmetic utilities are supplied in interval_utils.py. Exact primitives are
embedded; no data files or previous results are inputs.
"""
from pathlib import Path
import json
import sympy as s
import interval_utils as ir
from certificate_checks import mandatory, decimal_bound

ROOT=Path(__file__).resolve().parent
t,z,r,aa,k,y=s.symbols('t z r aa k y',real=True)
aT,bT,kap,D=s.Rational(1,20),s.Rational(1,10),s.Rational(1,50),s.Rational(1,2)
B=s.Rational(248617,240000)
checks=[]

def require(label,condition):
    mandatory(bool(condition), label)
    checks.append(label)

def R(v):return 3*v*v+2*v**3
def Phi(v):return v**3+v**4/2
def Rp(v):return 6*v+6*v*v
def integral(expr,lo=r,hi=bT):return s.factor(s.integrate(expr,(t,lo,hi))/(bT-aT))
def ev(expr,boxes):
    num,den=s.fraction(s.cancel(expr))
    return ir.evaluate(num,boxes)/ir.evaluate(den,boxes)
def brief(box):return {k:v for k,v in box.output().items() if k.startswith('decimal_')}

def main():
    Q=s.Rational(11,10)-kap+(t+kap)*R(t+z)-Phi(t+z)
    pricing=s.factor(integral(Q,aT,bT)-B)
    zbox=ir.bisect(pricing,z,s.Rational(393,1000),s.Rational(394,1000),{},75)
    require('unique actuarial pricing root',
        ir.evaluate(s.diff(pricing,z),{z:ir.Box(ir.frac('393/1000'),ir.frac('394/1000'))}).hi<0)
    rbox=ir.bisect(Q.subs(t,r)-B,r,s.Rational(76,1000),s.Rational(77,1000),{z:zbox},50)
    require('unique zero-margin root on borrower interval',
        ir.bernstein(s.diff(Q,t),t,aT,bT,{z:zbox}).lo>0)
    floor=ir.evaluate(integral(Q-B),{z:zbox,r:rbox})
    yr=r+z;yb=bT+z
    # Numerators use the common positive denominator R(yr).
    F1=R(yr)*(t+z-D)+(D-yr)*R(t+z)
    F2=R(yr)*Rp(t+z)-Rp(yr)*R(t+z)
    C1N=integral((t+kap)*s.diff(F1,t)-F1)
    C2N=integral((t+kap)*s.diff(F2,t)-F2)
    C1=s.cancel(C1N/R(yr));C2=s.cancel(C2N/R(yr))
    optimum=s.cancel(-C1N/C2N)
    coefficient=s.cancel(-C1N**2/(2*C2N*R(yr)))
    dual=s.cancel((bT+kap)*R(yb)/(bT-aT)-2*integral(R(t+z)))
    Dw_r=(D-yr)/R(yr)-(D-yb)/R(yb)
    Dr_r=Rp(yr)/R(yr)-Rp(yb)/R(yb)
    area1=2*integral((D-t-z)-R(t+z)*(D-yb)/R(yb))
    area2=-integral(Rp(t+z)-R(t+z)*Rp(yb)/R(yb))
    require('global endpoint dual linear coefficient identity',s.simplify(C1-dual*Dw_r-area1)==0)
    require('global endpoint dual quadratic coefficient identity',s.simplify(C2+dual*Dr_r-2*area2)==0)
    boxes={z:zbox,r:rbox}
    c1box,c2box,abox,cbox=[ev(expr,boxes) for expr in[C1,C2,optimum,coefficient]]
    dualbox=ev(dual,boxes)
    require('positive global root-feasibility dual multiplier',dualbox.lo>0)
    require('linear tangent profit coefficient positive',c1box.lo>0)
    require('quadratic tangent profit coefficient negative',c2box.hi<0)
    require('unique positive contrast optimizer',abox.lo>0)
    require('positive optimized tangent coefficient',cbox.lo>0)
    require('completion-of-square identity exact',
        s.simplify(C1*aa+C2*aa**2/2-coefficient-C2*(aa-optimum)**2/2)==0)
    # The endpoint objective is increasing in endpoint height.
    require('strict positive endpoint-height derivative lower bound',
        ev(2*r-bT+kap,boxes).lo>0)
    # The critical ratio a_c decreases on the entire active response interval.
    critical=(3-4*y*y)/(3*y*(3+4*y+2*y*y))
    n,d=s.fraction(critical)
    decreasing=s.expand(s.diff(n,y)*d-n*s.diff(d,y))
    require('critical ratio strictly decreases globally on borrower interval',
        ir.bernstein(decreasing,y,s.Rational(443,1000),s.Rational(494,1000),{}).hi<0)
    critbox=ev(critical.subs(y,yb),{z:zbox})
    require('optimizer below critical ratio throughout borrower interval',abox.hi<critbox.lo)
    gainN=aa*F1+aa**2*F2/2
    require('tangent gain exactly zero at profit cutoff',s.simplify(gainN.subs(t,r))==0)
    norm_der=s.expand(s.diff(gainN,t)*R(t+z)-gainN*Rp(t+z))
    require('normalized tangent gain strictly increasing over entire continuum',
        ir.bernstein(norm_der,t,aT,bT,{**boxes,aa:abox}).lo>0)
    # A buffered extremal law keeps every finite-witness threshold in its +1 piece.
    # psi=0 on[0,.2), -1/2 on[.2,.4), +1 on[.4,.5].
    require('buffered law derivative has zero integral',
        -s.Rational(1,2)*(s.Rational(2,5)-s.Rational(1,5))
          +(D-s.Rational(2,5))==0)
    eta=s.Rational(1,1000)
    require('perturbed densities nonnegative on compensated interval',
        6*s.Rational(1,5)+6*s.Rational(1,5)**2-eta/2>0)
    require('perturbed densities nonnegative on positive interval',
        6*s.Rational(2,5)+6*s.Rational(2,5)**2-eta>0)
    # h=a*eta is exactly represented by the algebraic coefficient formula.
    # Introduce aa as a root-isolated nuisance interval for polynomial evaluation.
    h=eta*aa
    gain=(Phi(t+z+h)+Phi(t+z-k)-2*Phi(t+z)
          +eta*((h+k)*(t+z-D)+(h*h-k*k)/2))/2
    tie=s.expand(gain.subs(t,r))
    kbox=ir.bisect(tie,k,s.Rational(132,10**6),s.Rational(134,10**6),
                  {**boxes,aa:abox},38)
    witness_boxes={**boxes,aa:abox,k:kbox}
    require('all finite-witness thresholds above buffered law boundary',
        ir.evaluate(aT+z-k-s.Rational(2,5),witness_boxes).lo>0)
    require('all finite-witness thresholds below upper collateral endpoint',
        ir.evaluate(D-bT-z-h,witness_boxes).lo>0)
    require('all finite-witness payments nonnegative',
        ir.evaluate(s.Rational(43,32)-z-h,witness_boxes).lo>0)
    actualnorm=s.expand(s.diff(gain,t)*R(t+z)-gain*Rp(t+z))
    require('actual normalized gain strictly increasing over whole type interval',
        ir.bernstein(actualnorm,t,aT,bT,witness_boxes).lo>0)
    # Direct paired receipts use the actual same law. For this linear-A region,
    # delta R pairing equals eta*(h+k), not an arbitrary optimistic recovery.
    paired_R=(R(t+z+h)+R(t+z-k)-2*R(t+z)+eta*(h+k))/2
    direct=(t+kap)*paired_R-gain
    require('actual cash-accounting derivative identity',
        s.simplify(direct-((t+kap)*s.diff(gain,t)-gain))==0)
    excess=ir.evaluate(integral(direct),witness_boxes)
    require('actual finite-radius witness exceeds positive-margin floor',excess.lo>0)
    actual=floor+excess
    remainder=excess-cbox*eta**2
    # This is a one-point interval, not a uniform remainder theorem.
    scaled=remainder/(eta**3)
    allowance=ir.Box.point(s.Rational(1800000,10**8))/cbox
    require('upper remainder allowance at eta0 is approximately 5.17',allowance.lo>ir.frac('5.16') and allowance.hi<ir.frac('5.18'))
    output={'scope':'exact tangent auxiliary coefficient, and one actual-law deterministic witness; full prelimit limsup/remainder must be proved separately',
      'checks':checks,'number_of_checks':len(checks),'incumbent_z':zbox.output(),
      'zero_margin_type':rbox.output(),'C1':c1box.output(),'C2':c2box.output(),
      'optimal_contrast':abox.output(),'optimal_tangent_coefficient':cbox.output(),
      'global_root_feasibility_dual_multiplier':dualbox.output(),
      'relative_upper_remainder_allowance_at_eta0':allowance.output(),
      'critical_contrast_at_upper_type':critbox.output(),'floor':floor.output(),
      'finite_witness':{'eta':str(eta),'law_direction':'0 on[0,1/5), -1/2 on[1/5,2/5), +1 on[2/5,1/2]',
         'first_state_cut':(abox*eta).output(),'second_state_raise':kbox.output(),
         'excess_profit':excess.output(),'actual_profit':actual.output(),
         'excess_minus_coefficient_eta_squared':remainder.output(),
         'one_point_cubic_scaled_remainder':scaled.output()},
      'uniform_upper_remainder_status':'Not certified by this script. Requires the separate full-menu localization/support-tangent proof.'}
    (ROOT/'stage3_coefficient_certificates.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'checks_passed':len(checks),'C1':brief(c1box),'C2':brief(c2box),
        'optimal_contrast':brief(abox),'optimal_tangent_coefficient':brief(cbox),
        'global_dual_multiplier':brief(dualbox),
        'finite_witness':{'eta':str(eta),'actual_profit':brief(actual),'excess':brief(excess),
           'one_point_cubic_scaled_remainder':brief(scaled)}},indent=2))

if __name__=='__main__':main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
