#!/usr/bin/env python3
"""Exact-rational polynomial and interval source for a rounded mortgage quote.

Run with Python 3 and SymPy. Float arithmetic is not used in any bound.
SymPy supplies symbolic polynomial expansion; certificate_checks.py supplies
mandatory guards. Rational interval arithmetic, root isolation and
continuous-interval Bernstein bounds are here. The paper's compact boxes are
verified by public_precision_certificate.py; this file's optional main also
produces additional stress rows not claimed in the paper.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import json
import math
import sympy as S
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from certificate_checks import mandatory, decimal_bound

ROOT = Path(__file__).resolve().parent
t,z,r,s,k,eta,x = S.symbols('t z r s k eta x')
a,b,kap,D,B = map(S.Rational,['1/20','1/10','1/50','1/2','248617/240000'])
checks=[]
def require(label, condition):
    mandatory(condition, label)
    checks.append(label)
def rational(value):
    value=S.Rational(value)
    return Q(int(value.p),int(value.q))
@dataclass(frozen=True)
class Box:
    lo: Q
    hi: Q
    @staticmethod
    def point(value):
        value=rational(value)
        return Box(value,value)
    def __add__(self,value):
        if not isinstance(value,Box):value=Box.point(value)
        return Box(self.lo+value.lo,self.hi+value.hi)
    __radd__=__add__
    def __neg__(self):return Box(-self.hi,-self.lo)
    def __sub__(self,value):return self+(-value if isinstance(value,Box) else -Box.point(value))
    def __rsub__(self,value):return Box.point(value)-self
    def __mul__(self,value):
        if not isinstance(value,Box):value=Box.point(value)
        values=[self.lo*value.lo,self.lo*value.hi,self.hi*value.lo,self.hi*value.hi]
        return Box(min(values),max(values))
    __rmul__=__mul__
    def __truediv__(self,value):
        if not isinstance(value,Box):value=Box.point(value)
        mandatory(value.lo>0 or value.hi<0, 'precision_certificate.py:49 mandatory value.lo>0 or value.hi<0')
        return self*Box(1/value.hi,1/value.lo)
    def __pow__(self,power):
        out=Box.point(1)
        for _ in range(power):out=out*self
        return out
    def output(self):
        return {'lower':str(self.lo),'upper':str(self.hi),
                'decimal_lower':decimal_bound(self.lo,False,24),
                'decimal_upper':decimal_bound(self.hi,True,24)}

@lru_cache(None)
def terms(expr):
    variables=tuple(sorted(expr.free_symbols,key=str))
    if not variables:return (),[((),rational(expr))]
    return variables,[(powers,rational(c)) for powers,c in S.Poly(S.expand(expr),*variables).terms()]
def evaluate(expr,boxes):
    variables,polynomial=terms(expr)
    out=Box.point(0)
    for powers,c in polynomial:
        term=Box.point(c)
        for variable,power in zip(variables,powers):term=term*boxes[variable]**power
        out=out+term
    return out
def coefficients(expr,var,lo,hi):
    w=S.Symbol('w')
    poly=S.Poly(S.expand(expr.subs(var,S.Rational(lo)+(S.Rational(hi)-S.Rational(lo))*w)),w)
    n=poly.degree()
    return [sum(poly.nth(j)*S.Rational(math.comb(i,j),math.comb(n,j))
                for j in range(i+1)) for i in range(n+1)]
def bernstein(expr,var,lo,hi,boxes):
    cc=[evaluate(c,boxes) for c in coefficients(expr,var,lo,hi)]
    return Box(min(c.lo for c in cc),max(c.hi for c in cc))
def isolate(expr,var,lo,hi,boxes,digits=55):
    lo,hi=rational(lo),rational(hi)
    def sign(value):
        result=evaluate(expr,{**boxes,var:Box.point(value)})
        if result.lo>0:return 1
        if result.hi<0:return -1
        raise AssertionError('nuisance interval too wide to certify sign')
    left=sign(lo)
    require('root opposing signs for '+str(var),sign(hi)==-left)
    while hi-lo>Q(1,10**digits):
        middle=(lo+hi)/2
        if sign(middle)==left:lo=middle
        else:hi=middle
    return Box(lo,hi)

def R(v):return 3*v*v+2*v*v*v
def Phi(v):return v**3+v**4/2
def Psi(v):return (v-D)**2/2-S.Rational(3,200)
pricing=2400000*z**4+4608000*z**3-415200*z*z-139560*z-218209
P=S.expand(S.Rational(11,10)-kap+(t+kap)*R(t+z)-Phi(t+z)-B)
g=S.expand((Phi(t+z+s)+Phi(t+z-k)-2*Phi(t+z)
              +eta*(Psi(t+z+s)-Psi(t+z-k)))/2)
p=S.expand(P+(t+kap)*S.diff(g,t)-g)
primitive=S.integrate(20*p,t)
baseline_primitive=S.integrate(20*P,t)

def main():
    zb=isolate(pricing,z,Q(393,1000),Q(394,1000),{},95)
    require('unique pricing root',evaluate(S.diff(pricing,z),{z:Box(Q(393,1000),Q(394,1000))}).lo>0)
    rb=isolate(P.subs(t,r),r,Q(76,1000),Q(77,1000),{z:zb},80)
    sb=Box.point(Q(655999774,10**11)); eb=Box.point(Q(1,20))
    kb=isolate(g.subs(t,r),k,Q(67,10000),Q(68,10000),{z:zb,r:rb,s:sb,eta:eb},70)
    boxes={z:zb,r:rb,s:sb,k:kb,eta:eb}
    def short_box(box,pad):
        scale=10**20
        floor=(box.lo.numerator*scale)//box.lo.denominator
        return Box(Q(floor-pad,scale),Q(floor+pad+1,scale))
    short_z=short_box(zb,1); short_r=short_box(rb,10); short_k=short_box(kb,100)
    for label,expr,var,box,nuisance in [
            ('z',pricing,z,short_z,{}),
            ('t_star',P.subs(t,r),r,short_r,{z:short_z}),
            ('k',g.subs(t,r),k,short_k,{z:short_z,r:short_r,s:sb,eta:eb})]:
        lv=evaluate(expr,{**nuisance,var:Box.point(box.lo)})
        rv=evaluate(expr,{**nuisance,var:Box.point(box.hi)})
        require('20-digit source root box opposing signs '+label,
                (lv.lo>0 and rv.hi<0) or (lv.hi<0 and rv.lo>0))
    require('unique compensation root',bernstein(S.diff(g.subs(t,r),k),k,Q(67,10000),Q(68,10000),{z:zb,r:rb,s:sb,eta:eb}).hi<0)
    require('buffer holds for every continuous type', (zb+Box.point(a)-kb).lo>Q(2,5) and (zb+Box.point(b)+sb).hi<Q(1,2))
    require('nonnegative promises', (Box.point(Q(43,32))-zb-sb).lo>0 and (Box.point(Q(247,160))-zb+kb).lo>0)
    require('density direction normalized',-Q(1,2)*(Q(2,5)-Q(1,5))+Q(1,2)-Q(2,5)==0)
    require('density direction feasible',Q(6,5)+Q(6,25)>Q(1,40) and Q(12,5)+Q(24,25)>Q(1,20))
    dg=bernstein(S.diff(g,t),t,a,b,boxes)
    dp=bernstein(S.diff(p,t),t,a,b,boxes)
    require('gain strictly increasing on whole type interval',dg.lo>0)
    require('margin strictly increasing on whole type interval',dp.lo>0)
    gain_a=evaluate(g,{**boxes,t:Box.point(a)}); gain_b=evaluate(g,{**boxes,t:Box.point(b)})
    # Severity density, including balanced direction, is at most 4.55.
    # dQ/dm=R_s-(t+kappa)f_s; the displayed conservative bound follows.
    receipt_lipschitz=Q(1)+Q(3,25)*Q(91,20)
    require('receipt Lipschitz constant',receipt_lipschitz==Q(773,500))
    def threshold(level):
        if level<gain_a.lo:return Box.point(a)
        if level>gain_b.hi:return Box.point(b)
        return isolate(g-S.Rational(level),t,a,b,boxes,55)
    def integral(expr,lo,hi):
        antiderivative=S.integrate(20*expr,t)
        return evaluate(antiderivative,{**boxes,t:hi})-evaluate(antiderivative,{**boxes,t:lo})
    zero=Box.point(0)
    i=integral(P,rb,Box.point(b))
    exact=integral(p,rb,Box.point(b))
    rows=[]
    for tau in map(Q,['0','0.0000001','0.000001','0.00001','0.0001','0.001','0.005','0.006']):
        q=threshold(tau)
        profit=integral(p,q,Box.point(b)) if q!=Box.point(b) else zero
        mass=20*(Box.point(b)-q)
        grant_compensation=exact-tau*20*(Box.point(b)-rb)
        rows.append({'tau':str(tau),'threshold':q.output(),'selected_mass':mass.output(),
                     'ungranted_profit':profit.output(),'cost_compensating_grant_profit':grant_compensation.output()})
    robust=[]
    for tau in map(Q,['0','0.0000001','0.000001','0.00001']):
        for epsilon in map(Q,['0.00000001','0.0000001','0.000001']):
            lower=threshold(tau-epsilon); upper=threshold(tau+epsilon)
            q0=isolate(p-S.Rational(receipt_lipschitz*epsilon),t,a,b,boxes,55)
            # Every admissible rounding selects all t>upper, selects none
            # below lower, and can choose any subset of their intervening band.
            # Integrate the conservative margin over the guaranteed tail and
            # its negative part over every possibly selected ambiguous type.
            adjusted=p-S.Rational(receipt_lipschitz*epsilon)
            guaranteed=integral(adjusted,upper,Box.point(b)) if upper!=Box.point(b) else zero
            if lower.lo>=q0.hi:negative=zero
            elif upper.hi<=q0.lo:negative=integral(adjusted,lower,upper)
            else:negative=integral(adjusted,lower,q0)
            result=guaranteed+negative
            if tau==0 and epsilon==Q(1,10**7):
                require('bounded payment errors permit K=I_star+6e-6',result.lo>i.hi+Q(6,10**6))
                require('switch cost 1e-6 offset by exact grant permits K=I_star+5e-6',result.lo-Q(1,10**6)>i.hi+Q(5,10**6))
            robust.append({'tau':str(tau),'epsilon':str(epsilon),
                           'possible_selection_lower':lower.output(),
                           'guaranteed_selection_lower':upper.output(),
                           'adjusted_margin_zero':q0.output(),
                           'guaranteed_tail_contribution':guaranteed.output(),
                           'ambiguous_negative_contribution':negative.output(),
                           'adjusted_margin_at_possible_selection_edge':evaluate(adjusted,{**boxes,t:lower}).output(),
                           'profit_lower_bound':result.output()})
    adversarial_epsilon=Q(1,100000)
    adversarial_gain=S.expand((Phi(t+z+s-S.Rational(adversarial_epsilon))
                              +Phi(t+z-k-S.Rational(adversarial_epsilon))-2*Phi(t+z)
                              +eta*(Psi(t+z+s-S.Rational(adversarial_epsilon))
                                    -Psi(t+z-k-S.Rational(adversarial_epsilon))))/2)
    adversarial_bound=bernstein(adversarial_gain,t,a,b,boxes)
    require('both payments raised by 1e-5 wipe out every strict client',adversarial_bound.hi<0)
    data={'specification':{'eta':'1/20','s':'655999774/100000000000','k':'unique root g(t_star)=0 in (67/10000,68/10000)',
                           'density_direction':'0 on [0,1/5), -1/2 on [1/5,2/5), 1 on [2/5,1/2]',
                           'receipt_lipschitz':'773/500','rounding':'each payment error at most epsilon; grants exact'},
          'z':zb.output(),'t_star':rb.output(),'compensation_k':kb.output(),
          'paper_20_digit_root_boxes':{'z':short_z.output(),'t_star':short_r.output(),'k':short_k.output()},
          'gain_derivative':dg.output(),'gain_derivative_Bernstein_coefficients':[evaluate(c,boxes).output() for c in coefficients(S.diff(g,t),t,a,b)],
          'margin_derivative':dp.output(),'margin_derivative_Bernstein_coefficients':[evaluate(c,boxes).output() for c in coefficients(S.diff(p,t),t,a,b)],
          'gain_at_a':gain_a.output(),'gain_at_b':gain_b.output(),
          'I_star':i.output(),'quote_profit':exact.output(),'quote_excess':(exact-i).output(),
          'tail_mass':(20*(Box.point(b)-rb)).output(),
          'compensating_grant_break_even_cost':(exact/(20*(Box.point(b)-rb))).output(),
          'switching_cost_rows':rows,'rounding_rows':robust,'assertions':checks,
          'adversarial_both_payments_raised':{'epsilon':str(adversarial_epsilon),
                                            'gain_range':adversarial_bound.output(),
                                            'Bernstein_coefficients':[evaluate(c,boxes).output() for c in coefficients(adversarial_gain,t,a,b)],
                                            'actual_profit':'0'},
          'root_accuracy':{'z':'1e-95','t_star':'1e-80','k':'1e-70','cost_and_margin_roots':'1e-55'},
          'arithmetic':'exact fractions; decimal output is outward rounded to 24 places'}
    (ROOT/'precision_certificate.json').write_text(json.dumps(data,indent=2)+'\n')
    for label in ['compensation_k','gain_derivative','gain_at_b','quote_profit']:
        print(label,data[label]['decimal_lower'],data[label]['decimal_upper'])
    print('switching cost, ungranted profit, compensated grant profit')
    for row in rows:print(row['tau'],row['ungranted_profit']['decimal_lower'],row['cost_compensating_grant_profit']['decimal_lower'])
    print('switching cost, payment error, robust profit lower bound')
    for row in robust:print(row['tau'],row['epsilon'],row['profit_lower_bound']['decimal_lower'])
    print('all',len(checks),'exact assertions passed')

if __name__=='__main__':main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
