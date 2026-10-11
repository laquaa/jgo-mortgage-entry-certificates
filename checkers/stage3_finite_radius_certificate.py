#!/usr/bin/env python3
"""Exact-rational coefficient and small-radius witness certificate.

Requires Python's standard library only. It uses the uniform mathematical
constants proved in Theorem 4.1 and Appendix B of the manuscript. The supplied
certificate_checks.py helper is required. No search,
simulation, observed data or assumptions about strict choices from mere envelope
distance enter this script. The actual offer's root and selection proof are the
root-compensation and normalized-slope inequalities in that theorem.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from certificate_checks import mandatory, decimal_bound

sys.set_int_max_str_digits(0)


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    @staticmethod
    def point(v):
        return I(F(v), F(v))

    def __add__(self, v):
        v = v if isinstance(v, I) else I.point(v)
        return I(self.lo + v.lo, self.hi + v.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, v):
        return self + (-v if isinstance(v, I) else -I.point(v))

    def __rsub__(self, v):
        return I.point(v) - self

    def __mul__(self, v):
        v = v if isinstance(v, I) else I.point(v)
        products = [self.lo*v.lo, self.lo*v.hi,
                    self.hi*v.lo, self.hi*v.hi]
        return I(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, v):
        v = v if isinstance(v, I) else I.point(v)
        mandatory(v.lo > 0 or v.hi < 0, 'stage3_finite_radius_certificate.py:54 mandatory v.lo > 0 or v.hi < 0')
        return self * I(1/v.hi, 1/v.lo)

    def __pow__(self, n):
        result = I.point(1)
        for _ in range(n):
            result = result * self
        return result

    def record(self):
        return dict(lower_fraction=str(self.lo), upper_fraction=str(self.hi),
                    lower_decimal=decimal_bound(self.lo,False,90),
                    upper_decimal=decimal_bound(self.hi,True,90))


def isolate(function, low, high, digits, derivative_floor=None):
    low, high = F(low), F(high)
    def sign(x):
        v = function(I.point(x))
        if v.lo > 0:
            return 1
        if v.hi < 0:
            return -1
        raise AssertionError('Refine nuisance root intervals before isolation.')
    left = sign(low)
    mandatory(sign(high) == -left, 'stage3_finite_radius_certificate.py:85 mandatory sign(high) == -left')
    tolerance = F(1, 10**digits)
    while high-low > tolerance:
        mid = (low+high)/2
        value = function(I.point(mid))
        if value.lo <= 0 <= value.hi:
            # Nuisance enclosures can cover zero before bisection finishes.
            # A proved absolute derivative floor gives a rigorous root box.
            mandatory(derivative_floor is not None, 'stage3_finite_radius_certificate.py:93 mandatory derivative_floor is not None')
            radius = max(abs(value.lo), abs(value.hi))/derivative_floor
            enclosure = I(max(low, mid-radius), min(high, mid+radius))
            mandatory(enclosure.hi-enclosure.lo <= tolerance, (
                'Nuisance precision too coarse', float(enclosure.hi-enclosure.lo),
                float(tolerance)))
            return enclosure
        if (1 if value.lo > 0 else -1) == left:
            low = mid
        else:
            high = mid
    return I(low, high)


D, a, b, kap = F(1, 2), F(1, 20), F(1, 10), F(1, 50)
B = F(248617, 240000)
width = b-a


def R(y):
    return 3*y**2+2*y**3


def Rp(y):
    return 6*y+6*y**2


def Rpp(y):
    return 6+12*y


def Phi(y):
    return y**3+y**4/2


def pricing(z):
    return 2400000*z**4+4608000*z**3-415200*z**2-139560*z-218209


def main():
    checks = []
    def need(label, condition):
        mandatory(bool(condition), label)
        checks.append(label)

    z = isolate(pricing, F(393, 1000), F(394, 1000), 105)
    def margin(t):
        y = t+z
        return F(11, 10)-kap-B+(t+kap)*R(y)-Phi(y)
    root = isolate(margin, F(76, 1000), F(77, 1000), 80)
    yr, Y = root+z, b+z
    Rr, RB = R(yr), R(Y)
    wr, rr = (yr-D)/Rr, Rp(yr)/Rr
    integral_R = Phi(Y)-Phi(yr)
    integral_A = ((Y-D)**2-(yr-D)**2)/2
    C1 = ((b+kap)*(Y-D-wr*RB)
          -2*(integral_A-wr*integral_R))/width
    C2 = ((b+kap)*(Rp(Y)-rr*RB)
          -2*(RB-Rr-rr*integral_R))/width
    astar = -C1/C2
    cstar = -(C1**2)/(2*C2)
    lam = ((b+kap)*RB-2*integral_R)/width
    def primitive_P(y):
        return (F(11,10)-kap-B)*y+y**4/2+F(3,10)*y**5+(kap-z)*Phi(y)
    floor = (primitive_P(Y)-primitive_P(yr))/width
    need('pricing bracket below .3932', z.hi < F(3932, 10000))
    need('positive C1', C1.lo > 0)
    need('negative C2', C2.hi < 0)
    need('positive global dual multiplier', lam.lo > 0)
    need('.13 < optimal contrast < .14', F(13,100) < astar.lo < astar.hi < F(14,100))
    def critical(y):
        return (3-4*y**2)/(3*y*(3+4*y+2*y**2))
    need('monotonicity threshold exceeds .25', critical(Y).lo > F(1,4))
    need('strict normalized leading slope lower .05',
         F(13,100)*(1-F(14,100)/F(1,4)) > F(1,20))

    eta, epsilon = F(1,10**12), F(2,10**30)
    need('within certified finite-radius range', 0 < eta <= F(1,10**8))
    s = astar*eta
    def gain(y, k):
        h = (s-k)/2
        q2 = (s**2+k**2)/4
        q3 = (s**3-k**3)/12
        q4 = (s**4+k**4)/48
        return (h*R(y)+q2*Rp(y)+q3*Rpp(y)+12*q4
                +eta*(s+k)*(y-D)/2+eta*(s**2-k**2)/4)
    # Along this whole bracket, -partial_k g=R_2(yr-k)/2>3/8.
    need('compensation derivative exceeds3/8 throughout bracket',
         R(yr-(s.hi+F(14,100)*eta**2)).lo > F(3,4))
    k = isolate(lambda value: gain(yr, value), s.lo,
                s.hi+F(14,100)*eta**2, 70, F(3,8))
    need('compensating raise is positive', k.lo > 0)
    need('compensation within proved .14 eta squared bracket',
         k.lo >= s.lo and k.hi <= s.hi+F(14,100)*eta**2)
    need('all entrant indices inside tail perturbation buffer',
         a+z.lo-k.hi > F(2,5) and b+z.hi+s.hi < D)
    need('both promised payments nonnegative',
         F(43,32)-z.hi-s.hi > 0 and F(247,160)-z.hi+k.lo > 0)
    need('balanced direction densities nonnegative', eta < F(36,25))

    h = (s-k)/2
    q2, q3, q4 = (s**2+k**2)/4, (s**3-k**3)/12, (s**4+k**4)/48
    def gain_primitive(y):
        return (h*Phi(y)+q2*R(y)+q3*Rp(y)+q4*Rpp(y)
                +eta*(s+k)*(y**2/2-D*y)/2
                +eta*(s**2-k**2)*y/4)
    actual_excess = ((b+kap)*gain(Y,k)
                    -2*(gain_primitive(Y)-gain_primitive(yr)))/width
    actual = floor+actual_excess
    upper = floor+cstar*eta**2+1800000*eta**3
    theorem_lower = floor+cstar*eta**2-2*eta**3
    gap = upper-actual
    need('actual offer positive excess', actual_excess.lo > 0)
    need('actual profit consistent with certified lower theorem',
         actual.lo > theorem_lower.hi)
    need('global certificate gap below requested epsilon', gap.hi < epsilon)
    need('normalized coefficient within cubic upper remainder',
         (upper-floor).hi/eta**2-cstar.lo <= F(1800001)*eta)

    # These validate the rational constant assembly, not sampled outcomes.
    need('quadratic corner defect below400', 6+6/(4*F(6,1000)) < 400)
    need('derivative localization below220000',
         (400*400+6*17+4*20)/F(3,4) < 220000)
    need('common mean below32000',
         (23000+17+F(9,4)*400+F(1,2)*20)/F(3,4) < 32000)
    need('state means within1e-5 supporting radius',
         21*F(1,10**8) < F(1,10**5))
    need('mean Taylor remainder below44800',
         F(32,6)*(F(20)+32000*F(1,10**8))*400 < 44800)
    need('law Taylor remainder below534', 400/F(3,4) < 534)
    need('uniform normalized expansion below46000', 44800+534 < 46000)
    need('root positive-part constant below630000',
         8*220000**2*5/(F(1,4)*20) < 630000**2)
    need('assembled upper remainder below1800000',
         F(12,5)*(630000+46000)+2*46000 < 1800000)
    need('lower utility remainder below.245',
         F(15,1000)*F(21,2)+F(5,10000)*28+F(2,100000)*28
         +F(7,100)*F(3,5)+F(1,100)*3 < F(245,1000))
    need('lower derivative remainder below.83',
         F(15,1000)*39+F(5,10000)*84+F(2,100000)*72
         +F(7,100)*2+F(1,100)*6 < F(83,100))

    output = dict(scope='Synthetic exact balanced-density robustness supremum; full menus/lotteries upper; actual one-item lower.',
                  eta=str(eta), eta0='1/100000000', epsilon=str(epsilon),
                  uniform_upper_cubic_constant=1800000,
                  uniform_lower_cubic_constant=2, checks=checks, count=len(checks),
                  z=z.record(), zero_margin_type=root.record(),
                  C1=C1.record(), C2=C2.record(), optimal_contrast=astar.record(),
                  coefficient=cstar.record(), dual_multiplier=lam.record(),
                  floor=floor.record(), first_state_cut=s.record(),
                  second_state_raise=k.record(), actual_excess=actual_excess.record(),
                  actual_lower_profit=actual.record(), all_class_upper_profit=upper.record(),
                  gap=gap.record(),
                  warning='Constants are proved analytic inequalities, not fitted numerical remainders. This does not infer the actual law or provide statistical coverage.')
    path=Path(__file__).resolve().with_name('stage3_finite_radius_certificate.json')
    path.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks), eta=str(eta),
                         coefficient=cstar.record()['lower_decimal'],
                         actual_excess=actual_excess.record()['lower_decimal'],
                         global_gap=gap.record()['upper_decimal'],
                         requested_epsilon=str(epsilon), output=str(path)),indent=2))


if __name__ == '__main__':
    main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
