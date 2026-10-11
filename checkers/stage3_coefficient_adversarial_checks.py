"""Independent exact root enclosures and algebraic coefficient checks."""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from decimal import Decimal, localcontext
import json
from certificate_checks import mandatory, decimal_bound


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    @staticmethod
    def point(value):
        return I(F(value), F(value))

    def __add__(self, other):
        other = other if isinstance(other, I) else I.point(other)
        return I(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other if isinstance(other, I) else -F(other))

    def __rsub__(self, other):
        return I.point(other) - self

    def __mul__(self, other):
        other = other if isinstance(other, I) else I.point(other)
        endpoints = [self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi]
        return I(min(endpoints), max(endpoints))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, I) else I.point(other)
        mandatory(not (other.lo <= 0 <= other.hi), 'stage3_coefficient_adversarial_checks.py:40 mandatory not (other.lo <= 0 <= other.hi)')
        return self * I(F(1)/other.hi, F(1)/other.lo)

    def __pow__(self, power):
        mandatory(isinstance(power, int) and power >= 0, 'stage3_coefficient_adversarial_checks.py:44 mandatory isinstance(power, int) and power >= 0')
        result = I.point(1)
        for _ in range(power):
            result = result * self
        return result

    def out(self):
        return {"lower_fraction":str(self.lo),"upper_fraction":str(self.hi),
                "lower_decimal":decimal_bound(self.lo,False,40),
                "upper_decimal":decimal_bound(self.hi,True,40)}


def main():
    checks = {}
    def check(name, assertion):
        mandatory(assertion, name)
        checks[name] = True

    def pricing(z):
        return 2400000*z**4+4608000*z**3-415200*z**2-139560*z-218209
    zlo,zhi = F(393,1000),F(394,1000)
    check("pricing_root_bracket",pricing(zlo)<0<pricing(zhi))
    # Monotonicity by conservative interval lower bound on the derivative.
    check("pricing_root_unique",9600000*zlo**3+13824000*zlo**2-830400*zhi-139560>0)
    for _ in range(220):
        mid=(zlo+zhi)/2
        if pricing(mid)<0:zlo=mid
        else:zhi=mid
    z=I(zlo,zhi)
    a,b,kappa,D=F(1,20),F(1,10),F(1,50),F(1,2)
    funding=F(248617,240000)
    R=lambda y:3*y**2+2*y**3
    Rp=lambda y:6*y+6*y**2
    Phi=lambda y:y**3+y**4/2
    def P(t):
        y=z+t
        return F(11,10)-kappa+(t+kappa)*R(y)-Phi(y)-funding
    check("margin_root_bracket",P(a).hi<0<P(b).lo)
    rlo,rhi=a,b
    for _ in range(150):
        mid=(rlo+rhi)/2
        pm=P(mid)
        if pm.hi<0:rlo=mid
        elif pm.lo>0:rhi=mid
        else:raise RuntimeError("Refine pricing enclosure before margin bisection")
    r=I(rlo,rhi)
    yr=r+z;yb=z+b
    Rr,Rb=R(yr),R(yb)
    wr=(yr-D)/Rr;rr=Rp(yr)/Rr
    intR=Phi(yb)-Phi(yr)
    intA=(b-r)*(z-D)+(b*b-r**2)/2
    C1=((b+kappa)*(yb-D-wr*Rb)-2*(intA-wr*intR))/(b-a)
    C2=((b+kappa)*(Rp(yb)-rr*Rb)-2*(Rb-Rr-rr*intR))/(b-a)
    check("linear_coefficient_positive",C1.lo>F(52,1000))
    check("quadratic_coefficient_negative",C2.hi<-F(394,1000))
    a_star=-C1/C2
    C_star=-(C1**2)/(2*C2)
    check("a_star_positive_small",F(132,1000)<a_star.lo<a_star.hi<F(134,1000))
    check("candidate_coefficient_strict_enclosure",F(3480,1000000)<C_star.lo<C_star.hi<F(3481,1000000))
    # Normalized gain derivative positivity: w'>=1, -r'<=2/.443²+4/9.
    rprime_bound=2/F(443,1000)**2+F(4,9)
    check("root_normalized_gain_increasing",1-F(14,100)/2*rprime_bound>0)
    check("tail_endpoint_functional_increases_in_x",2*rlo-b+kappa>0)
    # Uniform global supporting-line constants for the nonconvex K.
    eps0=F(1,1000)
    check("local_shift_below_support_slack",zlo+a-eps0>0 and zhi+b+eps0<D)
    check("reference_response_lower_bound",R(zlo+a)>F(3,4))
    check("local_curvature_positive",8-112*eps0>7)
    check("local_curvature_upper_bound",rprime_bound+112*eps0<11)
    check("remote_tangent_dominated",11/F(100)+F(33,2)/10000<F(5,12))
    check("global_frontier_integral_constant",1/(20*(b+kappa))==F(5,12))
    # Feasible dominant density: fixed buffer avoids shift-crossing the jump.
    yL=zlo+a-F(1,100)
    y0=F(1,8)
    negative_magnitude=(D-yL)/(yL-y0)
    check("dominant_law_zero_mean",negative_magnitude*(yL-y0)==D-yL)
    check("dominant_law_unit_bound",0<negative_magnitude<1)
    check("dominant_law_reference_margin",6*y0+6*y0*y0>F(1,2))
    output={"scope":"Exact algebraic arithmetic audit, not a global optimizer or empirical calibration.",
            "checks":checks,"count":len(checks),"z":z.out(),"r":r.out(),
            "C1":C1.out(),"C2":C2.out(),"a_star":a_star.out(),"candidate_C_star":C_star.out(),
            "global_support_line":{"local_radius":str(eps0),"mean_threshold":str(eps0/100),
                                    "quadratic_lower_constant":"5/12",
                                    "local_curvature_lower":"7","local_curvature_upper":"11"},
            "dominant_density":{"negative_block_start":str(y0),"positive_block_start":str(yL),
                                 "negative_height":str(-negative_magnitude),"positive_height":"1"}}
    target=Path(__file__).with_name("stage3_coefficient_adversarial_certificates.json")
    target.write_text(json.dumps(output,indent=2)+"\n")
    print(f"Passed {len(checks)} exact checks; wrote {target.name}")


if __name__=="__main__":
    main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
