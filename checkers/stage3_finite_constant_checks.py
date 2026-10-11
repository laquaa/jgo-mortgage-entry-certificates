"""Exact rational audit of the manuscript's uniform cubic constants.

No simulation or global optimizer is used. The mathematical reductions and
branch arguments appear in Theorem 4.1 and Appendix B of the manuscript.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import sympy as sp
from certificate_checks import mandatory, decimal_bound


def main():
    checks={}
    def require(name, condition):
        mandatory(condition, name)
        checks[name]=True
    eta0=Q(1,100000000)
    a,b,kappa,D=Q(1,20),Q(1,10),Q(1,50),Q(1,2)
    ymin,ymax=Q(443,1000),Q(1233,2500)
    R=lambda y:3*y*y+2*y**3
    F=lambda y:6*y+6*y*y
    Phi=lambda y:y**3+y**4/2
    rho=lambda y:F(y)/R(y)
    pricing=lambda z:2400000*z**4+4608000*z**3-415200*z*z-139560*z-218209
    require("stronger_pricing_bracket",pricing(Q(393,1000))<0<pricing(Q(983,2500)))
    require("Rmin_three_quarters",R(ymin)>Q(3,4))
    require("reverse_hazard_min_four_and_half",rho(ymax)>Q(9,2))
    require("reverse_hazard_max_six",2/ymin+Q(2,3)<6)
    require("beta_bound_quarter",1-6*(b+kappa)>Q(1,4))
    require("grant_bound_three_tenths",(a+kappa)*Q(9,2)>Q(3,10))
    require("paired_K_below_one",1+(a+kappa)*(1-D*Q(9,2))<=1)
    require("reference_density_max",F(D)==Q(9,2))
    require("moment_localization_constant",(1/Q(1,20))**2==400)
    require("grant_localization_constant",1/(4*Q(1,20)*Q(3,10))<17)
    require("positive_gain_integral_constant",1/(4*Q(1,20)*Q(1,4))==20)

    # Explicit polynomial anchor defect at the interior and both corners.
    y,x=sp.symbols("y x")
    interior=(3*x*x*y+3*x*x+2*x*y*y+6*x*y+6*x+y**3+3*y*y+3*y)/(y*(2*y+3))
    derivative_num=sp.factor(sp.diff(interior,x)*y*(2*y+3))
    require("interior_ratio_increases_in_shifted_threshold",sp.expand(derivative_num)==6*x*y+6*x+2*y*y+6*y+6)
    N=lambda y:-4*y**3+32*y*y+45*y-15
    require("interior_anchor_ratio_below_six",N(ymin)>0)
    require("interior_max_numerator_monotone",-12*D**2+64*ymin+45>0)
    require("below_support_anchor_ratio_below_six",15-ymax*ymax>0)
    slack=Q(6,1000)
    require("above_support_anchor_ratio_below_400",6+6/(4*slack)==256 and 256<400)
    defect=400
    osc=(defect*400+(1+6*D)*20+6*17)/Q(3,4)
    require("normalized_gain_derivative_constant",osc==213576 and osc<220000)
    terminal=220000*(b-a)+20/Q(3,4)
    require("terminal_normalized_gain_constant",terminal<11027)
    require("whole_gain_sup_constant",11027+220000*(b-a)<23000)
    mean=(23000+17+(Q(9,2)/2)*400+D*20)/Q(3,4)
    require("common_mean_constant",mean<32000)
    state_mean=20+32000*eta0
    require("state_means_sharp_bound",state_mean<Q(201,10))
    require("state_means_coarse_bound",state_mean<21)
    require("state_means_within_support_tangent",21*eta0<Q(1,100000))
    require("selection_regret_constant",1/(4*Q(1,20))==5)
    require("margin_type_derivative_min",(a+kappa)*F(ymin)>Q(1,4))
    require("positive_window_not_full_lower_interval",Q(5,2)*(Q(3,40)-a)**2>5*eta0**2)
    require("root_normalized_gain_constant",2*(2*220000)**2<630000**2)

    # Reference third-order term can use sharper means or averaged moments.
    require("reference_K_third_derivative_constant",2*12/Q(3,4)==32)
    reference=Q(32,6)*Q(201,10)**3
    require("reference_Taylor_44800_using_sharp_means",reference<44800)
    require("loose_21_bound_is_insufficient",Q(32,6)*21**3==49392 and 49392>44800)
    require("reference_Taylor_44800_using_mean_second_moment",Q(32,6)*21*400==44800)
    law=Q(400)/Q(3,4)
    require("law_normalized_endpoint_remainder",law==Q(1600,3) and law<534)
    require("normalized_combined_remainder",44800+534<46000)
    require("root_H_constant",630000+46000==676000)
    require("dual_multiplier_constant",20*(b+kappa)==Q(12,5))
    require("integrated_upper_weight_bound",40*(b-a)==2)
    upper=Q(12,5)*676000+46000*2
    require("assembled_upper_constant",upper==1714400 and upper<1800000)

    # Lower exact root-compensated quote.
    funding=Q(248617,240000)
    margin=lambda t,z:Q(11,10)-kappa+(t+kappa)*R(t+z)-Phi(t+z)-funding
    require("stronger_margin_root_bracket",margin(Q(19,250),Q(393,1000))<0<margin(Q(77,1000),Q(983,2500)))
    require("root_compensation_start_positive",Q(1,2)*Q(13,100)**2*4-Q(14,100)*Q(31,1000)>0)
    require("root_compensation_start_upper",Q(1,2)*Q(14,100)**2*Q(9,2)+eta0**2*Q(14,100)**4*12/24<Q(5,100))
    require("root_compensation_v_cap",Q(3,8)*Q(14,100)>Q(5,100))
    max_shift=Q(14,100)*eta0+Q(14,100)*eta0**2
    require("lower_shifts_inside_saturated_buffer",ymin-max_shift>Q(2,5) and ymax+max_shift<D)
    require("lower_shifted_reference_response",R(ymin-max_shift)>Q(3,4))
    require("lower_densities_nonnegative",F(Q(3,10))>eta0)
    require("lower_state_promises_nonnegative",Q(43,32)-Q(983,2500)-max_shift>0)
    # Root-normalized profile bounds from R<=1, R>=3/4,
    # R'<=9/2, R''<=12 and |y_root-D|<=31/1000.
    require("lower_L2_profile_bounds",Q(9,2)+6==Q(21,2) and 12+6*Q(9,2)==39)
    require("lower_L3_profile_bounds",12+12/Q(3,4)==28 and 12+(12/Q(3,4))*Q(9,2)==84)
    require("lower_L4_profile_bounds",12+12/Q(3,4)==28 and (12/Q(3,4))*Q(9,2)==72)
    require("lower_L1_profile_bounds",D-ymin+Q(31,1000)/Q(3,4)<Q(6,10) and 1+Q(31,1000)*Q(9,2)/Q(3,4)<2)
    require("lower_L0_profile_bounds",1+1/Q(3,4)<3 and Q(9,2)/Q(3,4)==6)
    A,V=Q(14,100),Q(14,100)
    coeff2=A*V/2+V*V*eta0/4
    coeff3=A*A*V*eta0/4+A*V*V*eta0**2/4+V**3*eta0**3/12
    coeff4=(A**4+(A+V*eta0)**4)*eta0/48
    coeff1=V/2
    coeff0=A*V*eta0/2+V*V*eta0**2/4
    require("lower_L2_coefficient_error",coeff2<Q(15,1000))
    require("lower_L3_coefficient_error",coeff3<Q(5,10000))
    require("lower_L4_coefficient_error",coeff4<Q(2,100000))
    require("lower_L1_coefficient_error",coeff1==Q(7,100))
    require("lower_L0_coefficient_error",coeff0<Q(1,100))
    error=Q(15,1000)*Q(21,2)+Q(5,10000)*28+Q(2,100000)*28+Q(7,100)*Q(6,10)+Q(1,100)*3
    error_derivative=Q(15,1000)*39+Q(5,10000)*84+Q(2,100000)*72+Q(7,100)*2+Q(1,100)*6
    require("lower_gain_error",error==Q(12203,50000) and error<Q(245,1000))
    require("lower_gain_derivative_error",error_derivative==Q(20711,25000) and error_derivative<Q(83,100))
    critical=lambda y:(3-4*y*y)/(3*y*(3+4*y+2*y*y))
    require("critical_tilt_coarse_494_fails",critical(Q(247,500))<Q(1,4))
    require("critical_tilt_exact_frozen_bracket_passes",critical(ymax)>Q(1,4))
    critical_symbolic=(3-4*y*y)/(3*y*(3+4*y+2*y*y))
    critical_derivative=(2*y+3)*(4*y**3-6*y*y-6*y-3)/(3*y*y*(2*y*y+4*y+3)**2)
    require("critical_tilt_derivative_identity",sp.simplify(sp.diff(critical_symbolic,y)-critical_derivative)==0)
    require("critical_tilt_monotone",4*D**3-3<0)
    require("leading_normalized_gain_derivative",Q(13,100)*(1-Q(14,100)/Q(1,4))>Q(5,100))
    require("normalized_error_derivative",1/Q(3,4)+Q(9,2)/Q(3,4)**2==Q(28,3))
    require("actual_quote_normalized_gain_increases",Q(5,100)-Q(28,3)*eta0>0)
    require("lower_direct_profit_constant",1+b+kappa<Q(2))

    out={"scope":"Exact rational assertions supporting the assembled mathematical finite-radius proof, not simulations.",
         "eta0":str(eta0),"checks":checks,"count":len(checks),
         "upper_constants":{"normalized_gain_derivative":str(osc),"terminal":str(terminal),
                            "common_mean":str(mean),"state_mean_at_eta0":str(state_mean),
                            "reference_Taylor_using_20.1":str(reference),"law":str(law),
                            "aggregate_cubic":str(upper),"safe_cubic":"1800000"},
         "lower_constants":{"gain_error":str(error),"derivative_error":str(error_derivative),"safe_cubic":"2"},
         "required_documentation_corrections":[
             "44800 follows either from the sharpened 20.1 eta state means or from 21 eta means combined with their averaged squared moment <=400 eta^2; cubing 21 alone is insufficient.",
             "a_c>=1/4 uses y_b<0.4932 from z<0.3932; it fails at the coarse endpoint y=0.494."]}
    target=Path(__file__).with_name("stage3_finite_constant_certificates.json")
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(f"Passed {len(checks)} exact assertions; wrote {target.name}")


if __name__=="__main__":main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
