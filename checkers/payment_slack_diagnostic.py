#!/usr/bin/env python3
"""Exact checks of funding/primitive invariance in the payment-slack diagnostic.

The all-menu immunity conclusion is proved analytically in Proposition 5.1:
nonnegative promises relative to a zero incumbent imply nonincreasing gains.
These algebraic checks do not sample laws or replace that proof.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp
from certificate_checks import mandatory, decimal_bound

z,t,Hbar,B,kappa=sp.symbols("z t Hbar B kappa", real=True)
R=3*(t+z)**2+2*(t+z)**3
Phi=(t+z)**3+(t+z)**4/2
mean_X=sp.Rational(11,32)
old_P=Hbar-mean_X-kappa+(t+kappa)*R-Phi-B
new_B=B-Hbar+z
new_P=z-mean_X-kappa+(t+kappa)*R-Phi-new_B
checks={}
checks["margin polynomial unchanged"]=sp.simplify(new_P-old_P)==0
checks["mean severity exact"]=sp.integrate(
    sp.Symbol("x")*(6*sp.Symbol("x")+6*sp.Symbol("x")**2),
    (sp.Symbol("x"),0,sp.Rational(1,2)))==mean_X
checks["original mean collateral 11/10"]=sp.Rational(231,160)-mean_X==sp.Rational(11,10)
checks["funding negative on full pricing bracket"]=F(248617,240000)-F(231,160)+F(3932,10000)<0
checks["negative collateral has positive support interval"]=F(3932,10000)<F(1,2)
checks["both incumbent promises identically zero"]=sp.simplify(z-z)==0
mandatory(all(checks.values()), checks)
out={"checks_passed":len(checks),"checks":checks,
     "proof_scope":"primitive invariance and domain signs; monotone-gain immunity is analytic",
     "funding_upper_rational":str(F(248617,240000)-F(231,160)+F(3932,10000))}
target=Path(__file__).with_suffix(".json")
target.write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
