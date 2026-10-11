#!/usr/bin/env python3
"""Exact rational interval utilities for the benchmark coefficient checks.

Root signs and continuous-interval Bernstein coefficients are enclosed using
exact fractions. This is a supplied source dependency, with no data-file inputs.
"""
from dataclasses import dataclass
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import hashlib
import json
import math
import sympy as sp
from certificate_checks import mandatory, decimal_bound

ROOT = Path(__file__).resolve().parent
t, z, r, k, q, y = sp.symbols('t z r k q y', real=True)
a, b = sp.Rational(1, 20), sp.Rational(1, 10)
kap, B, theta = sp.Rational(1, 50), sp.Rational(248617, 240000), sp.Rational(1, 4)
D, lam = sp.Rational(1, 2), sp.Rational(1, 10)
checks = []

def require(label, value):
    mandatory(bool(value), label)
    checks.append(label)

def frac(value):
    value = sp.Rational(value)
    return Fraction(int(value.p), int(value.q))

@dataclass(frozen=True)
class Box:
    lo: Fraction
    hi: Fraction
    @staticmethod
    def point(x):
        x = frac(x)
        return Box(x, x)
    def __add__(self, x):
        if not isinstance(x, Box): x = Box.point(x)
        return Box(self.lo+x.lo, self.hi+x.hi)
    __radd__ = __add__
    def __neg__(self): return Box(-self.hi, -self.lo)
    def __sub__(self, x): return self + (-x if isinstance(x, Box) else -Box.point(x))
    def __rsub__(self, x): return Box.point(x)-self
    def __mul__(self, x):
        if not isinstance(x, Box): x = Box.point(x)
        v = [self.lo*x.lo, self.lo*x.hi, self.hi*x.lo, self.hi*x.hi]
        return Box(min(v), max(v))
    __rmul__ = __mul__
    def __truediv__(self, x):
        if not isinstance(x, Box): x = Box.point(x)
        mandatory(x.lo > 0 or x.hi < 0, 'interval_utils.py:54 mandatory x.lo > 0 or x.hi < 0')
        return self*Box(1/x.hi, 1/x.lo)
    def __pow__(self, degree):
        v = Box.point(1)
        for _ in range(degree): v = v*self
        return v
    def output(self):
        return {'rational_lower':str(self.lo), 'rational_upper':str(self.hi),
                'decimal_lower':decimal_bound(self.lo),
                'decimal_upper':decimal_bound(self.hi,True)}

def evaluate(expr, boxes):
    expr=sp.expand(expr)
    variables=tuple(sorted(expr.free_symbols,key=str))
    if not variables: return Box.point(expr)
    value=Box.point(0)
    for powers, coef in sp.Poly(expr,*variables).terms():
        term=Box.point(coef)
        for var,power in zip(variables,powers): term=term*boxes[var]**power
        value=value+term
    return value

def bisect(expr, var, lo, hi, boxes, digits):
    lo,hi=frac(lo),frac(hi)
    def sign(x):
        v=evaluate(expr,{**boxes,var:Box(x,x)})
        if v.lo>0:return 1
        if v.hi<0:return -1
        raise AssertionError('nuisance box wider than requested root isolation')
    left=sign(lo)
    require(str(var)+' root endpoints oppose',sign(hi)==-left)
    while hi-lo>Fraction(1,10**digits):
        mid=(lo+hi)/2
        if sign(mid)==left:lo=mid
        else:hi=mid
    return Box(lo,hi)

def bernstein(expr, var, lo, hi, boxes):
    w=sp.Symbol('w')
    p=sp.Poly(sp.expand(expr.subs(var,lo+(hi-lo)*w)),w)
    n=p.degree()
    coeff=[sum(p.nth(j)*sp.Rational(math.comb(i,j),math.comb(n,j))
               for j in range(i+1)) for i in range(n+1)]
    v=[evaluate(c,boxes) for c in coeff]
    return Box(min(c.lo for c in v),max(c.hi for c in v))

def response(x, th=theta):return (1-th)*x**2/D**2+th*x**3/D**3
def phi(x, th=theta):return (1-th)*x**3/(3*D**2)+th*x**4/(4*D**3)
def meanQ(x):return sp.Rational(11,10)-kap+(t+kap)*response(t+x)-phi(t+x)
def integral(expr, lo=a, hi=b):return sp.factor(sp.integrate(expr,(t,lo,hi))/(b-a))


# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
