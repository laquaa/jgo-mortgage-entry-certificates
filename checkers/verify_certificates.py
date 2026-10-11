#!/usr/bin/env python3
"""Exact-rational certificates for the five manuscript witnesses.

Requires only Python 3's standard library. Run: python3 verify_certificates.py
Writes certificates.json next to this file; imports the supplied mandatory-check
helper certificate_checks.py and requires no data files or previous results.
This checks numerical inequalities. The manuscript proves all-law/menu coverage.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import json
from certificate_checks import mandatory, decimal_bound, exact_fraction

DEN = 10**65
CHECKS = []
DIVISOR_LOWER = None
def require(label, value):
    mandatory(value, label)
    CHECKS.append(label)
def down(v): return Q((v.numerator*DEN)//v.denominator, DEN)
def up(v): return -down(-v)

@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q
    def __post_init__(self):
        mandatory(not isinstance(self.lo,float) and not isinstance(self.hi,float) and self.lo<=self.hi,'valid rational interval endpoints')
    @staticmethod
    def p(v):
        v=exact_fraction(v)
        return I(v,v)
    def __add__(self, v):
        if not isinstance(v,I): v=I.p(v)
        return I(down(self.lo+v.lo),up(self.hi+v.hi))
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,v): return self+(-v if isinstance(v,I) else -I.p(v))
    def __rsub__(self,v): return I.p(v)-self
    def __mul__(self,v):
        if not isinstance(v,I): v=I.p(v)
        vals=[self.lo*v.lo,self.lo*v.hi,self.hi*v.lo,self.hi*v.hi]
        return I(down(min(vals)),up(max(vals)))
    __rmul__=__mul__
    def __truediv__(self,v):
        if not isinstance(v,I): v=I.p(v)
        mandatory(v.lo>0 or v.hi<0, 'division interval contains zero')
        return self*I(down(1/v.hi),up(1/v.lo))
    def __pow__(self,n):
        mandatory(isinstance(n,int) and n>=0,'nonnegative integral interval power')
        result=I.p(1)
        for _ in range(n): result=result*self
        return result
    def minimum(self,v): return I(min(self.lo,v.lo),min(self.hi,v.hi))
    def output(self):
        return {'rational_lower':str(self.lo),'rational_upper':str(self.hi),
                'decimal_lower':decimal(self.lo,35,False),
                'decimal_upper':decimal(self.hi,35,True)}

def decimal(v, digits, upward):
    return decimal_bound(v,upward,digits)

a,b,kap,D,B=map(Q,['1/20','1/10','1/50','1/2','248617/240000'])
H1,H2=Q(43,32),Q(247,160)
C0=Q(11,10)-kap-B
RADII=list(map(Q,['.0002','.002','.01','.02','.05']))
CUTS=list(map(Q,['.00002656585','.00026552461','.00132472248',
                 '.00264253147','.00655999774']))
TABLE_LOW=list(map(Q,['1.39209387e-10','1.39209623e-8','3.48038362e-7',
                      '1.39233232e-6','8.70992255e-6']))
TABLE_UP=list(map(Q,['1.47844077e-10','1.48121748e-8','3.73454566e-7',
                     '1.51018877e-6','9.77402061e-6']))
TABLE_WIDTH=list(map(Q,['8.635e-12','8.913e-10','2.542e-8','1.179e-7','1.065e-6']))
TABLE_DERIVATIVE=list(map(Q,['2.9113e-9','2.9113e-7','7.2786e-6','2.9118e-5','1.8213e-4']))
TABLE_M=list(map(Q,['.0078441568348087','.0078590288527393','.0079261818593246','.0080126478606944','.0082907495921007']))
TABLE_J=list(map(Q,['.0036961018020530','.0037030338034969','.0037342943558823','.0037744446697578','.0039027345910477']))

def R(y): return y*y*(3+2*y)
def Rp(y): return 6*y+6*y*y
def Rpp(y): return 6+12*y
def Phi(y): return y**3+y**4/2
def price(z): return 2400000*z**4+4608000*z**3-415200*z*z-139560*z-218209
def price_prime(z): return 9600000*z**3+13824000*z*z-830400*z-139560
def P(t,z):
    y=t+z
    return C0+(t+kap)*R(y)-Phi(y)
def Pprimitive(y,z):
    ell=kap-z
    return Q(3,10)*y**5+(1+ell)*y**4/2+ell*y**3+C0*y

def scalar_root(fun,left,right,width):
    """Opposing exact signs and a separately checked derivative give the root."""
    left,right=Q(left),Q(right)
    lval=fun(left)
    require('scalar root has opposing signs',lval*fun(right)<0)
    while right-left>width:
        mid=(left+right)/2
        if (fun(mid)>0)==(lval>0): left=mid
        else: right=mid
    return I(left,right)

def nuisance_root(fun,left,right,width,label):
    left,right=Q(left),Q(right)
    def sign(v):
        val=fun(I.p(v))
        if val.lo>0: return 1
        if val.hi<0: return -1
        raise AssertionError(label+': nuisance interval too wide')
    leftsign=sign(left)
    require(label+' has opposing signs',sign(right)==-leftsign)
    while right-left>width:
        mid=(left+right)/2
        if sign(mid)==leftsign: left=mid
        else: right=mid
    return I(left,right)

def root_decimal_hull(box,digits):
    """A short paper box containing a more accurate isolated root."""
    den=10**digits
    left=Q((box.lo.numerator*den)//box.lo.denominator,den)
    right=Q(-((-box.hi.numerator*den)//box.hi.denominator),den)
    return I(left,right)

def twice_gain(y,s,k,eta):
    # Exact quartic, with cancellation performed symbolically in this formula.
    return ((s-k)*R(y)+(s*s+k*k)*(3*y+3*y*y)
            +(s**3-k**3)*(1+2*y)+(s**4+k**4)/2
            +eta*((s+k)*(y-D)+(s*s-k*k)/2))

def gain_coefficients(s,k,eta):
    h=(s-k)/2
    q=(s*s+k*k)/4
    r=(s**3-k**3)/12
    u=(s**4+k**4)/48
    return [6*r+12*u-eta*D*(s+k)/2+eta*(s*s-k*k)/4,
            6*q+12*r+eta*(s+k)/2,3*h+6*q,2*h]

def eval_poly(coefs,y):
    result=I.p(0)
    for c in reversed(coefs): result=result*y+c
    return result
def gain_prime(y,s,k,eta):
    h=(s-k)/2
    return (h*Rp(y)+(s*s+k*k)*Rpp(y)/4
            +(s**3-k**3)+eta*(s+k)/2)
def gain_second(y,s,k): return (s-k)*Rpp(y)/2+3*(s*s+k*k)

def Cfinite(t,z,eta):
    global DIVISOR_LOWER
    y=t+z
    r=6*(1+y)/(y*(3+2*y))
    k1=(3-4*y*y)/(y*(3+2*y))
    A0=3*(2*y*y+4*y+3)/(3+2*y)
    A1=2*(4*y*y+6*y+3)/(y*(3+2*y))
    fy=(y*y+3*y+3)/(3+2*y)
    H=eta*k1/(fy-eta*r/2)
    mandatory(0<H.lo<=H.hi<y.lo, 'verify_certificates.py:153 mandatory 0<H.lo<=H.hi<y.lo')
    for _ in range(5):
        m=A0-eta*r/2-A1*H+r*H*H/2
        mandatory(m.lo>2, 'verify_certificates.py:156 mandatory m.lo>2')
        H=H.minimum(eta*k1/m)
    m=A0-eta*r/2-A1*H+r*H*H/2
    mandatory(m.lo>2, 'verify_certificates.py:159 mandatory m.lo>2')
    DIVISOR_LOWER=m.lo if DIVISOR_LOWER is None else min(DIVISOR_LOWER,m.lo)
    result=(t+kap)*k1*k1/(4*m)
    mandatory(result.lo>0, 'verify_certificates.py:161 mandatory result.lo>0')
    return result,H

def main():
    z=scalar_root(price,Q('.393'),Q('.394'),Q(1,10**80))
    require('pricing root is unique',price_prime(I(Q('.393'),Q('.394'))).lo>0)
    root=nuisance_root(lambda t:P(t,z),Q('.076'),Q('.077'),Q(1,10**55),'type root')
    require('margin slope exceeds 1/4',(I(a,b)+kap).__mul__(Rp(I(a,b)+z)).lo>Q(1,4))
    zpaper=root_decimal_hull(z,60)
    tpaper=root_decimal_hull(root,50)
    # These signs certify the finite printed boxes as well as the bisection roots.
    require('paper pricing root box signs',price(zpaper.lo)<0<price(zpaper.hi))
    require('paper type root box signs',P(I.p(tpaper.lo),zpaper).hi<0<P(I.p(tpaper.hi),zpaper).lo)
    require('generated pricing box equals printed Appendix C input',zpaper.lo==Q('.393103890949911708598915717576988616660911028234782729295715') and zpaper.hi-zpaper.lo==Q(1,10**60))
    require('generated type box equals printed Appendix C input',tpaper.lo==Q('.07637758394939579289694434827763576112329960259796') and tpaper.hi-tpaper.lo==Q(1,10**50))
    istar=20*(Pprimitive(zpaper+b,zpaper)-Pprimitive(zpaper+tpaper,zpaper))
    require('positive-radius limiting opportunity',istar.lo>0)
    require('buffered direction has zero mass',-Q(1,2)*Q(1,5)+Q(1,10)==0)
    require('buffered law remains nonnegative',Rp(Q(1,5))>Q(1,20)/2)
    require('upper localization H0 below support',Q(7,85)<Q(443,1000))
    require('uniform bootstrap denominators exceed two',3-Q(3,20)-10*Q(7,85)>2)
    rows=[]
    for eta,s,printed_low,printed_up,printed_width in zip(RADII,CUTS,TABLE_LOW,TABLE_UP,TABLE_WIDTH):
        index=len(rows)
        tag=str(eta)
        yr=z+root
        rootk=nuisance_root(lambda k:twice_gain(yr,s,k,eta),0,Q(1,100),
                            Q(1,10**50),'compensation '+tag)
        kpaper=root_decimal_hull(rootk,30)
        yrpaper=zpaper+tpaper
        tie_left=twice_gain(yrpaper,s,I.p(kpaper.lo),eta)
        tie_right=twice_gain(yrpaper,s,I.p(kpaper.hi),eta)
        require(tag+' printed compensation root box signs',
                tie_left.lo>Q('3.9e-31') and tie_right.hi<-Q('3.2e-31'))
        derivative_k=-(R(yrpaper-kpaper)-eta*(yrpaper-kpaper-D))
        require(tag+' compensation derivative less than -3/4',derivative_k.hi<-Q(3,4))
        # All following witness bounds use the printed paper boxes, not hidden roots.
        kk=kpaper
        require(tag+' all shifted types in buffer',
                (zpaper+a-kk).lo>Q(2,5) and (zpaper+b+s).hi<D)
        require(tag+' both promises positive',
                (H1-zpaper-s).lo>0 and (H2-zpaper+kk).lo>0)
        require(tag+' printed shifted-index and promise margins',(zpaper+a-kk).lo>Q('.4363') and (zpaper+b+s).hi<Q('.4997') and (H1-zpaper-s).lo>Q('.944') and (H2-zpaper+kk).lo>Q('1.150'))
        require(tag+' common densities remain positive',Rp(Q(1,5))>eta/2 and Rp(Q(2,5))>eta)
        ya=zpaper+a
        db0=gain_prime(ya,s,kk,eta)
        db1=db0+(b-a)*gain_second(ya,s,kk)/2
        db2=gain_prime(zpaper+b,s,kk,eta)
        bernstein=[db0,db1,db2]
        require(tag+' all quadratic Bernstein coefficients positive',min(v.lo for v in bernstein)>0)
        require(tag+' printed Bernstein lower margin',min(v.lo for v in bernstein)>=TABLE_DERIVATIVE[index])
        coefs=gain_coefficients(s,kk,eta)
        integral=I.p(0)
        for i,c in enumerate(coefs):
            integral=integral+c*((zpaper+b)**(i+1)-(zpaper+tpaper)**(i+1))/(i+1)
        extra=20*((b+kap)*eval_poly(coefs,zpaper+b)-2*integral)
        require(tag+' printed lower excess is certified',extra.lo>=printed_low)
        # Direct full-class upper: 1024 whole-cell interval evaluations.
        n=1024
        Cmax=Q(0)
        Hmax=Q(0)
        for j in range(n):
            lo=a+(b-a)*j/n
            hi=a+(b-a)*(j+1)/n
            c,hh=Cfinite(I(lo,hi),zpaper,eta)
            Cmax=max(Cmax,c.hi)
            Hmax=max(Hmax,hh.hi)
        coefficient=I.p(0)
        start=tpaper.lo
        for j in range(n):
            lo=start+(b-start)*j/n
            hi=start+(b-start)*(j+1)/n
            c,_=Cfinite(I(lo,hi),zpaper,eta)
            coefficient=coefficient+20*c*(hi-lo)
        boundary=40*Cmax*Cmax*eta**4
        upper_extra=eta*eta*coefficient+boundary
        paper_M=Q(decimal(Cmax,16,True))
        paper_J=Q(decimal(coefficient.hi,16,True))
        paper_upper=eta**2*paper_J+40*paper_M**2*eta**4
        require(tag+' printed Appendix C M/J ingredients',Cmax<=TABLE_M[index] and coefficient.hi<=TABLE_J[index])
        require(tag+' rounded paper upper ingredients certify table',paper_upper<=printed_up)
        require(tag+' printed upper excess is certified',upper_extra.hi<=printed_up)
        require(tag+' printed width is certified',upper_extra.hi-extra.lo<=printed_width)
        rows.append({'eta':str(eta),'first_cut_exact':str(s),
            'second_raise_paper_root_box':kpaper.output(),
            'twice_gain_left_endpoint':tie_left.output(),
            'twice_gain_right_endpoint':tie_right.output(),
            'twice_gain_compensation_derivative':derivative_k.output(),
            'gain_derivative_Bernstein_coefficients':[v.output() for v in bernstein],
            'actual_excess':extra.output(),'actual_profit':(istar+extra).output(),
            'Cmax':I.p(Cmax).output(),'Hmax':I.p(Hmax).output(),
            'paper_Cmax_upper':str(paper_M),'paper_integral_upper':str(paper_J),
            'paper_upper_excess':str(paper_upper),
            'reported_targets':{'lower_excess':str(printed_low),'upper_excess':str(printed_up),'width':str(printed_width),'derivative_lower':str(TABLE_DERIVATIVE[index]),'Cmax_upper':str(TABLE_M[index]),'integral_upper':str(TABLE_J[index])},
            'coefficient_integral':coefficient.output(),'boundary_budget':I.p(boundary).output(),
            'global_upper_excess':upper_extra.output(),'global_upper_profit':(istar+upper_extra).output(),
            'upper_minus_lower':(upper_extra-extra).output()})
    output={'status':'pass','arithmetic':'outward rational intervals after every addition, multiplication and reciprocal, denominator 10^65; minimum by endpoint minima; positive-divisor guards; recursion m0,H1,...,m4,H5,m5,C5; rational k1 expression',
            'dependencies':'Python 3 standard library and supplied certificate_checks.py; no data inputs',
            'z_paper_box':zpaper.output(),'type_root_paper_box':tpaper.output(),
            'pricing_root_endpoint_values':[str(price(zpaper.lo)),str(price(zpaper.hi))],
            'pricing_derivative':price_prime(I(Q('.393'),Q('.394'))).output(),
            'type_root_endpoint_values':[P(I.p(tpaper.lo),zpaper).output(),P(I.p(tpaper.hi),zpaper).output()],
            'type_margin_derivative':((I(a,b)+kap)*Rp(I(a,b)+zpaper)).output(),
            'minimum_computed_m5_lower':str(DIVISOR_LOWER),'continuous_type_cells_per_partition':1024,'iterations':5,'positive_integral_lower_endpoint':str(tpaper.lo),
            'I_star':istar.output(),'rows':rows,'checks':CHECKS,
            'coverage':'all-law/menu coverage is the analytic theorem; this script verifies its numerical bounds and five actual quotes'}
    target=Path(__file__).resolve().with_name('certificates.json')
    target.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'checks':len(CHECKS),'output':str(target),'rows':[
        {'eta':v['eta'],'lower_excess':v['actual_excess']['decimal_lower'],
         'upper_excess':v['global_upper_excess']['decimal_upper']} for v in rows]},indent=2))

if __name__=='__main__': main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
