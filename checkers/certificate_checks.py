"""Mandatory checks and lossless display helpers; standard library only.

These checks remain active under Python -O. A failure raises an uncaught
CertificateFailure, giving the certificate process a nonzero exit status.
"""
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING


class CertificateFailure(RuntimeError):
    pass


def mandatory(condition, label):
    if not bool(condition):
        raise CertificateFailure(str(label))


def exact_fraction(value):
    if isinstance(value, float):
        raise CertificateFailure('binary floating input is forbidden')
    return Fraction(value)


def decimal_bound(value, upward=False, digits=25):
    value = exact_fraction(value)
    mandatory(isinstance(digits,int) and digits>=0,'nonnegative decimal display precision')
    integer_digits=max(1,len(str(abs(value.numerator)))-len(str(value.denominator))+1)
    with localcontext() as ctx:
        ctx.prec = max(110, digits + 45, integer_digits + digits + 10)
        ctx.rounding = ROUND_CEILING if upward else ROUND_FLOOR
        return str((Decimal(value.numerator) / Decimal(value.denominator)).quantize(
            Decimal(1).scaleb(-digits),
            rounding=ctx.rounding))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
