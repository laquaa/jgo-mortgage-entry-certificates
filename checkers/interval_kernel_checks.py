#!/usr/bin/env python3
"""Small exact inclusion and mandatory-failure checks for the shipped kernel."""
from fractions import Fraction as Q
from pathlib import Path
import ast
import json
import subprocess
import sys
import tempfile
import shutil
import verify_certificates as v
from certificate_checks import mandatory, CertificateFailure, decimal_bound

checks = []


def check(label, condition):
    mandatory(condition, label)
    checks.append(label)


def includes(box, values):
    return all(box.lo <= value <= box.hi for value in values)


def stable_failure_output(output, command):
    """Keep legacy traceback frames, removing interpreter display annotations.

    Python 3.13 adds range underlines and the source of the -c command even
    under -O. The raw trace is retained in the reproduction console log.
    """
    lines = []
    after_string_frame = False
    for line in output.splitlines(keepends=True):
        stripped = line.strip()
        if stripped and set(stripped) <= {'~', '^'}:
            continue
        if after_string_frame and stripped == command:
            after_string_frame = False
            continue
        after_string_frame = line.startswith('  File "<string>",')
        lines.append(line)
    return ''.join(lines)


a = v.I(Q(-3), Q(-2))
b = v.I(Q(-5), Q(1))
check('signed multiplication encloses all endpoint products',
      includes(a*b, [x*y for x in (a.lo,a.hi) for y in (b.lo,b.hi)]))
check('positive reciprocal', (v.I.p(1)/v.I(Q(2),Q(4))).lo==Q(1,4)
      and (v.I.p(1)/v.I(Q(2),Q(4))).hi==Q(1,2))
check('negative reciprocal', (v.I.p(1)/v.I(Q(-4),Q(-2))).lo==Q(-1,2)
      and (v.I.p(1)/v.I(Q(-4),Q(-2))).hi==Q(-1,4))
check('minimum endpoint enclosure', includes(a.minimum(b),
      [min(x,y) for x in (a.lo,a.hi) for y in (b.lo,b.hi)]))
check('subtraction and endpoint reversal', (a-b).lo==Q(-4) and (a-b).hi==Q(3))
check('repeated multiplication powers include signed point powers',
      includes(v.I(Q(-2),Q(3))**4,[Q(-2)**4,Q(0)**4,Q(3)**4]))
unit=Q(1,v.DEN)
check('positive lattice rounding',v.down(Q(2,3)*unit)==0 and v.up(Q(2,3)*unit)==unit)
check('negative lattice rounding',v.down(-Q(2,3)*unit)==-unit and v.up(-Q(2,3)*unit)==0)
check('lattice boundary unchanged',v.down(-2*unit)==-2*unit and v.up(-2*unit)==-2*unit)
for value in [Q('1.01')-Q(1,10**150),Q('1.01')+Q(1,10**150),
              -Q('1.01')-Q(1,10**150),-Q('1.01')+Q(1,10**150),
              Q(10**200)+Q(1,3)]:
    check('directed decimal division and quantization preserve inclusion: '+str(value),
          Q(decimal_bound(value,False,2))<=value<=Q(decimal_bound(value,True,2)))
near=v.I(Q('1.00001')-Q(1,10**150),Q('1.00001')+Q(1,10**150)).output()
check('benchmark decimal serializer encloses near-boundary endpoints',
      Q(near['decimal_lower'])<=Q(near['rational_lower']) and
      Q(near['rational_upper'])<=Q(near['decimal_upper']))
for label,operation in [('zero-containing divisor',lambda:v.I.p(1)/v.I(Q(-1),Q(1))),
                        ('negative power',lambda:v.I.p(1)**-1),
                        ('unordered endpoints',lambda:v.I(Q(2),Q(1))),
                        ('binary float input',lambda:v.I.p(0.05))]:
    try:
        operation()
    except CertificateFailure:
        check(label+' fails explicitly',True)
    else:
        raise CertificateFailure(label+' incorrectly succeeded')

root=Path(__file__).resolve().parent
for path in sorted(root.rglob('*.py')):
    check('no removable mandatory assert in '+str(path.relative_to(root)),
          not any(isinstance(node,ast.Assert) for node in ast.walk(ast.parse(path.read_text()))))

# Exercise a real printed-bound failure with optimization enabled in a fresh
# no-data directory, rather than only testing a constant false predicate.
with tempfile.TemporaryDirectory(prefix='jgo-mandatory-failure-') as folder:
    dest=Path(folder)
    for name in ('verify_certificates.py','certificate_checks.py'):
        shutil.copy2(root/name,dest/name)
    command=[sys.executable,'-O','-B','-c',
             'import verify_certificates as v; v.TABLE_UP[0]=v.Q(0); v.main()']
    run=subprocess.run(command,cwd=dest,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    check('contradictory printed bound exits nonzero under -O',run.returncode!=0)
    check('failure identifies violated upper certificate',
          run.stdout.splitlines()[-1:] == ['certificate_checks.CertificateFailure: 1/5000 rounded paper upper ingredients certify table'])
    check('failed run produces no successful certificate',not (dest/'certificates.json').exists())
    failure={'command':['python3','-O','-B','-c',command[-1]],
             'injected_change':'comparison target TABLE_UP[0]=0, temporary copy only',
             'exit_code':run.returncode,
             'output':stable_failure_output(run.stdout,command[-1]).replace(str(dest),'{temporary_checker_directory}'),
             'success_output_absent':True}
    print('Raw optimized-Python failure probe:', flush=True)
    print(run.stdout.replace(str(dest),'{temporary_checker_directory}'), end='', flush=True)

output={'status':'pass','checks':checks,'count':len(checks),
        'arithmetic':'exact fractions; all decisions rational',
        'optimized_python_failure_probe':failure}
root.joinpath('interval_kernel_checks.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':len(checks),'optimization':sys.flags.optimize,
                  'real_bound_failure_exit':failure['exit_code']},indent=2))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
