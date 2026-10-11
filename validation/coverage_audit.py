#!/usr/bin/env python3
"""Audit generated report coverage and lossless outward decimal summaries.

This post-run audit reads reports. The mathematical verifiers themselves start
from embedded exact primitives and never read these reports as inputs.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
import sys
sys.set_int_max_str_digits(0)  # Shipped lossless fractions can exceed 4,300 digits.
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checkers'))
from certificate_checks import mandatory

root=Path(__file__).resolve().parents[1]/'checkers'
checks=[];decimal_pairs=0
def check(label,condition):
    mandatory(condition,label);checks.append(label)
def displays(obj):
    global decimal_pairs
    if isinstance(obj,float):raise RuntimeError('binary float in mathematical report')
    if isinstance(obj,list):
        for value in obj:displays(value)
    if not isinstance(obj,dict):return
    for lower,upper,dl,du in [
        ('rational_lower','rational_upper','decimal_lower','decimal_upper'),
        ('lower','upper','decimal_lower','decimal_upper'),
        ('lower_fraction','upper_fraction','lower_decimal','upper_decimal')]:
        if all(key in obj for key in (lower,upper,dl,du)):
            mandatory(Q(obj[dl])<=Q(obj[lower])<=Q(obj[upper])<=Q(obj[du]),'outward report display')
            decimal_pairs+=1
    if all(key in obj for key in ('rational','decimal_lower','decimal_upper')):
        mandatory(Q(obj['decimal_lower'])<=Q(obj['rational'])<=Q(obj['decimal_upper']),'outward scalar report display')
        decimal_pairs+=1
    for value in obj.values():displays(value)
reports={}
for path in sorted(root.glob('*.json')):
    obj=json.loads(path.read_text());displays(obj)
    reports[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
check('ten advertised suite outputs exist',len(reports)==10)
family=json.loads((root/'power_family_certificates.json').read_text())
check('three discrete exponents',[r['nu'] for r in family['results']]==[2,3,4])
for row in family['results']:
    cells=row['cells'];nu=row['nu']
    check('256 unique original joint cells nu '+str(nu),len(cells)==256 and
          {(c['z_index'],c['lambda_index']) for c in cells}=={(i,j) for i in range(16) for j in range(16)})
    for c in cells:
        i,j=c['z_index'],c['lambda_index']
        check('full parameter bounds '+str((nu,i,j)),
              Q(c['z_box']['lower'])==Q(13,20)+Q(1,10)*i/16 and
              Q(c['z_box']['upper'])==Q(13,20)+Q(1,10)*(i+1)/16 and
              Q(c['lambda_box']['lower'])==-Q(1,4)+Q(3,4)*j/16 and
              Q(c['lambda_box']['upper'])==-Q(1,4)+Q(3,4)*(j+1)/16 and c['status']=='pass')
    for name,aggregate in row['bounds'].items():
        check('exact recorded hull '+str((nu,name)),
              Q(aggregate['lower'])==min(Q(c['bounds'][name]['lower']) for c in cells) and
              Q(aggregate['upper'])==max(Q(c['bounds'][name]['upper']) for c in cells))
precision=json.loads((root/'public_precision_certificate.json').read_text())
check('all twelve defining root sign records',len(precision['root_checks'])==12)
check('all seven gain/margin Bernstein coefficients',
      len(precision['Bernstein_coefficients']['gain_derivative'])==3 and
      len(precision['Bernstein_coefficients']['margin_derivative'])==4)
check('four adverse coefficients',len(precision['adverse_Bernstein_coefficients'])==4)
check('three signed-band rows',len(precision['rows'])==3)
result={'status':'pass','checks':len(checks),'check_labels':checks,
        'outward_decimal_intervals_checked':decimal_pairs,'verified_joint_cells':768,
        'report_sha256':reports,'scope':'report completeness and rational display inclusion, in addition to the verifiers mathematical obligations'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('status','checks','outward_decimal_intervals_checked','verified_joint_cells')}))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
