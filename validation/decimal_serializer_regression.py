#!/usr/bin/env python3
"""Signed near-integer regressions for the three existing interval serializers."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
import json
import sys

root=Path(__file__).resolve().parents[1]/'checkers'
sys.path.insert(0,str(root))
import verify_certificates as v
import interval_utils as u
from certificate_checks import mandatory
spec=importlib.util.spec_from_file_location('serializer_precision',root/'numeric'/'precision_certificate.py')
p=importlib.util.module_from_spec(spec);sys.modules[spec.name]=p;spec.loader.exec_module(p)
serializers=[('verify.decimal',lambda q:{'decimal_lower':v.decimal(q,35,False),'decimal_upper':v.decimal(q,35,True)}),
             ('interval_utils.Box.output',lambda q:u.Box(q,q).output()),
             ('precision.Box.output',lambda q:p.Box(q,q).output())]
records=[]
for name,serialize in serializers:
    for sign in (-1,1):
        for offset in (-1,0,1):
            q=sign*Q(10**250+offset,10**250)
            result=serialize(q)
            mandatory(Q(result['decimal_lower'])<=q<=Q(result['decimal_upper']),name+' signed near-integer enclosure')
            records.append({'serializer':name,'exact':str(q),'decimal_lower':result['decimal_lower'],
                            'decimal_upper':result['decimal_upper'],'status':'pass'})
result={'status':'pass','checks':len(records),'records':records,
        'scope':'all three former nearest-quotient serializers; signed exact and near-integer fractions'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':len(records),'optimization':sys.flags.optimize}))

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
