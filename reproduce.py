#!/usr/bin/env python3
"""Reproduce ten certificate suites in fresh source-only directories.

Default: mathematical checks under normal Python and Python -O.
--article-only rebuilds the accompanying source tree without repeating math.
All paths resolve from this file; no data or earlier certificate is an input.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time

SUPPLEMENT=Path(__file__).resolve().parent
ROOT=SUPPLEMENT.parent
CHECKERS=SUPPLEMENT/'checkers'
LOGS=SUPPLEMENT/'execution'
JOBS=[
 ('coefficient','stage3_coefficient_checks.py','stage3_coefficient_certificates.json'),
 ('uniform_constants','stage3_finite_constant_checks.py','stage3_finite_constant_certificates.json'),
 ('coefficient_adversarial','stage3_coefficient_adversarial_checks.py','stage3_coefficient_adversarial_certificates.json'),
 ('tiny_quote','stage3_finite_radius_certificate.py','stage3_finite_radius_certificate.json'),
 ('payment_slack','payment_slack_diagnostic.py','payment_slack_diagnostic.json'),
 ('disclosed_quotes_and_direct_upper','verify_certificates.py','certificates.json'),
 ('continuous_power_family','power_family_checks.py','power_family_certificates.json'),
 ('power_direct_upper','power_family_bounds.py','power_family_bounds.json'),
 ('paper_precision_boxes','public_precision_certificate.py','public_precision_certificate.json'),
 ('interval_kernel_and_failure','interval_kernel_checks.py','interval_kernel_checks.json'),
]

def timestamp():return dt.datetime.now(dt.timezone.utc).isoformat()
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def canonical_digest(path):
    obj=json.loads(path.read_text())
    raw=json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()
def no_float(value):
    if isinstance(value,float):return False
    if isinstance(value,dict):return all(no_float(v) for v in value.values())
    if isinstance(value,list):return all(no_float(v) for v in value)
    return True
def versions():
    packages={}
    for name in ['sympy','mpmath','matplotlib']:
        try:packages[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:packages[name]=None
    return {'python':sys.version,'python_executable':sys.executable,
            'platform':platform.platform(),'packages':packages}
def relative(path):return str(path.relative_to(ROOT))
def validate_sources():
    expected=json.loads((SUPPLEMENT/'SOURCE_MANIFEST.json').read_text())
    actual={str(p.relative_to(CHECKERS)):digest(p) for p in CHECKERS.rglob('*.py')}
    if actual!=expected:raise RuntimeError('checker source manifest differs; no run certified')
    return expected
def mathematical_summary(name,output):
    d=json.loads(output.read_text())
    if name=='continuous_power_family':
        return {'exponents':[x['nu'] for x in d['results']],
                'all_cells_verified':all(x['all_cells_verified'] for x in d['results']),
                'verified_cells':sum(len(x['cells']) for x in d['results']),
                'bad_cells':sum(x['bad_cells'] for x in d['results'])}
    if name=='power_direct_upper':
        return {'response_exponents':[x['nu'] for x in d['rows']],
                'alpha':d['uniform']['alpha']['rational'],
                'defect':d['uniform']['defect']['rational'],
                'direct_Cmax':d['uniform']['direct_Cmax']['rational'],
                'boundary_constant':d['uniform']['boundary_constant']['rational']}
    count=d.get('count',d.get('number_of_checks',d.get('checks_passed')))
    if count is None:count=len(d.get('checks',[]))
    summary={'stated_check_count':count}
    if name=='disclosed_quotes_and_direct_upper':summary['certified_radii']=[x['eta'] for x in d['rows']]
    if name=='paper_precision_boxes':summary.update(robust_lower_cases=len(d['rows']),root_sign_records=len(d['root_checks']))
    return summary

def run_math(modes):
    validate_sources()
    expected=json.loads((SUPPLEMENT/'EXPECTED_OUTPUTS.json').read_text())
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
    env['PYTHONDONTWRITEBYTECODE']='1'
    ledger={'schema':2,'scope':'ten suites under normal and/or optimized Python',
            'status':'running','started_utc':timestamp(),'environment':versions(),
            'invoked_from':str(Path.cwd()),'modes':modes,
            'source_manifest_sha256':digest(SUPPLEMENT/'SOURCE_MANIFEST.json'),
            'expected_outputs_manifest_sha256':digest(SUPPLEMENT/'EXPECTED_OUTPUTS.json'),
            'checks':'exception-based mandatory guards, including under -O',
            'inputs':'fresh temporary copies of shipped Python sources only',
            'normalization':'canonical JSON, no output fields excluded',
            'article_build_performed':False,'jobs':[]}
    path=SUPPLEMENT/'execution_ledger.json'
    path.write_text(json.dumps(ledger,indent=2)+'\n')
    for _,_,output in JOBS:(CHECKERS/output).unlink(missing_ok=True)
    try:
        for optimization in modes:
            mode='optimized' if optimization else 'normal'
            mode_logs=LOGS/mode;mode_logs.mkdir(parents=True,exist_ok=True)
            with tempfile.TemporaryDirectory(prefix='jgo-source-only-') as folder:
                work=Path(folder)
                for source in CHECKERS.rglob('*.py'):
                    dest=work/source.relative_to(CHECKERS)
                    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
                for name,script,output_name in JOBS:
                    output=work/output_name;output.unlink(missing_ok=True)
                    command=[sys.executable]+(['-O'] if optimization else [])+['-B',script]
                    start=time.monotonic();begun=timestamp()
                    result=subprocess.run(command,cwd=work,env=env,text=True,
                                          encoding='utf-8',errors='replace',
                                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                    log=mode_logs/(name+'.txt')
                    log.write_text(result.stdout.replace(str(work),'{checker_directory}'))
                    actual=canonical_digest(output) if output.exists() else None
                    lossless=output.exists() and no_float(json.loads(output.read_text()))
                    expected_hash=expected[output_name]['canonical_sha256']
                    passed=result.returncode==0 and actual==expected_hash and lossless
                    entry={'job':name,'optimization':optimization,'started_utc':begun,
                           'finished_utc':timestamp(),'elapsed_seconds':time.monotonic()-start,
                           'actual_command':command,'actual_cwd':str(work),
                           'portable_command':['python3']+(['-O'] if optimization else [])+['-B',script],
                           'portable_cwd':'checkers/ (source-only copy)',
                           'exit_code':result.returncode,'status':'pass' if passed else 'fail',
                           'log':relative(log),'log_sha256':digest(log),
                           'output':relative(CHECKERS/output_name),
                           'output_sha256':digest(output) if output.exists() else None,
                           'expected_canonical_sha256':expected_hash,
                           'actual_canonical_sha256':actual,'expected_output_match':actual==expected_hash,
                           'no_binary_float_fields':lossless,'script_sha256':digest(work/script)}
                    if result.returncode==0 and output.exists():entry['mathematical_summary']=mathematical_summary(name,output)
                    ledger['jobs'].append(entry);path.write_text(json.dumps(ledger,indent=2)+'\n')
                    print(mode+'/'+name+': '+('PASS' if passed else 'FAILED'),flush=True)
                    if not passed:raise RuntimeError('certificate failure: '+mode+'/'+name+'\n'+result.stdout[-5000:])
                    shutil.copy2(output,CHECKERS/output_name)
        validate_sources()
        ledger.update(status='pass',finished_utc=timestamp(),all_requested_suites_passed=True,
                      all_expected_outputs_match=True,source_hashes_unchanged=True)
        path.write_text(json.dumps(ledger,indent=2)+'\n')
    except BaseException as exc:
        for _,_,output in JOBS:(CHECKERS/output).unlink(missing_ok=True)
        ledger.update(status='failed',finished_utc=timestamp(),error=str(exc),all_requested_suites_passed=False)
        path.write_text(json.dumps(ledger,indent=2)+'\n');raise
    return ledger

def run_article():
    LOGS.mkdir(exist_ok=True)
    script=ROOT/'source'/'build.sh'
    if not script.exists():raise RuntimeError('accompanying manuscript source tree is required')
    start=time.monotonic();began=timestamp()
    result=subprocess.run(['bash',str(script)],cwd=ROOT,text=True,encoding='utf-8',errors='replace',
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    log=LOGS/'article_build.txt';log.write_text(result.stdout)
    record={'scope':'figure and PDF build only','status':'pass' if result.returncode==0 else 'failed',
            'started_utc':began,'finished_utc':timestamp(),'elapsed_seconds':time.monotonic()-start,
            'environment':versions(),'portable_command':'python3 reproduce.py --article-only',
            'cwd':'.','exit_code':result.returncode,'log':relative(log),'log_sha256':digest(log),
            'pdf_sha256':digest(ROOT/'source'/'main.pdf') if result.returncode==0 else None}
    (SUPPLEMENT/'article_execution_ledger.json').write_text(json.dumps(record,indent=2)+'\n')
    print('article_build: '+('PASS' if result.returncode==0 else 'FAILED'),flush=True)
    if result.returncode:raise SystemExit(result.returncode)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--math-only',action='store_true',help='run ten suites (default)')
    mode.add_argument('--article-only',action='store_true',help='build accompanying figure/PDF')
    mode.add_argument('--with-article',action='store_true',help='run math then build figure/PDF')
    parser.add_argument('--optimization',choices=['both','normal','optimized'],default='both')
    args=parser.parse_args()
    if not args.article_only:run_math({'both':[0,1],'normal':[0],'optimized':[1]}[args.optimization])
    if args.article_only or args.with_article:run_article()

if __name__=='__main__':main()

# Publication metadata
# Article: Global certificates for endogenous mortgage entry under joint collateral-law constraints
# Journal: Journal of Global Optimization
# Author and corresponding author: Letong Zhou
# Affiliation: EDHEC Business School, Paris, France
# Corresponding email: letong.zhou@edhec.com
# Resource: Online Resource 2 (ESM_2.zip)
