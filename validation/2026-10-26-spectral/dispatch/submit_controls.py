"""Idempotent fixed-scope dispatch; intent is durable before each sbatch."""
import datetime,hashlib,json,pathlib,subprocess,tomllib
R=pathlib.Path('/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct');B=R/'outputs/spectral-v5';D=R/'spectral-dispatch-v1';M=B/'dispatch';M.mkdir(exist_ok=True)
NUM=R/'spectral-20261026-v5-source'
def write(p,x):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
# Passed source/config preflight and all four completed control cases are gates.
pre=json.loads((R/'outputs/spectral-control-pairing-preflight.json').read_text());assert pre['status']=='passed'
reused=[(1,1),(1,2),(16,1),(16,2)]
for s,rep in reused:
 p=B/f'control-{s}-{rep}'
 assert tomllib.loads((p/'complete.toml').read_text())['status'] in ('completed','completed_with_failures')
 assert tomllib.loads((p/'independent-check.toml').read_text())['status']=='passed'
# Numerical sources remain byte-identical to the passed v5 check/case source.
manifest=json.loads((NUM/'source.json').read_text())
for f,digest in manifest['files'].items():assert hashlib.sha256((NUM/f).read_bytes()).hexdigest()==digest
plan_path=M/'controls-plan.json'
if not plan_path.exists():
 tasks=[]
 for s in range(1,17):
  for start in range(1,101,10):
   cases=[('control',s,rep) for rep in range(start,start+10) if (s,rep) not in reused]
   task=M/f'control-worker-{len(tasks)+1:03d}.tsv';assert not task.exists()
   task.write_text(''.join(f'{k}\t{i}\t{rep}\n' for k,i,rep in cases));tasks.append(dict(task=str(task),datasets=cases))
 pools=[]
 for i in range(0,len(tasks),4):
  index=i//4+1;tag=f'dfgm-spectral-controls-v5-{index:03d}';tsv=M/f'control-pool-{index:03d}.tsv'
  tsv.write_text(''.join(t['task']+'\n' for t in tasks[i:i+4]))
  pools.append(dict(index=index,tag=tag,tsv=str(tsv),workers=tasks[i:i+4]))
 plan=dict(pools=pools,expected_datasets=1600,reused_cases=reused,new_datasets=1596,numerical_source=str(NUM),
  orchestration_sha256={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ('spectral_batch.jl','control-worker.sh','controls.sbatch')})
 assert sum(len(t['datasets']) for p in pools for t in p['workers'])==1596
 write(plan_path,plan)
plan=json.loads(plan_path.read_text());journal_path=M/'control-submission-journal.json'
journal=json.loads(journal_path.read_text()) if journal_path.exists() else dict(pools={},source=str(NUM),scope='16x100 paired controls, same160/40 split; fixed6-score menu; no baseline fits')
for pool in plan['pools']:
 key=str(pool['index']);record_path=M/f'control-submission-{key}.json'
 if record_path.exists():
  record=json.loads(record_path.read_text())
  if str(record.get('job_id','')).isdigit():journal['pools'][key]=record;write(journal_path,journal);continue
  raise RuntimeError(f'Uncertain previous submission; reconcile tag before retry: {record_path}')
 cmd=['sbatch','--parsable',f'--job-name={pool["tag"]}',f'--comment={pool["tag"]}',
      '-o',str(M/f'control-pool-{pool["index"]:03d}-%j.out'),'-e',str(M/f'control-pool-{pool["index"]:03d}-%j.err'),str(D/'controls.sbatch'),pool['tsv']]
 record={**pool,'state':'submitting','submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd}
 write(record_path,record)
 reply=subprocess.run(cmd,capture_output=True,text=True);record.update(stdout=reply.stdout,stderr=reply.stderr,returncode=reply.returncode)
 job=reply.stdout.strip().split(';')[0]
 if reply.returncode or not job.isdigit():
  record['state']='uncertain';write(record_path,record);raise RuntimeError('Ambiguous submission reply retained; do not blindly retry')
 record.update(state='submitted',job_id=job);write(record_path,record);journal['pools'][key]=record;write(journal_path,journal)
 print(json.dumps(dict(pool=key,job_id=job,datasets=sum(len(t['datasets']) for t in pool['workers']))),flush=True)
print(json.dumps(dict(status='submitted',pools=len(journal['pools']),datasets=1596,reused=4)),flush=True)
