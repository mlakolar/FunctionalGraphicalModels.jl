"""Idempotent fixed-scope dispatch; intent is durable before each sbatch."""
import datetime,hashlib,json,pathlib,subprocess,tomllib
R=pathlib.Path('/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct');B=R/'outputs/spectral-v5';D=R/'spectral-dispatch-v1';M=B/'dispatch';M.mkdir(exist_ok=True)
NUM=R/'spectral-20261026-v5-source'
def write(p,x):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
# Full valid p50 canaries are primary outputs, not refitted.
reused=[(1,1),(7,1)]
for setting,rep in reused:
 source=R/f'outputs/spectral-v5-canary/paper-{setting}-{rep}'
 assert tomllib.loads((source/'complete.toml').read_text())['status'] in ('completed','completed_with_failures')
 assert tomllib.loads((source/'independent-check.toml').read_text())['status']=='passed'
 dest=B/f'paper-{setting}-{rep}'
 if dest.is_symlink():assert dest.resolve()==source.resolve()
 else:assert not dest.exists();dest.symlink_to(source)
# Numerical sources remain byte-identical to the passed v5 check/case source.
manifest=json.loads((NUM/'source.json').read_text())
for f,digest in manifest['files'].items():assert hashlib.sha256((NUM/f).read_bytes()).hexdigest()==digest
plan_path=M/'paper50-plan.json'
if not plan_path.exists():
 tasks=[]
 for s in (1,2,3,7,8,9):
  for start in range(1,21,2):
   cases=[('paper',s,rep) for rep in range(start,start+2) if (s,rep) not in reused]
   task=M/f'paper50-worker-{len(tasks)+1:03d}.tsv';assert not task.exists()
   task.write_text(''.join(f'{k}\t{i}\t{rep}\n' for k,i,rep in cases));tasks.append(dict(task=str(task),datasets=cases))
 pools=[]
 for i in range(0,len(tasks),4):
  index=i//4+1;tag=f'dfgm-spectral-paper50-v5-{index:03d}';tsv=M/f'paper50-pool-{index:03d}.tsv'
  tsv.write_text(''.join(t['task']+'\n' for t in tasks[i:i+4]))
  pools.append(dict(index=index,tag=tag,tsv=str(tsv),workers=tasks[i:i+4]))
 plan=dict(pools=pools,expected_datasets=120,reused_cases=reused,new_datasets=118,numerical_source=str(NUM),
  orchestration_sha256={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ('spectral_batch.jl','paper-worker.sh','paper.sbatch')})
 assert sum(len(t['datasets']) for p in pools for t in p['workers'])==118
 write(plan_path,plan)
plan=json.loads(plan_path.read_text());journal_path=M/'paper50-submission-journal.json'
journal=json.loads(journal_path.read_text()) if journal_path.exists() else dict(pools={},source=str(NUM),scope='six p50 paper settings x20 paired seeds; six-score exact-grid protocol; no baseline fits')
for pool in plan['pools']:
 key=str(pool['index']);record_path=M/f'paper50-submission-{key}.json'
 if record_path.exists():
  record=json.loads(record_path.read_text())
  if str(record.get('job_id','')).isdigit():journal['pools'][key]=record;write(journal_path,journal);continue
  raise RuntimeError(f'Uncertain previous submission; reconcile tag before retry: {record_path}')
 cmd=['sbatch','--parsable',f'--job-name={pool["tag"]}',f'--comment={pool["tag"]}',
      '-o',str(M/f'paper50-pool-{pool["index"]:03d}-%j.out'),'-e',str(M/f'paper50-pool-{pool["index"]:03d}-%j.err'),str(D/'paper.sbatch'),pool['tsv']]
 record={**pool,'state':'submitting','submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd}
 write(record_path,record)
 reply=subprocess.run(cmd,capture_output=True,text=True);record.update(stdout=reply.stdout,stderr=reply.stderr,returncode=reply.returncode)
 job=reply.stdout.strip().split(';')[0]
 if reply.returncode or not job.isdigit():
  record['state']='uncertain';write(record_path,record);raise RuntimeError('Ambiguous submission reply retained; do not blindly retry')
 record.update(state='submitted',job_id=job);write(record_path,record);journal['pools'][key]=record;write(journal_path,journal)
 print(json.dumps(dict(pool=key,job_id=job,datasets=sum(len(t['datasets']) for t in pool['workers']))),flush=True)
print(json.dumps(dict(status='submitted',pools=len(journal['pools']),datasets=118,reused=2)),flush=True)
