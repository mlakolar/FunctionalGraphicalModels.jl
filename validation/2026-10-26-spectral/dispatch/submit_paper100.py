"""Idempotent fixed-scope dispatch; intent is durable before each sbatch."""
import datetime,hashlib,json,pathlib,subprocess,tomllib
R=pathlib.Path('/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct');B=R/'outputs/spectral-v5';D=R/'spectral-dispatch-v1';M=B/'dispatch';M.mkdir(exist_ok=True)
NUM=R/'spectral-20261026-v5-source'
def write(p,x):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
# Full p100 canaries must be complete and independently checked. A high
# numerical-failure rate blocks expansion; certified infeasibility is preserved.
import math,csv
reused=[(6,1),(12,1)];elapsed=[];cost=[]
for setting,rep in reused:
 source=R/f'outputs/spectral-v4-canary/paper-{setting}-{rep}'
 receipt=tomllib.loads((source/'complete.toml').read_text())
 assert receipt['status'] in ('completed','completed_with_failures')
 assert tomllib.loads((source/'independent-check.toml').read_text())['status']=='passed'
 with (source/'fits.csv').open() as stream:fits=list(csv.DictReader(stream))
 numerical=sum(x['status']=='numerical_failure' for x in fits)
 assert numerical<=.05*len(fits), 'Numerical failures exceed prespecified5% expansion gate; diagnose before submitting'
 elapsed.append(receipt['elapsed_seconds']);cost.append(dict(setting=setting,replicate=rep,elapsed_seconds=receipt['elapsed_seconds'],candidate_times=len(fits),numerical_failures=numerical))
 dest=B/f'paper-{setting}-{rep}'
 if dest.is_symlink():assert dest.resolve()==source.resolve()
 else:assert not dest.exists();dest.symlink_to(source)
wall_hours=max(1,min(8,math.ceil((1.5*max(elapsed)+600)/3600)))
write(D/'p100-cost-gate.json',dict(status='passed',cases=cost,wall_hours_per_dataset=wall_hours,formula='ceil(1.5*max observed full-case seconds+600)/3600; min1h,max8h; scheduling estimate, not a statistical guarantee'))
# Numerical sources remain byte-identical to the passed v5 check/case source.
manifest=json.loads((NUM/'source.json').read_text())
for f,digest in manifest['files'].items():assert hashlib.sha256((NUM/f).read_bytes()).hexdigest()==digest
plan_path=M/'paper100-plan.json'
if not plan_path.exists():
 tasks=[]
 for s in (4,5,6,10,11,12):
  for start in range(1,21):
   cases=[('paper',s,rep) for rep in range(start,start+1) if (s,rep) not in reused]
   task=M/f'paper100-worker-{len(tasks)+1:03d}.tsv';assert not task.exists()
   task.write_text(''.join(f'{k}\t{i}\t{rep}\n' for k,i,rep in cases));tasks.append(dict(task=str(task),datasets=cases))
 pools=[]
 for i in range(0,len(tasks),4):
  index=i//4+1;tag=f'dfgm-spectral-paper100-v5-{index:03d}';tsv=M/f'paper100-pool-{index:03d}.tsv'
  tsv.write_text(''.join(t['task']+'\n' for t in tasks[i:i+4]))
  pools.append(dict(index=index,tag=tag,tsv=str(tsv),workers=tasks[i:i+4]))
 plan=dict(pools=pools,expected_datasets=120,reused_cases=reused,new_datasets=118,wall_hours=wall_hours,numerical_source=str(NUM),
  orchestration_sha256={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ('spectral_batch.jl','paper-worker.sh','paper100.sbatch')})
 assert sum(len(t['datasets']) for p in pools for t in p['workers'])==118
 write(plan_path,plan)
plan=json.loads(plan_path.read_text());journal_path=M/'paper100-submission-journal.json'
journal=json.loads(journal_path.read_text()) if journal_path.exists() else dict(pools={},source=str(NUM),scope='six p100 paper settings x20 paired seeds; six-score exact-grid protocol; no baseline fits')
for pool in plan['pools']:
 key=str(pool['index']);record_path=M/f'paper100-submission-{key}.json'
 if record_path.exists():
  record=json.loads(record_path.read_text())
  if str(record.get('job_id','')).isdigit():journal['pools'][key]=record;write(journal_path,journal);continue
  raise RuntimeError(f'Uncertain previous submission; reconcile tag before retry: {record_path}')
 cmd=['sbatch','--parsable',f'--time={wall_hours}:00:00',f'--ntasks={len(pool["workers"])}',f'--mem={8*len(pool["workers"])}G',f'--job-name={pool["tag"]}',f'--comment={pool["tag"]}',
      '-o',str(M/f'paper100-pool-{pool["index"]:03d}-%j.out'),'-e',str(M/f'paper100-pool-{pool["index"]:03d}-%j.err'),str(D/'paper100.sbatch'),pool['tsv']]
 record={**pool,'state':'submitting','submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd}
 write(record_path,record)
 reply=subprocess.run(cmd,capture_output=True,text=True);record.update(stdout=reply.stdout,stderr=reply.stderr,returncode=reply.returncode)
 job=reply.stdout.strip().split(';')[0]
 if reply.returncode or not job.isdigit():
  record['state']='uncertain';write(record_path,record);raise RuntimeError('Ambiguous submission reply retained; do not blindly retry')
 record.update(state='submitted',job_id=job);write(record_path,record);journal['pools'][key]=record;write(journal_path,journal)
 print(json.dumps(dict(pool=key,job_id=job,datasets=sum(len(t['datasets']) for t in pool['workers']))),flush=True)
print(json.dumps(dict(status='submitted',pools=len(journal['pools']),datasets=118,reused=2)),flush=True)
