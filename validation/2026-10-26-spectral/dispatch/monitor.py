import pathlib,json,subprocess,collections,datetime,tomllib,csv
R=pathlib.Path('/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct');B=R/'outputs/spectral-v5';D=B/'dispatch';ids=[];pools=[]
for f in ('control-submission-journal.json','paper50-submission-journal.json','paper100-submission-journal.json'):
 if (D/f).exists():
  j=json.loads((D/f).read_text());ids += [p['job_id'] for p in j['pools'].values()];pools+=list(j['pools'].values())
q=subprocess.check_output(['squeue','-h','-j',','.join(ids),'-o','%i|%T|%C|%r|%M|%l'],text=True).splitlines();rows=[x.split('|') for x in q]
a=subprocess.check_output(['sacct','-X','-j',','.join(ids),'--format=JobIDRaw,State,Elapsed,ExitCode','--parsable2'],text=True)
account=list(csv.DictReader(a.splitlines(),delimiter='|'));coverage=collections.defaultdict(collections.Counter);pipeline=[]
for kind,ns,nr in [('paper',12,20),('control',16,100)]:
 for s in range(1,ns+1):
  for rep in range(1,nr+1):
   p=B/f'{kind}-{s}-{rep}'
   if kind=='paper' and s in (6,12) and rep==1 and not p.exists():p=R/f'outputs/spectral-v4-canary/paper-{s}-{rep}'
   key=f'{kind}-{s}'
   if (p/'pipeline-failure.txt').exists():pipeline.append(f'{kind}:{s}:{rep}')
   if not (p/'complete.toml').exists():coverage[key]['unfinished']+=1;continue
   receipt=tomllib.loads((p/'complete.toml').read_text());coverage[key]['complete']+=1
   coverage[key]['candidate_time_failures']+=receipt['candidate_time_failures']
   if (p/'independent-check.toml').exists() and tomllib.loads((p/'independent-check.toml').read_text())['status']=='passed':coverage[key]['verified']+=1
sample_ids=[x[0] for x in rows if x[1]=='RUNNING'][:2]
resource=subprocess.check_output(['sstat','--allsteps','--jobs='+','.join(sample_ids),'--format=JobID,MaxRSS,AveCPU,AveRSS','--parsable2'],text=True) if sample_ids else ''
r=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),job_ids=ids,queue=rows,accounting=account,
 states=dict(collections.Counter(x[1] for x in rows)),running_cpus=sum(int(x[2]) for x in rows if x[1]=='RUNNING'),
 coverage={k:dict(v) for k,v in coverage.items()},pipeline_failures=pipeline,sampled_resources=resource,
 scope='Fixed known output keys and exact owned IDs; valid completed_with_failures outputs are not refit.')
(D/'latest-monitor.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:v for k,v in r.items() if k not in ('job_ids','queue','accounting')},indent=2))
