"""Supplement the PD-selected campaign with default SZ returned-matrix risk.
Read saved metrics only. Independent replicates are the sampling units.
"""
import collections,csv,json,math,pathlib,statistics,sys,tomllib
root,out=map(pathlib.Path,sys.argv[1:]);out.mkdir(exist_ok=False)
parent=json.loads((root/'outputs/paper-v4/latest-report.json').read_text());assert parent['exit_code']==0
(out/'validated-parent.json').write_text(json.dumps(parent,indent=2)+'\n')
quantities=('operator','frobenius','normalized_squared_frobenius','fpr','tpr')
baselines=('scaled_lasso','scaled_lasso_tuned','DFGM_C','DFGM_G','aclime','l2_aclime','DGM_C','DGM_G','GGL','Naive_C','Naive_G')
values=collections.defaultdict(dict);coverage=[];applicable={}
def read(path):
 with path.open(newline='') as f:yield from csv.DictReader(f)
def save(name,data):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
for s in range(1,13):
 for rep in range(1,101):
  folder=root/f'outputs/paper-v4/setting-{s}-rep-{rep}'
  receipt=tomllib.loads((folder/'complete.toml').read_text());config=tomllib.loads((folder/'settings.toml').read_text())
  assert receipt['status'] in ('completed','completed_with_failures')
  assert (config['setting_index'],config['replicate'],config['seed'])==(s,rep,2026102500+10000*s+rep)
  if rep==1:applicable[s]=set(config['expected_candidates'])
  q=len(config['evaluation_grid']);fit=next(r for r in read(folder/'fits.csv') if r['method']=='scaled_lasso' and r['candidate']=='1')
  with (folder/'selected.csv').open(newline='') as f:choices={(m,t):k for m,t,k,status in csv.reader(f) if status=='selected'}
  records=collections.defaultdict(list)
  for r in read(folder/'metrics.csv'):
   m=r['method']
   if m=='scaled_lasso' and r['candidate']=='1':records['scaled_lasso_returned',r['target']].append(r)
   if m in baselines and (rep<=20 or m not in ('aclime','l2_aclime')) and choices.get((m,r['time_index']))==r['candidate']:records[m,r['target']].append(r)
  returned=records['scaled_lasso_returned','actual'];assert len(returned)==q if fit['status']=='ok' else not returned
  coverage.append(dict(setting=s,replicate=rep,source=str(folder.resolve()),seed=config['seed'],expected_times=q,
   default_status=fit['status'],returned_times=len(returned),positive_definite_times=int(fit['pd_slices']),
   finite_BIC_times=sum(math.isfinite(float(r['bic'])) for r in returned),default_fit_seconds=fit['seconds'],
   complete_returned_curve=len(returned)==q,complete_PD_curve=int(fit['pd_slices'])==q))
  for (m,target),rs in records.items():
   for quantity in quantities:
    xs=[float(r[quantity]) for r in rs]
    if len(xs)==q and all(map(math.isfinite,xs)):values[s,m,target,quantity][rep]=statistics.mean(xs)
summary=[];paired=[];repvalues=[]
for s in range(1,13):
 for m in ('scaled_lasso_returned',*(b for b in baselines if b in applicable[s])):
  intended=20 if m in ('aclime','l2_aclime') else 100
  for target in ('actual','nominal'):
   for quantity in quantities:
    vs=values[s,m,target,quantity];xs=list(vs.values())
    summary.append(dict(setting=s,method=m,target=target,quantity=quantity,intended_replicates=intended,valid_replicates=len(xs),missing_replicates=intended-len(xs),mean=statistics.mean(xs) if xs else math.nan,mcse=statistics.stdev(xs)/math.sqrt(len(xs)) if len(xs)>1 else math.nan))
    for rep,value in sorted(vs.items()):repvalues.append(dict(setting=s,replicate=rep,method=m,target=target,quantity=quantity,value=value))
 for baseline in (b for b in baselines if b in applicable[s]):
  intended=20 if baseline in ('aclime','l2_aclime') else 100
  for target in ('actual','nominal'):
   for quantity in quantities:
    first=values[s,'scaled_lasso_returned',target,quantity];second=values[s,baseline,target,quantity]
    ids=sorted({i for i in first.keys()&second.keys() if i<=intended});diffs=[first[i]-second[i] for i in ids]
    paired.append(dict(setting=s,first='scaled_lasso_returned',second=baseline,target=target,quantity=quantity,intended_replicates=intended,first_valid=sum(i<=intended for i in first),second_valid=sum(i<=intended for i in second),paired_replicates=len(ids),mean_difference=statistics.mean(diffs) if diffs else math.nan,mcse=statistics.stdev(diffs)/math.sqrt(len(diffs)) if len(diffs)>1 else math.nan))
save('summary.csv',summary);save('paired-comparisons.csv',paired);save('coverage.csv',coverage);save('replicate-values.csv',repvalues)
assert len(coverage)==1200 and len({(r['setting'],r['replicate']) for r in coverage})==1200
(out/'receipt.json').write_text(json.dumps({'status':'completed','datasets':1200,'validated_original_report':parent,'applicable_methods':{s:sorted(ms) for s,ms in applicable.items()},'graph_scope':'Selected FPR/TPR for every applicable original comparator, with replicate-level MCSE and undefined rates retained as missing.','scope':'All100 default SZ returned matrices per12 paper settings, including finite indefinite matrices. Baseline BIC selections remain separately labeled; ACLIME/l2 use fixedfirst20. No fits or truth tuning.','actual_target_primary':True,'times_averaged_within_replicate':True,'paired_MCSE':'sample SD of independent replicate differences / sqrt(paired n)','limitations':['PD-selected baselines and returned default SZ differ in selection population; paired n remains explicit.','Undefined graph rates are excluded per quantity, never imputed.','Default saved fit duration is its original API wall time, not the six-penalty search or an algorithm-only speedup.'],'returned_complete':sum(r['complete_returned_curve'] for r in coverage),'PD_complete':sum(r['complete_PD_curve'] for r in coverage),'failed_default_fits':sum(r['default_status']!='ok' for r in coverage)},indent=2)+'\n')
