"""Paired new-method report from saved CSVs; never fits or resimulates.

Run under Slurm. Scope is exactly 12x20 paper and 16x100 control seeds.
Returned finite matrices and PD-screened selections remain separate populations.
"""
import collections,csv,json,math,pathlib,statistics,sys,tomllib

root,output=map(pathlib.Path,sys.argv[1:3]);assert not output.exists();output.mkdir()
new=root/'outputs/spectral-v5'
values=collections.defaultdict(dict);coverage=[];failures=[];runtimes=[];pdrows=[]
quantities=('operator','frobenius','normalized_squared_frobenius','fpr','tpr')
newmethods=('spectral_l2_aclime','spectral_l2_aclime_tuned','spectral_l2_aclime_PD')
expected={(kind,s):n for kind,ns,n in [('paper',12,20),('control',16,100)] for s in range(1,ns+1)}

def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def save(name,data,fields=None):
    with (output/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(data[0]));w.writeheader();w.writerows(data)
def finite(x):
    try:return math.isfinite(float(x))
    except (ValueError,TypeError):return False

def average_times(kind,s,rep,method,target,quantity,records,nt):
    xs=[float(r[quantity]) for r in records if r['target']==target]
    if len(xs)==nt and all(map(math.isfinite,xs)):
        values[kind,s,method,target,quantity][rep]=statistics.mean(xs)

for (kind,s),nreps in expected.items():
    for rep in range(1,nreps+1):
        folder=new/f'{kind}-{s}-{rep}'
        if not (folder/'complete.toml').exists():
            coverage.append(dict(kind=kind,setting=s,replicate=rep,status='incomplete',source='',independent_check='missing'));continue
        config=tomllib.loads((folder/'settings.toml').read_text());receipt=tomllib.loads((folder/'complete.toml').read_text())
        assert (config['kind'],config['setting_index'],config['replicate'])==(kind,s,rep)
        assert config['seed']==2026102500+10000*s+rep and not config['point_pilot']
        assert config['score_menu']==[.01,.03,.1,.3,1.,2.] and config['grid_mode']=='exact_theory_grid'
        independent=tomllib.loads((folder/'independent-check.toml').read_text())
        assert independent['status']=='passed'
        coverage.append(dict(kind=kind,setting=s,replicate=rep,status=receipt['status'],source=config['source_sha256']['src/spectral_direct.jl'],independent_check='passed'))
        nt=len(config['evaluation_grid']);metrics=rows(folder/'metrics.csv');fitrows=rows(folder/'fits.csv')
        assert len(metrics)==6*nt*2 and len(fitrows)==6*nt
        selected={int(t):int(k) for method,t,k,status in csv.reader((folder/'selected.csv').open()) if method=='spectral_l2_aclime_tuned' and status=='selected'}
        for method in newmethods:
            chosen=[r for r in metrics if r['status']=='optimal' and int(r['candidate'])==(selected.get(int(r['time_index']),0) if method=='spectral_l2_aclime_tuned' else 6)
                    and (method!='spectral_l2_aclime_PD' or r['pd']=='true')]
            for target in ('actual','nominal'):
                for q in quantities:average_times(kind,s,rep,method,target,q,chosen,nt)
            # PD means a property of the returned matrix, never an imputed fit.
            unique=[r for r in chosen if r['target']=='actual']
            pdrows.append(dict(kind=kind,setting=s,replicate=rep,method=method,intended_times=nt,
                returned_optimal_times=len(unique),positive_definite_times=sum(r['pd']=='true' for r in unique),
                all_times_returned=len(unique)==nt,all_times_positive_definite=len(unique)==nt and all(r['pd']=='true' for r in unique)))
        counts=collections.Counter(r['status'] for r in fitrows);colcounts=collections.Counter(r['status'] for r in rows(folder/'column-status.csv'))
        failures.append(dict(kind=kind,setting=s,replicate=rep,candidate_time_fits=len(fitrows),optimal=counts['optimal'],
            infeasible_fallback=counts['infeasible_zero_fallback'],numerical_failure=counts['numerical_failure'],
            infeasible_columns=colcounts['infeasible_zero_fallback'],numerical_failure_columns=colcounts['numerical_failure']))
        runtime=dict(kind=kind,setting=s,replicate=rep,candidates=6,**{k:config[k] for k in ('threads','time_tasks','solver_threads','hostname','cpu_model','path_reuse','warmup','covariance_seconds','calibration_seconds')},
            source=config['source_sha256']['src/spectral_direct.jl'],elapsed_seconds=receipt['elapsed_seconds'],fit_wall_seconds=receipt['fit_wall_seconds'],
            path_sum_seconds=receipt['fit_path_sum_seconds'],scoring_seconds=receipt['scoring_seconds'])
        for name in ('precision_seconds','model_setup_seconds','solve_seconds','residual_seconds','symmetry_seconds','cap_seconds'):
            runtime[name]=sum(float(r[name]) for r in fitrows)
        runtimes.append(runtime)
        if kind=='paper':
            baseline=pathlib.Path(config['identity']['baseline_folder']);old=rows(baseline/'metrics.csv')
            choices={(m,t):k for m,t,k,status in csv.reader((baseline/'selected.csv').open()) if status=='selected'}
            for method in sorted({r['method'] for r in old}):
                chosen=[r for r in old if choices.get((method,r['time_index']))==r['candidate']]
                for target in ('actual','nominal'):
                    for q in quantities:average_times(kind,s,rep,method,target,q,chosen,nt)
            chosen=[r for r in old if r['method']=='scaled_lasso' and r['candidate']=='1']
            for target in ('actual','nominal'):
                for q in quantities:average_times(kind,s,rep,'scaled_lasso_returned',target,q,chosen,nt)
        else:
            for r in rows(root/f'outputs/controls100/stress-{s}/stress.csv'):
                if int(r['replicate'])!=rep or r['status']!='ok' or r['selected']!='true':continue
                for q in quantities[:3]:
                    value=math.sqrt(50*float(r['normalized_squared_frobenius'])) if q=='frobenius' else float(r[q])
                    values[kind,s,r['method'],'actual',q][rep]=value
            for r in rows(root/f'outputs/scaled100/calibration-{s}/fits.csv'):
                if int(r['replicate'])!=rep or r['status']!='ok' or r['selected']!='true':continue
                for q in quantities[:3]:
                    value=math.sqrt(50*float(r['normalized_squared_frobenius'])) if q=='frobenius' else float(r[q])
                    values[kind,s,'scaled_lasso_tuned','actual',q][rep]=value

summary=[];paired=[]
for (kind,s),n in expected.items():
    baselines=sorted({key[2] for key in values if key[0]==kind and key[1]==s and key[2] not in newmethods})
    for m in (*newmethods,*baselines):
        for target in (('actual','nominal') if kind=='paper' else ('actual',)):
            for q in quantities:
                v=values[kind,s,m,target,q];x=list(v.values())
                summary.append(dict(kind=kind,setting=s,method=m,target=target,quantity=q,intended_replicates=n,valid_replicates=len(x),missing_replicates=n-len(x),
                    mean=statistics.mean(x) if x else math.nan,mcse=statistics.stdev(x)/math.sqrt(len(x)) if len(x)>1 else math.nan))
    for m in newmethods:
        for b in baselines:
            for target in (('actual','nominal') if kind=='paper' else ('actual',)):
                for q in quantities:
                    a=values[kind,s,m,target,q];bv=values[kind,s,b,target,q];ids=sorted(a.keys()&bv.keys());diff=[a[i]-bv[i] for i in ids]
                    paired.append(dict(kind=kind,setting=s,new_method=m,baseline=b,target=target,quantity=q,intended_replicates=n,
                        new_valid=len(a),baseline_valid=len(bv),paired_replicates=len(ids),new_minus_baseline=statistics.mean(diff) if diff else math.nan,
                        mcse=statistics.stdev(diff)/math.sqrt(len(diff)) if len(diff)>1 else math.nan))
save('coverage.csv',coverage);save('summary.csv',summary);save('paired-comparisons.csv',paired)
if failures:save('failures.csv',failures);save('runtime.csv',runtimes);save('pd-coverage.csv',pdrows)
(output/'receipt.json').write_text(json.dumps({'datasets':len(coverage),'intended':1840,'completed':sum(r['status']!='incomplete' for r in coverage),
    'objective':'Compare our implementations; actual precision primary, no published-table reproduction objective.',
    'scope':'Paper12x20, controls16x100; times averaged within independent seeds before MCSE. No baseline fits repeated.',
    'populations':{'spectral_l2_aclime':'fixed Kscore2 returned optimal finite matrices, regardless PD','spectral_l2_aclime_PD':'same fixed matrices after PD screening',
        'spectral_l2_aclime_tuned':'six-score BIC selection for paper, heldout residual for controls','scaled_lasso_returned':'fixed default SZ before PD screening, paper only',
        'other_paper_baselines':'existing stored BIC selections'},
    'limitations':['Control baseline graph and PD properties were not saved; marked unavailable, no refits.',
        'Paired complete cases can be selective. Empirical constants do not certify theorem assumptions.',
        'Fit timing includes first-call compilation and path reuse; shared menus/stages cannot be called algorithm-only speedups.']},indent=2)+'\n')
