"""Paired new-method report from saved CSVs; never fits or resimulates.

Run under Slurm. Scope is exactly 12x20 paper and 16x100 control seeds.
Returned finite matrices and PD-screened selections remain separate populations.
"""
import collections,csv,json,math,pathlib,statistics,sys,tomllib

root,output=map(pathlib.Path,sys.argv[1:3]);assert not output.exists();output.mkdir()
new=root/'outputs/spectral-v5'
values=collections.defaultdict(dict);coverage=[];failures=[];runtimes=[];pdrows=[];rocrows=[];baseline_times=[]
quantities=('operator','frobenius','normalized_squared_frobenius','fpr','tpr')
newmethods=('spectral_l2_aclime','spectral_l2_aclime_tuned','spectral_l2_aclime_PD')
expected={(kind,s):n for kind,ns,n in [('paper',12,20),('control',16,100)] for s in range(1,ns+1)}
interim='--interim-paper50-controls' in sys.argv
if interim:expected={key:n for key,n in expected.items() if key[0]=='control' or key[1] in (1,2,3,7,8,9)}

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
            coverage.append(dict(kind=kind,setting=s,replicate=rep,status='incomplete',source='',independent_check='missing',baseline_status='not_checked'));continue
        config=tomllib.loads((folder/'settings.toml').read_text());receipt=tomllib.loads((folder/'complete.toml').read_text())
        assert (config['kind'],config['setting_index'],config['replicate'])==(kind,s,rep)
        assert config['seed']==2026102500+10000*s+rep and not config['point_pilot']
        assert config['score_menu']==[.01,.03,.1,.3,1.,2.] and config['grid_mode']=='exact_theory_grid'
        assert config['Aconf']==1 and config['C_K']==1 and config['centering']=='known_zero'
        if kind=='paper':assert config['n']==200 and config['noise_variance']==.5 and config['L']==1 and config['kappa_bound']==10
        else:assert config['n']==160 and config['p']==50 and config['kappa_bound']==2
        independent=tomllib.loads((folder/'independent-check.toml').read_text())
        assert independent['status']=='passed'
        coverage.append(dict(kind=kind,setting=s,replicate=rep,status=receipt['status'],source=config['source_sha256']['src/spectral_direct.jl'],independent_check='passed',baseline_status='not_applicable' if kind=='control' else 'pending'))
        nt=len(config['evaluation_grid']);metrics=rows(folder/'metrics.csv');fitrows=rows(folder/'fits.csv')
        assert len(metrics)==6*nt*2 and len(fitrows)==6*nt
        assert len({(r['target'],r['time_index'],r['candidate']) for r in metrics})==len(metrics)
        for r in rows(folder/'auc.csv'):
            rocrows.append(dict(kind=kind,setting=s,replicate=rep,method=r['method'],target=r['target'],time_index=r['time_index'],
                observed_area=r['observed_area'],full_auc=r['full_auc'],minimum_fpr=r['minimum_fpr'],maximum_fpr=r['maximum_fpr'],
                successful_candidates=r['successful_candidates'],expected_candidates=r['expected_candidates']))
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
        counts=collections.Counter(r['status'] for r in fitrows);columns=rows(folder/'column-status.csv');colcounts=collections.Counter(r['status'] for r in columns)
        assert len(columns)==6*nt*config['p']
        for r in columns:
            if r['status'] in ('optimal','optimal_analytic_zero'):assert finite(r['max_relative_violation']) and float(r['max_relative_violation'])<=1e-5
            if r['status']=='infeasible_zero_fallback':assert r['termination']=='INFEASIBLE' and r['dual_status']=='INFEASIBILITY_CERTIFICATE'
        failures.append(dict(kind=kind,setting=s,replicate=rep,candidate_time_fits=len(fitrows),optimal=counts['optimal'],
            infeasible_fallback=counts['infeasible_zero_fallback'],numerical_failure=counts['numerical_failure'],
            infeasible_columns=colcounts['infeasible_zero_fallback'],numerical_failure_columns=colcounts['numerical_failure']))
        runtime=dict(kind=kind,setting=s,replicate=rep,candidates=6,**{k:config[k] for k in ('threads','time_tasks','solver_threads','hostname','cpu_model','path_reuse','warmup','covariance_seconds','calibration_seconds')},
            source=config['source_sha256']['src/spectral_direct.jl'],elapsed_seconds=receipt['elapsed_seconds'],
            recorded_fit_wall_seconds_UNRELIABLE=receipt['fit_wall_seconds'],recorded_path_sum_seconds_UNRELIABLE=receipt['fit_path_sum_seconds'],
            aggregate_precision_timing_status='valid' if config.get('timing_schema')=='private_point_wall_clocks_v2' else 'unavailable_shared_clock_bug_in_frozen_v4_v5; use valid column/model/solve records',
            precision_wall_seconds=receipt['fit_wall_seconds'] if config.get('timing_schema')=='private_point_wall_clocks_v2' else math.nan,
            scoring_seconds=receipt['scoring_seconds'])
        for name in ('precision_seconds','model_setup_seconds','solve_seconds','residual_seconds','symmetry_seconds','cap_seconds'):
            runtime[name]=sum(float(r[name]) for r in fitrows)
        runtimes.append(runtime)
        if kind=='paper':
            # Recovery can replace an incomplete original producer. Its cached
            # data provenance remains in identity.baseline_folder; comparisons
            # must read the current completed canonical baseline, never refit it.
            baseline=(root/f'outputs/paper-v4/setting-{s}-rep-{rep}').resolve()
            if not (baseline/'complete.toml').exists():
                coverage[-1]['baseline_status']='incomplete'
                assert '--require-complete' not in sys.argv, f'Baseline still incomplete: {baseline}'
                continue
            coverage[-1]['baseline_status']='complete'
            oldconfig=tomllib.loads((baseline/'settings.toml').read_text())
            assert oldconfig['seed']==config['seed'] and oldconfig['evaluation_grid']==config['evaluation_grid'] and oldconfig['n']==200 and oldconfig['p']==config['p']
            old=rows(baseline/'metrics.csv')
            choices={(m,t):k for m,t,k,status in csv.reader((baseline/'selected.csv').open()) if status=='selected'}
            for method in sorted({r['method'] for r in old}):
                chosen=[r for r in old if r['method']==method and choices.get((method,r['time_index']))==r['candidate']]
                intended_selected={t for (m,t) in choices if m==method}
                assert len(chosen)==2*len(intended_selected), f'Selected metric coverage differs for {method}, {baseline}'
                for target in ('actual','nominal'):
                    for q in quantities:average_times(kind,s,rep,method,target,q,chosen,nt)
            chosen=[r for r in old if r['method']=='scaled_lasso' and r['candidate']=='1']
            for target in ('actual','nominal'):
                for q in quantities:average_times(kind,s,rep,'scaled_lasso_returned',target,q,chosen,nt)
            for r in rows(baseline/'auc.csv'):
                rocrows.append(dict(kind=kind,setting=s,replicate=rep,method=r['method'],target=r['target'],time_index=r['time_index'],
                    observed_area=r['observed_area'],full_auc=r['full_auc'],minimum_fpr=r['minimum_fpr'],maximum_fpr=r['maximum_fpr'],
                    successful_candidates=r['successful_candidates'],expected_candidates=r['expected_candidates']))
            times=collections.defaultdict(list)
            for r in rows(baseline/'fits.csv'):times[r['method']].append(r)
            for m,rs in times.items():
                baseline_times.append(dict(kind=kind,setting=s,replicate=rep,method=m,candidates=len(rs),
                    saved_menu_seconds=sum(float(r['seconds']) for r in rs),failed_candidates=sum(r['status']=='failed' for r in rs),
                    source=str(baseline),timing_scope='legacy stored whole-path/candidate durations; includes path/pilot/symmetrization and source-specific reuse; not an algorithm-only controlled benchmark'))
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
    # Prespecify method coverage, including methods with no valid selections.
    # Inferring this list from successful values would conceal missingness.
    baselines=['aclime','l2_aclime','scaled_lasso','scaled_lasso_tuned']
    if kind=='paper':
        baselines+=['DFGM_C','DFGM_G','scaled_lasso_returned']
        if s<=6:baselines+=['DGM_C','DGM_G','GGL','Naive_C','Naive_G']
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
if rocrows:save('roc-coverage.csv',rocrows)
rocgroups=collections.defaultdict(list)
for r in rocrows:rocgroups[r['kind'],r['setting'],r['method'],r['target']].append(r)
rocmissing=[]
for (kind,s,method,target),rs in sorted(rocgroups.items()):
    nt=1 if kind=='control' else ({1:10,2:25,3:50,4:10,5:25,6:50}.get(s,21))
    planned=expected[kind,s]*nt
    rocmissing.append(dict(kind=kind,setting=s,method=method,target=target,intended_times=planned,recorded_times=len(rs),
        complete_candidate_path_times=sum(r['successful_candidates']==r['expected_candidates'] for r in rs),
        finite_observed_area_times=sum(finite(r['observed_area']) for r in rs),finite_full_auc_times=sum(finite(r['full_auc']) for r in rs),
        missing_full_auc_times=planned-sum(finite(r['full_auc']) for r in rs)))
if rocmissing:save('roc-missingness.csv',rocmissing)
if baseline_times:save('baseline-saved-fit-times.csv',baseline_times)
if '--require-complete' in sys.argv:
    assert not interim and len(coverage)==1840 and all(r['status']!='incomplete' and r['independent_check']=='passed' and r['baseline_status'] in ('complete','not_applicable') for r in coverage), 'Full1840 coverage not verified'
(output/'receipt.json').write_text(json.dumps({'datasets':len(coverage),'intended':sum(expected.values()),'overall_study_intended':1840,'interim_paper50_controls':interim,'cohorts':list(map(list,expected)),
    'completed':sum(r['status']!='incomplete' for r in coverage),
    'objective':'Compare our implementations; actual precision primary, no published-table reproduction objective.',
    'scope':'Paper12x20, controls16x100; times averaged within independent seeds before MCSE. No baseline fits repeated.',
    'populations':{'spectral_l2_aclime':'fixed Kscore2 returned optimal finite matrices, regardless PD','spectral_l2_aclime_PD':'same fixed matrices after PD screening',
        'spectral_l2_aclime_tuned':'six-score BIC selection for paper, heldout residual for controls','scaled_lasso_returned':'fixed default SZ before PD screening, paper only',
        'other_paper_baselines':'existing stored BIC selections'},
    'limitations':['Control baseline graph and PD properties were not saved; marked unavailable, no refits.',
        'Paired complete cases can be selective. Empirical constants do not certify theorem assumptions.',
        'Frozen v4/v5 parallel aggregate precision wall/path clocks are unreliable; their recorded values are explicitly labeled. Per-column model/solve/residual and covariance/scoring timings remain valid. Use the separate controlled timing experiment for wall-time comparisons.',
        'Shared menus/stages cannot be called algorithm-only speedups.']},indent=2)+'\n')
