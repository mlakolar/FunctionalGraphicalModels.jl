#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
source_dir=$root/scaled-original-v2-source
out=$root/outputs/scope-amendment-gate-$SLURM_JOB_ID
mkdir "$out"
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_GC_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_HEAP_SIZE_HINT=2500M
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$source_dir"
julia=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
"$julia" --project="$root/common-fast-v2-source" -e 'include("experiments/dfgm2020.jl"); for (s,d) in ((6,:common),(12,:irregular)); folder=joinpath(ARGS[1],"outputs/paper-v4/setting-$s-rep-21"); c=load_dfgm_covariance_cache(folder;setting_index=s,replicate=21,seed=2026102500+10000s+21,n=200,p=100,design=d,schedule=:T50,mode=:paper,integration_points=101,work_points=31); @assert size(c.covariance.direct.matrices,1)==100; rejected=false; try load_dfgm_covariance_cache(folder;setting_index=s,replicate=21,seed=0,n=200,p=100,design=d,schedule=:T50,mode=:paper,integration_points=101,work_points=31); catch; rejected=true; end; @assert rejected; end; println("Both frozen covariance variants accepted with exact provenance; wrong seed rejected")' "$root"
for setting in 1 7; do
    extra=()
    if [[ $setting == 7 ]]; then extra+=(--evaluation-points 3); fi
    code=0
    "$julia" --project="$root/common-fast-v2-source" experiments/dfgm2020.jl --setting "$setting" --replicate 1001 \
        --n 30 --p 10 --mode pilot --integration-points 31 --work-points 15 --nlambda 3 \
        --method-subset scaled_original --common-fast "$((setting == 1))" \
        --out "$out/canary-$setting" "${extra[@]}" || code=$?
    printf '%s\t%s\n' "$setting" "$code" >> "$out/worker-exits.tsv"
done
python3 - "$out" "$source_dir" "$root" <<'PY'
import copy,csv,hashlib,json,pathlib,shutil,subprocess,sys,tomllib
out,source,root=map(pathlib.Path,sys.argv[1:]);sys.path.insert(0,str(source/'experiments'))
from validate_dfgm_outputs import source_settings,validate
for setting,expected in ((1,28),(7,13)):
    folder=out/f'canary-{setting}'
    settings=tomllib.loads((folder/'settings.toml').read_text())
    assert settings['method_subset']=='scaled_original'
    assert not {'aclime','l2_aclime'}.intersection(settings['expected_candidates'])
    assert sum(settings['expected_candidates'].values())==expected
    view=out/f'view-{setting}';view.mkdir();(view/f'setting-{setting}-rep-1001').symlink_to(folder)
    (view/'tasks.tsv').write_text(f'{setting}\t1001\n')
    assert validate(view,allow_recorded_fit_failures=True)
locks=json.loads((source/'production-source-locks.json').read_text())
settings=tomllib.loads((root/'outputs/paper-v4/setting-9-rep-21/settings.toml').read_text())
new=next(v for v in locks['variants'] if v['name']=='scaled-original-v2-source')
settings['source_sha256']=new['source_sha256'];settings['method_subset']='scaled_original'
settings['expected_candidates'].pop('aclime');settings['expected_candidates'].pop('l2_aclime')
assert source_settings(settings,locks)==new['name']
bad=copy.deepcopy(settings);bad['replicate']=1;bad['seed']=2026102500+90000+1
try:source_settings(bad,locks)
except AssertionError:pass
else:raise AssertionError('first20 incorrectly accepted as subset')
view=out/'report-scope';view.mkdir();shutil.copyfile(source/'comparison-scope.json',view/'comparison-scope.json')
for rep in (1,21):
    prior=root/f'outputs/paper-v4/setting-9-rep-{rep}';dest=view/prior.name;dest.mkdir()
    for name in ('settings.toml','complete.toml','selected.csv','metrics.csv','auc.csv'):
        shutil.copyfile(prior/name,dest/name)
(view/'tasks.tsv').write_text('9\t1\n9\t21\n')
def report():
    for script in ('summarize_dfgm.py','paired_dfgm.py'):
        subprocess.run([sys.executable,str(source/'experiments'/script),str(view)],check=True)
    return {n:(view/n).read_bytes() for n in ('summary.csv','paired-comparisons.csv')}
before=report()
with (view/'summary.csv').open(newline='') as f:
    assert all(int(r['valid_replicates'])==1 for r in csv.DictReader(f) if r['method'] in ('aclime','l2_aclime'))
with (view/'paired-comparisons.csv').open(newline='') as f:
    rows=list(csv.DictReader(f))
    for r in rows:assert int(r['expected_replicates'])==(1 if {'aclime','l2_aclime'}.intersection((r['first'],r['second'])) else 2)
path=view/'setting-9-rep-21/metrics.csv'
with path.open(newline='') as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
for row in rows:
    if row['method'] in ('aclime','l2_aclime'):
        for key in ('operator','frobenius','normalized_squared_frobenius'):row[key]='1000000000'
with path.open('w',newline='') as f:writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
assert report()==before,'excluded seed changed primary reports'
(out/'gate-receipt.json').write_text(json.dumps({'status':'passed','scope':'New subset driver, two covariance provenance checks, source rejection guard and primary paired-report exclusion. Fresh reduced canaries only; no primary refits.'},indent=2)+'\n')
PY
