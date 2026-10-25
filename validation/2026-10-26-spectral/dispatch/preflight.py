import pathlib,tomllib,hashlib,json,datetime
r=pathlib.Path('/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct')
a=r/'production-v1-source/experiments/direct_stress.jl';b=r/'common-fast-v2-source/experiments/direct_stress.jl'
assert a.read_bytes()==b.read_bytes()
checks=[]
configs=[(d,m,c,v) for d in (2,10) for m in (8,30) for c in (0.,4.) for v in (.5,2.)]
for i,(d,m,c,v) in enumerate(configs,1):
 p=r/f'outputs/controls100/stress-{i}/settings.toml';meta=tomllib.loads(p.read_text())
 expected=dict(index=i,degree=d,m=m,curvature=c,noise_variance=v,p=50,n=200,replicates=100,seed_recipe='2026102500 + 10000*index + replicate')
 assert all(meta[k]==value for k,value in expected.items())
 checks.append(dict(index=i,settings_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
out=r/'outputs/spectral-v5';out.mkdir(exist_ok=True)
for i in (1,16):
 source=r/f'outputs/spectral-v4-canary/control-{i}-1';dest=out/f'control-{i}-1'
 assert tomllib.loads((source/'independent-check.toml').read_text())['status']=='passed'
 assert (source/'complete.toml').exists()
 if dest.is_symlink():assert dest.resolve()==source.resolve()
 else:assert not dest.exists();dest.symlink_to(source)
receipt=dict(status='passed',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 controls=checks,generating_function_sha256=hashlib.sha256(a.read_bytes()).hexdigest(),
 scope='Exact original control source/configuration/seed recipe, known split160/40. No baseline fits or simulations repeated.')
(r/'outputs/spectral-control-pairing-preflight.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
