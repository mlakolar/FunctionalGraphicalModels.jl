"""Plot saved paired reports. Run in the campaign plotting environment under Slurm.
Usage: python plot_spectral_comparison.py REPORT TIMING OUTPUT
"""
import csv, pathlib, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
report,timing,out=map(pathlib.Path,sys.argv[1:4]);out.mkdir(exist_ok=False)
def rows(name):
    with (report/name).open() as f:return list(csv.DictReader(f))
def save(fig,name):
    fig.savefig(out/(name+'.png'),dpi=170,bbox_inches='tight')
    fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
paired=rows('paired-comparisons.csv');summary=rows('summary.csv')
fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
for ax,quantity,label in zip(axes,['operator','normalized_squared_frobenius'],['Operator norm error','Squared Frobenius error / p']):
    rs=sorted([r for r in paired if r['kind']=='control' and r['new_method']=='spectral_l2_aclime_tuned' and r['baseline']=='scaled_lasso_tuned' and r['target']=='actual' and r['quantity']==quantity],key=lambda r:int(r['setting']))
    assert len(rs)==16 and all(int(r['paired_replicates'])==100 for r in rs)
    ax.errorbar([int(r['setting']) for r in rs],[float(r['new_minus_baseline']) for r in rs],yerr=[2*float(r['mcse']) for r in rs],fmt='o',capsize=3)
    ax.axhline(0,color='black',linewidth=.8);ax.set(xlabel='Control configuration',ylabel='Spectral tuned minus SZ tuned',title=label,xticks=range(1,17));ax.grid(alpha=.2)
fig.suptitle('Paired control differences: 100 independent replicates, bars = 2 MCSE')
save(fig,'control-paired-errors')
methods=[('spectral_l2_aclime','Spectral fixed'),('spectral_l2_aclime_tuned','Spectral tuned'),('scaled_lasso_returned','SZ default returned'),('scaled_lasso_tuned','SZ tuned (BIC)'),('l2_aclime','l2-ACLIME (BIC)'),('DFGM_C','DFGM_C (BIC)')]
fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained')
for method,label in methods:
    rs=sorted([r for r in summary if r['kind']=='paper' and r['method']==method and r['target']=='actual' and r['quantity']=='operator'],key=lambda r:int(r['setting']))
    axes[0].errorbar([int(r['setting']) for r in rs],[float(r['mean']) for r in rs],yerr=[2*float(r['mcse']) for r in rs],fmt='o-',label=label,alpha=.8,capsize=2)
    axes[1].plot([int(r['setting']) for r in rs],[int(r['valid_replicates']) for r in rs],'o-',label=label,alpha=.8)
axes[0].set(yscale='symlog',ylabel='Mean operator norm error',title='Actual precision target; time average within replicate, bars = 2 MCSE')
axes[0].legend(ncol=3,fontsize=9);axes[1].set(xlabel='Paper simulation setting',ylabel='Complete replicate curves',ylim=(-.5,20.5),yticks=[0,5,10,15,20])
for ax in axes:ax.set_xticks(range(1,13));ax.grid(alpha=.2)
pending={}
for r in rows('coverage.csv'):
    if r['kind']=='paper' and r['baseline_status']=='incomplete':pending[r['setting']]=pending.get(r['setting'],0)+1
qualification='; '.join(f"Setting {s}: {n} baseline datasets pending" for s,n in pending.items())
fig.suptitle('Paper comparison: first 20 fixed seeds per setting; missing curves are not imputed'+('\n'+qualification if qualification else ''))
save(fig,'paper-errors-and-coverage')
fails=rows('failures.csv');groups={}
for r in fails:
    key=(r['kind'],int(r['setting']));v=groups.setdefault(key,[0,0,0])
    for i,k in enumerate(['optimal','infeasible_fallback','numerical_failure']):v[i]+=int(r[k])
fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
for ax,kind in zip(axes,['paper','control']):
    ks=sorted(k for k in groups if k[0]==kind);bottom=np.zeros(len(ks))
    for i,label in enumerate(['Optimal','Certified infeasible fallback','Numerical failure']):
        counts=np.array([groups[k][i] for k in ks]);total=np.array([sum(groups[k]) for k in ks]);fraction=counts/total
        ax.bar([k[1] for k in ks],fraction,bottom=bottom,label=label);bottom+=fraction
    ax.set(xlabel='Setting',ylabel='Fraction of candidate-time fits',title=kind.capitalize(),ylim=(0,1),xticks=[k[1] for k in ks])
axes[1].legend(loc='lower center',fontsize=8)
fig.suptitle('All six scores retained: infeasibility and numerical failures remain distinct')
save(fig,'candidate-outcomes')
with (timing/'timings.csv').open() as f:tr=list(csv.DictReader(f))
fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
for ax,setting in zip(axes,['7','6']):
    rs=[r for r in tr if r['setting']==setting and r['sample']=='1'];names=[r['method'].replace('_',' ') for r in rs]
    ax.barh(names,[float(r['seconds'])/int(r['candidates']) for r in rs]);ax.set(xscale='log',xlabel='Seconds per candidate (warm repeat)',title='p = '+('50' if setting=='7' else '100'))
    for i,r in enumerate(rs):
        ax.text(float(r['seconds'])/int(r['candidates'])*1.1,i,f"{r['accepted']}/{r['candidates']} accepted",va='center',fontsize=8)
fig.suptitle('Two saved covariance inputs, one CPU; fixed spectral optimum is zero\nDifferent statistical calibrations; excludes covariance and selection scoring')
save(fig,'controlled-point-timing')
