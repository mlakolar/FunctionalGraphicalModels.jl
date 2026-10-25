"""Plot existing paired selected-fit contrasts; no fitting or resimulation."""
import csv,pathlib,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root,out=map(pathlib.Path,sys.argv[1:]);out.mkdir(exist_ok=False)
with (root/'paired-comparisons.csv').open() as f:rows=list(csv.DictReader(f))
fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
for ax,quantity,label in zip(axes,['operator','normalized_squared_frobenius'],['Operator norm error','Squared Frobenius error / p']):
 for method,offset in [('DFGM_C',-.10),('DFGM_G',.10)]:
  rs=sorted([r for r in rows if r['first']==method and r['second']=='scaled_lasso_tuned' and r['quantity']==quantity and r['target']=='actual'],key=lambda r:int(r['setting']))
  assert len(rs)==12
  ax.errorbar([int(r['setting'])+offset for r in rs],[-float(r['mean_difference']) for r in rs],yerr=[2*float(r['mcse']) for r in rs],fmt='o',capsize=3,label='minus '+method)
 ax.axhline(0,color='black',linewidth=.8);ax.set_yscale('symlog',linthresh=.01)
 ax.set(xticks=range(1,13),xticklabels=[f"{r['setting']}\n{r['paired_replicates']}" for r in rs],xlabel='Simulation setting; paired n on second line',ylabel='SZ tuned error minus original-method error',title=label);ax.grid(alpha=.2)
axes[0].legend();fig.suptitle('Actual precision target: paired means and 2 MCSE\nBoth methods use recorded BIC selections; time averaged within independent seeds')
for ext in ('png','pdf','svg'):fig.savefig(out/f'sz-tuned-vs-dfgm-paired.{ext}',dpi=170,bbox_inches='tight')
plt.close(fig)
