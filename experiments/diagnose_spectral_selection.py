"""Read-only explanation of PD/BIC missingness from recorded column statuses.
No refits, eigenvalue computations, truth losses or calibration changes.
"""
import collections,csv,json,pathlib,sys
root,out=map(pathlib.Path,sys.argv[1:]);out.mkdir(exist_ok=False)
counts=[];times=[];replicates=[]
def rows(path):
 with path.open(newline='') as f:yield from csv.DictReader(f)
def save(name,data):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
for s in range(1,13):
 for rep in range(1,21):
  folder=root/f'outputs/spectral-v5/paper-{s}-{rep}'
  zeros=collections.Counter();total=collections.Counter()
  for row in rows(folder/'column-status.csv'):
   key=(int(row['time_index']),int(row['candidate']));total[key]+=1
   zeros[key]+=row['status']=='optimal_analytic_zero'
  fits=list(rows(folder/'fits.csv'));bytime=collections.defaultdict(list)
  for r in fits:
   key=(int(r['time_index']),int(r['candidate']));nzero=zeros[key]
   category='rejected_'+r['status'] if r['status']!='optimal' else ('optimal_all_analytic_zero' if nzero==total[key] else 'optimal_some_analytic_zero' if nzero else 'optimal_no_analytic_zero')
   counts.append(dict(setting=s,replicate=rep,time_index=key[0],candidate=key[1],category=category,analytic_zero_columns=nzero,columns=total[key],positive_definite=r['positive_definite']))
   bytime[key[0]].append((r,nzero))
  with (folder/'selected.csv').open(newline='') as f:
   selections={int(t):status for method,t,c,status in csv.reader(f) if method=='spectral_l2_aclime_tuned'}
  for t,rs in bytime.items():
   optimal=[(r,z) for r,z in rs if r['status']=='optimal']
   pd=sum(r['positive_definite']=='true' for r,z in optimal)
   analytic_blocked=sum(z>0 for r,z in optimal)
   category='selected' if selections[t]=='selected' else 'all_optimal_candidates_have_analytic_zero_column' if analytic_blocked==len(optimal) else 'nonPD_without_analytic_zero_column_exists' if pd==0 else 'PD_candidate_but_no_BIC_selection'
   times.append(dict(setting=s,replicate=rep,time_index=t,optimal_candidates=len(optimal),optimal_candidates_with_analytic_zero_column=analytic_blocked,optimal_PD_candidates=pd,selection_status=selections[t],category=category))
  rs=[r for r in times if r['setting']==s and r['replicate']==rep]
  replicates.append(dict(setting=s,replicate=rep,evaluation_times=len(rs),selected_times=sum(r['category']=='selected' for r in rs),analytic_zero_blocked_times=sum(r['category']=='all_optimal_candidates_have_analytic_zero_column' for r in rs),other_no_selection_times=sum(r['category'] not in ('selected','all_optimal_candidates_have_analytic_zero_column') for r in rs)))
summary=[]
for s in range(1,13):
 rs=[r for r in times if r['setting']==s];cs=[r for r in counts if r['setting']==s]
 summary.append(dict(setting=s,replicates=20,evaluation_times=len(rs),selected_times=sum(r['category']=='selected' for r in rs),analytic_zero_blocked_times=sum(r['category']=='all_optimal_candidates_have_analytic_zero_column' for r in rs),other_no_selection_times=sum(r['category'] not in ('selected','all_optimal_candidates_have_analytic_zero_column') for r in rs),complete_selection_replicates=sum(r['selected_times']==r['evaluation_times'] for r in replicates if r['setting']==s),reported_PD_with_analytic_zero_columns=sum(r['positive_definite']=='true' and r['analytic_zero_columns']>0 and r['category'].startswith('optimal_') for r in cs)))
save('summary.csv',summary);save('candidate-reasons.csv',counts);save('time-reasons.csv',times);save('replicate-reasons.csv',replicates)
(out/'receipt.json').write_text(json.dumps({'status':'completed','datasets':240,'candidate_times':len(counts),'evaluation_times':len(times),'scope':'Paper fixed first20 seeds per12 settings. Saved statuses only; no truth tuning or refits. An analytic zero precision column implies a zero diagonal retained by symmetrization and therefore prevents strict PD in exact arithmetic. Other nonPD candidates require a separate explanation.','reported_PD_zero_column_anomalies':sum(r['reported_PD_with_analytic_zero_columns'] for r in summary)},indent=2)+'\n')
