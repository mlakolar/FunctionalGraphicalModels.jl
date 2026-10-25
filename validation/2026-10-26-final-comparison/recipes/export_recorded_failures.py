import csv,json,pathlib,sys
root=pathlib.Path('/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct')
report=root/'outputs/paper-v4/reports/report-recovery-13447246-0'
validation=json.loads((report/'validation.json').read_text());selected=[r for r in validation if r['fit_failures']]
writer=None
for check in selected:
 folder=root/f"outputs/paper-v4/setting-{check['setting']}-rep-{check['replicate']}"
 with (folder/'fits.csv').open(newline='') as f:
  for row in csv.DictReader(f):
   if row['status']!='failed':continue
   record={'setting':check['setting'],'replicate':check['replicate'],'source_variant':check['source_variant'],**row}
   if writer is None:writer=csv.DictWriter(sys.stdout,fieldnames=list(record));writer.writeheader()
   writer.writerow(record)
