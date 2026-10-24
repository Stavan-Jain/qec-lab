import pathlib,json,re,os
stage=pathlib.Path(os.environ['QEC360_STAGE']).resolve()
manifest=pathlib.Path(os.environ.get('QEC360_RESEARCH_MANIFEST', str(pathlib.Path(__file__).resolve().parent.parent/'research_manifest.json')))
j=json.loads(manifest.read_text());names=[]
for step in j['steps']:
 if step.get('intermediate_only') or 'mark_research_noncomputable' in step['script']:continue
 names += step['outputs']
if 'pending_batched_cap2_step'in j:names+=j['pending_batched_cap2_step']['outputs']
# Computable packed operations are useful generic checkers; all their data is in OrbitBatchData.
names=[x for x in dict.fromkeys(names)if x!='OrbitBatchCheck.lean']
for name in names:
 p=stage/name
 if p.exists():
  s=p.read_text();p.write_text(re.sub(r'(?m)^def ', 'noncomputable def ',s))
print('Marked explicit noncomputable data in',len(names),'generated source paths')
