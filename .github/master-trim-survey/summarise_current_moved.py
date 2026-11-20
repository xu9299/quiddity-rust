from pathlib import Path
import hashlib, json, os, sys
from summarise_moved_order import differences
p=Path(sys.argv[1])
report={'source':os.environ['SOURCE'],'tree':os.environ['TREE'],'orders':{}}
records={}
for key in ['mass-first','recognition-first']:
 d=p/key
 docs=[json.loads((d/name).read_bytes()) for name in ['moved.json','read-back.json']]
 delta=list(differences(*docs)); issues=json.loads((d/'original-test-problems.json').read_bytes())
 assert bool(delta)==any('recognises otherwise' in v for v in issues)
 report['orders'][key]={'differences':len(delta),'first_differences':delta[:5],'original_test_problems':issues}
 records[key]=docs
report['query_order_differences']=[len(list(differences(records['mass-first'][i],records['recognition-first'][i]))) for i in (0,1)]
(p/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(report)
