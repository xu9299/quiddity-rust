from pathlib import Path
import hashlib
import json
import sys
from summarise_moved_order import differences


out = Path(sys.argv[1])
result = {}
read = lambda p: json.loads(p.read_bytes())
hashfile = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for role in ['baseline','candidate']:
    data = {}
    docs = {}
    for order in ['recognition-first','mass-first']:
        folder = out / role / order
        docs[order] = [read(folder / name) for name in ['moved.json', 'read-back.json']]
        diffs = list(differences(*docs[order]))
        issues = read(folder / 'original-test-problems.json')
        assert bool(diffs) == any('recognises otherwise' in text for text in issues), (role, order)
        data[order] = {'original_test': issues, 'outside_recognition_tolerance': len(diffs),
                       'differences': diffs[:64], 'document_sha256': [hashfile(folder/name) for name in ['moved.json','read-back.json']]}
    for i, name in enumerate(['moved','read_back']):
        changes = list(differences(docs['recognition-first'][i], docs['mass-first'][i]))
        data[name + '_order_dependence'] = {'count':len(changes),'differences':changes[:64]}
    data['commit'],data['tree'] = (out / role / 'identity.txt').read_text().splitlines()
    data['fixture_sha256'] = (out / role / 'fixture.sha256').read_text().split()[0]
    result[role] = data
assert result['baseline']['fixture_sha256'] == result['candidate']['fixture_sha256']
(out / 'roundtrip-summary.json').write_text(json.dumps(result, indent=2)+'\n',encoding='utf8')
print(json.dumps({k:{o:v[o]['outside_recognition_tolerance'] for o in ['mass-first','recognition-first']} for k,v in result.items()}))
