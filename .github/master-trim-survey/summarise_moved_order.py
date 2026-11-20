from pathlib import Path
import hashlib
import json
import sys


def differences(a, b, path='', tolerance=1e-6):
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(a.keys() | b.keys()):
            if key not in a or key not in b:
                yield {'path': f'{path}/{key}', 'present': [key in a, key in b]}
            else:
                yield from differences(a[key], b[key], f'{path}/{key}', tolerance)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield {'path': f'{path}/length', 'before': len(a), 'after': len(b)}
        for i, (x, y) in enumerate(zip(a, b)):
            yield from differences(x, y, f'{path}/{i}', tolerance)
    elif type(a) in (int, float) and type(b) in (int, float):
        if abs(a - b) > tolerance * max(abs(a), abs(b), 1):
            yield {'path': path, 'before': a, 'after': b}
    elif type(a) is not type(b) or a != b:
        yield {'path': path, 'before': a, 'after': b}


def summarize(cold, warm):
    result = {}
    for role, doc in [('recognition-first', cold), ('mass-first', warm)]:
        violations = list(differences(doc['moved'], doc['read_back']))
        result[role] = {'outside_recognition_tolerance': len(violations), 'differences': violations[:30]}
    assert len(warm['area_pairs']) == warm['faces'] and len(warm['volume_pairs']) == warm['solids']
    for label in ['area_pairs', 'volume_pairs']:
        vals = warm[label]
        result[label] = {'pairs': len(vals), 'unequal_availability': sum((a is None) != (b is None) for a,b in vals),
                         'outside_1e_9': sum(abs(a-b) > 1e-9 * max(1,abs(a),abs(b)) for a,b in vals if a is not None and b is not None)}
    for key in ['moved', 'read_back']:
        diffs = list(differences(cold[key], warm[key]))
        result[key + '_order_dependence'] = {'count': len(diffs), 'differences': diffs[:30]}
    return result


if __name__ == '__main__':
    cold, warm, output = map(Path, sys.argv[1:])
    report = summarize(json.loads(cold.read_bytes()), json.loads(warm.read_bytes()))
    for k, f in [('recognition-first',cold),('mass-first',warm)]:
        report[k]['raw_sha256'] = hashlib.sha256(f.read_bytes()).hexdigest()
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k: v.get('outside_recognition_tolerance',v.get('count',v)) for k,v in report.items()}))
