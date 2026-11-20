from pathlib import Path
import collections, gzip, hashlib, json, os, re, subprocess, sys
root, corpus, out, exe = [Path(x).resolve() for x in sys.argv[1:]]
rows = {}
for ident, axes, source in [
    ('mfcadpp_10163', ['x'], corpus/'mfcadpp/10163.step'),
    ('section_0dc2847e',['z'],root/'tests/fixtures/captured/section_geometry/parts/0dc2847ef983c5bb.step.gz'),
    ('recess_722424d8',['z','x'],root/'tests/fixtures/captured/section_recess_geometry/parts/722424d84733ed7a.step.gz')
]:
    data=gzip.decompress(source.read_bytes()) if source.suffix=='.gz' else source.read_bytes()
    modified,n=re.subn(rb'\.PCURVE_S[12]\.',b'.CURVE_3D.',data)
    assert n and modified != data
    outputs={}
    for label,raw in [('declared',data),('spatial',modified)]:
        file=out/(ident+'-'+label+'.step')
        file.write_bytes(raw)
        for axis in axes:
            stem=ident+'-'+label+'-'+axis
            result=subprocess.run([exe,file,axis],capture_output=True,check=False)
            (out/(stem+'.log')).write_bytes(result.stderr)
            result.check_returncode()
            vals=json.loads(result.stdout)
            assert [q['placement'] for q in vals]==['identity','quarter_turn']
            outputs[label+'_'+axis]=vals
        file.unlink()
    rows[ident]={'source_sha256':hashlib.sha256(data).hexdigest(),'changed_selectors':n,'comparisons':outputs}
(out/'measurements.json').write_text(json.dumps(rows,indent=2)+'\n')
