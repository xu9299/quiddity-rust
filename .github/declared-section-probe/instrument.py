from pathlib import Path
import hashlib, json, shutil, sys
root, evidence = map(Path, sys.argv[1:])
tools = Path(__file__).resolve().parent
records = {}
for stem, method, key, fields in [
    ('cylindrical_pockets','fn proofs_or_none(','floor','floor'),
    ('cylindrical_passages','pub fn prove(','planar','(planar,cylinder,walls)'),
]:
    file = root / 'src/features' / (stem + '.rs')
    original = file.read_bytes()
    content = original.decode().splitlines(keepends=True)
    first = next(i for i,line in enumerate(content) if line.startswith(method))
    last = next(i for i in range(first+1,len(content)) if content[i].startswith('}'))
    changes = {}
    for i in range(first,last):
        stripped = content[i].strip()
        if stripped not in ('return None;','continue;','let part = ctx.part;') and not stripped.startswith(('let Some(cell_volume)','let cell_volume = volume','results.push(CylindricalPocketProof','Some(CylindricalPassageProof')): continue
        label = ('reject' if stripped in ('return None;','continue;') else 'enter' if stripped.startswith('let part') else 'reached')
        trace = 'eprintln!("TRACE '+stem+' '+str(i+1)+' '+label+' {:?}", '+fields+');'
        content[i] = content[i][:-len(content[i].lstrip())] + trace + '\n' + content[i]
        changes[i+1]=stripped
    instrumented = ''.join(content).encode()
    file.write_bytes(instrumented)
    records[stem]={'original_sha256':hashlib.sha256(original).hexdigest(),'instrumented_sha256':hashlib.sha256(instrumented).hexdigest(),'lines':changes}
for rel, previous, traced in [
    ('src/features/entry_treatments.rs', '    Some(common_volume(ctx.solid_classifier(owner), &probe)? / whole)', '    let common = common_volume(ctx.solid_classifier(owner), &probe);\n    eprintln!("VOLUME owner={owner} common={common:?} whole={whole:?}");\n    Some(common? / whole)'),
    ('src/features/support_patches.rs','    covered(patch, supports)\n}', '    let answer = covered(patch, supports);\n    eprintln!("COVER face={} supports={:?} ok={answer}", patch.1, supports.iter().map(|s| s.1).collect::<Vec<_>>());\n    answer\n}'),
]:
    file=root/rel
    original=file.read_bytes()
    source=original.decode().replace('\r\n','\n')
    assert source.count(previous)==1,rel
    changed=source.replace(previous,traced).encode()
    file.write_bytes(changed)
    records[rel]={'original_sha256':hashlib.sha256(original).hexdigest(),'instrumented_sha256':hashlib.sha256(changed).hexdigest()}
dst = root/'examples/section_diagnostic.rs'
assert not dst.exists()
shutil.copyfile(tools / dst.name, dst)
(evidence/'instrumentation.json').write_text(json.dumps(records,indent=2)+'\n')
