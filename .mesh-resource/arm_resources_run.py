"""在同一 ARM Mac 顺序构建和测量五个固定版本；硬件计数不可用时明确保留缺失。"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import time

STAGES=('base','pr1','pr2','pr3','pr4')
LABELS={'instructions_retired':'instructions retired','maximum_resident_set_size':'maximum resident set size',
        'peak_memory_footprint':'peak memory footprint','cycles_elapsed':'cycles elapsed'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')


def run(*args,cwd=None):
    return subprocess.check_output(args,cwd=cwd).decode().strip()


def parse_time(text):
    result={}
    for key,label in LABELS.items():
        found=re.findall(r'^\s*(\d+)\s+'+re.escape(label)+r'\s*$',text,re.MULTILINE)
        assert len(found)<=1,(key,found)
        result[key]=int(found[0]) if found else None
    found=re.findall(r'^\s*([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys\s*$',text,re.MULTILINE)
    assert len(found)==1,text
    result.update(zip(('real_seconds','user_seconds','system_seconds'),map(float,found[0])))
    assert result['maximum_resident_set_size'] and result['maximum_resident_set_size']>0
    return result


def orders():
    forward=[tuple(STAGES[(i+offset)%5] for i in (0,1,4,2,3)) for offset in range(5)]
    return forward+[tuple(reversed(order)) for order in forward]


def capability(out):
    assert platform.system()=='Darwin' and platform.machine()=='arm64'
    source=out/'counter-probe.c'
    source.write_text('#include <stdint.h>\n#include <stdio.h>\nint main(void) { volatile uint64_t x=7; for (uint64_t i=0;i<60000000;i++) { x=x*2862933555777941757ULL+3037000493ULL; } printf("%llu\\n",(unsigned long long)x); return 0; }\n',encoding='utf-8')
    subprocess.run(['clang','-O2',str(source),'-o',str(out/'counter-probe')],check=True)
    counters=out/'counter-probe.time.txt'
    process=subprocess.run(['/usr/bin/time','-l','-o',str(counters),str(out/'counter-probe')],capture_output=True,check=True)
    assert process.stdout.strip().isdigit()
    stats=parse_time(counters.read_text())
    info={'utc':datetime.now(timezone.utc).isoformat(),'platform':platform.platform(),'machine':platform.machine(),
        'sw_vers':run('sw_vers'),'cpu':run('sysctl','-n','machdep.cpu.brand_string'),
        'memory_bytes':run('sysctl','-n','hw.memsize'),'measurement_tool':'/usr/bin/time -l','probe':stats,
        'instructions_available':bool(stats['instructions_retired'] and stats['instructions_retired']>1000000),
        'runner_image':{key:os.environ.get(key) for key in ('ImageOS','ImageVersion','RUNNER_ARCH','RUNNER_OS')},
        'time_binary_sha256':digest(Path('/usr/bin/time'))}
    virtual=subprocess.run(['sysctl','-n','kern.hv_vmm_present'],capture_output=True,text=True)
    info['hypervisor_present']=virtual.stdout.strip() if virtual.returncode==0 else None
    write(out/'capability.json',info)
    print(json.dumps({'capability':info}),flush=True)


def build(source,corpus,out,plan,probe):
    assert run('git','rev-parse','HEAD',cwd=corpus)==plan['corpus_commit']
    assert not run('git','status','--porcelain','--untracked-files=no',cwd=source)
    assert run('rustc','-V').startswith('rustc 1.99.1 ')
    target=Path(os.environ['CARGO_TARGET_DIR'])
    binaries={}
    for stage in STAGES:
        identity=plan['sources'][stage]
        subprocess.run(['git','checkout','--detach',identity['source']],cwd=source,check=True)
        assert run('git','rev-parse','HEAD^{tree}',cwd=source)==identity['tree']
        # Ensure the package is reconsidered even with a restored incremental cache.
        (source/'crates/haecceity/src/lib.rs').touch()
        argv=['cargo','build','--release','--locked','-p','haecceity','--message-format=json']
        with (out/f'build-{stage}.jsonl').open('wb') as stdout,(out/f'build-{stage}.log').open('wb') as stderr:
            subprocess.run(argv,cwd=source,stdout=stdout,stderr=stderr,check=True)
        assert not run('git','status','--porcelain','--untracked-files=no',cwd=source)
        events=[json.loads(line) for line in (out/f'build-{stage}.jsonl').read_text().splitlines()]
        libs=[Path(name) for event in events if event.get('reason')=='compiler-artifact' and event['target']['name']=='haecceity'
              for name in event['filenames'] if name.endswith('.rlib')]
        assert len(libs)==1
        library=out/f'libhaecceity-{stage}.rlib'
        shutil.copy2(libs[0],library)
        binary=out/f'mesh-resource-{stage}'
        command=['rustc','--edition=2024','-C','opt-level=3','-L','dependency='+str(target/'release/deps'),
                 '--extern','haecceity='+str(library),str(probe),'-o',str(binary)]
        subprocess.run(command,check=True)
        binaries[stage]={'path':str(binary),'sha256':digest(binary),'library_sha256':digest(library),**identity,
                         'cargo_command':argv,'rustc_command':command}
        print(json.dumps({'built':stage,'source':identity['source']}),flush=True)
    write(out/'build-identities.json',{'binaries':binaries,'rustc':run('rustc','-Vv'),'cargo':run('cargo','-V'),
        'uv':run('uv','--version'),'python':platform.python_version(),'controller_sha256':digest(Path(__file__)),
        'probe_sha256':digest(probe),'plan_sha256':digest(out/'plan.json'),
        'environment':{key:os.environ.get(key) for key in ('RUSTFLAGS','RUSTUP_TOOLCHAIN','CARGO_BUILD_JOBS','CARGO_TARGET_DIR')}})
    return binaries


def background_work():
    names=run('ps','-axo','comm=').splitlines()
    return [name for name in names if Path(name.strip()).name in ('cargo','rustc','clang','clang++','ld','lld','swiftc')]


def measure(corpus,out,plan,binaries):
    assert not (out/'samples.jsonl').exists(), 'Do not overwrite or combine independent measurement runs'
    assert not background_work()
    expected={}
    count=0
    begin=time.monotonic()
    for round_id,order in enumerate(orders()):
        cases=plan['cases'][round_id%len(plan['cases']):]+plan['cases'][:round_id%len(plan['cases'])]
        for case,path in cases:
            for calls in ((1,10) if round_id%2==0 else (10,1)):
                for stage in order:
                    assert not background_work()
                    stem=f'{round_id:02}-{case}-{calls:02}-{stage}'
                    timer=out/(stem+'.time.txt')
                    stdout=out/(stem+'.stdout.json')
                    stderr=out/(stem+'.stderr.log')
                    command=['/usr/bin/time','-l','-o',str(timer),binaries[stage]['path'],str(corpus/path),str(calls)]
                    with stdout.open('wb') as output,stderr.open('wb') as errors:
                        result=subprocess.run(command,stdout=output,stderr=errors,timeout=900)
                    assert result.returncode==0,(stem,result.returncode)
                    values=json.loads(stdout.read_text())
                    assert values['calls']==calls and len(values['mesh_seconds'])==calls
                    shape=tuple(values[key] for key in ('faces','points','triangles','refusals'))
                    assert expected.setdefault((case,stage),shape)==shape,(case,stage,shape)
                    metrics=parse_time(timer.read_text())
                    entry={'round':round_id,'case':case,'stage':stage,'calls':calls,'command':command,**values,**metrics,
                           'time_sha256':digest(timer),'stdout_sha256':digest(stdout),'stderr_sha256':digest(stderr)}
                    with (out/'samples.jsonl').open('a',encoding='utf-8') as stream:
                        stream.write(json.dumps(entry)+'\n')
                    count+=1
                print(json.dumps({'round':round_id,'case':case,'calls':calls,'completed_processes':count,'elapsed_seconds':round(time.monotonic()-begin,1)}),flush=True)
    assert count==len(orders())*len(plan['cases'])*2*len(STAGES)
    rows=[json.loads(line) for line in (out/'samples.jsonl').read_text().splitlines()]
    available=all(row['instructions_retired'] and row['instructions_retired']>0 for row in rows)
    for stage in STAGES:
        assert digest(Path(binaries[stage]['path']))==binaries[stage]['sha256']
    write(out/'completion.json',{'all_processes_succeeded':True,'processes':count,'mesh_calls':sum(row['calls'] for row in rows),
        'instructions_available':available,'required_metrics_complete':bool(available),'peak_memory_measured':True,
        'samples_sha256':digest(out/'samples.jsonl'),'binary_hashes_unchanged':True})


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('phase',choices=('capability','measure'))
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--corpus',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--probe',type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if args.phase=='capability':
        capability(args.out)
        return
    assert (args.out/'capability.json').exists()
    plan=json.loads(args.plan.read_text())
    assert plan['controller_sha256']==digest(Path(__file__)) and plan['probe_sha256']==digest(args.probe)
    for case,path in plan['cases']:
        assert digest(args.corpus/path)==plan['compressed_input_sha256'][case]
    shutil.copy2(args.plan,args.out/'plan.json')
    binaries=build(args.source,args.corpus,args.out,plan,args.probe)
    measure(args.corpus,args.out,plan,binaries)


if __name__=='__main__':
    main()
