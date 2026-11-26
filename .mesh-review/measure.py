"""Measure a fixed, calibrated, user-mode instruction event on every version in one job."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from ablation import BASE, HEAD, GROUPS, generate, git_bytes, identity

HERE=Path(__file__).resolve().parent
def save(path,value):
    path.write_text(json.dumps(value,indent=2)+"\n",encoding="utf8")
def digest(p):
    with p.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()
def main(source,corpus,evidence):
    evidence.mkdir(exist_ok=True)
    builddir=Path(os.environ["RUNNER_TEMP"])/"mesh-split-workers"
    builddir.mkdir(exist_ok=True)
    countobj=builddir/"counter.o"
    calibration=builddir/"calibrate"
    cfile=HERE/"counter.c"
    subprocess.run(["cc","-O2","-c",cfile,"-o",countobj],check=True)
    subprocess.run(["cc","-O2","-DSTUDY_CALIBRATION",cfile,"-pthread","-o",calibration],check=True)
    trials=[]
    prefix=None
    for cpu in sorted(os.sched_getaffinity(0)):
        cmd=["sudo","-n","--","taskset","-c",str(cpu)]
        vals=[]
        for n in (2000000,20000000):
            a=subprocess.run([*cmd,calibration,str(n)],capture_output=True)
            vals.append([n,a.returncode,a.stderr.decode().strip()])
            if a.returncode:break
        trials.append([cpu,vals])
        if len(vals)==2 and vals[-1][1]==0:
            lo,hi=[json.loads(row[2]) for row in vals]
            if 5.5<(hi["instructions"]-lo["instructions"])/18000000<8.5 and all(r["running"]>=0.999*r["enabled"]>0 for r in (lo,hi)):
                prefix=cmd
                break
    save(evidence/"calibration.json",{"trials":trials,"command":prefix,"counter_sha256":digest(cfile)})
    cases={"cgb243":"cadgenbench_inputs/cgb243.step.gz","cgb242":"cadgenbench_inputs/cgb242.step.gz","string_post":"gramel/string_post.step"}
    stages=["pr349","stripped",*GROUPS]
    rounds=len(stages)
    event_kind="hardware user-mode retired instructions"
    if prefix is None:
        subprocess.run(["sudo","-n","apt-get","update","-qq"],check=True)
        subprocess.run(["sudo","-n","apt-get","install","-y","-qq","valgrind"],check=True)
        subprocess.run(["cc","-O2","-c",HERE/"callgrind_counter.c","-o",countobj],check=True)
        cases={"cgb243":cases["cgb243"]}
        rounds=2
        event_kind="Callgrind architectural guest instructions (Ir) on the worker"
        save(evidence/"callgrind.json",{"version":subprocess.check_output(["valgrind","--version"]).decode().strip(),"counter_sha256":digest(HERE/"callgrind_counter.c")})
    save(evidence/"plan.json",{"stages":stages,"groups":GROUPS,"rounds":rounds,"event":event_kind,"inputs":{k:digest(corpus/v) for k,v in cases.items()},
        "repo":HEAD,"parent":BASE,"corpus":git_bytes(corpus,"rev-parse","HEAD").decode().strip(),"rustc":subprocess.check_output(["rustc","-Vv"]).decode(),
        "kernel":subprocess.check_output(["uname","-a"]).decode(),"processor":subprocess.check_output(["lscpu"]).decode(),
        "worker_sha256":digest(HERE/"worker.rs"),"controller_sha256":digest(Path(__file__))})
    binaries={}
    for pos,stage in enumerate(stages):
        if pos in (0,1):
            subprocess.run(["git","checkout","--detach",BASE if pos==0 else HEAD],cwd=source,check=True)
        files=generate(source,GROUPS[:max(pos-1,0)]) if pos else {}
        for name,content in files.items():
            p=source/name
            if p.read_bytes()!=content:p.write_bytes(content)
        if pos==len(stages)-1 or pos==0:
            assert not git_bytes(source,"status","--porcelain","--untracked-files=no")
        with (evidence/(stage+"-build.log")).open("wb") as log:
            output=subprocess.check_output(["cargo","build","--release","--locked","-p","haecceity","--message-format=json"],cwd=source,stderr=log)
        ev=[json.loads(line) for line in output.splitlines()]
        event=next(r for r in ev if r.get("reason")=="compiler-artifact" and r["target"]["name"]=="haecceity")
        lib=next(Path(p) for p in event["filenames"] if p.endswith(".rlib"))
        binary=builddir/stage
        deps=Path(os.environ["CARGO_TARGET_DIR"])/"release/deps"
        subprocess.run(["rustc","--edition=2024","-C","opt-level=3","-C","link-arg="+str(countobj),
                "-L","dependency="+str(deps),"--extern","haecceity="+str(lib),HERE/"worker.rs","-o",binary],check=True)
        binaries[stage]=binary
        save(evidence/(stage+"-identity.json"),{"revision":BASE if pos==0 else HEAD,"enabled":[] if pos==0 else GROUPS[:pos-1],
                "files":identity(files),"library_sha256":digest(lib),"binary_sha256":digest(binary),"features":event["features"],"profile":event["profile"]})
        print("Built",stage,flush=True)
    # Every stage occupies every position once for each input in the rotated schedule.
    with (evidence/"samples.jsonl").open("w",encoding="utf8") as stream:
        for rnd in range(rounds):
            names=list(cases)
            for ic,name in enumerate(names[rnd%len(names):]+names[:rnd%len(names)]):
                rotation=(rnd+2*names.index(name))%len(stages)
                order=stages[rotation:]+stages[:rotation] if prefix else (stages if rnd==0 else list(reversed(stages)))
                for position,stage in enumerate(order):
                    args=[binaries[stage],corpus/cases[name],"1"]
                    if prefix:
                        result=subprocess.run([*prefix,*args],capture_output=True,check=True)
                        event=json.loads(result.stderr)
                        assert event["instructions"]>0 and event["running"]>=event["enabled"]*0.999>0
                    else:
                        pattern=f"cg-{rnd}-{stage}"
                        base=builddir/(pattern+".%p")
                        result=subprocess.run(["valgrind","--tool=callgrind","--cache-sim=no","--branch-sim=no","--instr-atstart=no","--separate-threads=yes","--callgrind-out-file="+str(base),*args],capture_output=True,check=True)
                        records=[]
                        for output in builddir.glob(pattern+".*"):
                            with output.open(errors="replace") as cg:
                                header=[line.strip() for _,line in zip(range(75),cg)]
                            if any("Trigger: Client Request: mesh-workload" in line for line in header):
                                assert "events: Ir" in header
                                records.append(int(next(line.split()[1] for line in header if line.startswith("summary:"))))
                        assert len(records)==1 and records[0]>0,(header,records)
                        event={"instructions":records[0],"event":"Callgrind Ir"}
                    row=dict(round=rnd,case=name,stage=stage,position=position,**json.loads(result.stdout),**event)
                    stream.write(json.dumps(row)+"\n");stream.flush()
            print("Finished round",rnd,flush=True)
    save(evidence/"completed.json",{"rounds":rounds,"measurements":rounds*len(stages)*len(cases),"event":event_kind,"samples_sha256":digest(evidence/"samples.jsonl")})
if __name__=="__main__":main(*(Path(s).resolve() for s in sys.argv[1:]))
