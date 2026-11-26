"""Measure one public STEP read and one whole-part mesh per fresh process on a Linux runner."""
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE, CORPUS, OUT = (Path(a).resolve() for a in sys.argv[1:4])
MODE = sys.argv[4]
PLAN = json.loads((HERE/"plan.json").read_bytes())

def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def write(name, data):
    (OUT/name).write_text(json.dumps(data, indent=2)+"\n", encoding="utf8")

def capture(*args, cwd=None):
    return subprocess.check_output(args,cwd=cwd,text=True).strip()

def calibration():
    calibrate=OUT/"calibrate"
    subprocess.run(["cc","-std=gnu17","-O2","-Wall","-Wextra","-Werror","-DSTUDY_CALIBRATION",HERE/"counter.c","-pthread","-o",calibrate],check=True)
    attempts=[]
    for cpu in sorted(os.sched_getaffinity(0)):
        prefix=["sudo","-n","--","taskset","-c",str(cpu)]
        values=[]
        for n in (2000000,20000000):
            r=subprocess.run([*prefix,calibrate,str(n)],capture_output=True,text=True)
            values.append({"iterations":n,"exit":r.returncode,"stderr":r.stderr.strip()})
            if r.returncode:
                break
        attempts.append({"cpu":cpu,"values":values})
        if len(values)==2 and all(v["exit"]==0 for v in values):
            a,b=[json.loads(v["stderr"]) for v in values]
            slope=(b["instructions"]-a["instructions"])/18000000
            # At -O2 GCC uses seven retired instructions per iteration on x86-64;
            # on AArch64 the multiply-add loop can use one fewer instruction.
            if 4.5<slope<8.5 and all(0.999<=v["running"]/v["enabled"]<=1.0 for v in (a,b)):
                write("calibration.json",{"hardware":True,"prefix":prefix,"slope":slope,"trials":attempts})
                return prefix
    write("calibration.json",{"hardware":False,"trials":attempts})
    return None

def build(instrumented):
    probe=HERE/"probe.rs"
    opts=[]
    if instrumented:
        text=probe.read_text("utf8")
        begin="    let start = Instant::now();"
        end="    drop(part);"
        assert end in text and begin in text
        text=text.replace(begin,"    let fd = unsafe { mesh_pmu_begin() };\n    assert!(fd >= 0);\n"+begin,1)
        text=text.replace(end,end+"\n    assert_eq!(unsafe { mesh_pmu_end(fd) }, 0);",1)
        text+='\nunsafe extern "C" { fn mesh_pmu_begin() -> i32; fn mesh_pmu_end(fd: i32) -> i32; }\n'
        probe=OUT/"counter_probe.rs"
        probe.write_text(text,"utf8")
        obj=OUT/"counter.o"
        subprocess.run(["cc","-std=gnu17","-O2","-c",HERE/"counter.c","-o",obj],check=True)
        opts=["-C","link-arg="+str(obj)]
    builds={}
    for stage,src in PLAN["sources"].items():
        subprocess.run(["git","checkout","--quiet","--detach",src["source"]],cwd=SOURCE,check=True)
        assert capture("git","rev-parse","HEAD^{tree}",cwd=SOURCE)==src["tree"]
        wrapper=OUT/(stage+"-workspace-rustc")
        wrapper.write_text('#!/bin/sh\nexec "$@"\n',"utf8")
        wrapper.chmod(0o755)
        manifest=OUT/("cargo-"+stage+".jsonl")
        with manifest.open("wb") as log:
            subprocess.run(["cargo","build","-p","haecceity","--lib","--release","--locked","--message-format=json"],cwd=SOURCE,check=True,env=dict(os.environ,RUSTC_WORKSPACE_WRAPPER=str(wrapper)),stdout=log)
        entries=[json.loads(l) for l in manifest.read_bytes().splitlines()]
        assert entries[-1]=={"reason":"build-finished","success":True}
        crate=next(r for r in entries if r.get("target",{}).get("name")=="haecceity")
        lib,=[Path(p) for p in crate["filenames"] if p.endswith(".rlib")]
        exe=OUT/("mesh-"+stage)
        deps=Path(os.environ["CARGO_TARGET_DIR"])/"release/deps"
        subprocess.run(["rustc","--crate-name","mesh_resource","--edition=2024","-C","opt-level=3",*opts,"-L","dependency="+str(deps),"--extern","haecceity="+str(lib),probe,"-o",exe],check=True)
        builds[stage]=dict(**src,rlib=digest(lib),binary=digest(exe),features=crate["features"],profile=crate["profile"],lock=digest(SOURCE/"Cargo.lock"))
    assert len({v["rlib"] for v in builds.values()})==len(builds)
    write("builds.json",builds)
    return builds

def main():
    OUT.mkdir()
    assert capture("git","rev-parse","HEAD",cwd=CORPUS)==PLAN["corpus"]
    assert capture("rustc","-V").startswith("rustc 1.99.1 ")
    for record in PLAN["inputs"].values():
        assert digest(CORPUS/record["path"])==record["sha256"]
    write("plan.json",PLAN)
    write("environment.json",dict(platform=platform.platform(),cpu=capture("lscpu"),rustc=capture("rustc","-vV"),paranoid=Path("/proc/sys/kernel/perf_event_paranoid").read_text().strip()))
    prefix=calibration()
    hardware=prefix is not None
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"],"a") as output:
            print("hardware="+str(hardware).lower(),file=output)
    if MODE=="hardware-only" and not hardware:
        print("This runner does not expose a calibrated hardware event.")
        return
    if not prefix:
        prefix=["taskset","-c",str(min(os.sched_getaffinity(0)))]
    builds=build(hardware)
    stages=list(builds)
    orders=[[stages[(i+k)%4] for k in (0,1,3,2)] for i in range(4)]
    orders+= [list(reversed(o)) for o in orders]
    inputs=list(PLAN["inputs"].items())
    rows=[]
    shapes={}
    with (OUT/"samples.jsonl").open("w",encoding="utf8") as stream:
        for round_id,order in enumerate(orders):
            for name,item in inputs[round_id%len(inputs):]+inputs[:round_id%len(inputs)]:
                for stage in order:
                    timer=OUT/"rss.txt"
                    cmd=[*prefix,"/usr/bin/time","-f","%M","-o",str(timer),str(OUT/("mesh-"+stage)),str(CORPUS/item["path"]),"1"]
                    result=subprocess.run(cmd,capture_output=True,text=True,check=True)
                    measured=json.loads(result.stdout)
                    rss=int(timer.read_text().strip())*1024
                    assert rss>0 and measured["calls"]==1
                    event=json.loads(result.stderr) if hardware else {}
                    if hardware:
                        assert event["instructions"]>0 and 0.999<=event["running"]/event["enabled"]<=1.0
                    shape=[measured[k] for k in ("faces","points","triangles","refusals")]
                    assert shapes.setdefault((stage,name),shape)==shape
                    row={"case":name,"round":round_id,"stage":stage,"rss_bytes":rss,**event,**measured}
                    print(json.dumps(row),file=stream,flush=True)
                    rows.append(row)
            print("Completed paired round",round_id+1,flush=True)
    metrics=["rss_bytes"]+(["instructions"] if hardware else [])
    summary={}
    for name,_ in inputs:
        summary[name]={}
        for index,stage in enumerate(stages):
            def samples(s,key):
                return [row[key] for row in rows if row["case"]==name and row["stage"]==s]
            fields={}
            for key in metrics:
                values=samples(stage,key)
                ratios=lambda base:statistics.median(x/y for x,y in zip(values,samples(base,key),strict=True))
                fields[key]={"median":statistics.median(values),"vs_baseline":ratios(stages[0]),"vs_previous":ratios(stages[max(0,index-1)])}
            summary[name][stage]={"shape":shapes[stage,name],**fields}
    write("summary.json",{"hardware":hardware,"pairs":len(orders),"processes":len(rows),"samples_sha256":digest(OUT/"samples.jsonl"),"results":summary})

if __name__=="__main__":
    main()

