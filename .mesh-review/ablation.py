"""Create successive measured configurations from the exact PR head without changing the public APIs."""
from pathlib import Path
import hashlib
import subprocess

BASE = "4a0ae2f93dc6fae6a30b834fdb290c6781c8c70f"
HEAD = "2e5878594a455c9564af19d3ff5024dc45f91c66"
GROUPS = ("operands","basis","surface","qualification","chords","crossings","uv")
SOURCES = ("mass/declared.rs","nurbs.rs","mass.rs","mass/declared/intersections.rs","uv.rs")
PREFIX = "crates/haecceity/src/"

def git_bytes(repo, *args):
    return subprocess.check_output(["git","-C",str(repo),*args])

def generate(repo, enabled):
    files={n:git_bytes(repo,"show",HEAD+":"+PREFIX+n).decode() for n in SOURCES}
    def replace(n, old, new, times=1):
        if files[n].count(old)!=times:
            raise ValueError((n,old,files[n].count(old),times))
        files[n]=files[n].replace(old,new)
    decl="mass/declared.rs"
    if "operands" not in enabled:
        # Build each occurrence in its own units; take a fresh carrier instead of its shared Arc.
        replace(decl,"let source = if let Some(edge) = edge {","let source = if let Some(edge) = edge.filter(|_| false) {")
        replace(decl,"self.cache[f]\n                .carrier\n                .get_or_init(|| Arc::new(self.faces[f].surface.clone()))\n                .clone()","Arc::new(self.faces[f].surface.clone())")
    if "basis" not in enabled:
        replace("nurbs.rs","if t.is_finite() {\n        match degree {","if false && t.is_finite() {\n        match degree {",2)
    if "surface" not in enabled:
        replace("mass.rs","if reuse => Some(Box::new(std::cell::RefCell::new(","if false && reuse => Some(Box::new(std::cell::RefCell::new(",2)
    if "qualification" not in enabled:
        replace(decl,"if mesh\n        && !trim.is_transferred()\n        && EdgeTrim::midpoint_bound", "if false && mesh\n        && !trim.is_transferred()\n        && EdgeTrim::midpoint_bound")
        replace(decl,"if mesh\n        && let Some(size) = EdgeTrim::affine_spatial_size", "if false && mesh\n        && let Some(size) = EdgeTrim::affine_spatial_size")
    if "chords" not in enabled:
        replace(decl,"if self.mesh\n                    && !self.uniform_mesh_stations","if false && self.mesh\n                    && !self.uniform_mesh_stations")
        replace(decl,"let bound = self\n                    .mesh\n                    .then(", "let bound = (false && self.mesh)\n                    .then(")
    if "crossings" not in enabled:
        replace("mass/declared/intersections.rs","            sweep_axes,\n","            sweep_axes: false && sweep_axes,\n")
    if "uv" not in enabled:
        name="uv.rs"
        start=files[name].index("        let key: Vec<_> = self",files[name].index("fn mesh_uv_loops"))
        end=files[name].index("        self.compute_uv_loops(face, None, Some(selected))\n            .map(std::borrow::Cow::Owned)",start)
        files[name]=files[name][:start]+files[name][end:]
        replace(name,".filter(|edge| edge.face == face)\n                        .and(own.native_uv.as_deref())",".filter(|edge| false && edge.face == face)\n                        .and(own.native_uv.as_deref())")
        replace(decl,"if mesh && let Some(Some(loops)) = self.cache[face].uv_loops.get()", "if false && mesh && let Some(Some(loops)) = self.cache[face].uv_loops.get()")
    return {PREFIX+n:s.encode() for n,s in files.items()}

def identity(files):
    return {n:hashlib.sha256(b).hexdigest() for n,b in files.items()}

if __name__=="__main__":
    import json,sys
    repo=Path(sys.argv[1])
    for i in range(len(GROUPS)+1):
        enabled=GROUPS[:i]
        files=generate(repo,enabled)
        if i==len(GROUPS):
            assert all(b==git_bytes(repo,"show",HEAD+":"+n) for n,b in files.items())
        print(i,enabled,hashlib.sha256(json.dumps(identity(files),sort_keys=True).encode()).hexdigest())
