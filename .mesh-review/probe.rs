use haecceity::{geom::{self, V3}, mesh::{BoundaryPoint, FaceMesh, MeshRefusal}, read_step_file, Part};
use serde_json::{json, Value};
use std::{collections::{BTreeMap, BTreeSet}, fs, hash::{Hash, Hasher}, path::Path, sync::Barrier};

type Mesh = Vec<Result<FaceMesh, MeshRefusal>>;
type Bits = [u64; 3];

fn hash<T: Hash>(value: T) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}

fn signature(m: &FaceMesh) -> (u64, u64) {
    let mut rim = m.boundary.iter().zip(&m.points)
        .filter_map(|(k,p)| k.map(|k| (k,p.map(f64::to_bits)))).collect::<Vec<_>>();
    rim.sort();
    let signature = hash((m.points.iter().map(|p| p.map(f64::to_bits)).collect::<Vec<_>>(),
        m.uv.iter().map(|(x,y)| (x.to_bits(),y.to_bits())).collect::<Vec<_>>(), &m.triangles,
        &m.boundary, m.folded_dropped.0, m.folded_dropped.1.to_bits(), m.converged));
    (signature,hash(rim))
}

fn counts(part: &Part, mesh: &Mesh) -> Value {
    let rows: Vec<_> = mesh.iter().enumerate().map(|(i,r)| {
        let id=part.face_source(i).map(|x|x.entity);
        let kind=format!("{:?}",part.faces[i].surface.kind());
        match r {
            Err(e)=>json!({"face":i,"step":id,"kind":kind,"refusal":e.reason}),
            Ok(m)=>{let (s,b)=signature(m);json!({"face":i,"step":id,"kind":kind,"triangles":m.triangles.len(),"points":m.points.len(),"boundary_n":m.boundary.iter().flatten().count(),"boundary_hash":b,"mesh_hash":s,"area":m.area(),"folded":m.folded_dropped.0,"converged":m.converged})}
        }
    }).collect();
    json!(rows)
}

fn segment_distance(p: V3, a: V3, b: V3) -> f64 {
    let ab=geom::sub(b,a);
    let denom=geom::dot(ab,ab);
    let t=if denom>0. { (geom::dot(geom::sub(p,a),ab)/denom).clamp(0.,1.) } else {0.};
    geom::dist(p,geom::add(a,geom::scale(ab,t)))
}

fn triangle_distance(p: V3, a: V3, b: V3, c: V3) -> f64 {
    let ab=geom::sub(b,a); let ac=geom::sub(c,a); let ap=geom::sub(p,a);
    let n=geom::cross(ab,ac); let nn=geom::dot(n,n);
    if nn>0. {
        let u=geom::dot(geom::cross(ap,ac),n)/nn;
        let v=geom::dot(geom::cross(ab,ap),n)/nn;
        if u>=0. && v>=0. && u+v<=1. {return geom::dot(ap,n).abs()/nn.sqrt();}
    }
    segment_distance(p,a,b).min(segment_distance(p,b,c)).min(segment_distance(p,c,a))
}

fn quality(part: &Part,i:usize,m:&FaceMesh,d:f64,probes:&[[f64;2]])->Value {
    let surf=&part.faces[i].surface;
    let mut max_deviation=0.0_f64;
    let mut inward=0; let mut unresolved=0; let mut unresolved_normals=0; let mut samples=0; let mut slivers=0;
    for tri in &m.triangles {
        let p=tri.map(|j|m.points[j]); let uv=tri.map(|j|m.uv[j]);
        let cross=geom::cross(geom::sub(p[1],p[0]),geom::sub(p[2],p[0]));
        let longest=(0..3).map(|j|geom::dist(p[j],p[(j+1)%3])).fold(0.,f64::max);
        for (k,w) in [[1./3.;3],[0.5,0.5,0.],[0.,0.5,0.5],[0.5,0.,0.5]].iter().enumerate() {
            let q=std::array::from_fn(|a| (0..3).map(|j|p[j][a]*w[j]).sum());
            let seed=((0..3).map(|j|uv[j].0*w[j]).sum(),(0..3).map(|j|uv[j].1*w[j]).sum());
            if let Some((u,v))=surf.parameters(q,Some(seed)) {
                samples+=1;
                max_deviation=max_deviation.max(geom::dist(q,surf.value(u,v)));
                if k==0 {
                    if geom::norm(cross)<=1e-6*d*longest {slivers+=1;}
                    else if let Some(n)=part.face_normal(i,u,v){if geom::dot(n,cross)<=0. {inward+=1;}}
                    else {unresolved_normals+=1;}
                }
            }else{unresolved+=1; if k==0 {unresolved_normals+=1;}}
        }
    }
    let witness = probes.iter().map(|&[u,v]| {
        let p=surf.value(u,v);
        let nearest=m.triangles.iter().map(|t|triangle_distance(p,m.points[t[0]],m.points[t[1]],m.points[t[2]])).fold(f64::INFINITY,f64::min);
        json!({"uv":[u,v],"point":p,"distance_to_mesh":nearest})
    }).collect::<Vec<_>>();
    let area=part.face_moments(i).map(|x|x.area);
    json!({"face":i,"step":part.face_source(i).map(|x|x.entity),"surface":format!("{:?}",surf.kind()),"area":m.area(),"reference":area,"area_error_pct":area.filter(|a|*a>0.).map(|a|100.*(m.area()/a-1.)),"triangles":m.triangles.len(),"converged":m.converged,"folds":m.folded_dropped.0,"inward":inward,"slivers":slivers,"unresolved_normals":unresolved_normals,"max_mesh_to_surface":max_deviation,"samples":samples,"projection_failures":unresolved,"witnesses":witness,"mesh_hash":signature(m).0})
}

fn shared(part:&Part, mesh:&Mesh, focus:&BTreeSet<usize>)->Value {
    let mut edges=BTreeMap::<usize,Vec<usize>>::new();
    for f in 0..part.faces.len(){ for e in part.face_edges(f) {edges.entry(e).or_default().push(f);} }
    let mut checked=0; let mut skipped=0; let mut unobserved=0; let mut mismatches=Vec::new(); let mut samples=0;
    for (e,mut fs) in edges {fs.sort();fs.dedup();
        if fs.len()<2 || !fs.iter().any(|f|focus.contains(f)){continue;}
        if fs.iter().any(|&f|mesh[f].is_err()){skipped+=1;continue;}
        let edge=&part.edges[e];
        let boundary=|f:usize|{
            let m=mesh[f].as_ref().unwrap();
            m.boundary.iter().zip(&m.points).filter_map(|(k,p)|match *k {
                Some(BoundaryPoint::Edge{edge,index}) if edge==e =>Some(((1,index),p.map(f64::to_bits))),
                Some(BoundaryPoint::Vertex(v)) if *p==edge.start || *p==edge.end =>Some(((0,v),p.map(f64::to_bits))),
                _=>None
            }).collect::<BTreeMap<_,Bits>>()
        };
        let first=boundary(fs[0]); let all:Vec<_>=fs.iter().skip(1).map(|&f|boundary(f)).collect();
        if first.is_empty() && all.iter().all(|m|m.is_empty()){unobserved+=1;continue;}
        checked+=1; samples+=first.len();
        if all.iter().any(|b|*b!=first){mismatches.push(json!({"edge":e,"faces":fs,"points_first":first.len(),"points_others":all.iter().map(|s|s.len()).collect::<Vec<_>>() }));}
    }
    json!({"checked":checked,"first_side_points":samples,"refused_pairs":skipped,"unobserved":unobserved,"mismatches":mismatches})
}

fn geometry(path:&Path,d:f64,focus:&BTreeSet<usize>,probes:&[[f64;2]])->Value {
    let original=read_step_file(path).unwrap();
    let mesh=original.mesh(d,0.3);
    let mut alone1=BTreeMap::new();
    for &i in focus {alone1.insert(i,original.triangulate(i,d,0.3));}
    let other=read_step_file(path).unwrap();
    let mut alone2=BTreeMap::new();
    for &i in focus.iter().rev(){alone2.insert(i,other.triangulate(i,d,0.3));}
    let reversed=other.mesh(d,0.3);
    let rows=focus.iter().map(|&i|match &mesh[i]{
        Ok(m)=>quality(&original,i,m,d,probes),
        Err(e)=>json!({"face":i,"refusal":e.reason})
    }).collect::<Vec<_>>();
    json!({"d":d,"results":rows,"shared":shared(&original,&mesh,focus),"whole_order_equal":mesh==reversed,"single_order_equal":alone1==alone2,
        "whole_repeat_equal":mesh==original.mesh(d,0.3)})
}

fn orders(path:&Path,d:f64,trials:usize)->Value {
    let reference=read_step_file(path).unwrap(); let cold=reference.mesh(d,0.3);
    let warm=read_step_file(path).unwrap();
    for f in (0..warm.faces.len()).rev(){std::hint::black_box(warm.uv_loops(f));}
    let warm_mesh=warm.mesh(d,0.3);
    let mut different=0;
    for _ in 0..trials{
        let part=read_step_file(path).unwrap(); let barrier=Barrier::new(2);
        let result=std::thread::scope(|scope|{
            scope.spawn(||{barrier.wait(); for f in 0..part.faces.len(){std::hint::black_box(part.uv_loops(f));}});
            barrier.wait();part.mesh(d,0.3)
        });
        if result!=cold{different+=1;}
    }
    json!({"cold_warm_equal":cold==warm_mesh,"races":trials,"race_different_from_cold":different,"cold":counts(&reference,&cold),"warm":counts(&warm,&warm_mesh)})
}

fn main(){
    let args:Vec<_>=std::env::args().collect();
    let path=Path::new(&args[2]); let part=read_step_file(path).unwrap();
    let d=if args[3]=="default"{part.bounds().diagonal()*0.001}else{args[3].parse().unwrap()};
    let row=match args[1].as_str(){
        "counts"=>counts(&part,&part.mesh(d,0.3)),
        "order"=>orders(path,d,args[5].parse().unwrap()),
        "case"=>{
            let focus=if args[5]=="all"{(0..part.faces.len()).collect()}else{args[5].split(',').map(|n|n.parse().unwrap()).collect()};
            let probes:Vec<[f64;2]>=args.get(6).map(|s|serde_json::from_str(s).unwrap()).unwrap_or_default();
            geometry(path,d,&focus,&probes)
        }
        _=>panic!("unknown mode")
    };
    fs::write(&args[4],serde_json::to_vec(&row).unwrap()).unwrap();
}
