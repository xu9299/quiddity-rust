//! Public caller query capture; no private inputs or geometry are embedded.
use super::*;
use crate::brep::RawCurve2d;
use crate::geom::Curve;
use serde_json::{Value, json};

fn kind(c: &Curve) -> &'static str {
    match c {
        Curve::Line { .. } => "line",
        Curve::Circle { .. } => "circle",
        Curve::Ellipse { .. } => "ellipse",
        Curve::Nurbs(_) => "spline",
    }
}


fn raw(c: &RawCurve2d) -> &'static str {
    match c {
        RawCurve2d::Unresolved => "unresolved",
        RawCurve2d::Line { .. } => "line",
        RawCurve2d::Circle { .. } => "circle",
        RawCurve2d::Ellipse { .. } => "ellipse",
        RawCurve2d::BSpline { .. } => "spline",
    }
}


fn surface(s: &Surface) -> String {
    match s {
        Surface::Plane { .. } => "plane".into(),
        Surface::Cylinder { .. } => "cylinder".into(),
        Surface::Freeform { kind, .. } | Surface::Other { kind } => (*kind).into(),
        _ => "other".into(),
    }
}


fn metadata((part, face): FaceRef<'_>) -> Value {
    let rows: Vec<_> = part.face_edges(face).into_iter().map(|e| {
        let edge = &part.edges[e];
        let alternatives = representations::alternatives((part, face),e);
        let primary: Vec<_> = part.edge_faces()[e].iter().flat_map(|&f|
            part.primary_pcurves(f).filter(move |(i,_)|*i==e).map(move |(_,c)|
                json!({"face":f,"surface":surface(&part.faces[f].surface),"kind":raw(c),"debug":format!("{c:?}")}))).collect();
        let pcurves: Vec<_> = part.file_pcurves(face).iter().filter(|(i,_)|*i==e).map(|(_,c)|raw(c)).collect();
        let feet: Vec<_> = alternatives.as_ref().into_iter().flatten().filter_map(|c| {
            let a = representations::straight_foot(c,edge.start)?;
            let b = representations::straight_foot(c,edge.end)?;
            Some(json!([geom::dist(a,edge.start),geom::dist(b,edge.end)]))
        }).collect();
        json!({"edge":e,"curve":kind(&edge.curve),"alternatives_available":alternatives.is_some(),
            "primary":primary,"pcurves":pcurves,"straight_endpoint_distances":feet})
    }).collect();
    json!({"face":face,"surface":surface(&part.faces[face].surface),
        "file_uncertainty":part.file_uncertainty(),"file_units":part.file_units(),"edges":rows})
}


fn explain(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Value {
    json!({"linear":linear::research_eligibility(patch,supports),
        "corresponding":corresponding::research_eligibility(patch,supports),
        "metadata":std::iter::once(patch).chain(supports.iter().copied())
            .filter(|s|same_surface(&patch.0.faces[patch.1].surface,&s.0.faces[s.1].surface))
            .map(metadata).collect::<Vec<_>>()})
}


pub(super) fn record(patch: FaceRef<'_>, supports: &[FaceRef<'_>], decision: Coverage) {
    let Some(path) = std::env::var_os("QUIDDITY_COVER_QUERY_TRACE") else {
        return;
    };
    static TRACE_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _guard = TRACE_LOCK.lock().unwrap();
    let geometry = |(part, face): FaceRef<'_>| {
        json!({
            "index":face,
            "face":format!("{:?}",part.faces[face]),
            "pcurves":format!("{:?}",part.file_pcurves(face)),
            "edges":part.face_edges(face).into_iter().map(|e| json!({"edge":e,"geometry":format!("{:?}",part.edges[e])})).collect::<Vec<_>>()
        })
    };
    let row = json!({"decision":format!("{decision:?}"),"detail":explain(patch,supports),
        "patch":geometry(patch),"supports":supports.iter().copied().map(geometry).collect::<Vec<_>>()});
    use std::io::Write;
    let mut out = std::fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(path)
        .unwrap();
    serde_json::to_writer(&mut out, &row).unwrap();
    out.write_all(b"\n").unwrap();
}

