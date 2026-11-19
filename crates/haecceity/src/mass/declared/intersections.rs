//! Interior crossings of selected pcurves in a common surface parameter chart.
use super::{DeclaredBoundary, EdgeTrim, on_carrier};
use crate::brep::Part;
use crate::geom::{self, Surface, V3};

struct Span<'a> {
    edge: usize,
    trim: &'a EdgeTrim,
    t: [f64; 2],
    uv: [V3; 2],
}

impl Span<'_> {
    fn left(&self) -> f64 {
        self.uv[0][0].min(self.uv[1][0])
    }
    fn right(&self) -> f64 {
        self.uv[0][0].max(self.uv[1][0])
    }
}

fn cross(a: V3, b: V3) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

/// Chord intersection supplies candidate parameters. Only a transverse
/// intersection of the actual finite pcurves, away from their STEP vertices
/// and surface singularities, can invalidate the trim.
fn crossing(a: &Span<'_>, b: &Span<'_>, part: &Part, surface: &Surface, tolerance: f64) -> bool {
    if a.edge == b.edge && (a.t[0] == b.t[1] || a.t[1] == b.t[0]) {
        return false;
    }
    let da = geom::sub(a.uv[1], a.uv[0]);
    let db = geom::sub(b.uv[1], b.uv[0]);
    let offset = geom::sub(b.uv[0], a.uv[0]);
    let det = cross(da, db);
    if det == 0. {
        return false;
    }
    let fa = cross(offset, db) / det;
    let fb = cross(offset, da) / det;
    if !(0. ..=1.).contains(&fa) || !(0. ..=1.).contains(&fb) {
        return false;
    }
    let mut s = a.t[0] + fa * (a.t[1] - a.t[0]);
    let mut t = b.t[0] + fb * (b.t[1] - b.t[0]);
    for _ in 0..12 {
        let pa = a.trim.curve.source.value(s);
        let pb = b.trim.curve.source.value(t);
        let d = geom::sub(pa, pb);
        let u = a.trim.curve.source.derivative(s);
        let v = b.trim.curve.source.derivative(t);
        let jacobian = cross(u, v);
        if jacobian.abs() <= 1e-10 * geom::norm(u) * geom::norm(v) || !jacobian.is_finite() {
            return false;
        }
        let scale = geom::norm(pa)
            .max(geom::norm(pb))
            .max(geom::norm(da))
            .max(geom::norm(db));
        if geom::norm(d) <= 64. * f64::EPSILON * scale {
            let q = a.trim.curve.carrier_uv(s);
            let Some((x, y)) = on_carrier(surface, a.trim.spans, q) else {
                return false;
            };
            let (su, sv) = surface.partials(x, y);
            if geom::norm(geom::cross(su, sv)) <= 1e-10 * geom::norm(su) * geom::norm(sv) {
                return false;
            }
            let point = surface.value(x, y);
            return [a.edge, b.edge].into_iter().all(|index| {
                let edge = &part.edges[index];
                geom::dist(point, edge.start) > tolerance && geom::dist(point, edge.end) > tolerance
            });
        }
        s -= cross(d, v) / jacobian;
        t -= cross(d, u) / jacobian;
        if !(a.t[0].min(a.t[1])..=a.t[0].max(a.t[1])).contains(&s)
            || !(b.t[0].min(b.t[1])..=b.t[0].max(b.t[1])).contains(&t)
        {
            return false;
        }
    }
    false
}

pub(super) fn crosses(
    boundary: &DeclaredBoundary,
    part: &Part,
    surface: &Surface,
    tolerance: f64,
) -> bool {
    let mut spans = Vec::new();
    for (&edge, trim) in &boundary.edges {
        if trim.curve.transfer.is_some() {
            continue;
        }
        // These coordinates precede the common unit and swept-surface mapping.
        // They describe a continuous chart even when a finite edge crosses a seam.
        let nodes = if trim.curve.line().is_some() {
            vec![trim.interval.0, trim.interval.1]
        } else if let Some(nodes) = trim.nodes() {
            nodes
        } else {
            continue;
        };
        let mut points = nodes.iter().map(|&t| trim.curve.source.value(t));
        let Some(mut prev) = points.next() else {
            continue;
        };
        for (params, point) in nodes.windows(2).zip(points) {
            spans.push(Span {
                edge,
                trim,
                t: [params[0], params[1]],
                uv: [prev, point],
            });
            prev = point;
        }
    }
    spans.sort_unstable_by(|a, b| a.left().total_cmp(&b.left()));
    for (i, a) in spans.iter().enumerate() {
        let lo = a.uv[0][1].min(a.uv[1][1]);
        let hi = a.uv[0][1].max(a.uv[1][1]);
        for b in &spans[i + 1..] {
            if b.left() > a.right() {
                break;
            }
            if b.uv[0][1].min(b.uv[1][1]) <= hi
                && b.uv[0][1].max(b.uv[1][1]) >= lo
                && crossing(a, b, part, surface, tolerance)
            {
                return true;
            }
        }
    }
    false
}
