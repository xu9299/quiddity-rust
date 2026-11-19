//! STEP selects the master representation per edge. Each usable selection keeps its
//! own parameter interval. The 3D edge and a pcurve need not use the same parameter.
use super::*;
use crate::brep::RawCurve2d;
use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock};
mod parameters;
use parameters::Axis;

#[derive(Debug)]
pub(crate) struct ParamCurve {
    source: Curve,
    axes: [Axis; 2],
    transfer: Option<Box<Transfer>>,
}

#[derive(Debug)]
struct Transfer {
    carrier: Arc<Surface>,
    face: Arc<Surface>,
    carrier_spans: (f64, f64),
    face_spans: (f64, f64),
}

/// The bounded evaluator clamps its inputs. A trim outside a carrier's native
/// patch has no point or tangent there. Crossing a genuinely closed seam is
/// different and is handled by wrapping before checking the domain.
fn on_carrier(surface: &Surface, spans: (f64, f64), q: V3) -> Option<(f64, f64)> {
    if !q[0].is_finite() || !q[1].is_finite() {
        return None;
    }
    let (u, v) = wrap_closed(surface, spans, (q[0], q[1]));
    if let Surface::Freeform { surface, .. } = surface {
        let (u0, u1, v0, v1) = surface.domain();
        let bounded = |x: f64, lo: f64, hi: f64| {
            let roundoff = 64. * f64::EPSILON * lo.abs().max(hi.abs()).max(hi - lo);
            (x >= lo - roundoff && x <= hi + roundoff).then(|| x.clamp(lo, hi))
        };
        return Some((bounded(u, u0, u1)?, bounded(v, v0, v1)?));
    }
    Some((u, v))
}

impl ParamCurve {
    fn line(&self) -> Option<(V3, V3)> {
        match &self.source {
            Curve::Line { origin, dir } => Some((*origin, *dir)),
            Curve::Nurbs(c) if c.degree() == 1 => {
                let [a, b] = c.control_points() else {
                    return None;
                };
                let dir = geom::sub(*b, *a);
                (geom::norm(dir) > 0.).then_some((*a, dir))
            }
            _ => None,
        }
    }

    fn carrier_uv(&self, t: f64) -> V3 {
        let p = self.source.value(t);
        [self.axes[0].at(p[0]).0, self.axes[1].at(p[1]).0, 0.]
    }
    fn lifted(&self, t: f64, tr: &Transfer) -> Option<((f64, f64), V3)> {
        let uv = on_carrier(&tr.carrier, tr.carrier_spans, self.carrier_uv(t))?;
        Some((uv, tr.carrier.value(uv.0, uv.1)))
    }
    pub(crate) fn value(&self, t: f64) -> V3 {
        if let Some(tr) = &self.transfer {
            self.lifted(t, tr)
                .and_then(|(_, p)| tr.face.parameters(p, None))
                .map_or([f64::NAN; 3], |(u, v)| [u, v, 0.])
        } else {
            self.carrier_uv(t)
        }
    }
    pub(crate) fn derivative(&self, t: f64) -> V3 {
        let p = self.source.value(t);
        let d = self.source.derivative(t);
        let velocity = [
            self.axes[0].at(p[0]).1 * d[0],
            self.axes[1].at(p[1]).1 * d[1],
            0.,
        ];
        let Some(tr) = &self.transfer else {
            return velocity;
        };
        let Some(((u, v), xyz)) = self.lifted(t, tr) else {
            return [f64::NAN; 3];
        };
        let (su, sv) = tr.carrier.partials(u, v);
        let tangent = geom::add(geom::scale(su, velocity[0]), geom::scale(sv, velocity[1]));
        let Some((u, v)) = tr.face.parameters(xyz, None) else {
            return [f64::NAN; 3];
        };
        let (du, dv) = uv_velocity(&tr.face, tr.face_spans, u, v, xyz, tangent);
        [du, dv, 0.]
    }
    fn parameter(&self, q: V3) -> f64 {
        let q = if let Some(tr) = &self.transfer {
            let (u, v) = wrap_closed(&tr.face, tr.face_spans, (q[0], q[1]));
            let p = tr.face.value(u, v);
            let Some((u, v)) = tr.carrier.parameters(p, None) else {
                return f64::NAN;
            };
            [u, v, 0.]
        } else {
            q
        };
        self.source
            .parameter([self.axes[0].inverse(q[0]), self.axes[1].inverse(q[1]), 0.])
    }
}

#[derive(Debug)]
pub(crate) struct EdgeTrim {
    pub(crate) curve: ParamCurve,
    spans: (f64, f64),
    pub(crate) interval: (f64, f64),
    pub(crate) deviation: f64,
    samples: OnceLock<Option<TrimSamples>>,
    nodes: OnceLock<Option<Vec<f64>>>,
    midpoint: OnceLock<Option<(V3, V3)>>,
    span_lengths: OnceLock<Option<Vec<f64>>>,
    spatial_curve: OnceLock<Option<Curve>>,
}

#[derive(Debug)]
pub(crate) struct TrimSamples {
    pub(crate) ts: Vec<f64>,
    pub(crate) points: Vec<V3>,
}

#[derive(Debug, Default)]
pub(crate) struct DeclaredBoundary {
    pub(crate) edges: BTreeMap<usize, EdgeTrim>,
    pub(crate) loops: OnceLock<Option<Vec<UvLoop>>>,
    pub(crate) domain: OnceLock<Option<crate::uv::FaceDomain>>,
}

fn stationary_parameters(ts: &mut Vec<f64>, dirs: &[V3], velocity: impl Fn(f64) -> V3) {
    let spans = ts.clone();
    for dir in dirs {
        for span in spans.windows(2) {
            let (mut lo, mut hi) = (span[0], span[1]);
            let left = geom::dot(velocity(lo), *dir);
            let right = geom::dot(velocity(hi), *dir);
            if left * right >= 0. {
                continue;
            }
            for _ in 0..52 {
                let mid = 0.5 * (lo + hi);
                if geom::dot(velocity(mid), *dir).is_sign_negative() == left.is_sign_negative() {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            ts.push(0.5 * (lo + hi));
        }
    }
}

impl EdgeTrim {
    /// Lift a planar curve by its affine frame. Analytic isoparameters on a
    /// quadric are also exact lines or circles. Other surface curves are
    /// evaluated on their carrier and have no standalone spatial primitive.
    fn spatial_curve(&self, face: &Surface) -> Option<&Curve> {
        self.spatial_curve
            .get_or_init(|| {
                let carrier = self.curve.transfer.as_ref().map_or(face, |tr| &tr.carrier);
                if let Surface::Plane { frame } = carrier {
                    let mut curve = self.curve.line().map_or_else(
                        || self.curve.source.clone(),
                        |(origin, dir)| Curve::Line { origin, dir },
                    );
                    match &mut curve {
                        Curve::Line { origin, dir } => {
                            *origin = frame.to_world(*origin);
                            *dir = frame.dir_to_world(*dir);
                        }
                        Curve::Circle { frame: f, .. } | Curve::Ellipse { frame: f, .. } => {
                            f.origin = frame.to_world(f.origin);
                            f.x = frame.dir_to_world(f.x);
                            f.y = frame.dir_to_world(f.y);
                            f.z = geom::cross(f.x, f.y);
                        }
                        Curve::Nurbs(n) => n.transform_points(|p| frame.to_world(p)),
                    }
                    return Some(curve);
                }
                let (origin, dir) = self.curve.line()?;
                let [u, v, _] = origin;
                let circle = |origin, x, y, radius: f64| {
                    (radius.abs() > 0.).then(|| Curve::Circle {
                        frame: geom::Frame {
                            origin,
                            x: geom::scale(x, radius.signum()),
                            y: geom::scale(y, radius.signum()),
                            z: geom::cross(x, y),
                        },
                        radius: radius.abs(),
                    })
                };
                let parallel = |frame: &geom::Frame, z, r| {
                    circle(frame.to_world([0., 0., z]), frame.x, frame.y, r)
                };
                match carrier {
                    Surface::Cylinder { frame, radius } if dir[1] == 0. => {
                        parallel(frame, v, *radius)
                    }
                    Surface::Cone {
                        frame,
                        radius,
                        semi_angle,
                    } if dir[1] == 0. => {
                        parallel(frame, v * semi_angle.cos(), radius + v * semi_angle.sin())
                    }
                    Surface::Sphere { frame, radius } if dir[1] == 0. => {
                        parallel(frame, radius * v.sin(), radius * v.cos())
                    }
                    Surface::Torus {
                        frame,
                        major,
                        minor,
                    } if dir[1] == 0. => parallel(frame, minor * v.sin(), major + minor * v.cos()),
                    Surface::Cylinder { .. } | Surface::Cone { .. } if dir[0] == 0. => {
                        Some(Curve::Line {
                            origin: carrier.value(u, 0.),
                            dir: carrier.partials(u, 0.).1,
                        })
                    }
                    Surface::Sphere { frame, radius } if dir[0] == 0. => circle(
                        frame.origin,
                        frame.dir_to_world([u.cos(), u.sin(), 0.]),
                        frame.z,
                        *radius,
                    ),
                    Surface::Torus {
                        frame,
                        major,
                        minor,
                    } if dir[0] == 0. => {
                        let radial = frame.dir_to_world([u.cos(), u.sin(), 0.]);
                        circle(
                            geom::add(frame.origin, geom::scale(radial, *major)),
                            radial,
                            frame.z,
                            *minor,
                        )
                    }
                    _ => None,
                }
            })
            .as_ref()
    }

    // The planar affine lift preserves the trim parameter. The elementary
    // isoparameters above instead run in the varying STEP surface coordinate.
    fn primitive_interval(&self, face: &Surface) -> (f64, f64) {
        let carrier = self.curve.transfer.as_ref().map_or(face, |tr| &tr.carrier);
        let (a, b) = self.interval;
        if matches!(carrier, Surface::Plane { .. }) {
            // A two-pole rational spline can trace a line at nonuniform speed.
            if matches!(self.curve.source, Curve::Nurbs(_))
                && let Some(line @ Curve::Line { .. }) = self.spatial_curve(face)
            {
                return (
                    line.parameter(self.point(face, a)),
                    line.parameter(self.point(face, b)),
                );
            }
            return (a, b);
        }
        let axis = if self.curve.line().is_some_and(|(_, d)| d[1] == 0.) {
            0
        } else {
            1
        };
        (
            self.curve.carrier_uv(a)[axis],
            self.curve.carrier_uv(b)[axis],
        )
    }

    pub(crate) fn point(&self, s: &Surface, t: f64) -> V3 {
        if let Some(tr) = &self.curve.transfer {
            return self.curve.lifted(t, tr).map_or([f64::NAN; 3], |(_, p)| p);
        }
        let (u, v) = self.surface_uv(s, t);
        s.value(u, v)
    }

    fn spatial_velocity(&self, s: &Surface, t: f64) -> V3 {
        let carrier = self.curve.transfer.as_ref().map_or(s, |tr| &tr.carrier);
        let spans = self
            .curve
            .transfer
            .as_ref()
            .map_or(self.spans, |tr| tr.carrier_spans);
        let uv = self.curve.carrier_uv(t);
        let (u, v) = wrap_closed(carrier, spans, (uv[0], uv[1]));
        let raw = self.curve.source.value(t);
        let d = self.curve.source.derivative(t);
        let (su, sv) = carrier.partials(u, v);
        geom::add(
            geom::scale(su, self.curve.axes[0].at(raw[0]).1 * d[0]),
            geom::scale(sv, self.curve.axes[1].at(raw[1]).1 * d[1]),
        )
    }

    fn span_lengths(&self, s: &Surface) -> Option<&[f64]> {
        self.span_lengths
            .get_or_init(|| {
                let nodes = self.nodes()?;
                nodes
                    .windows(2)
                    .map(|w| {
                        let value =
                            integrate(w, |t| [geom::norm(self.spatial_velocity(s, t))])[0].abs();
                        value.is_finite().then_some(value)
                    })
                    .collect()
            })
            .as_deref()
    }

    fn length(&self, s: &Surface) -> Option<f64> {
        if let Some(c @ (Curve::Line { .. } | Curve::Circle { .. })) = self.spatial_curve(s) {
            let (a, b) = self.primitive_interval(s);
            return Some(geom::norm(c.derivative(a)) * (b - a).abs());
        }
        Some(self.span_lengths(s)?.iter().sum())
    }

    fn point_and_tangent(&self, s: &Surface, fraction: f64) -> Option<(V3, V3)> {
        if fraction == 0.5 {
            return *self.midpoint.get_or_init(|| self.at_length(s, fraction));
        }
        self.at_length(s, fraction)
    }

    fn at_length(&self, s: &Surface, fraction: f64) -> Option<(V3, V3)> {
        if !(0.0..=1.0).contains(&fraction) {
            return None;
        }
        if let Some(curve @ (Curve::Line { .. } | Curve::Circle { .. })) = self.spatial_curve(s) {
            let (a, b) = self.primitive_interval(s);
            let t = a + fraction * (b - a);
            return Some((
                curve.value(t),
                geom::scale(curve.derivative(t), (b - a).signum()),
            ));
        }
        let nodes = self.nodes()?;
        let speed = |t| [geom::norm(self.spatial_velocity(s, t))];
        let lengths = self.span_lengths(s)?;
        let mut remaining = fraction * lengths.iter().sum::<f64>();
        for (i, (w, length)) in nodes.windows(2).zip(lengths.iter().copied()).enumerate() {
            if remaining > length && i + 2 < nodes.len() {
                remaining -= length;
                continue;
            }
            let (mut a, mut b) = (w[0], w[1]);
            for _ in 0..48 {
                let m = 0.5 * (a + b);
                if integrate(&[w[0], m], speed)[0].abs() < remaining {
                    a = m;
                } else {
                    b = m;
                }
            }
            let t = if fraction == 0. {
                self.interval.0
            } else if fraction == 1. {
                self.interval.1
            } else {
                0.5 * (a + b)
            };
            return Some((
                self.point(s, t),
                geom::scale(
                    self.spatial_velocity(s, t),
                    (self.interval.1 - self.interval.0).signum(),
                ),
            ));
        }
        None
    }

    pub(crate) fn nearest(
        &self,
        s: &Surface,
        edge: &crate::brep::Edge,
        q: V3,
    ) -> Option<(V3, bool)> {
        if let Some(curve @ (Curve::Line { .. } | Curve::Circle { .. })) = self.spatial_curve(s) {
            let (a, b) = self.primitive_interval(s);
            let (lo, hi) = (a.min(b), a.max(b));
            let t = curve.parameter(q);
            let t = curve
                .period()
                .map_or(t, |period| lo + (t - lo).rem_euclid(period));
            if edge.is_closed() || (lo..hi).contains(&t) && t != lo {
                return Some((curve.value(t), true));
            }
            let (x, y) = (curve.value(a), curve.value(b));
            return Some((
                if geom::dist(x, q) <= geom::dist(y, q) {
                    x
                } else {
                    y
                },
                false,
            ));
        }
        let samples = self.samples(s, edge)?;
        let mut best = (self.interval.0, self.point(s, self.interval.0));
        let end = (self.interval.1, self.point(s, self.interval.1));
        if geom::dist(q, end.1) < geom::dist(q, best.1) {
            best = end;
        }
        // Chords select candidate intervals only. Minimise on the lifted curve
        // and compare its actual endpoints so a vertex remains a vertex.
        for (w, p) in samples.ts.windows(2).zip(samples.points.windows(2)) {
            let chord = geom::sub(p[1], p[0]);
            let fraction =
                (geom::dot(geom::sub(q, p[0]), chord) / geom::dot(chord, chord)).clamp(0., 1.);
            let foot = geom::add(p[0], geom::scale(chord, fraction));
            if geom::dist(q, foot) > geom::dist(q, best.1) + crate::sampling::CHORD_TOLERANCE {
                continue;
            }
            let (mut a, mut b) = (w[0], w[1]);
            for _ in 0..56 {
                let l = (2. * a + b) / 3.;
                let r = (a + 2. * b) / 3.;
                if geom::dist(q, self.point(s, l)) <= geom::dist(q, self.point(s, r)) {
                    b = r;
                } else {
                    a = l;
                }
            }
            let t = 0.5 * (a + b);
            let p = self.point(s, t);
            if geom::dist(q, p) < geom::dist(q, best.1) {
                best = (t, p);
            }
        }
        Some((
            best.1,
            edge.is_closed() || (best.0 != self.interval.0 && best.0 != self.interval.1),
        ))
    }
    pub(crate) fn samples(&self, s: &Surface, edge: &crate::brep::Edge) -> Option<&TrimSamples> {
        self.samples
            .get_or_init(|| {
                let ts = if matches!(self.spatial_curve(s), Some(Curve::Line { .. })) {
                    vec![self.interval.0, self.interval.1]
                } else {
                    self.nodes()?
                };
                let eval = |t| self.point(s, t);
                let mut samples = TrimSamples {
                    ts: vec![*ts.first()?],
                    points: vec![edge.start],
                };
                for w in ts.windows(2) {
                    let mut stack = vec![(w[0], eval(w[0]), w[1], eval(w[1]), 0)];
                    while let Some((a, p, b, q, depth)) = stack.pop() {
                        if samples.points.len() >= MAX_PANELS {
                            return None;
                        }
                        let t = 0.5 * (a + b);
                        let m = eval(t);
                        if depth < 14
                            && geom::dist(m, geom::scale(geom::add(p, q), 0.5))
                                > crate::sampling::CHORD_TOLERANCE
                        {
                            stack.push((t, m, b, q, depth + 1));
                            stack.push((a, p, t, m, depth + 1));
                        } else {
                            samples.ts.push(b);
                            samples.points.push(q);
                        }
                    }
                }
                *samples.points.last_mut()? = edge.end;
                Some(samples)
            })
            .as_ref()
    }
    pub(crate) fn surface_uv(&self, s: &Surface, t: f64) -> (f64, f64) {
        let p = self.curve.value(t);
        wrap_closed(s, self.spans, (p[0], p[1]))
    }
    pub(super) fn checked_uv(&self, s: &Surface, t: f64) -> Option<V3> {
        let p = self.curve.value(t);
        on_carrier(s, self.spans, p)?;
        Some(p)
    }
    pub(crate) fn uv_extrema(&self) -> Vec<f64> {
        if self.curve.transfer.is_some() {
            let mut ts = self.nodes().unwrap_or_default();
            stationary_parameters(&mut ts, &[[1., 0., 0.], [0., 1., 0.]], |t| {
                self.curve.derivative(t)
            });
            return ts;
        }
        let mut ts = crate::sampling::extreme_parameters(
            &self.curve.source,
            self.interval,
            &[[1., 0., 0.], [0., 1., 0.]],
        );
        ts.extend([self.interval.0, self.interval.1]);
        ts
    }

    /// Bracket stationary coordinates of the lifted trim on each knot span.
    /// On a plane the lifting map is affine, so use the ordinary curve extrema.
    pub(crate) fn world_extrema(&self, surface: &Surface, dirs: &[V3]) -> Option<Vec<V3>> {
        let mut ts = self.nodes()?;
        if let Surface::Plane { frame } = surface
            && self.curve.transfer.is_none()
        {
            let local: Vec<_> = dirs.iter().map(|d| frame.dir_to_local(*d)).collect();
            ts.extend(crate::sampling::extreme_parameters(
                &self.curve.source,
                self.interval,
                &local,
            ));
        } else {
            let velocity = |t| {
                let (u, v) = self.surface_uv(surface, t);
                let d = self.curve.derivative(t);
                let (su, sv) = surface.partials(u, v);
                geom::add(geom::scale(su, d[0]), geom::scale(sv, d[1]))
            };
            stationary_parameters(&mut ts, dirs, velocity);
        }
        Some(
            ts.into_iter()
                .map(|t| {
                    let (u, v) = self.surface_uv(surface, t);
                    surface.value(u, v)
                })
                .collect(),
        )
    }
    pub(crate) fn oriented(&self, forward: bool) -> (f64, f64) {
        if forward {
            self.interval
        } else {
            (self.interval.1, self.interval.0)
        }
    }

    pub(crate) fn nodes(&self) -> Option<Vec<f64>> {
        self.nodes.get_or_init(|| self.compute_nodes()).clone()
    }

    fn compute_nodes(&self) -> Option<Vec<f64>> {
        let (a, b) = self.interval;
        let breaks = match &self.curve.source {
            Curve::Nurbs(c) => c.breaks(),
            _ => Vec::new(),
        };
        let mut ts = cuts(a, b, 16, &breaks)?;
        // A foreign spline's partials can jump at a knot even when this
        // pcurve and the receiving face are smooth. Split the path there
        // before budgeting the higher-moment integrals on the receiving face.
        if let Some(tr) = &self.curve.transfer
            && let Surface::Freeform { surface, .. } = &*tr.carrier
        {
            let spans = [tr.carrier_spans.0, tr.carrier_spans.1];
            let knots = [
                crate::nurbs::breaks(&surface.knots_u, surface.degree_u),
                crate::nurbs::breaks(&surface.knots_v, surface.degree_v),
            ];
            let grid = ts.clone();
            for w in grid.windows(2) {
                let from = self.curve.carrier_uv(w[0]);
                let to = self.curve.carrier_uv(w[1]);
                for axis in 0..2 {
                    let (lo, hi) = (from[axis].min(to[axis]), from[axis].max(to[axis]));
                    if hi <= lo {
                        continue;
                    }
                    for &knot in &knots[axis] {
                        let span = spans[axis];
                        let mut crossing = if span.is_finite() {
                            knot + ((lo - knot) / span).ceil() * span
                        } else {
                            knot
                        };
                        while crossing < hi {
                            if crossing > lo {
                                let (mut left, mut right) = (w[0], w[1]);
                                for _ in 0..52 {
                                    let mid = 0.5 * (left + right);
                                    if (self.curve.carrier_uv(mid)[axis] < crossing)
                                        == (from[axis] < crossing)
                                    {
                                        left = mid;
                                    } else {
                                        right = mid;
                                    }
                                }
                                ts.push(0.5 * (left + right));
                                if ts.len() >= MAX_PANELS {
                                    return None;
                                }
                            }
                            if !span.is_finite() {
                                break;
                            }
                            let next = crossing + span;
                            if next <= crossing {
                                return None;
                            }
                            crossing = next;
                        }
                    }
                }
            }
            ts.sort_by(|x, y| {
                if a <= b {
                    x.total_cmp(y)
                } else {
                    y.total_cmp(x)
                }
            });
            ts.dedup();
        }
        Some(ts)
    }
}

impl DeclaredBoundary {
    /// Keep #281's physical disagreement guard on the total trimmed area,
    /// including the perimeter of holes and both uses of a seam.
    pub(super) fn within_strip(&self, part: &Part, face: usize, area: f64) -> Option<bool> {
        let deviation = self.edges.values().map(|c| c.deviation).fold(0.0, f64::max);
        if deviation == 0. {
            return Some(true);
        }
        let f = &part.faces[face];
        let mut perimeter = 0.;
        for &(edge, _) in f.loops.iter().flat_map(|lp| &lp.edges) {
            let ed = &part.edges[edge];
            let trim = self.edges.get(&edge);
            let cut = if let Some(c) = trim {
                c.nodes()?
            } else {
                let (a, b) =
                    edge_interval(&ed.curve, ed.start, ed.end, ed.same_sense, ed.is_closed());
                curve_cuts(&f.surface, &ed.curve, a, b)?
            };
            let speed = |t| {
                [if let Some(c) = trim {
                    let (u, v) = c.surface_uv(&f.surface, t);
                    let d = c.curve.derivative(t);
                    let (su, sv) = f.surface.partials(u, v);
                    geom::norm(geom::add(geom::scale(su, d[0]), geom::scale(sv, d[1])))
                } else {
                    geom::norm(ed.curve.derivative(t))
                }]
            };
            let coarse = integrate(&cut, speed)[0].abs();
            let (sum, error) = preferred::integrate_inner(
                &cut,
                1e-11 * coarse,
                speed,
                IntegralRule {
                    components: Components::Prefix(1),
                    observe_mass: None,
                },
            )?;
            perimeter += sum[0].abs() + error;
        }
        let ratio = deviation * perimeter / area;
        Some(ratio.is_finite() && ratio <= 0.001)
    }
}

impl Part {
    pub(crate) fn spatial_trim(&self, e: usize) -> Option<(&Surface, &EdgeTrim)> {
        if !self.has_master_pcurves() {
            return None;
        }
        let owner = self.primary_owner(e).map(|(f, _)| f);
        if owner.is_none() && !self.detached_pcurves.contains_key(&e) {
            return None;
        }
        owner
            .into_iter()
            .chain(
                self.edge_faces()[e]
                    .iter()
                    .copied()
                    .filter(|&f| Some(f) != owner),
            )
            .find_map(|f| {
                Some((
                    &self.faces[f].surface,
                    self.measured_boundary(f)?.edges.get(&e)?,
                ))
            })
    }

    /// One physical discretisation per topological edge for mesh joins and ray candidates.
    pub fn boundary_edge_samples(&self, e: usize) -> &[V3] {
        self.spatial_trim(e)
            .and_then(|(s, c)| c.samples(s, &self.edges[e]))
            .map_or(&self.edges[e].samples, |samples| &samples.points)
    }

    /// The standalone curve of the selected boundary when one exists. A general
    /// surface curve is still available through its samples and closest-point query.
    pub fn boundary_edge_curve(&self, e: usize) -> Option<&Curve> {
        match self.spatial_trim(e) {
            Some((s, c)) => c.spatial_curve(s),
            None => Some(&self.edges[e].curve),
        }
    }

    /// Start and end parameters on the selected standalone curve, in vertex order.
    pub fn boundary_edge_interval(&self, e: usize) -> Option<(f64, f64)> {
        if let Some((s, c)) = self.spatial_trim(e) {
            c.spatial_curve(s)?;
            return Some(c.primitive_interval(s));
        }
        let edge = &self.edges[e];
        Some(edge_interval(
            &edge.curve,
            edge.start,
            edge.end,
            edge.same_sense,
            edge.is_closed(),
        ))
    }

    /// Directional extrema of an edge whose master is a surface curve.
    pub fn declared_edge_extrema(&self, e: usize, dirs: &[V3]) -> Option<Vec<V3>> {
        let (s, c) = self.spatial_trim(e)?;
        if let Some(curve) = c.spatial_curve(s) {
            let interval = c.primitive_interval(s);
            let mut points = vec![curve.value(interval.0), curve.value(interval.1)];
            points.extend(crate::sampling::extremes_along(curve, interval, dirs));
            return Some(points);
        }
        let mut ts = c.nodes()?;
        stationary_parameters(&mut ts, dirs, |t| c.spatial_velocity(s, t));
        Some(ts.into_iter().map(|t| c.point(s, t)).collect())
    }

    /// Point and tangent at a fraction of the selected surface curve's length
    /// from the start vertex. A spatial edge uses the caller's ordinary curve.
    pub fn declared_edge_at_length(&self, e: usize, fraction: f64) -> Option<(V3, V3)> {
        let (s, c) = self.spatial_trim(e)?;
        c.point_and_tangent(s, fraction)
    }

    /// Physical length of a file-selected surface curve, in millimetres.
    pub fn declared_edge_length(&self, edge: usize) -> Option<f64> {
        let (surface, trim) = self.spatial_trim(edge)?;
        trim.length(surface)
    }

    /// Upper bound in mm on the deviation from the segment between the STEP
    /// vertices, for a straight UV trim on a cylinder. In affine UV parameter
    /// on [0,1] the axial component is linear and the second derivative has
    /// norm r·Δu². Linear interpolation therefore contributes at most r·Δu²/8.
    pub fn declared_cylinder_chord_error(&self, e: usize) -> Option<f64> {
        let (face, trim) = self.spatial_trim(e)?;
        let carrier = trim.curve.transfer.as_ref().map_or(face, |tr| &tr.carrier);
        let Surface::Cylinder { radius, .. } = carrier else {
            return None;
        };
        trim.curve.line()?;
        let edge = &self.edges[e];
        if edge.is_closed() {
            return None;
        }
        let (a, b) = trim.interval;
        let turn = trim.curve.carrier_uv(b)[0] - trim.curve.carrier_uv(a)[0];
        let ends = geom::dist(trim.point(face, a), edge.start)
            .max(geom::dist(trim.point(face, b), edge.end));
        Some(ends + radius.abs() * turn * turn / 8.)
    }

    /// The arc-length midpoint of the selected spatial boundary of edge `e`.
    pub fn edge_midpoint(&self, e: usize) -> V3 {
        self.declared_edge_at_length(e, 0.5)
            .map(|(p, _)| p)
            .unwrap_or_else(|| self.edges[e].midpoint())
    }
    /// A physical scale independent of trimmed-face queries (which themselves
    /// use this choice). A rigid placement or a different query order cannot
    /// introduce a recursive dependency on face bounds.
    pub(super) fn trim_reference_length(&self, face: usize) -> f64 {
        let mut bounds = geom::Bounds::empty();
        for lp in &self.faces[face].loops {
            for &(e, _) in &lp.edges {
                for &p in &self.edges[e].samples {
                    bounds.add(p);
                }
            }
            if let Some(p) = lp.vertex {
                bounds.add(p);
            }
        }
        bounds.diagonal()
    }

    pub(crate) fn declared_boundary(&self, face: usize) -> Option<&DeclaredBoundary> {
        self.cache[face]
            .declared
            .get_or_init(|| self.read_declared_boundary(face).map(Box::new))
            .as_deref()
    }

    fn read_declared_boundary(&self, face: usize) -> Option<DeclaredBoundary> {
        if !self.has_master_pcurves() {
            return None;
        }
        if !self.faces[face]
            .loops
            .iter()
            .flat_map(|lp| &lp.edges)
            .any(|(e, _)| self.primary_owner(*e).is_some() || self.detached_pcurves.contains_key(e))
        {
            return None;
        }
        let entries = self.primary_pcurves(face);
        let surface = &self.faces[face].surface;
        let units = self.file_units()?;
        let make_curve = |raw: &RawCurve2d, carrier: &Surface, owner: Option<usize>| {
            let (scales, axes) = match carrier.file_parameter_scales(units) {
                Some(scales) => (scales, [Axis::ID; 2]),
                None => ((1., 1.), parameters::swept_axes(self, owner?)?),
            };
            Some(ParamCurve {
                source: scaled_pcurve(raw, scales.0, scales.1)?,
                axes,
                transfer: None,
            })
        };
        let length = self.trim_reference_length(face);
        let tolerance = self.file_uncertainty().unwrap_or(0.).max(0.) + geom::COORD_FLOOR;
        let mut result = DeclaredBoundary::default();
        let mut previous = BTreeMap::new();
        for (index, raw) in entries {
            if let Some(other) = previous.insert(*index, raw) {
                // A seam lists the same slots once per use. A second distinct
                // selection for a single edge does not identify a unique curve.
                let occurrences = self.faces[face]
                    .loops
                    .iter()
                    .flat_map(|lp| &lp.edges)
                    .filter(|(e, _)| e == index)
                    .count();
                if other != raw || occurrences < 2 {
                    return None;
                }
                continue;
            }
            let curve = make_curve(raw, surface, Some(face))?;
            let edge = &self.edges[*index];
            result
                .edges
                .insert(*index, resolve(surface, edge, curve, tolerance, length)?);
        }
        let shared = |f: usize| {
            self.cache[f]
                .carrier
                .get_or_init(|| Arc::new(self.faces[f].surface.clone()))
                .clone()
        };
        for &(e, _) in self.faces[face].loops.iter().flat_map(|lp| &lp.edges) {
            if result.edges.contains_key(&e) {
                continue;
            }
            let (carrier, raw, owner) = if let Some((f, raw)) = self.primary_owner(e) {
                (shared(f), raw, Some(f))
            } else if let Some((s, raw)) = self.detached_pcurves.get(&e) {
                (Arc::new(s.placed(&self.edge_placement(e))), raw, None)
            } else {
                continue;
            };
            let mut curve = make_curve(raw, &carrier, owner)?;
            curve.transfer = Some(Box::new(Transfer {
                carrier_spans: closed_spans(&carrier),
                face_spans: closed_spans(surface),
                carrier,
                face: shared(face),
            }));
            result.edges.insert(
                e,
                resolve(surface, &self.edges[e], curve, tolerance, length)?,
            );
        }
        (!result.edges.is_empty()).then_some(result)
    }
}

fn scaled_pcurve(raw: &RawCurve2d, su: f64, sv: f64) -> Option<Curve> {
    // A plane scales both axes equally. Preserve the circle or ellipse and
    // its angular parameter so downstream conic calculations stay analytic.
    if let RawCurve2d::Conic {
        origin,
        reference,
        major,
        minor,
    } = raw
        && su == sv
        && su.is_finite()
        && su > 0.
    {
        let [x, y] = reference.unwrap_or([1., 0.]);
        let x = geom::unit([x, y, 0.])?;
        let z = [0., 0., 1.];
        let frame = geom::Frame {
            origin: [origin[0] * su, origin[1] * sv, 0.],
            x,
            y: geom::cross(z, x),
            z,
        };
        let a = major * su;
        let b = minor.unwrap_or(*major) * sv;
        if !(a.is_finite() && b.is_finite() && a > 0. && b > 0.) {
            return None;
        }
        return Some(if a == b {
            Curve::Circle { frame, radius: a }
        } else {
            Curve::Ellipse {
                frame,
                major: a,
                minor: b,
            }
        });
    }
    preferred::raw_param_curve(&raw.scaled(su, sv)?)
}

fn resolve(
    surface: &Surface,
    edge: &crate::brep::Edge,
    curve: ParamCurve,
    tolerance: f64,
    length: f64,
) -> Option<EdgeTrim> {
    let spatial = edge_interval(
        &edge.curve,
        edge.start,
        edge.end,
        edge.same_sense,
        edge.is_closed(),
    );
    let spans = closed_spans(surface);
    let lift = |t| {
        let q = curve.value(t);
        let (u, v) = wrap_closed(surface, spans, (q[0], q[1]));
        surface.value(u, v)
    };
    if !(length.is_finite() && length > 0.) {
        return None;
    }
    let natural = match &curve.source {
        Curve::Nurbs(c) => Some(c.domain()),
        _ => None,
    };
    let period = curve.source.period();
    let valid = |(a, b): (f64, f64)| {
        if !a.is_finite() || !b.is_finite() || a == b {
            return false;
        }
        let in_domain = natural.is_none_or(|(lo, hi)| {
            let eps = 64. * f64::EPSILON * lo.abs().max(hi.abs()).max(hi - lo);
            if period.is_some() {
                (b - a).abs() <= hi - lo + eps
            } else {
                a.min(b) >= lo - eps && a.max(b) <= hi + eps
            }
        });
        // A bounded NURBS evaluator clamps outside parameters to its endpoint.
        // Such a match must not suppress inversion of the pcurve's own endpoints.
        let on_surface = |t| on_carrier(surface, spans, curve.value(t)).is_some();
        in_domain
            && on_surface(a)
            && on_surface(b)
            && geom::dist(lift(a), edge.start) <= tolerance
            && geom::dist(lift(b), edge.end) <= tolerance
    };
    let (pu, pv) = surface.periodic();
    let periods = [
        if pu { std::f64::consts::TAU } else { spans.0 },
        if pv { std::f64::consts::TAU } else { spans.1 },
    ];
    let guess = |p, near: V3| {
        let (u, v) = surface.parameters(p, Some((near[0], near[1])))?;
        let mut q = [u, v, 0.];
        for k in 0..2 {
            if periods[k].is_finite() {
                q[k] += periods[k] * ((near[k] - q[k]) / periods[k]).round();
            }
        }
        Some(curve.parameter(q))
    };
    let guessed = || {
        let a = guess(edge.start, curve.value(natural.map_or(0., |ab| ab.0)))?;
        let mut b = guess(edge.end, curve.value(a))?;
        if let Some(period) = period {
            if edge.is_closed() {
                b = a + if edge.same_sense { period } else { -period };
            } else if (b > a) != edge.same_sense {
                b += if edge.same_sense { period } else { -period };
            }
        }
        Some((a, b))
    };
    let interval = [
        Some(spatial),
        natural.filter(|_| edge.same_sense),
        natural.map(|(a, b)| (b, a)).filter(|_| !edge.same_sense),
    ]
    .into_iter()
    .flatten()
    .find(|&ab| valid(ab))
    .or_else(|| guessed().filter(|&ab| valid(ab)))?;
    let mut trim = EdgeTrim {
        curve,
        spans,
        interval,
        deviation: 0.,
        samples: OnceLock::new(),
        nodes: OnceLock::new(),
        midpoint: OnceLock::new(),
        span_lengths: OnceLock::new(),
        spatial_curve: OnceLock::new(),
    };
    let (a, b) = interval;
    let knots = match &trim.curve.source {
        Curve::Nurbs(c) => c.breaks(),
        _ => Vec::new(),
    };
    let mut breaks = cuts(a, b, 1, &knots)?;
    if let Curve::Nurbs(c) = &edge.curve {
        for k in c.breaks() {
            let x = (k - spatial.0) / (spatial.1 - spatial.0);
            if (0.0..1.).contains(&x) {
                breaks.push(a + (b - a) * x);
            }
        }
        breaks.sort_by(|x, y| {
            if b >= a {
                x.total_cmp(y)
            } else {
                y.total_cmp(x)
            }
        });
        breaks.dedup();
    }
    // Compare geometry in physical units. Corresponding vertices do not imply
    // that independent curve parameters advance at the same speed (e.g. the
    // rational-quadratic and angular parametrisations of exactly the same arc).
    let bound = 0.001 * length;
    for span in breaks.windows(2) {
        for j in 0..=32 {
            let t = span[0] + (span[1] - span[0]) * (j as f64 / 32.);
            let (u, v) = on_carrier(surface, spans, trim.curve.value(t))?;
            let p = surface.value(u, v);
            if let Some(tr) = &trim.curve.transfer {
                let (_, xyz) = trim.curve.lifted(t, tr)?;
                let residual = geom::dist(p, xyz);
                if !residual.is_finite() || residual > tolerance {
                    return None;
                }
                trim.deviation = trim.deviation.max(residual);
            }
            let t3 = spatial.0 + (spatial.1 - spatial.0) * ((t - a) / (b - a));
            let mut d = geom::dist(p, edge.curve.value(t3));
            if d > geom::COORD_FLOOR {
                let mut foot = edge.curve.parameter(p);
                if let Some(period) = edge.curve.period() {
                    foot += period * ((t3 - foot) / period).round();
                }
                let foot = foot.clamp(spatial.0.min(spatial.1), spatial.0.max(spatial.1));
                d = d.min(geom::dist(p, edge.curve.value(foot)));
            }
            if !d.is_finite() || d > bound {
                return None;
            }
            trim.deviation = trim.deviation.max(d);
        }
    }
    Some(trim)
}

#[cfg(test)]
mod tests;
