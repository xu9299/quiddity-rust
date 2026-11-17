//! Coverage established by a bounded correspondence of complete finite boundaries.
//!
//! This path proves only acceptance: each target loop must have a distinct support loop,
//! with every finite edge paired in cyclic order (either orientation). A pointwise curve
//! displacement d confines the changed parity region to a tube of area at most 2*d*L+pi*d²,
//! where L bounds the reference curve's length. Independent pcurves and endpoint joins are
//! charged separately. A coarse polygon supplies only a conservative lower area bound for
//! the original relative budget; it is never the coverage evidence.

use super::{FaceRef, conic_bounds, linear, representations, same_surface};
use crate::brep::Edge;
use crate::geom::{self, Curve, Surface, V3};
use crate::sampling::edge_interval;

type Loops = Vec<Vec<(f64, f64)>>;

#[derive(Clone)]
struct Trim {
    curve: Curve,
    a: f64,
    b: f64,
    start: V3,
    end: V3,
}

impl Trim {
    fn new(edge: &Edge, curve: Curve, forward: bool, limit: f64) -> Option<Self> {
        let (a, b) = edge_interval(
            &curve,
            edge.start,
            edge.end,
            edge.same_sense,
            edge.vertices.0 == edge.vertices.1,
        );
        if ![a, b].into_iter().all(f64::is_finite)
            || geom::dist(
                edge.start,
                representations::straight_foot(&curve, edge.start)
                    .unwrap_or_else(|| curve.value(a)),
            ) > limit
            || geom::dist(
                edge.end,
                representations::straight_foot(&curve, edge.end).unwrap_or_else(|| curve.value(b)),
            ) > limit
        {
            return None;
        }
        let trim = Self {
            curve,
            a,
            b,
            start: edge.start,
            end: edge.end,
        };
        Some(if forward { trim } else { trim.reversed() })
    }

    fn reversed(&self) -> Self {
        Self {
            curve: self.curve.clone(),
            a: self.b,
            b: self.a,
            start: self.end,
            end: self.start,
        }
    }

    fn split_near(&self, point: V3) -> Option<(Self, Self)> {
        let mut t = self.curve.parameter(point);
        if conic(&self.curve).is_some() {
            t = geom::nearest_turn(t, 0.5 * (self.a + self.b));
        }
        let fraction = (t - self.a) / (self.b - self.a);
        if !fraction.is_finite() || fraction <= 0.0 || fraction >= 1.0 {
            return None;
        }
        let middle = self.curve.value(t);
        let mut first = self.clone();
        first.b = t;
        first.end = middle;
        let mut last = self.clone();
        last.a = t;
        last.start = middle;
        Some((first, last))
    }

    fn length(&self) -> Option<f64> {
        if let Some(a) = representations::straight_foot(&self.curve, self.start) {
            return Some(geom::dist(
                a,
                representations::straight_foot(&self.curve, self.end)?,
            ));
        }
        Some(match &self.curve {
            Curve::Line { .. } => geom::dist(self.curve.value(self.a), self.curve.value(self.b)),
            _ => {
                let (_, x, y) = conic(&self.curve)?;
                // The largest singular value bounds every unit combination of the two
                // conic axes, including the derivative's (-sin(t), cos(t)) coefficients.
                let (xx, yy, xy) = (geom::dot(x, x), geom::dot(y, y), geom::dot(x, y));
                let maximum = 0.5 * (xx + yy + ((xx - yy).powi(2) + 4.0 * xy * xy).sqrt());
                maximum.sqrt() * (self.b - self.a).abs() * (1.0 + 64.0 * f64::EPSILON)
            }
        })
    }
}

struct Boundary {
    loops: Vec<Vec<Trim>>,
    error: f64,
}

pub(super) fn covered(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<bool> {
    let surface = &patch.0.faces[patch.1].surface;
    let mine = boundary(patch, surface)?;
    let lower_area = lower_area(&mine, patch)?;
    for &support in supports {
        if !same_surface(surface, &support.0.faces[support.1].surface) {
            continue;
        }
        let Some(theirs) = boundary(support, surface) else {
            continue;
        };
        if let Some(error) = correspondence(&mine, &theirs, surface)
            && error <= 1e-9 * lower_area
        {
            return Some(true);
        }
    }
    cylindrical_union(patch, supports)
}

/// Removing one loop changes a parity region only inside that loop. A complete one-loop
/// support matching that finite boundary fills the change, up to this charged area bound.
/// This helper never recursively removes loops or treats an unmatched loop as filled.
pub(super) fn loop_match_error(
    owner: FaceRef<'_>,
    loop_index: usize,
    supports: &[FaceRef<'_>],
) -> Option<f64> {
    let surface = &owner.0.faces[owner.1].surface;
    let mine = boundary(owner, surface)?;
    let contour = mine.loops.get(loop_index)?;
    supports
        .iter()
        .copied()
        .filter(|s| !(std::ptr::eq(s.0, owner.0) && s.1 == owner.1))
        .filter(|s| {
            s.0.faces[s.1].loops.len() == 1 && same_surface(surface, &s.0.faces[s.1].surface)
        })
        .filter_map(|s| {
            let other = boundary(s, surface)?;
            let error = mine.error + other.error + loop_error(contour, &other.loops[0], surface)?;
            error.is_finite().then_some(error)
        })
        .min_by(f64::total_cmp)
}

fn size(curve: &Curve) -> f64 {
    let magnitude = |p: V3| p.into_iter().map(f64::abs).fold(1.0, f64::max);
    match curve {
        Curve::Line { origin, dir } => magnitude(*origin).max(magnitude(*dir)),
        Curve::Circle { frame, radius } => magnitude(frame.origin).max(radius.abs()),
        Curve::Ellipse {
            frame,
            major,
            minor,
        } => magnitude(frame.origin).max(major.abs()).max(minor.abs()),
        Curve::Nurbs(c) => c
            .control_points()
            .iter()
            .copied()
            .map(magnitude)
            .fold(1.0, f64::max),
    }
}

fn boundary((part, face): FaceRef<'_>, surface: &Surface) -> Option<Boundary> {
    let mut loops = Vec::new();
    let mut error = 0.0;
    let uncertainty = part
        .file_uncertainty()
        .filter(|x| x.is_finite() && *x >= 0.0)
        .unwrap_or(f64::INFINITY);
    for lp in &part.faces[face].loops {
        if lp.edges.is_empty() || lp.edges.len() > 512 {
            return None;
        }
        let mut trims = Vec::new();
        for &(index, forward) in &lp.edges {
            if part.unresolved_edges().contains(&index) {
                return None;
            }
            let edge = &part.edges[index];
            let pad = 128.0
                * f64::EPSILON
                * size(&edge.curve)
                    .max(edge.start.into_iter().map(f64::abs).fold(1.0, f64::max))
                    .max(edge.end.into_iter().map(f64::abs).fold(1.0, f64::max));
            let trim = Trim::new(edge, edge.curve.clone(), forward, uncertainty + pad)?;
            if trims.last().is_some_and(|t: &Trim| t.end != trim.start) {
                return None;
            }
            error += tube(&trim, surface, pad)?;
            for (v, t) in [(trim.start, trim.a), (trim.end, trim.b)] {
                let d = 4.0 * (geom::dist(v, trim.curve.value(t)) + pad);
                error += 4.0 * d * d;
            }
            for curve in representations::alternatives((part, face), index)? {
                let alternate = Trim::new(edge, curve, forward, uncertainty + pad)?;
                error += tube(&trim, surface, displacement(&trim, &alternate, surface)?)?;
            }
            trims.push(trim);
        }
        if trims.last()?.end != trims.first()?.start {
            return None;
        }
        loops.push(trims);
    }
    (!loops.is_empty() && error.is_finite()).then_some(Boundary { loops, error })
}

fn conic(curve: &Curve) -> Option<(V3, V3, V3)> {
    match curve {
        Curve::Circle { frame, radius } => Some((
            frame.origin,
            geom::scale(frame.x, *radius),
            geom::scale(frame.y, *radius),
        )),
        Curve::Ellipse {
            frame,
            major,
            minor,
        } => Some((
            frame.origin,
            geom::scale(frame.x, *major),
            geom::scale(frame.y, *minor),
        )),
        _ => None,
    }
}

/// Compare equal fractions of the finite intervals, rather than untrimmed conic boxes.
fn displacement(a: &Trim, b: &Trim, surface: &Surface) -> Option<f64> {
    let pad = 128.0 * f64::EPSILON * size(&a.curve).max(size(&b.curve));
    if matches!(b.curve, Curve::Nurbs(_)) && conic(&a.curve).is_some() {
        return if let Surface::Plane { frame } = surface {
            conic_bounds::planar_deviation(&a.curve, (a.a, a.b), &b.curve, (b.a, b.b), frame.z)
        } else {
            conic_bounds::deviation(&a.curve, (a.a, a.b), &b.curve, (b.a, b.b))
        };
    }
    let normal = match surface {
        Surface::Plane { frame } => Some(geom::unit(frame.z)?),
        _ => None,
    };
    let length =
        |p| geom::norm(normal.map_or(p, |n| geom::sub(p, geom::scale(n, geom::dot(p, n)))));
    let distance = |p, q| length(geom::sub(p, q));
    if let (Some(a0), Some(b0)) = (
        representations::straight_foot(&a.curve, a.start),
        representations::straight_foot(&b.curve, b.start),
    ) {
        return Some(
            distance(a0, b0).max(distance(
                representations::straight_foot(&a.curve, a.end)?,
                representations::straight_foot(&b.curve, b.end)?,
            )) + pad,
        );
    }
    let mut aligned = b.clone();
    let whole_conic = |t: &Trim| {
        t.start == t.end
            && conic(&t.curve).is_some()
            && ((t.b - t.a).abs() - std::f64::consts::TAU).abs()
                <= 256.0 * f64::EPSILON * std::f64::consts::TAU
    };
    if whole_conic(a) && whole_conic(b) {
        let sense = (b.b - b.a).signum();
        aligned.a = b.curve.parameter(a.curve.value(a.a));
        aligned.b = aligned.a + sense * std::f64::consts::TAU;
    }
    let b = &aligned;
    let normalized = |t: &Trim| {
        let (o, x, y) = conic(&t.curve)?;
        let p = geom::add(geom::scale(x, t.a.cos()), geom::scale(y, t.a.sin()));
        let q = geom::scale(
            geom::add(geom::scale(x, -t.a.sin()), geom::scale(y, t.a.cos())),
            (t.b - t.a).signum(),
        );
        Some((o, p, q, (t.b - t.a).abs()))
    };
    let (o, x, y, span) = normalized(a)?;
    let (p, u, v, other_span) = normalized(b)?;
    let d = distance(o, p)
        + distance(x, u)
        + distance(y, v)
        + (length(u) + length(v)) * (span - other_span).abs()
        + pad;
    d.is_finite().then_some(d)
}

/// A lower radial bound keeps the cylindrical parameter map away from its axis.
fn radial_lower(trim: &Trim, surface: &Surface) -> Option<f64> {
    let Surface::Cylinder { frame, .. } = surface else {
        return Some(f64::INFINITY);
    };
    let across = |p| geom::sub(p, geom::scale(frame.z, geom::dot(p, frame.z)));
    if let Some(p) = representations::straight_foot(&trim.curve, trim.start) {
        let a = across(geom::sub(p, frame.origin));
        let b = across(geom::sub(
            representations::straight_foot(&trim.curve, trim.end)?,
            frame.origin,
        ));
        let d = geom::sub(b, a);
        let t = if geom::dot(d, d) > 0.0 {
            (-geom::dot(a, d) / geom::dot(d, d)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        return Some(geom::norm(geom::add(a, geom::scale(d, t))));
    }
    let (o, a, b) = conic(&trim.curve)?;
    let (c, a, b) = (across(geom::sub(o, frame.origin)), across(a), across(b));
    let square = geom::dot(c, c) + 0.5 * (geom::dot(a, a) + geom::dot(b, b))
        - 2.0 * geom::dot(c, a).abs()
        - 2.0 * geom::dot(c, b).abs()
        - 0.5 * (geom::dot(a, a) - geom::dot(b, b)).abs()
        - geom::dot(a, b).abs();
    Some(square.max(0.0).sqrt())
}

fn tube(trim: &Trim, surface: &Surface, distance: f64) -> Option<f64> {
    let (d, length) = match surface {
        Surface::Plane { .. } => (distance, trim.length()?),
        Surface::Cylinder { radius, .. } if *radius > 0.0 && distance < 0.25 * radius => {
            let radial = radial_lower(trim, surface)? - 128.0 * f64::EPSILON * size(&trim.curve);
            if radial < 0.75 * radius {
                return None;
            }
            // In physical parameters (R*u,v), the map's Jacobian has orthogonal rows
            // with norms R/rho and 1. A straight homotopy of displacement d stays at
            // rho >= radial-d; the reference curve itself stays at rho >= radial.
            let padding = 1.0 + 128.0 * f64::EPSILON;
            let displacement_factor = (radius / (radial - distance)).max(1.0) * padding;
            let length_factor = (radius / radial).max(1.0) * padding;
            (
                displacement_factor * distance,
                length_factor * trim.length()?,
            )
        }
        _ => return None,
    };
    let area = 2.0 * d * length + std::f64::consts::PI * d * d;
    (d >= 0.0 && area.is_finite()).then_some(area)
}

fn correspondence(a: &Boundary, b: &Boundary, surface: &Surface) -> Option<f64> {
    if a.loops.len() != b.loops.len() {
        return None;
    }
    let mut used = vec![false; b.loops.len()];
    let mut error = a.error + b.error;
    for mine in &a.loops {
        let best = b
            .loops
            .iter()
            .enumerate()
            .filter(|(i, _)| !used[*i])
            .filter_map(|(i, other)| loop_error(mine, other, surface).map(|e| (i, e)))
            .min_by(|a, b| a.1.total_cmp(&b.1))?;
        used[best.0] = true;
        error += best.1;
    }
    error.is_finite().then_some(error)
}

fn loop_error(a: &[Trim], b: &[Trim], surface: &Surface) -> Option<f64> {
    let mut best: Option<f64> = None;
    for reverse in [false, true] {
        for shift in 0..b.len() {
            let ordered: Vec<_> = (0..b.len())
                .map(|i| {
                    let index = if reverse {
                        (shift + b.len() - i) % b.len()
                    } else {
                        (shift + i) % b.len()
                    };
                    if reverse {
                        b[index].reversed()
                    } else {
                        b[index].clone()
                    }
                })
                .collect();
            let cost = if a.len() == b.len() {
                a.iter().zip(&ordered).try_fold(0.0, |sum, (mine, other)| {
                    Some(sum + tube(mine, surface, displacement(mine, other, surface)?)?)
                })
            } else if a.len() + b.len() <= 64 {
                subdivided_loop_error(a, &ordered, surface)
            } else {
                None
            };
            if let Some(cost) = cost {
                best = Some(best.map_or(cost, |x| x.min(cost)));
            }
        }
    }
    best
}

/// Splits only reparameterize the original finite curves. Every step consumes at least
/// one whole remaining edge; the final correspondence must consume both closed loops.
/// Closest parameters propose splits, while the complete homotopy tubes prove the bound.
fn subdivided_loop_error(a: &[Trim], b: &[Trim], surface: &Surface) -> Option<f64> {
    let key = |i, j, mine: &Trim, other: &Trim| (i + j, i, j, mine.a.to_bits(), other.a.to_bits());
    let mine = a.first()?.clone();
    let other = b.first()?.clone();
    let mut paths = std::collections::BTreeMap::new();
    paths.insert(key(0, 0, &mine, &other), (0.0, mine, other));
    let mut best: Option<f64> = None;
    let mut visited = 0;
    while let Some(((_, i, j, _, _), (error, mine, other))) = paths.pop_first() {
        visited += 1;
        if visited + paths.len() > 4096 {
            return None;
        }
        if best.is_some_and(|bound| error >= bound) {
            continue;
        }
        let cost = |x: &Trim, y: &Trim| tube(x, surface, displacement(x, y, surface)?);
        let mut choices = Vec::new();
        if let Some(c) = cost(&mine, &other) {
            choices.push((c, None, None));
        }
        if let Some((first, last)) = mine.split_near(other.curve.value(other.b))
            && let Some(c) = cost(&first, &other)
        {
            choices.push((c, Some(last), None));
        }
        if let Some((first, last)) = other.split_near(mine.curve.value(mine.b))
            && let Some(c) = cost(&mine, &first)
        {
            choices.push((c, None, Some(last)));
        }
        for (bound, rest_a, rest_b) in choices {
            let total = error + bound;
            if !total.is_finite() || best.is_some_and(|cost| total >= cost) {
                continue;
            }
            let next_i = i + usize::from(rest_a.is_none());
            let next_j = j + usize::from(rest_b.is_none());
            if next_i == a.len() || next_j == b.len() {
                if next_i == a.len() && next_j == b.len() {
                    best = Some(best.map_or(total, |bound| bound.min(total)));
                }
                continue;
            }
            let next_a = rest_a.unwrap_or_else(|| a[next_i].clone());
            let next_b = rest_b.unwrap_or_else(|| b[next_j].clone());
            let state = key(next_i, next_j, &next_a, &next_b);
            if paths.get(&state).is_none_or(|(cost, _, _)| total < *cost) {
                paths.insert(state, (total, next_a, next_b));
            }
        }
    }
    best
}

/// Coarse chord/control error is subtracted from the denominator, never added to coverage.
fn parameter_polygon(
    boundary: &Boundary,
    surface: &Surface,
    plane: Option<&crate::geom::Frame>,
) -> Option<(Loops, f64)> {
    let mut error = boundary.error;
    let mut loops = Vec::new();
    for lp in &boundary.loops {
        let mut points: Vec<(f64, f64)> = Vec::new();
        for trim in lp {
            match surface {
                Surface::Plane { .. } => {
                    let frame = plane?;
                    let project = |p| {
                        let q = geom::sub(p, frame.origin);
                        (geom::dot(q, frame.x), geom::dot(q, frame.y))
                    };
                    points.push(project(trim.start));
                    let n = if conic(&trim.curve).is_some() {
                        ((trim.b - trim.a).abs() / 0.1).ceil().max(1.0) as usize
                    } else {
                        1
                    };
                    if n > 128 {
                        return None;
                    }
                    for i in 0..=n {
                        points.push(project(
                            trim.curve
                                .value(trim.a + (trim.b - trim.a) * i as f64 / n as f64),
                        ));
                    }
                    points.push(project(trim.end));
                    if let Some((_, a, b)) = conic(&trim.curve) {
                        let step = (trim.b - trim.a).abs() / n as f64;
                        let d = (geom::norm(a) + geom::norm(b)) * step * step / 8.0;
                        error += tube(trim, surface, d)?;
                    }
                }
                Surface::Cylinder { frame, radius } => {
                    let uv = |p| {
                        let (u, v) = surface.parameters(p, None)?;
                        let u = points
                            .last()
                            .map_or(u, |q| geom::nearest_turn(u, q.0 / radius));
                        Some((u * radius, v))
                    };
                    let vertex = uv(trim.start)?;
                    points.push(vertex);
                    let start = trim.curve.value(trim.a);
                    let end = trim.curve.value(trim.b);
                    if representations::straight_origin(&trim.curve).is_some() {
                        let a = surface.parameters(start, None)?;
                        let b = surface.parameters(end, None)?;
                        let ideal_a = surface.value(a.0, a.1);
                        let ideal_b = surface.value(a.0, b.1);
                        error += tube(
                            trim,
                            surface,
                            geom::dist(start, ideal_a).max(geom::dist(end, ideal_b)),
                        )?;
                        let u = geom::nearest_turn(a.0, vertex.0 / radius) * radius;
                        points.extend([(u, a.1), (u, b.1)]);
                    } else {
                        let (o, a, b) = conic(&trim.curve)?;
                        let h = geom::dot(geom::sub(o, frame.origin), frame.z);
                        let origin = frame.to_world([0.0, 0.0, h]);
                        let x =
                            geom::unit(geom::sub(a, geom::scale(frame.z, geom::dot(a, frame.z))))?;
                        let sense = geom::dot(geom::cross(a, b), frame.z).signum();
                        if sense == 0.0 {
                            return None;
                        }
                        let y = geom::scale(geom::cross(frame.z, x), sense);
                        let d = geom::dist(o, origin)
                            + geom::dist(a, geom::scale(x, *radius))
                            + geom::dist(b, geom::scale(y, *radius));
                        error += tube(trim, surface, d)?;
                        let p = geom::add(
                            origin,
                            geom::scale(
                                geom::add(
                                    geom::scale(x, trim.a.cos()),
                                    geom::scale(y, trim.a.sin()),
                                ),
                                *radius,
                            ),
                        );
                        let u =
                            geom::nearest_turn(surface.parameters(p, None)?.0, vertex.0 / radius)
                                * radius;
                        points.extend([(u, h), (u + sense * (trim.b - trim.a) * radius, h)]);
                    }
                    let (u, v) = surface.parameters(trim.end, None)?;
                    let u = geom::nearest_turn(u, points.last()?.0 / radius) * radius;
                    points.push((u, v));
                }
                _ => return None,
            }
        }
        let first = *points.first()?;
        let last = *points.last()?;
        let extent = points
            .iter()
            .flat_map(|p| [p.0.abs(), p.1.abs()])
            .fold(1.0, f64::max);
        if (last.0 - first.0).hypot(last.1 - first.1) > 256.0 * f64::EPSILON * extent {
            return None;
        }
        points.push(first);
        loops.push(points);
    }
    Some((loops, error))
}

fn lower_area(boundary: &Boundary, patch: FaceRef<'_>) -> Option<f64> {
    let surface = &patch.0.faces[patch.1].surface;
    let plane = linear::canonical_frame(patch);
    let (mut loops, mut error) = parameter_polygon(boundary, surface, plane.as_ref())?;
    let count = loops.iter().map(Vec::len).sum::<usize>() as f64;
    let low = loops
        .iter()
        .flatten()
        .fold((f64::INFINITY, f64::INFINITY), |a, p| {
            (a.0.min(p.0), a.1.min(p.1))
        });
    let mut extent = 1.0_f64;
    for p in loops.iter_mut().flatten() {
        p.0 -= low.0;
        p.1 -= low.1;
        extent = extent.max(p.0.abs()).max(p.1.abs());
    }
    error += 256.0 * f64::EPSILON * count * count * extent * extent;
    let lower = linear::polygon_area(&loops)? - error;
    (lower.is_finite() && lower > 0.0).then_some(lower)
}

/// A finite cylindrical parameter cut is an isometry in (R*u,v). Each retained boundary
/// has its independent representation and projection error charged before polygon union.
/// Only single-loop regions narrower than a turn are eligible, so the cut is injective;
/// copied support turns establish coverage but never erase a hole or prove a rejection.
fn cylindrical_union(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<bool> {
    let surface @ Surface::Cylinder { radius, .. } = &patch.0.faces[patch.1].surface else {
        return None;
    };
    let period = std::f64::consts::TAU * radius;
    if !period.is_finite() || period <= 0.0 {
        return None;
    }
    let footprint = |face| {
        let b = boundary(face, surface)?;
        if b.loops.len() != 1 {
            return None;
        }
        let (mut loops, error) = parameter_polygon(&b, surface, None)?;
        for lp in &mut loops {
            lp.dedup();
        }
        let low = loops
            .iter()
            .flatten()
            .fold((f64::INFINITY, f64::INFINITY), |a, p| {
                (a.0.min(p.0), a.1.min(p.1))
            });
        let high = loops
            .iter()
            .flatten()
            .fold((f64::NEG_INFINITY, f64::NEG_INFINITY), |a, p| {
                (a.0.max(p.0), a.1.max(p.1))
            });
        (high.0 - low.0 < period * (1.0 - 256.0 * f64::EPSILON))
            .then_some((loops, error, low, high))
    };
    let (mut mine, mine_error, low, high) = footprint(patch)?;
    let centre = 0.5 * (low.0 + high.0);
    let mut theirs = Vec::new();
    let mut error = mine_error;
    for &support in supports {
        if !same_surface(surface, &support.0.faces[support.1].surface) {
            continue;
        }
        let Some((loops, bound, a, b)) = footprint(support) else {
            continue;
        };
        let turn = ((centre - 0.5 * (a.0 + b.0)) / period).round();
        if !turn.is_finite() || turn.abs() > 16.0 {
            continue;
        }
        for k in [-1.0, 0.0, 1.0] {
            let shift = (turn + k) * period;
            if b.0 + shift < low.0 || a.0 + shift > high.0 || b.1 < low.1 || a.1 > high.1 {
                continue;
            }
            let mut copy = loops.clone();
            for p in copy.iter_mut().flatten() {
                p.0 += shift;
            }
            theirs.push(copy);
            error += bound;
        }
    }
    let count = mine
        .iter()
        .chain(theirs.iter().flatten())
        .map(Vec::len)
        .sum::<usize>();
    if count > 512 {
        return None;
    }
    let mut extent = 1.0_f64;
    let mut magnitude = period.max(low.0.abs()).max(low.1.abs());
    for p in mine.iter_mut().chain(theirs.iter_mut().flatten()).flatten() {
        magnitude = magnitude.max(p.0.abs()).max(p.1.abs());
        p.0 -= low.0;
        p.1 -= low.1;
        extent = extent.max(p.0.abs()).max(p.1.abs());
    }
    let count = count as f64;
    let round =
        256.0 * f64::EPSILON * (count * count * extent * extent + count * extent * magnitude);
    let lower = linear::polygon_area(&mine)? - mine_error - round;
    let gap = linear::polygon_difference_area(&mine, &theirs)?;
    (lower > 0.0 && gap + error + round <= 1e-9 * lower).then_some(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{Face, Loop, Part, RawCurve2d};
    use crate::cover::{Coverage, assess};
    use crate::geom::Frame;

    fn frame(angle: f64) -> Frame {
        Frame {
            origin: [0.0; 3],
            x: [angle.cos(), angle.sin(), 0.0],
            y: [-angle.sin(), angle.cos(), 0.0],
            z: [0.0, 0.0, 1.0],
        }
    }

    fn part(surface: Surface, curves: Vec<(Curve, V3, V3, bool)>) -> Part {
        let edges: Vec<_> = curves
            .into_iter()
            .enumerate()
            .map(|(i, (curve, start, end, same_sense))| Edge {
                curve,
                start,
                end,
                vertices: if start == end { (i, i) } else { (i, i + 1) },
                same_sense,
                samples: vec![],
            })
            .collect();
        let face = Face {
            surface,
            reversed: false,
            loops: vec![Loop {
                edges: (0..edges.len()).map(|i| (i, true)).collect(),
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        };
        Part::new(vec![face], edges, vec![])
    }

    fn disk(radius: f64, angle: f64) -> Part {
        let frame = frame(angle);
        let curve = Curve::Circle { frame, radius };
        let vertex = curve.value(0.0);
        part(
            Surface::Plane { frame },
            vec![(curve, vertex, vertex, true)],
        )
    }

    #[test]
    fn complete_conics_do_not_depend_on_the_seam_or_parameter_axes() {
        let patch = disk(1.0, 0.0);
        for angle in [0.0, 0.37, 1.2, 3.7] {
            let support = disk(1.0 - 1e-12, angle);
            assert_eq!(assess((&patch, 0), &[(&support, 0)]), Coverage::Covered);
        }
        let smaller = disk(1.0 - 1e-4, 0.0);
        assert_ne!(assess((&patch, 0), &[(&smaller, 0)]), Coverage::Covered);
    }

    #[test]
    fn planar_correspondence_charges_in_plane_motion_and_keeps_3d_coherence() {
        let patch = disk(1.0, 0.0);
        for (axis, shift, expected) in [(2, 1e-8, true), (0, 1e-8, false), (2, 1e-5, false)] {
            let mut support = disk(1.0, 0.0);
            let Curve::Circle { frame, .. } = &mut support.edges[0].curve else {
                unreachable!()
            };
            frame.origin[axis] = shift;
            let poles = vec![
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
                [-1.0, 1.0],
                [-1.0, 0.0],
                [-1.0, -1.0],
                [0.0, -1.0],
                [1.0, -1.0],
                [1.0, 0.0],
            ];
            let support = support.with_file_pcurves(
                std::collections::BTreeMap::from([(
                    0,
                    vec![(
                        0,
                        RawCurve2d::BSpline {
                            degree: 2,
                            points: poles,
                            multiplicities: vec![3, 2, 2, 2, 3],
                            knots: vec![0.0, 0.25, 0.5, 0.75, 1.0],
                            weights: Some(
                                (0..9)
                                    .map(|i| {
                                        if i % 2 == 0 {
                                            1.0
                                        } else {
                                            std::f64::consts::FRAC_1_SQRT_2
                                        }
                                    })
                                    .collect(),
                            ),
                        },
                    )],
                )]),
                Default::default(),
                (1.0, 1.0),
                Some(1e-7),
            );
            assert_eq!(
                assess((&patch, 0), &[(&support, 0)]) == Coverage::Covered,
                expected,
                "axis={axis}, shift={shift}"
            );
        }
    }

    #[test]
    fn opposite_finite_arcs_are_not_the_same_boundary() {
        let a = [1.0, 0.0, 0.0];
        let b = [-1.0, 0.0, 0.0];
        let half = |sense| {
            part(
                Surface::Plane { frame: frame(0.0) },
                vec![
                    (
                        Curve::Circle {
                            frame: frame(0.0),
                            radius: 1.0,
                        },
                        a,
                        b,
                        sense,
                    ),
                    (
                        Curve::Line {
                            origin: b,
                            dir: [1.0, 0.0, 0.0],
                        },
                        b,
                        a,
                        true,
                    ),
                ],
            )
        };
        let upper = half(true);
        let lower = half(false);
        assert_eq!(assess((&upper, 0), &[(&upper, 0)]), Coverage::Covered);
        assert_ne!(assess((&upper, 0), &[(&lower, 0)]), Coverage::Covered);
    }

    #[test]
    fn an_auxiliary_or_unsupported_pcurve_cannot_be_omitted() {
        let patch = disk(1.0, 0.0);
        for curve in [
            RawCurve2d::Unresolved,
            RawCurve2d::Circle {
                point: [0.1, 0.0],
                reference: None,
                radius: 1.0,
            },
        ] {
            let p = disk(1.0, 0.0).with_file_pcurves(
                std::collections::BTreeMap::from([(0, vec![(0, curve)])]),
                Default::default(),
                (1.0, 1.0),
                Some(1e-7),
            );
            assert_ne!(assess((&patch, 0), &[(&p, 0)]), Coverage::Covered);
        }
    }

    fn cylinder(angle: f64, height: f64) -> Part {
        let radius = 2.0;
        let arc = |z| Curve::Circle {
            frame: Frame {
                origin: [0.0, 0.0, z],
                ..frame(0.0)
            },
            radius,
        };
        let a = [radius, 0.0, 0.0];
        let b = [0.0, radius, 0.0];
        let c = [0.0, radius, height];
        let d = [radius, 0.0, height];
        part(
            Surface::Cylinder {
                frame: Frame {
                    origin: [0.0, 0.0, 100.0],
                    ..frame(angle)
                },
                radius,
            },
            vec![
                (arc(0.0), a, b, true),
                (
                    Curve::Line {
                        origin: b,
                        dir: [0.0, 0.0, 1.0],
                    },
                    b,
                    c,
                    true,
                ),
                (arc(height), c, d, false),
                (
                    Curve::Line {
                        origin: d,
                        dir: [0.0, 0.0, -1.0],
                    },
                    d,
                    a,
                    true,
                ),
            ],
        )
    }

    #[test]
    fn cylindrical_strips_keep_their_finite_turn_and_axial_extent() {
        let patch = cylinder(0.0, 3.0);
        for angle in [0.0, 0.37, 1.2, 3.7] {
            let support = cylinder(angle, 3.0);
            assert_eq!(assess((&patch, 0), &[(&support, 0)]), Coverage::Covered);
        }
        let shorter = cylinder(0.0, 2.999);
        assert_ne!(assess((&patch, 0), &[(&shorter, 0)]), Coverage::Covered);
    }

    #[test]
    fn extra_holes_are_not_hidden_by_an_outer_boundary_match() {
        let patch = disk(1.0, 0.0);
        let mut ring = disk(1.0, 0.0);
        let hole = disk(0.25, 0.0).edges[0].clone();
        ring.edges.push(hole);
        ring.faces[0].loops.push(Loop {
            edges: vec![(1, false)],
            vertex: None,
        });
        let ring = Part::new(ring.faces.clone(), ring.edges.clone(), vec![]);
        assert_ne!(assess((&patch, 0), &[(&ring, 0)]), Coverage::Covered);
        assert_eq!(assess((&ring, 0), &[(&ring, 0)]), Coverage::Covered);
    }

    #[test]
    fn a_cylindrical_homotopy_cannot_approach_the_axis() {
        let curve = Curve::Circle {
            frame: frame(0.0),
            radius: 0.2,
        };
        let trim = Trim {
            start: curve.value(0.0),
            end: curve.value(0.0),
            curve,
            a: 0.0,
            b: std::f64::consts::TAU,
        };
        let surface = Surface::Cylinder {
            frame: frame(0.0),
            radius: 2.0,
        };
        assert!(tube(&trim, &surface, 1e-12).is_none());
    }
}
