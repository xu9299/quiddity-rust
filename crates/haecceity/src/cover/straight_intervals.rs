//! Finite straight images and affine pcurve construction enclosed before rounding.

use super::interval::{self, Interval as I};
use crate::brep::RawCurve2d;
use crate::geom::{Curve, Frame, Surface, V3};

pub(super) type Point = [I; 3];
pub(super) type Ends = [Point; 2];

pub(super) fn point(p: V3) -> Point {
    p.map(I::point)
}

pub(super) fn norm_upper(p: Point) -> f64 {
    let sum = p.into_iter().fold(I::point(0.0), |s, v| {
        let a = I::point(v.abs_upper());
        s + a * a
    });
    sum.hi.sqrt().next_up()
}

pub(super) fn world(frame: &Frame, uv: [I; 2]) -> Point {
    std::array::from_fn(|k| {
        I::point(frame.origin[k]) + I::point(frame.x[k]) * uv[0] + I::point(frame.y[k]) * uv[1]
    })
}

pub(super) fn project(frame: &Frame, p: Point) -> [I; 2] {
    let q = std::array::from_fn(|k| p[k] - I::point(frame.origin[k]));
    [
        interval::dot(q, point(frame.x)),
        interval::dot(q, point(frame.y)),
    ]
}

fn foot(origin: Point, direction: Point, p: V3, bounded: bool) -> Option<Point> {
    let denominator = interval::dot(direction, direction);
    if denominator.lo <= 0.0 || !denominator.hi.is_finite() {
        return None;
    }
    let delta = std::array::from_fn(|k| I::point(p[k]) - origin[k]);
    let mut t = interval::dot(delta, direction) / denominator;
    if bounded {
        t = I {
            lo: t.lo.clamp(0.0, 1.0),
            hi: t.hi.clamp(0.0, 1.0),
        };
    }
    let result = std::array::from_fn(|k| origin[k] + direction[k] * t);
    result
        .iter()
        .all(|v| v.lo.is_finite() && v.hi.is_finite())
        .then_some(result)
}

fn ends(origin: Point, direction: Point, vertices: [V3; 2], bounded: bool) -> Option<Ends> {
    Some([
        foot(origin, direction, vertices[0], bounded)?,
        foot(origin, direction, vertices[1], bounded)?,
    ])
}

pub(super) fn curve_ends(curve: &Curve, vertices: [V3; 2]) -> Option<Ends> {
    match curve {
        Curve::Line { origin, dir } => ends(point(*origin), point(*dir), vertices, false),
        Curve::Nurbs(c)
            if c.degree() == 1
                && c.control_points().len() == 2
                && c.weights().iter().all(|w| w.is_finite() && *w > 0.0) =>
        {
            let a = point(c.control_points()[0]);
            let b = point(c.control_points()[1]);
            ends(a, std::array::from_fn(|k| b[k] - a[k]), vertices, true)
        }
        _ => None,
    }
}

pub(super) fn pcurve_ends(
    raw: &RawCurve2d,
    surface: &Surface,
    (mm, _): (f64, f64),
    vertices: [V3; 2],
) -> Option<Ends> {
    let Surface::Plane { frame } = surface else {
        return None;
    };
    let mapped = |p: [f64; 2]| world(frame, p.map(|x| I::point(x) * I::point(mm)));
    match raw {
        RawCurve2d::Line {
            point,
            direction,
            magnitude,
        } => {
            let uv = direction.map(|v| I::point(v) * I::point(*magnitude) * I::point(mm));
            let dir = std::array::from_fn(|k| {
                I::point(frame.x[k]) * uv[0] + I::point(frame.y[k]) * uv[1]
            });
            ends(mapped(*point), dir, vertices, false)
        }
        RawCurve2d::BSpline {
            degree: 1,
            points,
            multiplicities,
            knots,
            weights,
        } if points.len() == 2
            && multiplicities == &[2, 2]
            && knots.len() == 2
            && knots.iter().all(|x| x.is_finite())
            && knots[0] < knots[1]
            && weights
                .as_ref()
                .is_none_or(|w| w.len() == 2 && w.iter().all(|v| v.is_finite() && *v > 0.0)) =>
        {
            let (a, b) = (mapped(points[0]), mapped(points[1]));
            ends(a, std::array::from_fn(|k| b[k] - a[k]), vertices, true)
        }
        _ => None,
    }
}
