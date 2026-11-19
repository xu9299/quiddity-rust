//! Circular-rim evidence from the selected boundary of a conical face.
use super::{EdgeTrim, Part};
use crate::geom::{self, Curve, Frame, Surface, V3};
use crate::nurbs::NurbsCurve;

impl Part {
    /// Radius and centre of a conical face's circular rim. A selected surface
    /// curve can provide this evidence when its entire locus agrees with a
    /// circle at coordinate precision. The measured boundary retains its own
    /// curve and parameterisation.
    pub fn cone_rim_circle(&self, face: usize, edge: usize) -> Option<(f64, V3)> {
        let Surface::Cone {
            frame,
            radius,
            semi_angle,
        } = &self.faces.get(face)?.surface
        else {
            return None;
        };
        if !self.edge_faces().get(edge)?.contains(&face) {
            return None;
        }
        let Some((surface, trim)) = self.spatial_trim(edge) else {
            return exact_circle(&self.edges[edge].curve);
        };
        let curve = trim.spatial_curve(surface);
        if let Some(circle) = curve.and_then(exact_circle) {
            return Some(circle);
        }
        if let Some(Curve::Nurbs(spline)) = curve {
            let height = frame.to_local(trim.point(surface, trim.interval.0))[2];
            let r = (radius + height * semi_angle.tan()).abs();
            if quadratic_rim(spline, frame, height, r) {
                return Some((r, frame.to_world([0., 0., height])));
            }
        }
        trim.parallel_circle(surface)
    }
}

fn exact_circle(curve: &Curve) -> Option<(f64, V3)> {
    match curve {
        Curve::Circle { frame, radius } => Some((*radius, frame.origin)),
        _ => None,
    }
}

impl EdgeTrim {
    /// For a linear pcurve on a cone or cylinder, v is physical distance on
    /// the generator. At every angle its distance to the middle parallel is
    /// at most half the endpoint difference in v.
    fn parallel_circle(&self, face: &Surface) -> Option<(f64, V3)> {
        self.curve.line()?;
        let carrier = self.curve.transfer.as_ref().map_or(face, |tr| &tr.carrier);
        let a = self.curve.carrier_uv(self.interval.0);
        let b = self.curve.carrier_uv(self.interval.1);
        let v = 0.5 * (a[1] + b[1]);
        let (frame, height, r) = match carrier {
            Surface::Cylinder { frame, radius } => (frame, v, *radius),
            Surface::Cone {
                frame,
                radius,
                semi_angle,
            } => (frame, v * semi_angle.cos(), radius + v * semi_angle.sin()),
            _ => return None,
        };
        let roundoff =
            16. * f64::EPSILON * (geom::norm(frame.origin) + r.abs() + a[1].abs() + b[1].abs());
        let error = 0.5 * (a[1] - b[1]).abs() + roundoff;
        let arc = (r * (a[0] - b[0])).abs();
        (error <= geom::COORD_FLOOR
            && r.abs() > geom::COORD_FLOOR
            && arc.is_finite()
            && arc > geom::COORD_FLOOR)
            .then(|| (r.abs(), frame.to_world([0., 0., height])))
    }
}

/// Bound the radial and axial errors of every point of a planar rational
/// quadratic. The cone supplies the centre and radius independently of the
/// edge's optional spatial representation. Interior knots of multiplicity
/// two separate quadratic Bezier spans. Other splines have no such proof.
fn quadratic_rim(curve: &NurbsCurve, frame: &Frame, height: f64, radius: f64) -> bool {
    if curve.degree() != 2 || !radius.is_finite() || radius <= geom::COORD_FLOOR {
        return false;
    }
    let knots = curve.knots();
    let n = curve.control_points().len();
    // Each active span [knots[i+2], knots[i+3]] uses poles [i,i+1,i+2].
    // Its repeated end knots give the Bezier basis, also at the domain ends.
    if n.is_multiple_of(2)
        || !(1..n)
            .step_by(2)
            .all(|i| knots[i] == knots[i + 1] && knots[i + 1] < knots[i + 2])
        || knots[n] != knots[n + 1]
    {
        return false;
    }
    let mut curved = false;
    for i in (0..n - 2).step_by(2) {
        let p = [i, i + 1, i + 2].map(|j| frame.to_local(curve.control_points()[j]));
        let weights = [i, i + 1, i + 2].map(|j| curve.weights()[j]);
        let max_w = weights[0].max(weights[1]).max(weights[2]);
        let w = weights.map(|v| v / max_w);
        let z = p.iter().map(|q| (q[2] - height).abs()).fold(0., f64::max);
        // Radial squared error has degree four after clearing the positive
        // rational denominator. Bernstein products contribute the binomial
        // factors 1, 1, (1/3,2/3), 1, 1 to the five coefficients.
        let residual = |a: usize, b: usize| {
            w[a] * w[b] * (p[a][0] * p[b][0] + p[a][1] * p[b][1] - radius * radius)
        };
        let h = [
            residual(0, 0),
            residual(0, 1),
            (residual(0, 2) + 2. * residual(1, 1)) / 3.,
            residual(1, 2),
            residual(2, 2),
        ];
        let min_w = w[0].min(w[1]).min(w[2]);
        let radial = h.into_iter().map(f64::abs).fold(0., f64::max) / (min_w * min_w * radius);
        // Include arithmetic roundoff in the coefficient bound. Extremely
        // ill-conditioned weights or large coordinates cannot certify a rim.
        let scale = p.iter().map(|q| q[0].hypot(q[1])).fold(radius, f64::max);
        let weight_ratio = 1. / min_w;
        let roundoff = 64.
            * f64::EPSILON
            * weight_ratio.powi(2)
            * (scale * scale / radius
                + height.abs()
                + geom::norm(frame.origin) * (scale / radius).max(1.));
        let error = z.hypot(radial) + roundoff;
        if !error.is_finite() || error > geom::COORD_FLOOR {
            return false;
        }
        let mid = curve.value(0.5 * (knots[i + 2] + knots[i + 3]));
        let ends = &curve.control_points()[i..=i + 2];
        let chord = geom::sub(ends[2], ends[0]);
        // An unresolved straight or point-like span cannot establish curvature.
        curved |= geom::norm(geom::cross(geom::sub(mid, ends[0]), chord))
            > geom::COORD_FLOOR * geom::norm(chord);
    }
    curved
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_rational_rim_is_checked_between_its_endpoints() {
        let axes = Frame {
            origin: [0.; 3],
            x: [1., 0., 0.],
            y: [0., 1., 0.],
            z: [0., 0., 1.],
        };
        let arc = |middle: V3, weight| {
            NurbsCurve::new(
                2,
                vec![[1., 0., 0.], middle, [0., 1., 0.]],
                vec![1., weight, 1.],
                vec![0., 0., 0., 1., 1., 1.],
            )
            .unwrap()
        };
        let w = std::f64::consts::FRAC_1_SQRT_2;
        assert!(quadratic_rim(&arc([1., 1., 0.], w), &axes, 0., 1.));
        assert!(!quadratic_rim(&arc([1., 1., 0.], 1.), &axes, 0., 1.));
        assert!(!quadratic_rim(&arc([1., 1., 1e-4], w), &axes, 0., 1.));
    }
}
