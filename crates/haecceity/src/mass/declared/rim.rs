//! Circular supports of selected edges incident on cylinders or cones.
use super::{EdgeTrim, Part};
use crate::geom::{self, Curve, Frame, Surface};
use crate::nurbs::NurbsCurve;

impl Part {
    /// The spatial support used to recognise geometric features. A selected
    /// pcurve can prove a circular rim at coordinate precision on an incident
    /// cylinder or cone. Use `boundary_edge_curve` and `boundary_edge_interval`
    /// for evaluation in the selected curve's own parameterisation.
    pub fn boundary_edge_support(&self, edge: usize) -> Option<&Curve> {
        let Some((carrier, trim)) = self.spatial_trim(edge) else {
            return Some(&self.edges[edge].curve);
        };
        let selected = trim.spatial_curve(carrier);
        if matches!(
            selected,
            Some(Curve::Line { .. } | Curve::Circle { .. } | Curve::Ellipse { .. })
        ) {
            return selected;
        }
        trim.rim
            .get_or_init(|| {
                let Some(Curve::Nurbs(spline)) = selected else {
                    return trim.parallel_circle(carrier);
                };
                let point = trim.point(carrier, trim.interval.0);
                self.edge_faces()[edge].iter().find_map(|&f| {
                    let (frame, radius, slope) = match &self.faces[f].surface {
                        Surface::Cylinder { frame, radius } => (frame, radius, 0.),
                        Surface::Cone {
                            frame,
                            radius,
                            semi_angle,
                        } => (frame, radius, semi_angle.tan()),
                        _ => return None,
                    };
                    let height = frame.to_local(point)[2];
                    let r = (radius + height * slope).abs();
                    quadratic_rim(spline, frame, height, r).then(|| circle(frame, height, r))
                })
            })
            .as_ref()
            .or(selected)
    }
}

fn circle(frame: &Frame, height: f64, radius: f64) -> Curve {
    Curve::Circle {
        frame: Frame {
            origin: frame.to_world([0., 0., height]),
            z: geom::cross(frame.x, frame.y),
            ..*frame
        },
        radius,
    }
}

impl EdgeTrim {
    /// For a linear pcurve on a cone or cylinder, v is physical distance on
    /// the generator. At every angle its distance to the middle parallel is
    /// at most half the endpoint difference in v.
    fn parallel_circle(&self, face: &Surface) -> Option<Curve> {
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
            && arc > geom::COORD_FLOOR
            && (a[0] - b[0]).abs() <= std::f64::consts::TAU + 64. * f64::EPSILON)
            .then(|| circle(frame, height, r.abs()))
    }
}

/// Bound the radial and axial errors of every point of a planar rational
/// quadratic. An incident cylinder or cone supplies the centre and radius
/// independently of the optional spatial curve. Interior knots of multiplicity
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
    let mut sweep = 0.0_f64;
    let mut curved = false;
    for i in (0..n - 2).step_by(2) {
        let p = [i, i + 1, i + 2].map(|j| frame.to_local(curve.control_points()[j]));
        // The angular speed has a positive denominator. Its quadratic
        // Bernstein coefficients have the signs of these three cross products.
        // Each span and the complete rim must advance in one direction.
        let cross = |a: usize, b: usize| p[a][0] * p[b][1] - p[a][1] * p[b][0];
        let scale = p.iter().map(|q| q[0].hypot(q[1])).fold(radius, f64::max);
        let angular_roundoff = 64. * f64::EPSILON * scale * (scale + geom::norm(frame.origin));
        let chord = cross(0, 2);
        let direction = chord.signum();
        let turn = chord.atan2(p[0][0] * p[2][0] + p[0][1] * p[2][1]);
        if !(chord.abs() > angular_roundoff
            && cross(0, 1) * direction > angular_roundoff
            && cross(1, 2) * direction > angular_roundoff)
            || i > 0 && direction != sweep.signum()
        {
            return false;
        }
        sweep += turn;
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
    curved && sweep.abs() <= std::f64::consts::TAU + 64. * f64::EPSILON * n as f64
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
        let arc = |middle: geom::V3, weight| {
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
        // An arc followed back to its start lies on the same circle throughout.
        let retraced = NurbsCurve::new(
            2,
            vec![
                [1., 0., 0.],
                [1., 1., 0.],
                [0., 1., 0.],
                [1., 1., 0.],
                [1., 0., 0.],
            ],
            vec![1., w, 1., w, 1.],
            vec![0., 0., 0., 1., 1., 2., 2., 2.],
        )
        .unwrap();
        assert!(!quadratic_rim(&retraced, &axes, 0., 1.));
        let mut poles = vec![[1., 0., 0.]];
        let mut weights = vec![1.];
        let mut knots = vec![0., 0., 0.];
        for span in 1..=8 {
            let [x, y, _] = *poles.last().unwrap();
            poles.extend([[x - y, x + y, 0.], [-y, x, 0.]]);
            weights.extend([w, 1.]);
            knots.extend([f64::from(span); 2]);
        }
        knots.push(8.);
        let two_turns = NurbsCurve::new(2, poles, weights, knots).unwrap();
        assert!(!quadratic_rim(&two_turns, &axes, 0., 1.));
    }
}
