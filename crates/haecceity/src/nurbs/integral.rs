//! Integral-local tensor contraction and centered span polynomials.
//!
//! With u varying, form Q_i = sum_j M_j(v) H_ij and R_i = sum_j M'_j(v) H_ij
//! from homogeneous controls H. Then H = sum_i N_i Q_i, H_u = sum_i N'_i Q_i,
//! and H_v = sum_i N_i R_i. The v-varying case is symmetric. Both partials are
//! required for the area density; an isocurve's tangent alone is insufficient.
//! The fixed-parameter bits and both spans key each contraction.
//!
//! Basis coefficients use x = (t - centre) / half within each knot interval.
//! Horner evaluates the basis and its derivative in x, then scales by 1 / half.
//! This changes floating-point association, not the geometry or quadrature rule.
//! See `docs/nurbs-integration.md` for limits and numerical validation requirements.

use super::{NurbsSurface, SurfaceEvaluation, V3, find_span};

#[cfg(test)]
mod tests;

/// A single integral owns both the specialized evaluator and its original fallback.
/// No prepared state is shared between surfaces, integral directions, or threads.
pub(crate) struct IntegralEvaluation<'a> {
    surface: &'a NurbsSurface,
    original: Option<SurfaceEvaluation<'a>>,
    prepared: Option<PreparedSurface<'a>>,
}

impl<'a> IntegralEvaluation<'a> {
    pub(crate) fn new(surface: &'a NurbsSurface, along_v: bool) -> Self {
        Self {
            surface,
            original: None,
            prepared: PreparedSurface::new(surface, along_v),
        }
    }

    pub(crate) fn value_and_partials(&mut self, u: f64, v: f64) -> (V3, V3, V3) {
        if let Some(value) = self.prepared.as_mut().and_then(|p| p.evaluate(u, v)) {
            value
        } else {
            self.original
                .get_or_insert_with(|| SurfaceEvaluation::new(self.surface))
                .value_and_partials(u, v)
        }
    }
}

#[derive(Clone, Copy)]
struct Basis {
    values: [f64; 8],
    len: usize,
}
impl std::ops::Deref for Basis {
    type Target = [f64];
    fn deref(&self) -> &[f64] {
        &self.values[..self.len]
    }
}
fn basis_and_derivative(knots: &[f64], degree: usize, span: usize, t: f64) -> (Basis, Basis) {
    debug_assert!(degree < SIZE);
    let (b, d) = super::basis_and_derivative(knots, degree, span, t);
    (
        Basis {
            values: std::array::from_fn(|i| b.values[i]),
            len: b.len,
        },
        Basis {
            values: std::array::from_fn(|i| d.values[i]),
            len: d.len,
        },
    )
}

const SIZE: usize = 8;
const AXIS_BYTES: usize = 2 * 1024 * 1024;
// Near-zero polynomial basis errors can be amplified by unequal rational weights.
// Keep a conservative active-net contrast limit; this is not a uniform error bound.
const MAX_WEIGHT_RATIO: f64 = 256.0;

#[derive(Clone)]
struct Polynomial {
    centre: f64,
    half: f64,
    coefficients: [[f64; SIZE]; SIZE],
}

impl Polynomial {
    fn new(knots: &[f64], degree: usize, span: usize) -> Option<Self> {
        if degree >= SIZE {
            return None;
        }
        let (lo, hi) = (knots[span], knots[span + 1]);
        let half = 0.5 * (hi - lo);
        let centre = lo + half;
        if !half.is_finite()
            || half <= 1e-150
            || half <= 32.0 * f64::EPSILON * lo.abs().max(hi.abs())
        {
            return None;
        }
        let mut n = [[0.0; SIZE]; SIZE];
        n[0][0] = 1.0;
        for j in 1..=degree {
            let mut saved = [0.0; SIZE];
            for r in 0..j {
                let right = knots[span + r + 1] - centre;
                let left = centre - knots[span + 1 - j + r];
                let denominator = knots[span + r + 1] - knots[span + 1 - j + r];
                if !denominator.is_finite() || denominator.abs() < 1e-300 {
                    return None;
                }
                let temp = n[r].map(|v| v / denominator);
                for k in 0..=j {
                    let shifted = if k > 0 { half * temp[k - 1] } else { 0.0 };
                    n[r][k] = saved[k] + right * temp[k] - shifted;
                    saved[k] = left * temp[k] + shifted;
                }
            }
            n[j] = saved;
        }
        n.iter().flatten().all(|v| v.is_finite()).then_some(Self {
            centre,
            half,
            coefficients: n,
        })
    }

    fn evaluate(&self, degree: usize, t: f64) -> (Basis, Basis) {
        let x = (t - self.centre) / self.half;
        let mut basis = Basis {
            values: [0.0; 8],
            len: degree + 1,
        };
        let mut derivative = basis;
        for i in 0..=degree {
            let coefficients = &self.coefficients[i];
            let (mut y, mut dy) = (coefficients[degree], 0.0);
            for k in (0..degree).rev() {
                dy = dy * x + y;
                y = y * x + coefficients[k];
            }
            basis.values[i] = y;
            derivative.values[i] = dy / self.half;
        }
        (basis, derivative)
    }
}

struct Entry {
    parameter: u64,
    span: usize,
    basis: Basis,
    derivative: Basis,
}

struct Axis {
    last: Option<Entry>,
    polynomials: Option<Vec<Option<Polynomial>>>,
}

impl Axis {
    fn new(knots: &[f64], degree: usize, n: usize) -> Self {
        let polynomials = (degree < SIZE
            && n <= AXIS_BYTES / std::mem::size_of::<Option<Polynomial>>())
        .then(|| {
            (0..n)
                .map(|span| {
                    (span >= degree)
                        .then(|| Polynomial::new(knots, degree, span))
                        .flatten()
                })
                .collect()
        });
        Self {
            last: None,
            polynomials,
        }
    }

    fn prepare(&mut self, t: f64, knots: &[f64], degree: usize, n: usize) -> Option<&Entry> {
        if self
            .last
            .as_ref()
            .is_none_or(|e| e.parameter != t.to_bits())
        {
            let span = find_span(knots, degree, n, t);
            let (basis, derivative) =
                if let Some(p) = self.polynomials.as_ref().and_then(|p| p[span].as_ref()) {
                    p.evaluate(degree, t)
                } else {
                    basis_and_derivative(knots, degree, span, t)
                };
            self.last = Some(Entry {
                parameter: t.to_bits(),
                span,
                basis,
                derivative,
            });
        }
        self.last.as_ref()
    }
}

struct PreparedSurface<'a> {
    surface: &'a NurbsSurface,
    u: Axis,
    v: Axis,
    along_v: bool,
    spans: Option<(usize, usize)>,
    fixed: Option<(u64, usize, usize)>,
    controls: Vec<[f64; 4]>,
    q: Vec<[f64; 4]>,
    r: Vec<[f64; 4]>,
    safe_controls: bool,
}

impl<'a> PreparedSurface<'a> {
    fn new(surface: &'a NurbsSurface, along_v: bool) -> Option<Self> {
        if surface.degree_u >= 8 || surface.degree_v >= 8 {
            return None;
        }
        let u = Axis::new(
            &surface.knots_u,
            surface.degree_u,
            surface.control_points.len(),
        );
        let v = Axis::new(
            &surface.knots_v,
            surface.degree_v,
            surface.control_points[0].len(),
        );
        Some(Self {
            surface,
            u,
            v,
            along_v,
            spans: None,
            fixed: None,
            controls: Vec::new(),
            q: Vec::new(),
            r: Vec::new(),
            safe_controls: true,
        })
    }

    fn evaluate(&mut self, u: f64, v: f64) -> Option<(V3, V3, V3)> {
        let s = self.surface;
        let (u0, u1, v0, v1) = s.domain();
        let (u, v) = (u.clamp(u0, u1), v.clamp(v0, v1));
        if !(u.is_finite() && v.is_finite()) {
            return None;
        }
        let pu = self
            .u
            .prepare(u, &s.knots_u, s.degree_u, s.control_points.len());
        let pv = self
            .v
            .prepare(v, &s.knots_v, s.degree_v, s.control_points[0].len());
        let (Some(pu), Some(pv)) = (pu, pv) else {
            return None;
        };
        if self.spans != Some((pu.span, pv.span)) {
            self.controls.clear();
            let (mut min_weight, mut max_weight) = (f64::INFINITY, 0.0_f64);
            for ci in pu.span - s.degree_u..=pu.span {
                for cj in pv.span - s.degree_v..=pv.span {
                    let c = s.control_points[ci][cj];
                    let w = s.weights[ci][cj];
                    min_weight = min_weight.min(w);
                    max_weight = max_weight.max(w);
                    self.controls.push([c[0] * w, c[1] * w, c[2] * w, w]);
                }
            }
            self.safe_controls = min_weight > 0.0
                && max_weight / min_weight <= MAX_WEIGHT_RATIO
                && self.controls.iter().flatten().all(|x| x.is_finite());
            self.spans = Some((pu.span, pv.span));
        }
        if !self.safe_controls {
            return None;
        }
        let (bu, du, bv, dv) = (&pu.basis, &pu.derivative, &pv.basis, &pv.derivative);
        let mut a = [0.0; 4];
        let (au, av);
        {
            let key = (
                if self.along_v {
                    pu.parameter
                } else {
                    pv.parameter
                },
                pu.span,
                pv.span,
            );
            if self.fixed != Some(key) {
                let len = if self.along_v { bv.len() } else { bu.len() };
                self.q.clear();
                self.q.resize(len, [0.0; 4]);
                self.r.clear();
                self.r.resize(len, [0.0; 4]);
                for (i, row) in self.controls.chunks_exact(bv.len()).enumerate() {
                    for (j, h) in row.iter().enumerate() {
                        let (slot, b, d) = if self.along_v {
                            (j, bu[i], du[i])
                        } else {
                            (i, bv[j], dv[j])
                        };
                        for (k, &h) in h.iter().enumerate() {
                            self.q[slot][k] += b * h;
                            self.r[slot][k] += d * h;
                        }
                    }
                }
                self.fixed = Some(key);
            }

            let (b, d) = if self.along_v { (bv, dv) } else { (bu, du) };
            let (mut varying, mut fixed) = ([0.0; 4], [0.0; 4]);
            for i in 0..b.len() {
                for k in 0..4 {
                    a[k] += b[i] * self.q[i][k];
                    varying[k] += d[i] * self.q[i][k];
                    fixed[k] += b[i] * self.r[i][k];
                }
            }
            (au, av) = if self.along_v {
                (fixed, varying)
            } else {
                (varying, fixed)
            };
        }
        let point = [a[0] / a[3], a[1] / a[3], a[2] / a[3]];
        let derive = |d: [f64; 4]| [0, 1, 2].map(|k| (d[k] - d[3] * point[k]) / a[3]);
        let result = (point, derive(au), derive(av));
        if [result.0, result.1, result.2]
            .iter()
            .flatten()
            .any(|x| !x.is_finite())
        {
            return None;
        }
        Some(result)
    }
}
