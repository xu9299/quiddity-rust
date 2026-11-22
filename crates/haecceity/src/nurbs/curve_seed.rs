//! Reuse inversion-table coefficients within one part, never geometric points.
use super::{NurbsCurve, V3, basis, find_span};
use std::collections::{HashMap, hash_map::Entry};
use std::sync::{Arc, Mutex, OnceLock};

type SharedSeed = Arc<OnceLock<Option<CurveSeed>>>;

/// Each curve registers only when its own inversion table is first requested.
/// A unique coefficient shape retains the ordinary table builder. The second
/// curve with the same units, degree, weights and knots initializes shared rows.
/// No curve, carrier, topology or point sample is retained here.
#[derive(Debug, Default)]
pub(crate) struct CurveInversionPool(Mutex<HashMap<Vec<u64>, SharedSeed>>);

impl CurveInversionPool {
    fn request(&self, curve: &NurbsCurve, units: [u64; 2]) -> Option<SharedSeed> {
        let (lo, hi) = curve.domain();
        if curve.degree > 3
            || ![lo, hi, hi - lo].iter().all(|x| x.is_finite())
            || !curve
                .control_points
                .iter()
                .flatten()
                .all(|x| x.abs() < 1e50)
        {
            return None;
        }
        let mut key = Vec::with_capacity(4 + curve.weights.len() + curve.knots.len());
        key.extend(units);
        key.extend([curve.degree as u64, curve.control_points.len() as u64]);
        key.extend(
            curve
                .weights
                .iter()
                .chain(&curve.knots)
                .map(|x| x.to_bits()),
        );
        let mut pool = self.0.lock().ok()?;
        match pool.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(Arc::default());
                None
            }
            Entry::Occupied(entry) => Some(entry.get().clone()),
        }
    }

    pub(super) fn samples(&self, curve: &NurbsCurve, units: [u64; 2]) -> Option<Vec<(f64, V3)>> {
        let shared = self.request(curve, units)?;
        // Neither coefficient construction nor curve evaluation holds the map lock.
        let seed = shared
            .get_or_init(|| CurveSeed::new_batched(curve))
            .as_ref()?;
        match curve.degree {
            0 => Some(seed.samples::<1>(&curve.control_points)),
            1 => Some(seed.samples::<2>(&curve.control_points)),
            2 => Some(seed.samples::<3>(&curve.control_points)),
            3 => Some(seed.samples::<4>(&curve.control_points)),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct CurveSeedRow {
    t: f64,
    first: usize,
    weights: [f64; 4],
    denominator: f64,
}

#[derive(Debug)]
struct CurveSeed(Vec<CurveSeedRow>);

impl CurveSeed {
    fn samples<const N: usize>(&self, points: &[V3]) -> Vec<(f64, V3)> {
        let mut out = Vec::with_capacity(self.0.len());
        let (chunks, remainder) = self.0.as_chunks::<4>();
        for rows in chunks {
            let mut acc = [[0.; 4]; 3];
            for j in 0..N {
                for (lane, row) in rows.iter().enumerate() {
                    let point = points[row.first + j];
                    let w = row.weights[j];
                    acc[0][lane] += point[0] * w;
                    acc[1][lane] += point[1] * w;
                    acc[2][lane] += point[2] * w;
                }
            }
            for (lane, row) in rows.iter().enumerate() {
                out.push((
                    row.t,
                    [
                        acc[0][lane] / row.denominator,
                        acc[1][lane] / row.denominator,
                        acc[2][lane] / row.denominator,
                    ],
                ));
            }
        }
        for row in remainder {
            let mut acc = [0.; 3];
            for j in 0..N {
                let point = points[row.first + j];
                let w = row.weights[j];
                acc[0] += point[0] * w;
                acc[1] += point[1] * w;
                acc[2] += point[2] * w;
            }
            out.push((row.t, acc.map(|x| x / row.denominator)));
        }
        out
    }
}
impl CurveSeed {
    fn four<const N: usize>(
        curve: &NurbsCurve,
        span: usize,
        ts: [f64; 4],
    ) -> ([[f64; 4]; 4], [f64; 4]) {
        let mut basis = [[0.; 4]; N];
        let mut left = [[0.; 4]; N];
        let mut right = [[0.; 4]; N];
        basis[0] = [1.; 4];
        for j in 1..N {
            for (lane, &t) in ts.iter().enumerate() {
                left[j][lane] = t - curve.knots[span + 1 - j];
                right[j][lane] = curve.knots[span + j] - t;
            }
            let mut saved = [0.; 4];
            for k in 0..j {
                for (lane, saved) in saved.iter_mut().enumerate() {
                    let denominator = right[k + 1][lane] + left[j - k][lane];
                    let temp = if denominator.abs() < 1e-300 {
                        0.
                    } else {
                        basis[k][lane] / denominator
                    };
                    basis[k][lane] = *saved + right[k + 1][lane] * temp;
                    *saved = left[j - k][lane] * temp;
                }
            }
            basis[j] = saved;
        }

        let mut weighted = [[0.; 4]; 4];
        let mut denominators = [0.; 4];
        for (j, values) in basis.iter().enumerate() {
            let i = span + 1 - N + j;
            for (lane, &value) in values.iter().enumerate() {
                let w = curve.weights[i] * value;
                weighted[j][lane] = w;
                denominators[lane] += w;
            }
        }
        (weighted, denominators)
    }

    fn checked_row(
        t: f64,
        first: usize,
        weights: [f64; 4],
        denominator: f64,
    ) -> Option<CurveSeedRow> {
        (weights.iter().all(|w| w.abs() < 1e50)
            && denominator.abs() > 1e-50
            && denominator.abs() < 1e50)
            .then_some(CurveSeedRow {
                t,
                first,
                weights,
                denominator,
            })
    }
    fn scalar_row(curve: &NurbsCurve, t: f64) -> Option<CurveSeedRow> {
        let (lo, hi) = curve.domain();
        let clipped = t.clamp(lo, hi);
        let span = find_span(
            &curve.knots,
            curve.degree,
            curve.control_points.len(),
            clipped,
        );
        let basis = basis(&curve.knots, curve.degree, span, clipped);
        let first = span - curve.degree;
        let mut weights = [0.; 4];
        let mut denominator = 0.;
        for (j, &b) in basis.iter().enumerate() {
            weights[j] = curve.weights[first + j] * b;
            denominator += weights[j];
        }
        Self::checked_row(t, first, weights, denominator)
    }
    fn new_batched(curve: &NurbsCurve) -> Option<Self> {
        let (lo, hi) = curve.domain();
        let t = |i| lo + (hi - lo) * i as f64 / 512_f64;
        let mut rows = Vec::with_capacity(513);
        for start in (0..512).step_by(4) {
            let ts: [f64; 4] = std::array::from_fn(|lane| t(start + lane));
            let clipped = ts.map(|t| t.clamp(lo, hi));
            let spans = clipped
                .map(|t| find_span(&curve.knots, curve.degree, curve.control_points.len(), t));
            if spans.iter().all(|&s| s == spans[0]) {
                let (weights, denominators) = match curve.degree {
                    0 => Self::four::<1>(curve, spans[0], clipped),
                    1 => Self::four::<2>(curve, spans[0], clipped),
                    2 => Self::four::<3>(curve, spans[0], clipped),
                    _ => Self::four::<4>(curve, spans[0], clipped),
                };
                for lane in 0..4 {
                    rows.push(Self::checked_row(
                        ts[lane],
                        spans[0] - curve.degree,
                        std::array::from_fn(|j| weights[j][lane]),
                        denominators[lane],
                    )?);
                }
            } else {
                for t in ts {
                    rows.push(Self::scalar_row(curve, t)?);
                }
            }
        }
        rows.push(Self::scalar_row(curve, t(512))?);
        Some(Self(rows))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve(shift: f64) -> NurbsCurve {
        NurbsCurve::new(
            2,
            vec![[shift, 0., 0.], [shift + 1., 2., 0.], [shift + 3., 0., 0.]],
            vec![1., 0.7, 1.],
            vec![0., 0., 0., 1., 1., 1.],
        )
        .unwrap()
    }

    fn check(curve: &NurbsCurve, pool: &CurveInversionPool, units: [u64; 2]) {
        curve.prepare_inversion(pool, units);
        let (lo, hi) = curve.domain();
        let table = curve.table.get().unwrap();
        assert_eq!(table.len(), 513);
        for (i, &(t, point)) in table.iter().enumerate() {
            let expected = lo + (hi - lo) * i as f64 / 512_f64;
            assert_eq!(t.to_bits(), expected.to_bits());
            assert_eq!(
                point.map(f64::to_bits),
                curve.value(expected).map(f64::to_bits)
            );
        }
    }

    #[test]
    fn lazy_coefficients_keep_points_units_and_scalar_bits_separate() {
        let pool = CurveInversionPool::default();
        let units = [1_f64.to_bits(); 2];
        let a = curve(0.);
        let b = curve(5.);
        check(&a, &pool, units);
        assert!(pool.0.lock().unwrap().values().all(|s| s.get().is_none()));
        check(&b, &pool, units);
        assert_eq!(pool.0.lock().unwrap().len(), 1);
        assert!(
            pool.0
                .lock()
                .unwrap()
                .values()
                .all(|s| s.get().unwrap().is_some())
        );
        assert_ne!(a.table.get().unwrap()[100].1, b.table.get().unwrap()[100].1);
        check(&curve(9.), &pool, [2_f64.to_bits(); 2]);
        let mut negative_zero = curve(10.);
        negative_zero.knots[0] = -0.;
        check(&negative_zero, &pool, units);
        let mut other_weight = curve(11.);
        other_weight.weights[1] = f64::from_bits(other_weight.weights[1].to_bits() + 1);
        check(&other_weight, &pool, units);
        assert_eq!(pool.0.lock().unwrap().len(), 4);
        let separate = CurveInversionPool::default();
        check(&curve(0.), &separate, units);
        assert!(
            separate
                .0
                .lock()
                .unwrap()
                .values()
                .all(|s| s.get().is_none())
        );
    }

    #[test]
    fn transformation_discards_points_and_reuses_only_coefficients() {
        let pool = CurveInversionPool::default();
        let units = [1_f64.to_bits(); 2];
        let mut c = curve(0.);
        check(&c, &pool, units);
        let before = c.table.get().unwrap()[200].1;
        c.transform_points(|p| [p[0] + 7., p[1] * 2., p[2]]);
        assert!(c.table.get().is_none());
        check(&c, &pool, units);
        assert_ne!(before, c.table.get().unwrap()[200].1);
        let ordinary = curve(0.);
        ordinary.invert([1., 1., 0.]);
        assert!(ordinary.table.get().is_some());
    }

    #[test]
    fn concurrent_queries_keep_each_curve_table_independent() {
        let pool = Arc::new(CurveInversionPool::default());
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..8)
                .map(|i| {
                    let pool = pool.clone();
                    scope.spawn(move || {
                        let c = curve(i as f64);
                        check(&c, &pool, [1_f64.to_bits(); 2]);
                        let original = curve(i as f64);
                        for p in [[0., 0., 0.], [1., 1., 0.], [3., 0., 0.]] {
                            assert_eq!(c.invert(p).to_bits(), original.invert(p).to_bits());
                        }
                    })
                })
                .collect();
            for handle in handles {
                handle.join().unwrap();
            }
        });
        assert_eq!(pool.0.lock().unwrap().len(), 1);
    }

    #[test]
    fn poisoned_map_preserves_the_original_builder() {
        let pool = CurveInversionPool::default();
        let _ = std::panic::catch_unwind(|| {
            let _lock = pool.0.lock().unwrap();
            panic!("poison only this test pool");
        });
        check(&curve(0.), &pool, [1_f64.to_bits(); 2]);
    }
}
