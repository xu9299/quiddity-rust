//! Area enclosures for literal finite polygons. Exact orientation signs select
//! possible crossing events; outward intervals retain every event's height band.
//! Vertex cuts preserve a fixed active-edge set in each band. Uncertain swaps pay
//! their complete pairwise strip area instead of the whole patch width.

use super::interval::Interval as I;
use super::{intervals, subtract};

type Point = (f64, f64);
type Segment = (Point, Point);
type Loops = Vec<Vec<Point>>;

struct Crossing {
    height: I,
    a: Segment,
    b: Segment,
}

fn at(segment: Segment, y: I) -> I {
    let (a, b) = segment;
    I::point(a.0)
        + (I::point(b.0) - I::point(a.0)) * (y - I::point(a.1)) / (I::point(b.1) - I::point(a.1))
}

fn difference(a: Point, b: Point) -> [I; 2] {
    [I::point(a.0) - I::point(b.0), I::point(a.1) - I::point(b.1)]
}

fn cross(a: [I; 2], b: [I; 2]) -> I {
    a[0] * b[1] - a[1] * b[0]
}

fn crossing(a: Segment, b: Segment) -> Option<I> {
    let orient = |a: Point, b: Point, c: Point| {
        robust::orient2d(
            robust::Coord { x: a.0, y: a.1 },
            robust::Coord { x: b.0, y: b.1 },
            robust::Coord { x: c.0, y: c.1 },
        )
    };
    let opposite = |x: f64, y: f64| (x < 0.0 && y > 0.0) || (x > 0.0 && y < 0.0);
    if !opposite(orient(a.0, a.1, b.0), orient(a.0, a.1, b.1))
        || !opposite(orient(b.0, b.1, a.0), orient(b.0, b.1, a.1))
    {
        return None;
    }
    let r = difference(a.1, a.0);
    let s = difference(b.1, b.0);
    let t = difference(b.0, a.0);
    let denominator = cross(r, s);
    let lo = a.0.1.min(a.1.1).max(b.0.1.min(b.1.1));
    let hi = a.0.1.max(a.1.1).min(b.0.1.max(b.1.1));
    if denominator.lo <= 0.0 && denominator.hi >= 0.0 {
        return Some(I { lo, hi });
    }
    let height = I::point(a.0.1) + (cross(t, s) / denominator) * r[1];
    Some(I {
        lo: height.lo.max(lo),
        hi: height.hi.min(hi),
    })
}

fn scan_error(all: &[&Loops], v: f64, exact_v: I, x: (f64, f64)) -> Option<f64> {
    let mut error = I::point(0.0);
    for loops in all {
        for edge in loops.iter().flat_map(|lp| lp.windows(2)) {
            let (a, b) = (edge[0], edge[1]);
            if a.1 == b.1 {
                continue;
            }
            let low = a.1.min(b.1);
            let high = a.1.max(b.1);
            if high <= exact_v.lo || low >= exact_v.hi {
                continue;
            }
            if low >= exact_v.lo || high <= exact_v.hi {
                return None;
            }
            let denominator = I::point(b.1) - I::point(a.1);
            if denominator.lo <= 0.0 && denominator.hi >= 0.0 {
                return None;
            }
            let at = I::point(a.0)
                + (I::point(b.0) - I::point(a.0)) * (exact_v - I::point(a.1)) / denominator;
            let at = I {
                lo: at.lo.clamp(x.0, x.1),
                hi: at.hi.clamp(x.0, x.1),
            };
            let nominal = (a.0 + (b.0 - a.0) * (v - a.1) / (b.1 - a.1)).clamp(x.0, x.1);
            error = error + I::point(at.radius_from(nominal));
        }
    }
    error.hi.is_finite().then_some(error.hi)
}

fn band(
    mine: &Loops,
    supports: &[Loops],
    lo: f64,
    hi: f64,
    x: (f64, f64),
    crossings: &[Crossing],
) -> Option<I> {
    let height = I::point(hi) - I::point(lo);
    let width = I::point(x.1) - I::point(x.0);
    let upper = (height * width).hi;
    let exact_v = I::point(lo) + I::point(0.5) * height;
    let v = lo + 0.5 * (hi - lo);
    let all: Vec<_> = std::iter::once(mine).chain(supports).collect();
    let Some(error) = scan_error(&all, v, exact_v, x) else {
        return Some(I { lo: 0.0, hi: upper });
    };
    let clamp = |(a, b): (f64, f64)| (a.clamp(x.0, x.1), b.clamp(x.0, x.1));
    let mut left: Vec<_> = intervals(mine, v).into_iter().map(clamp).collect();
    for support in supports {
        for cut in intervals(support, v).into_iter().map(clamp) {
            left = subtract(&left, cut);
        }
    }
    let mut length = I::point(0.0);
    for (a, b) in left {
        let span = I::point(b) - I::point(a);
        length = length
            + I {
                lo: span.lo.max(0.0),
                hi: span.hi.max(0.0),
            };
    }
    // Sorting the perturbed scan crossings minimizes their L1 matching distance.
    // Their parity-region symmetric difference is therefore at most the sum of
    // the individual crossing displacements, even when their interval boxes overlap.
    let length = I {
        lo: (length - I::point(error)).lo.max(0.0),
        hi: (length + I::point(error)).hi.min(width.hi),
    };
    let mut swap_error = I::point(0.0);
    for event in crossings {
        if event.height.lo >= hi || event.height.hi <= lo {
            continue;
        }
        // Each segment pair changes order at most once. Relative to the ordering
        // at the midpoint, a swap changes a Boolean/parity scan length by at most
        // twice the pair separation. That separation is affine, so its absolute
        // maximum on this vertex-free band is attained at one of its endpoints.
        // Summing these complete strips also covers overlapping crossing bands.
        let separation = [lo, hi]
            .map(|y| (at(event.a, I::point(y)) - at(event.b, I::point(y))).abs_upper())
            .into_iter()
            .fold(0.0_f64, f64::max);
        swap_error = swap_error + I::point(2.0) * I::point(separation) * height;
    }
    let area = height * length;
    (area.lo.is_finite() && area.hi.is_finite() && swap_error.hi.is_finite()).then_some(I {
        lo: (area - swap_error).lo.max(0.0),
        hi: (area + swap_error).hi.min(upper),
    })
}

pub(super) fn area(mine: &Loops, supports: &[Loops], bounds: (Point, Point)) -> Option<I> {
    let (low, high) = bounds;
    if low.0 >= high.0 || low.1 >= high.1 {
        return Some(I::point(0.0));
    }
    let mut cuts = vec![low.1, high.1];
    let mut crossings = Vec::new();
    let mut segments = Vec::new();
    let mut total = 0;
    for edge in std::iter::once(mine)
        .chain(supports)
        .flatten()
        .flat_map(|lp| lp.windows(2))
    {
        let (a, b) = (edge[0], edge[1]);
        total += 1;
        if total > 16384
            || [a.0, a.1, b.0, b.1]
                .into_iter()
                .any(|v| !v.is_finite() || v.abs() > 1e75 || (v != 0.0 && v.abs() < 1e-75))
        {
            return None;
        }
        for y in [a.1, b.1] {
            if low.1 < y && y < high.1 {
                cuts.push(y);
            }
        }
        if a.1 != b.1
            && a.0.max(b.0) >= low.0
            && a.0.min(b.0) <= high.0
            && a.1.max(b.1) >= low.1
            && a.1.min(b.1) <= high.1
        {
            segments.push((a, b));
            if segments.len() > 512 {
                return None;
            }
        }
    }
    for (i, &a) in segments.iter().enumerate() {
        for &b in &segments[..i] {
            if let Some(mut event) = crossing(a, b) {
                event.lo = event.lo.max(low.1);
                event.hi = event.hi.min(high.1);
                if !event.lo.is_finite() || !event.hi.is_finite() {
                    return None;
                }
                if event.lo <= event.hi {
                    cuts.extend([event.lo, event.hi]);
                    crossings.push(Crossing {
                        height: event,
                        a,
                        b,
                    });
                }
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    if cuts.len().saturating_mul(crossings.len()) > 2_000_000 {
        return None;
    }
    let mut result = I::point(0.0);
    for endpoints in cuts.windows(2) {
        if endpoints[0] < endpoints[1] {
            result = result
                + band(
                    mine,
                    supports,
                    endpoints[0],
                    endpoints[1],
                    (low.0, high.0),
                    &crossings,
                )?;
        }
    }
    (result.lo.is_finite() && result.hi.is_finite()).then_some(I {
        lo: result.lo.max(0.0),
        hi: result.hi,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn almost_coincident_crossings_pay_their_pair_strip_not_the_patch_width() {
        let patch = vec![vec![
            (0.125, 0.25),
            (1.125, 0.25),
            (1.875, 1.25),
            (0.875, 1.25),
            (0.125, 0.25),
        ]];
        let d = 2.0_f64.powi(-51);
        let support = vec![vec![
            (0.125 + d, 0.25),
            (2.125, 0.25),
            (2.875, 1.25),
            (0.875 - d, 1.25),
            (0.125 + d, 0.25),
        ]];
        let bound = area(&patch, &[support], ((0.125, 0.25), (1.875, 1.25))).unwrap();
        // The left edge crosses at the exact midpoint; its triangular deficit is d/4.
        assert!(bound.lo <= d / 4.0 && bound.hi >= d / 4.0);
        assert!(bound.hi - bound.lo < 1e-12, "{bound:?}");
    }

    #[test]
    fn exact_triangular_gap_and_crossing_events_stay_in_the_area_bracket() {
        let patch = vec![vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 0.0),
        ]];
        let half = vec![vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)]];
        let bound = area(&patch, &[half], ((0.0, 0.0), (1.0, 1.0))).unwrap();
        assert!(bound.lo <= 0.5 && bound.hi >= 0.5 && bound.hi - bound.lo < 1e-12);
        let a = ((0.0, 0.0), (1.0, 1.0));
        let b = ((0.0, 1.0), (1.0, 0.0));
        let crossing = crossing(a, b).unwrap();
        assert!(crossing.lo <= 0.5 && crossing.hi >= 0.5);
    }
}
