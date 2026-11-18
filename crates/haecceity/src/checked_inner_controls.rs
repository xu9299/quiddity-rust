//! Analytic inner-rule and canonical-mass controls, with no external geometry.
use super::*;

#[test]
fn smooth_inner_panels_keep_the_established_sum_bits() {
    let rule = IntegralRule {
        components: Components::Prefix(3),
        observe_mass: None,
    };
    for cuts in [vec![0., 0.1, 0.7, 1.], vec![1., 0.7, 0.1, 0.]] {
        let f = |t: f64| [1. + t * t, 1e6 * (t - 0.5), t.powi(5)];
        assert_eq!(
            checked_inner_integrate(&cuts, f, rule)
                .unwrap()
                .map(f64::to_bits),
            integrate(&cuts, f).map(f64::to_bits)
        );
    }
}

#[test]
fn a_fixed_inner_rule_cannot_hide_a_narrow_jacobian_minimum() {
    let a = 0.37;
    let epsilon = 0.001;
    let f = |v: f64| [(4. * (v - a).powi(2) + epsilon * epsilon).sqrt()];
    let primitive = |t: f64| {
        (t * (t * t + epsilon * epsilon).sqrt() + epsilon * epsilon * (t / epsilon).asinh()) / 4.
    };
    let exact = primitive(2. * (1. - a)) - primitive(-2. * a);
    assert!((integrate(&[0., 1.], f)[0] - exact).abs() > 1e-5);
    assert!(
        checked_inner_integrate(
            &[0., 1.],
            f,
            IntegralRule {
                components: Components::Prefix(1),
                observe_mass: None
            }
        )
        .is_none()
    );
}

fn regular_patch(reverse: bool) -> Part {
    let a = 0.5;
    let epsilon = 0.001;
    let row = |u: f64| {
        vec![
            [u, a * a, 0.],
            [u, a * a - a, epsilon / 2.],
            [u, (1. - a) * (1. - a), epsilon],
        ]
    };
    let surface = crate::nurbs::NurbsSurface::new(
        1,
        2,
        [0., 0.25, 0.5, 0.75, 1.].into_iter().map(row).collect(),
        vec![vec![1.; 3]; 5],
        vec![0., 0., 0.25, 0.5, 0.75, 1., 1.],
        vec![0., 0., 0., 1., 1., 1.],
    )
    .unwrap();
    let points = [
        surface.value(0., 0.),
        surface.value(1., 0.),
        surface.value(1., 1.),
        surface.value(0., 1.),
    ];
    let edges = (0..4)
        .map(|i| {
            let curve = if i % 2 == 0 {
                Curve::Line {
                    origin: points[i],
                    dir: geom::sub(points[(i + 1) % 4], points[i]),
                }
            } else {
                let mut poles = row(if i == 1 { 1. } else { 0. });
                if i == 3 {
                    poles.reverse();
                }
                Curve::Nurbs(
                    crate::nurbs::NurbsCurve::new(
                        2,
                        poles,
                        vec![1.; 3],
                        vec![0., 0., 0., 1., 1., 1.],
                    )
                    .unwrap(),
                )
            };
            crate::brep::Edge {
                curve,
                start: points[i],
                end: points[(i + 1) % 4],
                vertices: (i, (i + 1) % 4),
                same_sense: true,
                samples: vec![points[i], points[(i + 1) % 4]],
            }
        })
        .collect();
    let uses = if reverse {
        (0..4).rev().map(|i| (i, false)).collect()
    } else {
        (0..4).map(|i| (i, true)).collect()
    };
    Part::new(
        vec![crate::brep::Face {
            surface: Surface::Freeform {
                kind: "synthetic",
                surface: Box::new(surface),
            },
            reversed: false,
            loops: vec![crate::brep::Loop {
                edges: uses,
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        edges,
        vec![],
    )
}

#[test]
fn rejected_cost_axis_keeps_analytic_area_and_canonical_mass_independent_of_moments() {
    let a = 0.5;
    let epsilon = 0.001;
    let primitive = |t: f64| {
        (t * (t * t + epsilon * epsilon).sqrt() + epsilon * epsilon * (t / epsilon).asinh()) / 4.
    };
    let exact = primitive(2. * (1. - a)) - primitive(-2. * a);
    for reverse in [false, true] {
        let part = regular_patch(reverse);
        assert!(part.short_axis_walk_is_closed(0));
        let mass = part.face_terms(0).unwrap();
        assert!(
            (mass[0] - exact).abs() < 1e-9 * exact,
            "{mass:?}, exact {exact}"
        );
        let before = mass.map(f64::to_bits);
        let moments = part.face_moments(0).unwrap();
        assert!((moments.area - exact).abs() < 1e-9 * exact);
        assert_eq!(part.face_terms(0).unwrap().map(f64::to_bits), before);
        let other = regular_patch(reverse);
        other.face_moments(0).unwrap();
        assert_eq!(other.face_terms(0).unwrap().map(f64::to_bits), before);
    }
}

#[test]
fn a_material_known_kink_cannot_hide_behind_an_aliased_embedded_pair() {
    let c = 0.018719847821034433;
    let exact = (1.0_f64 - c).powi(2) / 2.;
    let coarse = kronrod((), 0., 1., |(), x| Some([(x - c).max(0.)])).unwrap();
    assert!(coarse.1 < 1e-14);
    assert!((coarse.0[0] - exact).abs() > 1e-5);
    for (a, b, sign) in [(0., 1., 1.), (1., 0., -1.)] {
        let floor = 1e-12;
        let mut eval = |(): (), a: f64, b: f64, locate: bool| {
            let panel = kronrod((), a, b, |(), x| Some([(x - c).max(0.)]))?;
            let kink = (locate && c > a.min(b) && c < a.max(b)).then_some((c, panel.2 > floor));
            Some((panel, kink))
        };
        let result =
            adaptive_with_kink_budget(&[a, b], (), Some(floor), false, true, &mut eval).unwrap();
        assert!((result.0[0] - sign * exact).abs() < 1e-13);
    }
}

#[test]
fn a_cost_walk_refuses_an_unvalidated_remainder_after_kink_budget_exhaustion() {
    let run = |require_budget| {
        let mut eval = |state: usize, a: f64, b: f64, locate: bool| {
            let panel = ([b - a], 0., (b - a).abs(), state + 1);
            Some((panel, locate.then_some((a + 0.01 * (b - a), true))))
        };
        adaptive_with_kink_budget(&[0., 1.], 0, Some(1e-12), false, require_budget, &mut eval)
    };
    assert!(run(true).is_none());
    assert!((run(false).unwrap().0[0] - 1.).abs() < 1e-12);
}
