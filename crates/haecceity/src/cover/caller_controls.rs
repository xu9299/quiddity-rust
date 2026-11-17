//! Independent cylindrical union controls with analytic axial and angular gaps.

use super::{Coverage, assess};
use crate::brep::{Edge, Face, Loop, Part};
use crate::geom::{Curve, Frame, Surface};

fn strip(start: f64, end: f64, low: f64, high: f64, phase: f64) -> Part {
    let frame = |z| Frame {
        origin: [0.0, 0.0, z],
        x: [1.0, 0.0, 0.0],
        y: [0.0, 1.0, 0.0],
        z: [0.0, 0.0, 1.0],
    };
    let arc = |z| Curve::Circle {
        frame: frame(z),
        radius: 2.0,
    };
    let points = [
        arc(low).value(start),
        arc(low).value(end),
        arc(high).value(end),
        arc(high).value(start),
    ];
    let edges = vec![
        Edge {
            curve: arc(low),
            start: points[0],
            end: points[1],
            vertices: (0, 1),
            same_sense: true,
            samples: vec![],
        },
        Edge {
            curve: Curve::Line {
                origin: points[1],
                dir: [0.0, 0.0, 1.0],
            },
            start: points[1],
            end: points[2],
            vertices: (1, 2),
            same_sense: true,
            samples: vec![],
        },
        Edge {
            curve: arc(high),
            start: points[3],
            end: points[2],
            vertices: (3, 2),
            same_sense: true,
            samples: vec![],
        },
        Edge {
            curve: Curve::Line {
                origin: points[0],
                dir: [0.0, 0.0, 1.0],
            },
            start: points[0],
            end: points[3],
            vertices: (0, 3),
            same_sense: true,
            samples: vec![],
        },
    ];
    let mut surface_frame = frame(0.0);
    surface_frame.x = [phase.cos(), phase.sin(), 0.0];
    surface_frame.y = [-phase.sin(), phase.cos(), 0.0];
    Part::new(
        vec![Face {
            surface: Surface::Cylinder {
                frame: surface_frame,
                radius: 2.0,
            },
            reversed: false,
            solid: None,
            pcurves: vec![],
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
        }],
        edges,
        vec![],
    )
}

#[test]
fn cylindrical_unions_preserve_seams_and_real_gaps() {
    for start in [0.1, 2.9, -3.1] {
        for phase in [0.0, 0.37, 1.2, 3.7] {
            let patch = strip(start, start + 0.9, 0.0, 3.0, phase);
            let lower = strip(start, start + 0.9, 0.0, 1.0, 0.0);
            let upper = strip(start, start + 0.9, 1.0, 3.0, 1.2);
            assert_eq!(
                assess((&patch, 0), &[(&lower, 0), (&upper, 0)]),
                Coverage::Covered
            );
            let missing = strip(start, start + 0.9, 1.0001, 3.0, 1.2);
            assert_ne!(
                assess((&patch, 0), &[(&lower, 0), (&missing, 0)]),
                Coverage::Covered
            );
            let left = strip(start, start + 0.4, 0.0, 3.0, 0.0);
            let right = strip(start + 0.4001, start + 0.9, 0.0, 3.0, 1.2);
            assert_ne!(
                assess((&patch, 0), &[(&left, 0), (&right, 0)]),
                Coverage::Covered
            );
        }
    }
}

#[test]
fn cylindrical_unions_do_not_omit_inner_loops() {
    let patch = strip(2.9, 3.8, 0.0, 3.0, 0.0);
    let outer = strip(2.9, 3.8, 0.0, 3.0, 0.0);
    let hole = strip(3.0, 3.1, 1.0, 2.0, 0.0);
    let mut faces = outer.faces.clone();
    let mut edges = outer.edges.clone();
    edges.extend(hole.edges);
    faces[0].loops.push(Loop {
        edges: vec![(4, true), (5, true), (6, false), (7, false)],
        vertex: None,
    });
    let support = Part::new(faces, edges, vec![]);
    assert_ne!(assess((&patch, 0), &[(&support, 0)]), Coverage::Covered);
}
