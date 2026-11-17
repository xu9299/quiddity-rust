//! Finite queries captured from the existing public section-recess fixtures.
//! These test coverage only; complete caller documents are checked separately.

use haecceity::brep::{Edge, Face, Loop, Part};
use haecceity::cover::{Coverage, assess};
use haecceity::geom::{Curve, Frame, Surface};

#[test]
fn public_split_cylindrical_seat() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Cylinder {
                frame: Frame {
                    origin: [0.0, 2.99999999999999, 11.100000000000652],
                    x: [
                        -0.5216708224207978,
                        1.4285387170004708e-15,
                        -0.8531468531471055,
                    ],
                    y: [
                        0.8531468531471055,
                        8.825462222138132e-16,
                        -0.5216708224207978,
                    ],
                    z: [0.0, -1.0, -1.6917684184796536e-15],
                },
                radius: 7.149999999998652,
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, false), (2, false), (3, true)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Circle {
                    frame: Frame {
                        origin: [0.0, 2.99999999999999, 11.100000000000652],
                        x: [
                            -0.5216708224207978,
                            1.4285387170004708e-15,
                            -0.8531468531471055,
                        ],
                        y: [
                            0.8531468531471055,
                            8.825462222138132e-16,
                            -0.5216708224207978,
                        ],
                        z: [0.0, -1.0, -1.6917684184796536e-15],
                    },
                    radius: 7.149999999998652,
                },
                start: [-3.729946380308, 3.0, 4.999999999999999],
                end: [3.729946380308, 3.0, 4.999999999999999],
                vertices: (2, 3),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [3.729946380308, -3.0, 5.0],
                    dir: [-0.0, 6.0, -1.3322676295500001e-15],
                },
                start: [3.729946380308, -3.0, 5.0],
                end: [3.729946380308, 3.0, 4.999999999999999],
                vertices: (1, 3),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Circle {
                    frame: Frame {
                        origin: [0.0, -3.00000000000001, 11.100000000000653],
                        x: [
                            -0.5216708224207978,
                            1.4285387170004708e-15,
                            -0.8531468531471055,
                        ],
                        y: [
                            0.8531468531471055,
                            8.825462222138132e-16,
                            -0.5216708224207978,
                        ],
                        z: [0.0, -1.0, -1.6917684184796536e-15],
                    },
                    radius: 7.149999999998652,
                },
                start: [-3.729946380308, -3.0, 5.0],
                end: [3.729946380308, -3.0, 5.0],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [-3.729946380308, -3.0, 5.0],
                    dir: [-0.0, 6.0, -1.3322676295500001e-15],
                },
                start: [-3.729946380308, -3.0, 5.0],
                end: [-3.729946380308, 3.0, 4.999999999999999],
                vertices: (0, 2),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_recesses/parts/00b0026a353bfcdc.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(
        assess((&patch, 0), &[(&source, 7), (&source, 8)]),
        Coverage::Covered
    );
}

#[test]
fn public_finite_horizontal_cylindrical_pcurves() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Cylinder {
                frame: Frame {
                    origin: [0.0, 30.0, 50.0],
                    x: [-1.0, 0.0, -1.9875321280002802e-16],
                    y: [1.98753212800028e-16, 0.0, -0.9999999999999999],
                    z: [0.0, -0.9999999999999999, 0.0],
                },
                radius: 71.5,
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, false), (2, false), (3, true)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Circle {
                    frame: Frame {
                        origin: [0.0, 30.0, 50.0],
                        x: [-1.0, 0.0, -1.9875321280002802e-16],
                        y: [1.98753212800028e-16, 0.0, -0.9999999999999999],
                        z: [0.0, -0.9999999999999999, 0.0],
                    },
                    radius: 71.5,
                },
                start: [-71.5, 30.0, 49.999999999999986],
                end: [71.5, 30.0, 49.999999999999986],
                vertices: (2, 3),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [71.5, -30.0, 50.0],
                    dir: [-0.0, 60.0, -1.33226762955e-14],
                },
                start: [71.5, -30.0, 50.0],
                end: [71.5, 30.0, 49.999999999999986],
                vertices: (1, 3),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Circle {
                    frame: Frame {
                        origin: [0.0, -30.0, 50.000000000000014],
                        x: [-1.0, 0.0, -1.9875321280002802e-16],
                        y: [1.98753212800028e-16, 0.0, -0.9999999999999999],
                        z: [0.0, -0.9999999999999999, 0.0],
                    },
                    radius: 71.5,
                },
                start: [-71.5, -30.0, 50.0],
                end: [71.5, -30.0, 50.0],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [-71.5, -30.0, 50.0],
                    dir: [-0.0, 60.0, -1.33226762955e-14],
                },
                start: [-71.5, -30.0, 50.0],
                end: [-71.5, 30.0, 49.999999999999986],
                vertices: (0, 2),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_recesses/parts/42a1f432aab027f5.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 7)]), Coverage::Covered);
}

#[test]
fn public_split_planar_end_arc() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [124.7177, -63.33610166537973, 91.1345],
                    x: [0.0, 1.0, 0.0],
                    y: [0.0, 0.0, -1.0],
                    z: [-1.0, 0.0, 0.0],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Circle {
                    frame: Frame {
                        origin: [124.7177, -63.7234, 91.2345],
                        x: [0.0, 0.9682458365517839, -0.2500000000002729],
                        y: [-0.0, -0.2500000000002729, -0.9682458365517839],
                        z: [-1.0, 0.0, 0.0],
                    },
                    radius: 0.4,
                },
                start: [124.7177, -63.33610166537973, 91.1345],
                end: [124.7177, -63.33610166537973, 91.33449999999999],
                vertices: (0, 1),
                same_sense: false,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [124.7177, -63.33610166537973, 91.33449999999999],
                    dir: [0.0, 6.212701665379733, 0.0],
                },
                start: [124.7177, -63.33610166537973, 91.33449999999999],
                end: [124.7177, -57.1234, 91.33449999999999],
                vertices: (1, 5),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [124.7177, -57.1234, 91.1345],
                    dir: [0.0, 0.0, 0.19999999999998863],
                },
                start: [124.7177, -57.1234, 91.1345],
                end: [124.7177, -57.1234, 91.33449999999999],
                vertices: (4, 5),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [124.7177, -63.33610166537973, 91.1345],
                    dir: [0.0, 6.212701665379733, 0.0],
                },
                start: [124.7177, -63.33610166537973, 91.1345],
                end: [124.7177, -57.1234, 91.1345],
                vertices: (0, 4),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_recesses/parts/332c6144380c21d7.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(
        assess((&patch, 0), &[(&source, 7), (&source, 8), (&source, 9)]),
        Coverage::Covered
    );
}
