//! Finite queries captured from the existing public section-recess fixtures.
//! These test coverage only; complete caller documents are checked separately.

use haecceity::brep::{Edge, Face, Loop, Part};
use haecceity::cover::{Coverage, assess};
use haecceity::geom::{Curve, Frame, Surface};

#[test]
fn public_circular_aperture_band_and_finite_wall_filler() {
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_recesses/parts/137798c6ae23078d.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    let Surface::Cylinder { frame, radius } = source.faces[11].surface else {
        panic!()
    };
    let rotated = Frame {
        x: frame.y,
        y: haecceity::geom::scale(frame.x, -1.0),
        ..frame
    };
    let a = Curve::Circle {
        frame: rotated,
        radius,
    };
    let b = Curve::Circle {
        frame: Frame {
            origin: haecceity::geom::add(frame.origin, haecceity::geom::scale(frame.z, 4.0)),
            ..rotated
        },
        radius,
    };
    let (p, q) = (a.value(0.0), b.value(0.0));
    let edge = |curve, start, end, vertices| Edge {
        curve,
        start,
        end,
        vertices,
        same_sense: true,
        samples: vec![],
    };
    let band = Part::new(
        vec![Face {
            surface: Surface::Cylinder {
                frame: rotated,
                radius,
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (2, true), (1, false), (2, false)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        vec![
            edge(a, p, p, (0, 0)),
            edge(b, q, q, (1, 1)),
            edge(
                Curve::Line {
                    origin: p,
                    dir: haecceity::geom::sub(q, p),
                },
                p,
                q,
                (0, 1),
            ),
        ],
        vec![],
    );
    assert_eq!(assess((&source, 11), &[(&band, 0)]), Coverage::Covered);
    assert_eq!(assess((&band, 0), &[(&source, 11)]), Coverage::Covered);
    let mut wall = source.faces[6].clone();
    wall.loops.truncate(1);
    wall.pcurves.clear();
    wall.solid = None;
    let wall = Part::new(vec![wall], source.edges.clone(), vec![]);
    let circle = source.edges[22].clone();
    let Curve::Circle { frame, .. } = circle.curve else {
        panic!()
    };
    let disk = Part::new(
        vec![Face {
            surface: Surface::Plane { frame },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        vec![circle],
        vec![],
    );
    assert_ne!(assess((&wall, 0), &[(&source, 6)]), Coverage::Covered);
    assert_eq!(
        assess((&wall, 0), &[(&source, 6), (&disk, 0)]),
        Coverage::Covered
    );
}

#[test]
fn public_planar_end_with_a_preferred_cylindrical_polynomial() {
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_geometry/parts/0883f49ece1bde7d.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    // The preferred adjacent cylinder pcurve must participate even on this planar end.
    assert_eq!(assess((&source, 3), &[(&source, 3)]), Coverage::Covered);
}

#[test]
fn public_mixed_pocket_cylindrical_generatrices() {
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
        "../../tests/fixtures/captured/section_recess_geometry/parts/0262e62c10a5adb5.step.gz",
    );
    let source = haecceity::step::read_step_file(&file).unwrap();
    let swept = haecceity::sweep::extrude_face(&source, 14, [0.0; 3], [0.0, 0.0, 19.0]).unwrap();
    let mut supports: Vec<_> = (6..=14).map(|face| (&source, face)).collect();
    supports.push((&swept, 1));
    for face in 0..swept.faces.len() {
        assert_eq!(
            assess((&swept, face), &supports),
            Coverage::Covered,
            "face {face}"
        );
    }
}

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

#[test]
fn captured_cylindrical_polynomial_0883f49ece1bde7d() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [11.0, 7.0, 19.0],
                    x: [0.0, 0.0, 1.0],
                    y: [0.49999999999997347, -0.8660254037844541, 0.0],
                    z: [0.8660254037844541, 0.49999999999997347, -0.0],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: Some(0),
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Line {
                    origin: [11.0, 7.0, 19.0],
                    dir: [-4.0, 6.928203230276001, 0.0],
                },
                start: [11.0, 7.0, 19.0],
                end: [7.0, 13.928203230276, 19.0],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [7.0, 13.928203230276, 19.0],
                    dir: [0.0, 0.0, 10.761663039293538],
                },
                start: [7.0, 13.928203230276, 19.0],
                end: [7.0, 13.928203230276, 29.761663039293538],
                vertices: (1, 7),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [11.0, 7.0, 11.000000000000002],
                        x: [
                            0.49999999999997335,
                            -0.8660254037844539,
                            4.4408920985002357e-16,
                        ],
                        y: [-2.22044604925e-16, 3.8459253727668584e-16, 1.0],
                        z: [-0.8660254037844539, -0.49999999999997335, 0.0],
                    },
                    major: 23.094010767584624,
                    minor: 20.0,
                },
                start: [11.0, 7.0, 31.0],
                end: [7.0, 13.928203230276, 29.761663039293538],
                vertices: (6, 7),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [11.0, 7.0, 19.0],
                    dir: [0.0, 0.0, 12.0],
                },
                start: [11.0, 7.0, 19.0],
                end: [11.0, 7.0, 31.0],
                vertices: (0, 6),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_geometry/parts/0883f49ece1bde7d.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 3)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_00eb6f057ad5d83c() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [123.76710448242248, -56.48285943842785, 90.8693194333309],
                    x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712],
                    y: [
                        0.049178299740918655,
                        -0.2447215247204723,
                        -0.968345429158975,
                    ],
                    z: [-0.7775521304173983, 0.599139657323494, -0.19090404790273907],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: Some(0),
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [123.28132399198653, -56.94439001073898, 91.39942022179841],
                        x: [-0.04917829974091864, 0.24472152472047226, 0.968345429158975],
                        y: [
                            -0.6268924782751484,
                            -0.7623273879120266,
                            0.16081907325087116,
                        ],
                        z: [0.7775521304173983, -0.599139657323494, 0.19090404790273904],
                    },
                    major: 1.1313708498985076,
                    minor: 0.8,
                },
                start: [123.76710448242248, -56.48285943842785, 90.8693194333309],
                end: [123.7828379746067, -56.3345281004093, 91.27076496319768],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [123.7828379746067, -56.3345281004093, 91.27076496319768],
                    dir: [0.7522709739301234, 0.9147928654943627, -0.19298288790102447],
                },
                start: [123.7828379746067, -56.3345281004093, 91.27076496319768],
                end: [124.53510894853682, -55.41973523491494, 91.07778207529665],
                vertices: (1, 5),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [124.55597353407724, -55.52356178469422, 90.66694790360307],
                    dir: [
                        -0.02086458554042281,
                        0.10382654977928496,
                        0.4108341716935797,
                    ],
                },
                start: [124.55597353407724, -55.52356178469422, 90.66694790360307],
                end: [124.53510894853682, -55.41973523491494, 91.07778207529665],
                vertices: (4, 5),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [123.76710448242248, -56.48285943842785, 90.8693194333309],
                    dir: [0.788869051654757, 0.9592976537336284, -0.20237152972782724],
                },
                start: [123.76710448242248, -56.48285943842785, 90.8693194333309],
                end: [124.55597353407724, -55.52356178469422, 90.66694790360307],
                vertices: (0, 4),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_recesses/parts/00eb6f057ad5d83c.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 14)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_022df9198c3a6396() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [13.510404482422526, -6.7594594384278, 3.7348194333308906],
                    x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712],
                    y: [
                        0.049178299740901196,
                        -0.24472152472045883,
                        -0.9683454291589796,
                    ],
                    z: [-0.7775521304173996, 0.5991396573234995, -0.1909040479027173],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: Some(0),
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [13.024623991986537, -7.220990010738985, 4.264920221798404],
                        x: [
                            -0.04917829974090107,
                            0.24472152472045883,
                            0.9683454291589793,
                        ],
                        y: [
                            -0.6268924782751484,
                            -0.7623273879120267,
                            0.16081907325087128,
                        ],
                        z: [0.7775521304173995, -0.5991396573234994, 0.19090404790271723],
                    },
                    major: 1.1313708498984818,
                    minor: 0.8,
                },
                start: [13.510404482422526, -6.7594594384278, 3.7348194333308906],
                end: [13.526137974606657, -6.611128100409364, 4.136264963197706],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [13.526137974606657, -6.611128100409364, 4.136264963197706],
                    dir: [0.7522709739301785, 0.9147928654944328, -0.19298288790104579],
                },
                start: [13.526137974606657, -6.611128100409364, 4.136264963197706],
                end: [14.278408948536836, -5.696335234914931, 3.9432820752966604],
                vertices: (1, 5),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [14.299273534077244, -5.800161784694222, 3.532447903603072],
                    dir: [
                        -0.020864585540408598,
                        0.10382654977929029,
                        0.4108341716935886,
                    ],
                },
                start: [14.299273534077244, -5.800161784694222, 3.532447903603072],
                end: [14.278408948536836, -5.696335234914931, 3.9432820752966604],
                vertices: (4, 5),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [13.510404482422526, -6.7594594384278, 3.7348194333308906],
                    dir: [0.7888690516547179, 0.9592976537335787, -0.2023715297278188],
                },
                start: [13.510404482422526, -6.7594594384278, 3.7348194333308906],
                end: [14.299273534077244, -5.800161784694222, 3.532447903603072],
                vertices: (0, 4),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_geometry/parts/022df9198c3a6396.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 14)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_11e9e0b548aeed6c() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [12.135444425482607, 11.093727383783524, 16.27115274453129],
                    x: [0.5150380749101773, -0.25061146493908626, 0.819713166317282],
                    y: [0.8351981583167706, 0.3618449018307332, -0.4141404391790088],
                    z: [-0.19282068804226837, 0.8979210213934143, 0.3956741356251556],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: Some(0),
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Line {
                    origin: [12.135444425482607, 11.093727383783524, 16.27115274453129],
                    dir: [-11.572845156741018, -5.013870035444908, 5.73849825701555],
                },
                start: [12.135444425482607, 11.093727383783524, 16.27115274453129],
                end: [0.5625992687415895, 6.079857348338615, 22.00965100154684],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [0.5625992687415895, 6.079857348338615, 22.00965100154684],
                    dir: [5.542666214589753, -2.6969961394581725, 8.82147688477897],
                },
                start: [0.5625992687415895, 6.079857348338615, 22.00965100154684],
                end: [6.105265483331342, 3.382861208880443, 30.83112788632581],
                vertices: (1, 4),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [8.015139826201194, 13.098619103296215, 9.713447413993034],
                        x: [
                            0.8351981583167707,
                            0.36184490183073315,
                            -0.41414043917900856,
                        ],
                        y: [0.5150380749101771, -0.25061146493908637, 0.8197131663172822],
                        z: [
                            0.19282068804226846,
                            -0.8979210213934147,
                            -0.3956741356251557,
                        ],
                    },
                    major: 39.999999999997875,
                    minor: 20.0,
                },
                start: [18.315901324404734, 8.086389804514488, 26.107710740338675],
                end: [6.105265483331342, 3.382861208880443, 30.83112788632581],
                vertices: (3, 4),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [12.135444425482607, 11.093727383783524, 16.27115274453129],
                    dir: [6.180456898922127, -3.007337579269036, 9.836557995807386],
                },
                start: [12.135444425482607, 11.093727383783524, 16.27115274453129],
                end: [18.315901324404734, 8.086389804514488, 26.107710740338675],
                vertices: (0, 3),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
        "../../tests/fixtures/captured/section_recess_geometry/parts/11e9e0b548aeed6c.step.gz",
    );
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 3)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_160a9f5d03a23b16() {
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                    x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712],
                    y: [
                        -0.1537427100993845,
                        -0.08131408698477327,
                        -0.9847594621780131,
                    ],
                    z: [-0.7637859646349486, 0.6420630599069372, 0.06622708909201423],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: Some(0),
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [122.44416614564037, -56.08988534512533, 92.18666734450109],
                        x: [0.1537427100993845, 0.08131408698477326, 0.9847594621780131],
                        y: [
                            -0.6268924782751486,
                            -0.7623273879120267,
                            0.16081907325087116,
                        ],
                        z: [0.763785964634949, -0.6420630599069375, -0.06622708909201425],
                    },
                    major: 9.237604307033733,
                    minor: 8.0,
                },
                start: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                end: [127.63665066041575, -49.95858754370141, 92.6285857012957],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [127.63665066041575, -49.95858754370141, 92.6285857012957],
                    dir: [7.611655235913062, 9.256090086845319, -1.9526464000855412],
                },
                start: [127.63665066041575, -49.95858754370141, 92.6285857012957],
                end: [135.24830589632882, -40.70249745685609, 90.67593930121016],
                vertices: (1, 4),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [134.44943534077242, -41.125017846942185, 85.55897903603072],
                    dir: [0.7988705555563911, 0.4225203900860919, 5.1169602651794435],
                },
                start: [134.44943534077242, -41.125017846942185, 85.55897903603072],
                end: [135.24830589632882, -40.70249745685609, 90.67593930121016],
                vertices: (3, 4),
                same_sense: true,
                samples: vec![],
            },
            Edge {
                curve: Curve::Line {
                    origin: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                    dir: [7.888690516547172, 9.592976537335787, -2.023715297278187],
                },
                start: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                end: [134.44943534077242, -41.125017846942185, 85.55897903603072],
                vertices: (0, 3),
                same_sense: true,
                samples: vec![],
            },
        ],
        vec![],
    );
    let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/captured/section_recesses/parts/160a9f5d03a23b16.step.gz");
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 12)]), Coverage::Covered);
}

#[test]
fn captured_triangular_boundary_starts_inside_the_other_ellipse() {
    // The second public 160a9f5d03a23b16 passage reaches this later support query.
    let patch = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [110.07574593702597, -78.10777387752127, 95.05629556399606],
                    x: [0.6268924782751484, 0.7623273879120268, -0.16081907325087094],
                    y: [
                        -0.7383294034776023,
                        0.5153858772186564,
                        -0.43502538951647896,
                    ],
                    z: [-0.24874788972461787, 0.3914515949677336, 0.8859402554093402],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, false), (1, true), (2, false), (3, false), (4, false)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [124.4692338543598, -58.156914654874846, 90.28233265549872],
                        x: [0.7383294034776022, -0.5153858772186565, 0.43502538951647884],
                        y: [0.6268924782751486, 0.7623273879120267, -0.1608190732508709],
                        z: [-0.24874788972461787, 0.3914515949677336, 0.8859402554093402],
                    },
                    major: 9.237604307034,
                    minor: 8.0,
                },
                start: [130.6742522234947, -53.05913150832528, 89.77208366779239],
                end: [129.484373680561, -52.05829555157864, 88.99578006949176],
                vertices: (20, 19),
                same_sense: true,
                samples: vec![
                    [130.6742522234947, -53.05913150832528, 89.77208366779239],
                    [130.6059314317566, -52.990878278726925, 89.72274344638242],
                    [130.53675783731995, -52.923342961279864, 89.67348098979394],
                    [130.4667410530796, -52.85653494121166, 89.62430314391618],
                    [130.39589080910645, -52.79046350267867, 89.57521674288087],
                    [130.3242169512954, -52.7251378274762, 89.52622860811188],
                    [130.251729439997, -52.660566993762565, 89.47734554737721],
                    [130.17843834863342, -52.59675997479749, 89.42857435384299],
                    [130.10435386229847, -52.53372563769514, 89.37992180512944],
                    [130.02948627634223, -52.471472742191835, 89.33139466236894],
                    [129.95384599494028, -52.4100099394288, 89.28299966926652],
                    [129.87744352964796, -52.34934577074986, 89.23474355116268],
                    [129.80028949793956, -52.28948866651453, 89.18663301409877],
                    [129.72239462173278, -52.23044694492642, 89.13867474388508],
                    [129.64376972589878, -52.17222881087732, 89.09087540517174],
                    [129.5644257367578, -52.114842354806946, 89.04324164052252],
                    [129.484373680561, -52.05829555157864, 88.99578006949176],
                ],
            },
            Edge {
                curve: Curve::Line {
                    origin: [110.07574593702597, -78.10777387752127, 95.05629556399606],
                    dir: [0.6268924782751486, 0.7623273879120266, -0.16081907325087094],
                },
                start: [130.6742522234947, -53.05913150832528, 89.77208366779239],
                end: [138.28590745940767, -43.803041421480074, 87.81943726770686],
                vertices: (20, 16),
                same_sense: true,
                samples: vec![
                    [130.6742522234947, -53.05913150832528, 89.77208366779239],
                    [138.28590745940767, -43.803041421480074, 87.81943726770686],
                ],
            },
            Edge {
                curve: Curve::Line {
                    origin: [134.62431624941107, -41.247092289708576, 85.66201927124419],
                    dir: [0.7383294034776025, -0.5153858772186565, 0.4350253895164788],
                },
                start: [134.44943534077242, -41.1250178469422, 85.55897903603073],
                end: [138.28590745940767, -43.803041421480074, 87.81943726770686],
                vertices: (18, 16),
                same_sense: true,
                samples: vec![
                    [134.44943534077242, -41.1250178469422, 85.55897903603073],
                    [138.28590745940767, -43.803041421480074, 87.81943726770686],
                ],
            },
            Edge {
                curve: Curve::Line {
                    origin: [106.23927381839076, -75.42975030298341, 92.79583733231993],
                    dir: [0.6268924782751486, 0.7623273879120266, -0.16081907325087094],
                },
                start: [126.56074482422548, -50.71799438427772, 87.58269433330886],
                end: [134.44943534077242, -41.1250178469422, 85.55897903603073],
                vertices: (23, 18),
                same_sense: true,
                samples: vec![
                    [126.56074482422548, -50.71799438427772, 87.58269433330886],
                    [134.44943534077242, -41.1250178469422, 85.55897903603073],
                ],
            },
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [124.4692338543598, -58.156914654874846, 90.28233265549872],
                        x: [0.7383294034776022, -0.5153858772186565, 0.43502538951647884],
                        y: [0.6268924782751486, 0.7623273879120267, -0.1608190732508709],
                        z: [-0.24874788972461787, 0.3914515949677336, 0.8859402554093402],
                    },
                    major: 9.237604307034,
                    minor: 8.0,
                },
                start: [129.484373680561, -52.05829555157864, 88.99578006949176],
                end: [126.56074482422548, -50.71799438427772, 87.58269433330886],
                vertices: (19, 23),
                same_sense: true,
                samples: vec![
                    [129.484373680561, -52.05829555157864, 88.99578006949176],
                    [129.40208451261736, -52.001546728803355, 88.94760111633947],
                    [129.31908355405778, -51.945686100995054, 88.89961475960864],
                    [129.2353827815885, -51.89072172861886, 88.85182792353856],
                    [129.15099427289618, -51.83666154281381, 88.80424750357855],
                    [129.065930204905, -51.78351334424845, 88.75688036539282],
                    [128.98020285201963, -51.73128480199521, 88.70973334386994],
                    [128.8938245843542, -51.6799834524238, 88.66281324213648],
                    [128.80680786594718, -51.629616698113736, 88.61612683057544],
                    [128.71916525296285, -51.580191806786154, 88.56968084584922],
                    [128.63090939187973, -51.53171591025516, 88.52348198992762],
                    [128.54205301766544, -51.484196003398694, 88.4775369291207],
                    [128.45260895193934, -51.4376389431492, 88.43185229311698],
                    [128.3625901011223, -51.39205144750423, 88.38643467402663],
                    [128.27200945457432, -51.34744009455703, 88.3412906254304],
                    [128.18088008272036, -51.30381132154735, 88.29642666143391],
                    [128.0892151351642, -51.26117142393261, 88.25184925572769],
                    [127.99702783879111, -51.21952655447945, 88.20756484065303],
                    [127.9043314958591, -51.178882722375924, 88.16357980627389],
                    [127.81113948207965, -51.13924579236442, 88.11990049945474],
                    [127.71746524468752, -51.10062148389536, 88.07653322294483],
                    [127.62332230050042, -51.06301537030194, 88.03348423446867],
                    [127.52872423396857, -51.0264328779959, 87.99075974582307],
                    [127.43368469521454, -50.99087928568455, 87.94836592198082],
                    [127.33821739806356, -50.956359723609005, 87.90630888020111],
                    [127.24233611806469, -50.922879172803995, 87.86459468914681],
                    [127.14605469050304, -50.89044246437905, 87.82322936800882],
                    [127.04938700840344, -50.85905427882143, 87.7822188856375],
                    [126.95234702052572, -50.82871914532075, 87.74156915968138],
                    [126.85494872935189, -50.799441441115405, 87.70128605573329],
                    [126.75720618906573, -50.771225390861, 87.66137538648398],
                    [126.65913350352477, -50.7440750660207, 87.62184291088337],
                    [126.56074482422548, -50.71799438427772, 87.58269433330886],
                ],
            },
        ],
        vec![],
    );
    let support = Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: [130.6742522234946, -53.05913150832541, 89.77208366779243],
                    x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712],
                    y: [0.7383294034775338, -0.5153858772185486, 0.43502538951672254],
                    z: [
                        0.24874788972482076,
                        -0.39145159496787546,
                        -0.8859402554092204,
                    ],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: vec![(0, true), (1, true), (2, false), (3, false)],
                vertex: None,
            }],
            solid: Some(0),
            pcurves: vec![],
        }],
        vec![
            Edge {
                curve: Curve::Ellipse {
                    frame: Frame {
                        origin: [124.4692338543596, -58.15691465487465, 90.2823326554989],
                        x: [
                            -0.7383294034775341,
                            0.5153858772185487,
                            -0.43502538951672276,
                        ],
                        y: [-0.6268924782751485, -0.7623273879120267, 0.1608190732508712],
                        z: [
                            -0.24874788972482095,
                            0.39145159496787574,
                            0.8859402554092208,
                        ],
                    },
                    major: 9.237604307033715,
                    minor: 8.0,
                },
                start: [130.6742522234946, -53.05913150832541, 89.77208366779243],
                end: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                vertices: (2, 0),
                same_sense: true,
                samples: vec![
                    [130.6742522234946, -53.05913150832541, 89.77208366779243],
                    [130.60506096104493, -52.990018732644664, 89.72211925255746],
                    [130.53499528948294, -52.921642284894645, 89.67223467263986],
                    [130.46406519378058, -52.85401190931813, 89.62243703702926],
                    [130.3922807820979, -52.78713724383608, 89.57273344232493],
                    [130.31965228434237, -52.72102781867413, 89.52313097172456],
                    [130.2461900507111, -52.6556930550045, 89.47363669401469],
                    [130.171904550216, -52.59114226360333, 89.42425766256349],
                    [130.0968063691915, -52.52738464352386, 89.37500091431546],
                    [130.02090620978638, -52.46442928078547, 89.32587346878877],
                    [129.94421488843815, -52.402285147078835, 89.27688232707472],
                    [129.86674333433197, -52.340961098487384, 89.22803447084017],
                    [129.7885025878429, -52.28046587422524, 89.17933686133254],
                    [129.70950379896274, -52.22080809539181, 89.13079643838779],
                    [129.62975822571093, -52.16199626374315, 89.08242011944138],
                    [129.5492772325302, -52.10403876048048, 89.03421479854255],
                    [129.4680722886671, -52.046943845055715, 88.9861873453718],
                    [129.38615496653748, -51.990719653994454, 88.93834460426193],
                    [129.30353694007732, -51.93537419973644, 88.89069339322263],
                    [129.22022998307918, -51.88091536949373, 88.84324050296892],
                    [129.13624596751416, -51.82735092412669, 88.79599269595333],
                    [129.05159686184015, -51.774688497038014, 88.74895670540226],
                    [128.96629472929627, -51.72293559308488, 88.7021392343564],
                    [128.88035172618362, -51.67209958750944, 88.65554695471546],
                    [128.79378010013303, -51.62218772488784, 88.60918650628747],
                    [128.70659218835956, -51.5732071180977, 88.56306449584244],
                    [128.61880041590442, -51.52516474730457, 88.51718749617089],
                    [128.5304172938642, -51.47806745896714, 88.47156204514712],
                    [128.44145541760804, -51.43192196486156, 88.42619464479763],
                    [128.35192746498265, -51.38673484112498, 88.38109176037436],
                    [128.26184619450547, -51.34251252731836, 88.33625981943342],
                    [128.1712244435467, -51.299261325508816, 88.29170521091916],
                    [128.08007512649965, -51.256987399371475, 88.24743428425356],
                    [127.98841123294054, -51.21569677331114, 88.2034533484315],
                    [127.8962458257771, -51.17539533160371, 88.15976867112158],
                    [127.80359203938724, -51.136088817557685, 88.11638647777303],
                    [127.71046307774718, -51.097782832695636, 88.07331295072841],
                    [127.61687221254974, -51.06048283595596, 88.03055422834262],
                    [127.52283278131307, -51.024194142914936, 87.98811640410821],
                    [127.42835818547991, -50.988921925029224, 87.94600552578686],
                    [127.33346188850776, -50.954671208898844, 87.90422759454763],
                    [127.2381574139502, -50.92144687555088, 87.8627885641117],
                    [127.14245834352975, -50.889253659743886, 87.82169433990389],
                    [127.04637831520226, -50.85809614929312, 87.78095077821115],
                    [126.9499310212134, -50.82797878441676, 87.74056368534792],
                    [126.85313020614747, -50.79890585710313, 87.70053881682875],
                    [126.75598966496857, -50.770881510499045, 87.66088187654799],
                    [126.65852324105477, -50.743909738319395, 87.62159851596704],
                    [126.56074482422525, -50.71799438427797, 87.58269433330891],
                ],
            },
            Edge {
                curve: Curve::Line {
                    origin: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                    dir: [7.888690516547172, 9.592976537335787, -2.023715297278187],
                },
                start: [126.56074482422525, -50.71799438427797, 87.58269433330891],
                end: [134.44943534077242, -41.125017846942185, 85.55897903603072],
                vertices: (0, 3),
                same_sense: true,
                samples: vec![
                    [126.56074482422525, -50.71799438427797, 87.58269433330891],
                    [134.44943534077242, -41.125017846942185, 85.55897903603072],
                ],
            },
            Edge {
                curve: Curve::Line {
                    origin: [138.28590745940767, -43.80304142148009, 87.81943726770689],
                    dir: [-3.8364721186352426, 2.678023574537903, -2.260458231676168],
                },
                start: [138.28590745940767, -43.80304142148009, 87.81943726770689],
                end: [134.44943534077242, -41.125017846942185, 85.55897903603072],
                vertices: (5, 3),
                same_sense: true,
                samples: vec![
                    [138.28590745940767, -43.80304142148009, 87.81943726770689],
                    [134.44943534077242, -41.125017846942185, 85.55897903603072],
                ],
            },
            Edge {
                curve: Curve::Line {
                    origin: [130.6742522234946, -53.05913150832541, 89.77208366779243],
                    dir: [7.6116552359130765, 9.256090086845319, -1.9526464000855412],
                },
                start: [130.6742522234946, -53.05913150832541, 89.77208366779243],
                end: [138.28590745940767, -43.80304142148009, 87.81943726770689],
                vertices: (2, 5),
                same_sense: true,
                samples: vec![
                    [130.6742522234946, -53.05913150832541, 89.77208366779243],
                    [138.28590745940767, -43.80304142148009, 87.81943726770689],
                ],
            },
        ],
        vec![],
    );
    assert_eq!(assess((&patch, 0), &[(&support, 0)]), Coverage::Covered);
    let Surface::Plane { frame } = patch.faces[0].surface else {
        panic!()
    };
    let shift = haecceity::geom::scale(frame.x, 0.001);
    let shifted = support
        .transformed(&[
            [1.0, 0.0, 0.0, shift[0]],
            [0.0, 1.0, 0.0, shift[1]],
            [0.0, 0.0, 1.0, shift[2]],
        ])
        .unwrap();
    assert_ne!(assess((&patch, 0), &[(&shifted, 0)]), Coverage::Covered);
}
