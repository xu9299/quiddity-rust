//! Finite queries captured from the existing public section-recess fixtures.
//! These test coverage only; complete caller documents are checked separately.

use haecceity::brep::{Edge, Face, Loop, Part};
use haecceity::cover::{Coverage, assess};
use haecceity::geom::{Curve, Frame, Surface};

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
    let patch = Part::new(vec![Face { surface: Surface::Plane { frame: Frame { origin: [11.0, 7.0, 19.0], x: [0.0, 0.0, 1.0], y: [0.49999999999997347, -0.8660254037844541, 0.0], z: [0.8660254037844541, 0.49999999999997347, -0.0] } }, reversed: false, loops: vec![Loop { edges: vec![(0, true), (1, true), (2, false), (3, false)], vertex: None }], solid: Some(0), pcurves: vec![] }], vec![Edge { curve: Curve::Line { origin: [11.0, 7.0, 19.0], dir: [-4.0, 6.928203230276001, 0.0] }, start: [11.0, 7.0, 19.0], end: [7.0, 13.928203230276, 19.0], vertices: (0, 1), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [7.0, 13.928203230276, 19.0], dir: [0.0, 0.0, 10.761663039293538] }, start: [7.0, 13.928203230276, 19.0], end: [7.0, 13.928203230276, 29.761663039293538], vertices: (1, 7), same_sense: true, samples: vec![] },
Edge { curve: Curve::Ellipse { frame: Frame { origin: [11.0, 7.0, 11.000000000000002], x: [0.49999999999997335, -0.8660254037844539, 4.4408920985002357e-16], y: [-2.22044604925e-16, 3.8459253727668584e-16, 1.0], z: [-0.8660254037844539, -0.49999999999997335, 0.0] }, major: 23.094010767584624, minor: 20.0 }, start: [11.0, 7.0, 31.0], end: [7.0, 13.928203230276, 29.761663039293538], vertices: (6, 7), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [11.0, 7.0, 19.0], dir: [0.0, 0.0, 12.0] }, start: [11.0, 7.0, 19.0], end: [11.0, 7.0, 31.0], vertices: (0, 6), same_sense: true, samples: vec![] }], vec![]);
    let file = std::path::Path::new(env!(CARGO_MANIFEST_DIR)).join(../../tests/fixtures/captured/section_geometry/parts/0883f49ece1bde7d.step.gz);
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 3)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_00eb6f057ad5d83c() {
    let patch = Part::new(vec![Face { surface: Surface::Plane { frame: Frame { origin: [123.76710448242248, -56.48285943842785, 90.8693194333309], x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712], y: [0.049178299740918655, -0.2447215247204723, -0.968345429158975], z: [-0.7775521304173983, 0.599139657323494, -0.19090404790273907] } }, reversed: false, loops: vec![Loop { edges: vec![(0, true), (1, true), (2, false), (3, false)], vertex: None }], solid: Some(0), pcurves: vec![] }], vec![Edge { curve: Curve::Ellipse { frame: Frame { origin: [123.28132399198653, -56.94439001073898, 91.39942022179841], x: [-0.04917829974091864, 0.24472152472047226, 0.968345429158975], y: [-0.6268924782751484, -0.7623273879120266, 0.16081907325087116], z: [0.7775521304173983, -0.599139657323494, 0.19090404790273904] }, major: 1.1313708498985076, minor: 0.8 }, start: [123.76710448242248, -56.48285943842785, 90.8693194333309], end: [123.7828379746067, -56.3345281004093, 91.27076496319768], vertices: (0, 1), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [123.7828379746067, -56.3345281004093, 91.27076496319768], dir: [0.7522709739301234, 0.9147928654943627, -0.19298288790102447] }, start: [123.7828379746067, -56.3345281004093, 91.27076496319768], end: [124.53510894853682, -55.41973523491494, 91.07778207529665], vertices: (1, 5), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [124.55597353407724, -55.52356178469422, 90.66694790360307], dir: [-0.02086458554042281, 0.10382654977928496, 0.4108341716935797] }, start: [124.55597353407724, -55.52356178469422, 90.66694790360307], end: [124.53510894853682, -55.41973523491494, 91.07778207529665], vertices: (4, 5), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [123.76710448242248, -56.48285943842785, 90.8693194333309], dir: [0.788869051654757, 0.9592976537336284, -0.20237152972782724] }, start: [123.76710448242248, -56.48285943842785, 90.8693194333309], end: [124.55597353407724, -55.52356178469422, 90.66694790360307], vertices: (0, 4), same_sense: true, samples: vec![] }], vec![]);
    let file = std::path::Path::new(env!(CARGO_MANIFEST_DIR)).join(../../tests/fixtures/captured/section_recesses/parts/00eb6f057ad5d83c.step.gz);
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 14)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_022df9198c3a6396() {
    let patch = Part::new(vec![Face { surface: Surface::Plane { frame: Frame { origin: [13.510404482422526, -6.7594594384278, 3.7348194333308906], x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712], y: [0.049178299740901196, -0.24472152472045883, -0.9683454291589796], z: [-0.7775521304173996, 0.5991396573234995, -0.1909040479027173] } }, reversed: false, loops: vec![Loop { edges: vec![(0, true), (1, true), (2, false), (3, false)], vertex: None }], solid: Some(0), pcurves: vec![] }], vec![Edge { curve: Curve::Ellipse { frame: Frame { origin: [13.024623991986537, -7.220990010738985, 4.264920221798404], x: [-0.04917829974090107, 0.24472152472045883, 0.9683454291589793], y: [-0.6268924782751484, -0.7623273879120267, 0.16081907325087128], z: [0.7775521304173995, -0.5991396573234994, 0.19090404790271723] }, major: 1.1313708498984818, minor: 0.8 }, start: [13.510404482422526, -6.7594594384278, 3.7348194333308906], end: [13.526137974606657, -6.611128100409364, 4.136264963197706], vertices: (0, 1), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [13.526137974606657, -6.611128100409364, 4.136264963197706], dir: [0.7522709739301785, 0.9147928654944328, -0.19298288790104579] }, start: [13.526137974606657, -6.611128100409364, 4.136264963197706], end: [14.278408948536836, -5.696335234914931, 3.9432820752966604], vertices: (1, 5), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [14.299273534077244, -5.800161784694222, 3.532447903603072], dir: [-0.020864585540408598, 0.10382654977929029, 0.4108341716935886] }, start: [14.299273534077244, -5.800161784694222, 3.532447903603072], end: [14.278408948536836, -5.696335234914931, 3.9432820752966604], vertices: (4, 5), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [13.510404482422526, -6.7594594384278, 3.7348194333308906], dir: [0.7888690516547179, 0.9592976537335787, -0.2023715297278188] }, start: [13.510404482422526, -6.7594594384278, 3.7348194333308906], end: [14.299273534077244, -5.800161784694222, 3.532447903603072], vertices: (0, 4), same_sense: true, samples: vec![] }], vec![]);
    let file = std::path::Path::new(env!(CARGO_MANIFEST_DIR)).join(../../tests/fixtures/captured/section_geometry/parts/022df9198c3a6396.step.gz);
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 14)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_11e9e0b548aeed6c() {
    let patch = Part::new(vec![Face { surface: Surface::Plane { frame: Frame { origin: [12.135444425482607, 11.093727383783524, 16.27115274453129], x: [0.5150380749101773, -0.25061146493908626, 0.819713166317282], y: [0.8351981583167706, 0.3618449018307332, -0.4141404391790088], z: [-0.19282068804226837, 0.8979210213934143, 0.3956741356251556] } }, reversed: false, loops: vec![Loop { edges: vec![(0, true), (1, true), (2, false), (3, false)], vertex: None }], solid: Some(0), pcurves: vec![] }], vec![Edge { curve: Curve::Line { origin: [12.135444425482607, 11.093727383783524, 16.27115274453129], dir: [-11.572845156741018, -5.013870035444908, 5.73849825701555] }, start: [12.135444425482607, 11.093727383783524, 16.27115274453129], end: [0.5625992687415895, 6.079857348338615, 22.00965100154684], vertices: (0, 1), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [0.5625992687415895, 6.079857348338615, 22.00965100154684], dir: [5.542666214589753, -2.6969961394581725, 8.82147688477897] }, start: [0.5625992687415895, 6.079857348338615, 22.00965100154684], end: [6.105265483331342, 3.382861208880443, 30.83112788632581], vertices: (1, 4), same_sense: true, samples: vec![] },
Edge { curve: Curve::Ellipse { frame: Frame { origin: [8.015139826201194, 13.098619103296215, 9.713447413993034], x: [0.8351981583167707, 0.36184490183073315, -0.41414043917900856], y: [0.5150380749101771, -0.25061146493908637, 0.8197131663172822], z: [0.19282068804226846, -0.8979210213934147, -0.3956741356251557] }, major: 39.999999999997875, minor: 20.0 }, start: [18.315901324404734, 8.086389804514488, 26.107710740338675], end: [6.105265483331342, 3.382861208880443, 30.83112788632581], vertices: (3, 4), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [12.135444425482607, 11.093727383783524, 16.27115274453129], dir: [6.180456898922127, -3.007337579269036, 9.836557995807386] }, start: [12.135444425482607, 11.093727383783524, 16.27115274453129], end: [18.315901324404734, 8.086389804514488, 26.107710740338675], vertices: (0, 3), same_sense: true, samples: vec![] }], vec![]);
    let file = std::path::Path::new(env!(CARGO_MANIFEST_DIR)).join(../../tests/fixtures/captured/section_recess_geometry/parts/11e9e0b548aeed6c.step.gz);
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 3)]), Coverage::Covered);
}

#[test]
fn captured_cylindrical_polynomial_160a9f5d03a23b16() {
    let patch = Part::new(vec![Face { surface: Surface::Plane { frame: Frame { origin: [126.56074482422525, -50.71799438427797, 87.58269433330891], x: [0.6268924782751485, 0.7623273879120267, -0.1608190732508712], y: [-0.1537427100993845, -0.08131408698477327, -0.9847594621780131], z: [-0.7637859646349486, 0.6420630599069372, 0.06622708909201423] } }, reversed: false, loops: vec![Loop { edges: vec![(0, true), (1, true), (2, false), (3, false)], vertex: None }], solid: Some(0), pcurves: vec![] }], vec![Edge { curve: Curve::Ellipse { frame: Frame { origin: [122.44416614564037, -56.08988534512533, 92.18666734450109], x: [0.1537427100993845, 0.08131408698477326, 0.9847594621780131], y: [-0.6268924782751486, -0.7623273879120267, 0.16081907325087116], z: [0.763785964634949, -0.6420630599069375, -0.06622708909201425] }, major: 9.237604307033733, minor: 8.0 }, start: [126.56074482422525, -50.71799438427797, 87.58269433330891], end: [127.63665066041575, -49.95858754370141, 92.6285857012957], vertices: (0, 1), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [127.63665066041575, -49.95858754370141, 92.6285857012957], dir: [7.611655235913062, 9.256090086845319, -1.9526464000855412] }, start: [127.63665066041575, -49.95858754370141, 92.6285857012957], end: [135.24830589632882, -40.70249745685609, 90.67593930121016], vertices: (1, 4), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [134.44943534077242, -41.125017846942185, 85.55897903603072], dir: [0.7988705555563911, 0.4225203900860919, 5.1169602651794435] }, start: [134.44943534077242, -41.125017846942185, 85.55897903603072], end: [135.24830589632882, -40.70249745685609, 90.67593930121016], vertices: (3, 4), same_sense: true, samples: vec![] },
Edge { curve: Curve::Line { origin: [126.56074482422525, -50.71799438427797, 87.58269433330891], dir: [7.888690516547172, 9.592976537335787, -2.023715297278187] }, start: [126.56074482422525, -50.71799438427797, 87.58269433330891], end: [134.44943534077242, -41.125017846942185, 85.55897903603072], vertices: (0, 3), same_sense: true, samples: vec![] }], vec![]);
    let file = std::path::Path::new(env!(CARGO_MANIFEST_DIR)).join(../../tests/fixtures/captured/section_recesses/parts/160a9f5d03a23b16.step.gz);
    let source = haecceity::step::read_step_file(&file).unwrap();
    assert_eq!(assess((&patch, 0), &[(&source, 12)]), Coverage::Covered);
}
