//! A planar face swept along a straight vector into a solid of its own — what the Python
//! implementation builds with `Solid.extrude` / `BRepPrimAPI_MakePrism` to ask a volume probe
//! about a region shaped like a face (`kernel::volume::Probe::Solid`).
//!
//! The prism is built exactly: the face itself at both ends and one side face per edge, a plane
//! along each straight edge and a cylinder along each circular arc whose axis runs with the
//! sweep. A closed circle (a hole's whole rim) sweeps to a seamless cylinder band bounded by its
//! two copies. A selected planar spline sweeps to its tensor product with a linear segment.
//! A boundary without a supported spatial curve cannot supply a swept solid.

use super::brep::{Edge, Face, Loop, Part, Solid};
use super::geom::{self, Curve, Frame, Surface, V3};
use super::sampling::sample_edge;

/// *face* of *part* moved by *offset* and swept along *sweep*, as a one-solid part, or `None`
/// when the face is not planar, the sweep runs in its plane, or an edge is not a straight line
/// or a circle about the sweep axis or a selected spline in the cap's plane.
pub fn extrude_face(part: &Part, face: usize, offset: V3, sweep: V3) -> Option<Part> {
    let f = &part.faces[face];
    let Surface::Plane { frame } = f.surface else {
        return None;
    };
    let up = geom::unit(sweep)?;
    if geom::dot(frame.z, up).abs() <= 1e-9 {
        return None;
    }
    let top_by = geom::add(offset, sweep);
    // Each source edge once: its copies at both ends.
    let source: Vec<usize> = part.face_edges(face);
    let vertex_ids = sweepable_vertices(part, &source, frame, up)?;
    let n_vertices = vertex_ids.len();
    // Edges: bottom copies, top copies, then one rising edge per vertex.
    let mut edges = Vec::with_capacity(2 * source.len() + n_vertices);
    for (level, by) in [offset, top_by].into_iter().enumerate() {
        for &e in &source {
            edges.push(copy_edge(part, e, by, level, &vertex_ids)?);
        }
    }
    let riser = push_risers(part, &source, &vertex_ids, offset, sweep, &mut edges);
    #[allow(
        clippy::expect_used,
        reason = "the caps and sides name only the face's own edges"
    )]
    let bottom = |e: usize| source.iter().position(|&x| x == e).expect("face edge");
    let top = |e: usize| bottom(e) + source.len();

    // The face's outward normal decides which end keeps its orientation: the end the sweep
    // leaves from faces against the sweep.
    let outward = if f.reversed {
        geom::scale(frame.z, -1.0)
    } else {
        frame.z
    };
    let along = geom::dot(outward, up) > 0.0;
    let mut faces = vec![
        cap(f, frame, offset, &bottom, along),
        cap(f, frame, top_by, &top, !along),
    ];
    // Sides: each edge of the far end (which faces along the sweep, its interior on the left of
    // its walk looking down the sweep) walked back, down, along the near end and up again.
    let n = source.len();
    let far_loops: Vec<Vec<(usize, bool)>> =
        faces[1].loops.iter().map(|lp| lp.edges.clone()).collect();
    for lp in &far_loops {
        for &(t, forward) in lp {
            faces.push(side_face(&edges, t, forward, n, n_vertices, &riser, up)?);
        }
    }
    let solid = Solid {
        faces: (0..faces.len()).collect(),
    };
    let mut swept = Part::new(faces, edges, vec![solid]);
    if let Some(&outer) = part.outer_bounds.get(&face) {
        swept.outer_bounds.insert(0, outer);
        swept.outer_bounds.insert(1, outer);
    }
    Some(swept)
}

/// The vertices of the *source* edges, each once in order of first use, or `None` when an edge
/// has no supported spatial curve for its side wall.
fn sweepable_vertices(part: &Part, source: &[usize], cap: Frame, up: V3) -> Option<Vec<usize>> {
    let mut vertex_ids: Vec<usize> = Vec::new();
    for &e in source {
        let edge = &part.edges[e];
        match part.boundary_edge_curve(e)? {
            Curve::Line { .. } if !edge.is_closed() => {}
            Curve::Nurbs(n)
                if part.spatial_trim(e).is_some()
                    && n.control_points().iter().all(|&p| {
                        let roundoff =
                            64. * f64::EPSILON * geom::norm(p).max(geom::norm(cap.origin)).max(1.);
                        geom::dot(geom::sub(p, cap.origin), cap.z).abs() <= roundoff
                    }) => {}
            Curve::Circle { frame: c, .. } if 1.0 - geom::dot(c.z, up).abs() <= 1e-9 => {}
            _ => return None,
        }
        for v in [edge.vertices.0, edge.vertices.1] {
            if !vertex_ids.contains(&v) {
                vertex_ids.push(v);
            }
        }
    }
    Some(vertex_ids)
}

/// *v*'s index among the swept vertices.
#[allow(
    clippy::expect_used,
    reason = "sweepable_vertices checks the curve type and collects both ends of each edge"
)]
fn vertex_index(vertex_ids: &[usize], v: usize) -> usize {
    vertex_ids.iter().position(|&x| x == v).expect("collected")
}

/// A copy of *edge* moved by *by*, its vertices those of end *level* (0 near, 1 far).
#[allow(
    clippy::unreachable,
    reason = "sweepable_vertices checks the curve type and collects both ends of each edge"
)]
fn copy_edge(part: &Part, e: usize, by: V3, level: usize, vertex_ids: &[usize]) -> Option<Edge> {
    let edge = &part.edges[e];
    let curve = part.boundary_edge_curve(e)?;
    let (a, b) = part.boundary_edge_interval(e)?;
    let n_vertices = vertex_ids.len();
    Some(Edge {
        curve: match curve {
            Curve::Line { origin, dir } => Curve::Line {
                origin: geom::add(*origin, by),
                dir: *dir,
            },
            Curve::Circle { frame, radius } => Curve::Circle {
                frame: Frame {
                    origin: geom::add(frame.origin, by),
                    ..*frame
                },
                radius: *radius,
            },
            Curve::Nurbs(n) => {
                let mut translated = n.clone();
                translated.transform_points(|p| geom::add(p, by));
                Curve::Nurbs(translated)
            }
            _ => unreachable!("refused above"),
        },
        start: geom::add(edge.start, by),
        end: geom::add(edge.end, by),
        vertices: (
            vertex_index(vertex_ids, edge.vertices.0) + level * n_vertices,
            vertex_index(vertex_ids, edge.vertices.1) + level * n_vertices,
        ),
        same_sense: if part.spatial_trim(e).is_some() && a != b {
            b > a
        } else {
            edge.same_sense
        },
        samples: part
            .boundary_edge_samples(e)
            .iter()
            .map(|&p| geom::add(p, by))
            .collect(),
    })
}

/// Pushes onto *edges* a rising edge for each vertex an open edge ends at (a closed circle's
/// vertex has none), and returns each vertex's rising edge.
fn push_risers(
    part: &Part,
    source: &[usize],
    vertex_ids: &[usize],
    offset: V3,
    sweep: V3,
    edges: &mut Vec<Edge>,
) -> Vec<Option<usize>> {
    let n_vertices = vertex_ids.len();
    let vid = |v: usize| vertex_index(vertex_ids, v);
    let mut corner: Vec<Option<V3>> = vec![None; n_vertices];
    for &e in source {
        let edge = &part.edges[e];
        corner[vid(edge.vertices.0)] = Some(edge.start);
        corner[vid(edge.vertices.1)] = Some(edge.end);
    }
    let mut riser: Vec<Option<usize>> = vec![None; n_vertices];
    for &e in source.iter().filter(|&&e| !part.edges[e].is_closed()) {
        let edge = &part.edges[e];
        riser[vid(edge.vertices.0)] = Some(0);
        riser[vid(edge.vertices.1)] = Some(0);
    }
    for (i, p) in corner.iter().enumerate() {
        if riser[i].is_none() {
            continue;
        }
        riser[i] = Some(edges.len());
        // Every swept vertex is an end of an open edge, so it has its corner.
        let Some(p) = p else {
            continue;
        };
        let start = geom::add(*p, offset);
        let end = geom::add(start, sweep);
        let curve = Curve::Line {
            origin: start,
            dir: sweep,
        };
        edges.push(Edge {
            samples: sample_edge(&curve, start, end, true, false),
            curve,
            start,
            end,
            vertices: (i, i + n_vertices),
            same_sense: true,
        });
    }
    riser
}

/// An end of the prism: face *f* (in the plane of *frame*) moved by *by*, its edges renumbered
/// by *edge_of*, and turned over when *flip*.
fn cap(f: &Face, frame: Frame, by: V3, edge_of: &dyn Fn(usize) -> usize, flip: bool) -> Face {
    Face {
        surface: Surface::Plane {
            frame: Frame {
                origin: geom::add(frame.origin, by),
                ..frame
            },
        },
        reversed: f.reversed != flip,
        loops: f
            .loops
            .iter()
            .map(|lp| {
                let mut edges: Vec<(usize, bool)> = lp
                    .edges
                    .iter()
                    .map(|&(e, fw)| (edge_of(e), fw != flip))
                    .collect();
                if flip {
                    edges.reverse();
                }
                Loop {
                    edges,
                    vertex: None,
                }
            })
            .collect(),
        solid: Some(0),
        pcurves: Vec::new(),
    }
}

/// The side swept by far-end edge *t* (used *forward* or not by the far end), whose near-end
/// copy is `t - n`; `None` when its surface's orientation cannot be found.
fn side_face(
    edges: &[Edge],
    t: usize,
    forward: bool,
    n: usize,
    n_vertices: usize,
    riser: &[Option<usize>],
    up: V3,
) -> Option<Face> {
    let edge = &edges[t];
    let (p, q) = if forward {
        (edge.vertices.0, edge.vertices.1)
    } else {
        (edge.vertices.1, edge.vertices.0)
    };
    let (p, q) = (p - n_vertices, q - n_vertices);
    // A closed circle's band is bounded by its two copies as separate loops.
    let side_loops = match (riser[p], riser[q]) {
        (Some(down), Some(up)) if !edge.is_closed() => {
            vec![vec![
                (t, !forward),
                (down, false),
                (t - n, forward),
                (up, true),
            ]]
        }
        _ => vec![vec![(t, !forward)], vec![(t - n, forward)]],
    };
    // The side's outward normal: the walk's tangent across the sweep.
    let samples = &edge.samples;
    let k = samples.len() / 2;
    let (a, c) = (
        samples[k.saturating_sub(1).min(samples.len() - 2)],
        samples[k.max(1)],
    );
    let walk = if forward {
        geom::sub(c, a)
    } else {
        geom::sub(a, c)
    };
    let mid = geom::scale(geom::add(a, c), 0.5);
    let out_normal = geom::unit(geom::cross(walk, up))?;
    let surface = match &edge.curve {
        Curve::Line { .. } => {
            let x = geom::unit(geom::sub(edge.end, edge.start))?;
            let z = out_normal;
            Surface::Plane {
                frame: Frame {
                    origin: edge.start,
                    x,
                    y: geom::cross(z, x),
                    z,
                },
            }
        }
        Curve::Circle { frame, radius } => Surface::Cylinder {
            frame: *frame,
            radius: *radius,
        },
        Curve::Nurbs(top) => {
            let Curve::Nurbs(bottom) = &edges[t - n].curve else {
                return None;
            };
            let poles = bottom
                .control_points()
                .iter()
                .zip(top.control_points())
                .map(|(&a, &b)| vec![a, b])
                .collect();
            let weights = top.weights().iter().map(|&w| vec![w, w]).collect();
            Surface::Freeform {
                kind: "EXTRUSION",
                surface: Box::new(crate::nurbs::NurbsSurface::new(
                    top.degree(),
                    1,
                    poles,
                    weights,
                    top.knots().to_vec(),
                    vec![0., 0., 1., 1.],
                )?),
            }
        }
        // sweepable_vertices refuses every other curve.
        _ => return None,
    };
    let (u, v) = surface.parameters(mid, None)?;
    let natural = surface.normal(u, v)?;
    Some(Face {
        surface,
        reversed: geom::dot(natural, out_normal) < 0.0,
        loops: side_loops
            .into_iter()
            .map(|edges| Loop {
                edges,
                vertex: None,
            })
            .collect(),
        solid: Some(0),
        pcurves: Vec::new(),
    })
}
