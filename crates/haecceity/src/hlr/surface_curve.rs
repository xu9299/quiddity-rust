//! Exact orthographic image of the surface curve selected by a STEP edge.
use super::View;
use crate::geom::{self, Curve, Surface};
use crate::mass::DrawingPcurve;
use std::sync::Arc;

/// A finite 2D stretch obtained by evaluating its pcurve on its named surface
/// and projecting that point. It owns its geometry and can outlive the part.
/// With the `json` feature the operands and parameter range can be serialized.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "json", derive(serde::Serialize, serde::Deserialize))]
pub struct ProjectedSurfaceCurve {
    geometry: Arc<DrawingPcurve>,
    axes: [[f64; 3]; 2],
    offset: [f64; 2],
    interval: (f64, f64),
}

impl ProjectedSurfaceCurve {
    pub(super) fn new(geometry: Arc<DrawingPcurve>, view: &View, interval: (f64, f64)) -> Self {
        Self {
            geometry,
            axes: [view.x, view.y],
            offset: [0.; 2],
            interval,
        }
    }

    /// Original parameter-space curve. The carrier parameters for evaluation
    /// include the STEP units and any swept-surface reparameterization.
    pub fn pcurve(&self) -> &Curve {
        &self.geometry.pcurve
    }

    /// The surface named by the selected pcurve, in file coordinates.
    pub fn carrier(&self) -> &Surface {
        &self.geometry.carrier
    }

    /// Parameters on the evaluated carrier at the selected curve's parameter.
    pub fn surface_parameters(&self, t: f64) -> (f64, f64) {
        self.geometry.uv(t)
    }

    /// The drawing parameter at the first and last points, in traversal order.
    pub fn range(&self) -> (f64, f64) {
        self.interval
    }

    pub(super) fn set_range(&mut self, ends: (f64, f64)) {
        self.interval = ends;
    }

    /// Finite parameter spans in drawing order, including the end parameters
    /// and interior breaks of the selected curve. A vector writer can begin
    /// its error-controlled subdivision at these parameters.
    pub fn span_parameters(&self) -> Vec<f64> {
        let (a, b) = self.range();
        let mut spans = vec![a, b];
        spans.extend(
            self.geometry
                .breaks
                .iter()
                .copied()
                .filter(|&t| a.min(b) < t && t < a.max(b)),
        );
        spans.sort_by(f64::total_cmp);
        spans.dedup();
        if b < a {
            spans.reverse();
        }
        spans
    }

    /// Evaluate the exact projected point, including any page placement.
    pub fn at(&self, t: f64) -> [f64; 2] {
        let (u, v) = self.surface_parameters(t);
        let p = self.carrier().value(u, v);
        std::array::from_fn(|i| geom::dot(p, self.axes[i]) + self.offset[i])
    }

    /// Return the same curve in page coordinates, scaled about `centre`
    /// and translated so that `centre` maps to `at`.
    pub fn placed(&self, scale: f64, centre: [f64; 2], at: [f64; 2]) -> Self {
        Self {
            geometry: Arc::clone(&self.geometry),
            axes: self.axes.map(|v| geom::scale(v, scale)),
            offset: std::array::from_fn(|i| at[i] + scale * (self.offset[i] - centre[i])),
            interval: self.interval,
        }
    }
}

/// A drawing cache retains the evaluated carrier and the STEP kind of a spline.
#[cfg(feature = "json")]
pub(crate) mod carrier_wire {
    use crate::geom::{Frame, Surface};
    use crate::nurbs::NurbsSurface;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    type SplineKind = &'static str;

    #[derive(Serialize, Deserialize)]
    #[serde(remote = "Surface")]
    enum Data {
        Plane {
            frame: Frame,
        },
        Cylinder {
            frame: Frame,
            radius: f64,
        },
        Cone {
            frame: Frame,
            radius: f64,
            semi_angle: f64,
        },
        Sphere {
            frame: Frame,
            radius: f64,
        },
        Torus {
            frame: Frame,
            major: f64,
            minor: f64,
        },
        Freeform {
            #[serde(deserialize_with = "spline_kind")]
            kind: SplineKind,
            surface: Box<NurbsSurface>,
        },
        #[serde(skip)]
        Other {
            #[serde(skip)]
            kind: SplineKind,
        },
    }

    fn spline_kind<'de, D: Deserializer<'de>>(d: D) -> Result<SplineKind, D::Error> {
        const KINDS: &[SplineKind] = &["BSPLINE", "BEZIER", "EXTRUSION", "REVOLUTION"];
        let name = String::deserialize(d)?;
        KINDS
            .iter()
            .copied()
            .find(|kind| *kind == name)
            .ok_or_else(|| serde::de::Error::unknown_variant(&name, KINDS))
    }

    pub(crate) fn serialize<S: Serializer>(s: &Surface, serializer: S) -> Result<S::Ok, S::Error> {
        Data::serialize(s, serializer)
    }
    pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
        d: D,
    ) -> Result<std::sync::Arc<Surface>, D::Error> {
        Data::deserialize(d).map(std::sync::Arc::new)
    }
}

#[cfg(all(test, feature = "json"))]
mod tests {
    use super::ProjectedSurfaceCurve;
    use serde_json::json;

    #[test]
    fn stored_spline_on_spline_retains_its_geometry_and_checks_knot_arrays() {
        // The bilinear patch (u,v,u*v), on the pcurve (t,t), projected along y.
        let axis = json!({"Linear":{"scale":1.,"offset":0.}});
        let data = json!({
            "geometry":{
                "pcurve":{"Nurbs":{"degree":1, "control_points":[[0.,0.,0.],[1.,1.,0.]],
                    "weights":[1.,1.], "knots":[0.,0.,1.,1.]}},
                "carrier":{"Freeform":{"kind":"BEZIER","surface":{"degree_u":1,"degree_v":1,
                    "control_points":[[[0.,0.,0.],[0.,1.,0.]],[[1.,0.,0.],[1.,1.,1.]]],
                    "weights":[[1.,1.],[1.,1.]],"knots_u":[0.,0.,1.,1.],"knots_v":[0.,0.,1.,1.]}}},
                "axes":[axis,axis], "spans":[null,null],"breaks":[0.,1.]
            },
            "axes":[[1.,0.,0.],[0.,0.,1.]], "offset":[0.,0.], "interval":[1.,0.]
        });
        let curve: ProjectedSurfaceCurve = serde_json::from_value(data.clone()).unwrap();
        let roundtrip: ProjectedSurfaceCurve =
            serde_json::from_str(&serde_json::to_string(&curve).unwrap()).unwrap();
        assert_eq!(roundtrip, curve);
        assert!(matches!(
            roundtrip.carrier(),
            crate::geom::Surface::Freeform { kind: "BEZIER", .. }
        ));
        assert_eq!(roundtrip.range(), (1., 0.));
        assert_eq!(roundtrip.span_parameters(), vec![1., 0.]);
        for t in [0., 0.25, 0.5, 0.75, 1.] {
            assert_eq!(roundtrip.at(t), [t, t * t]);
        }
        let mut invalid = data.clone();
        invalid["geometry"]["pcurve"]["Nurbs"]["knots"] = json!([0., 1.]);
        assert!(serde_json::from_value::<ProjectedSurfaceCurve>(invalid).is_err());
        let mut invalid = data;
        invalid["geometry"]["carrier"]["Freeform"]["surface"]["weights"] = json!([[1., 1.]]);
        assert!(serde_json::from_value::<ProjectedSurfaceCurve>(invalid).is_err());
    }
}
