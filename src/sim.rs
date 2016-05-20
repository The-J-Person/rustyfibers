//
// sim.rs
// Where the simulation itself takes place
//

use fiber::{Geometry3D,Torus,Cylinder,Cone};
use math::{Point,Vector,Plane};

pub struct Wire {
    n1: f64,
    n2: f64,
    layer_entry: Plane,
    entry: Cylinder,
    layer1: Plane,
    connector_a: Cone,
    layer2: Plane,
    layer3: Plane,
    connector_b: Cone,
    layer4: Plane,
    exit: Cylinder,
    layer_exit: Plane,
    center: Geometry3D,
}

impl Wire {
    /// Generates a new Wire
    /// A -     length of external tubes
    /// B -     length of internal tube
    /// C -     length of the interconnecting tube-cones
    /// n1 -    outer index of refraction
    /// n2 -    inner index of refraction
    /// angle - angle of the inner tube's bend, between 0 and 180 degrees.
    ///         At an angle of 0 degrees, the inner tube is a cylinder. Otherwise, it is a torus.
    pub fn generate(A: f64, B: f64, C: f64, n1: f64, n2: f64, angle: f64){

    }
}
