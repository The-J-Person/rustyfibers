//
// sim.rs
// Where the simulation itself takes place
//

use fiber::{Geometry3D,Torus,Cylinder,Cone};
use math::{Point,Vector,Plane};
use std::sync::Arc;
use std::f64;

#[allow(dead_code)]
pub struct Wire {
    n1: f64,
    n2: f64,
    layer_entry: Plane,
    entry: Cylinder,
    layer1: Plane,
    entry_connector: Cone,
    layer2: Plane,
    layer3: Plane,
    exit_connector: Cone,
    layer4: Plane,
    exit: Cylinder,
    layer_exit: Plane,
    center: Box<Geometry3D>,
}

impl Wire {
    /// Generates a new Wire
    /// A -     length of outer tubes
    /// B -     length of inner tube
    /// C -     length of the interconnecting tube-cones
    /// n1 -    outer index of refraction
    /// n2 -    inner index of refraction
    /// r1 -    outer tube radius
    /// r2 -    inner tube radius
    /// angle - angle of the inner tube's bend, between 0 and 180 degrees.
    ///         At an angle of 0 degrees, the inner tube is a cylinder. Otherwise, it is a torus.
    #[allow(non_snake_case)]
    pub fn generate(A: f64, B: f64, C: f64,
                    n1: f64, n2: f64, r1: f64, r2: f64, angle: f64) ->  Arc<Wire>{

        let center: Box<Geometry3D>;
        let entry_centerpoint: Point;
        let exit_centerpoint: Point;
        let entry_sidepoint: Point;
        let exit_sidepoint: Point;

        //
        // Centerpiece construction
        //

        if angle>180.+f64::EPSILON {
            panic!("This software does not handle angles over 180 degrees.");
        }
        else if angle<0. {
            panic!("This software does not handle negative angles.");
        }
        else if angle==0. {
            //central tube is a cylinder
            //we don't really care how we orient it, so it'll be the y-axis and centered on 0,0,0
            center = Box::new(Cylinder{c: Vector{p: Point{x: 0.,y: -B/2.,z: 0.},x: 0., y: 1., z: 0.}
                ,r: r2, length: B});
            entry_centerpoint = Point{x: 0., y: -B/2., z: 0.};
            entry_sidepoint = Point{x: r2, y: -B/2., z: 0.};
            exit_centerpoint = Point{x: 0., y: B/2., z: 0.};
            exit_sidepoint = Point{x: r2, y: B/2., z: 0.};
        }
        else {
            //central tube is a torus
            let phi_max = angle*f64::consts::PI/180.;//convert to radians
            //phi*R = length of tube, therefore
            let R = B/phi_max;
            center = Box::new(Torus{R: R, r: r2, phi_max: phi_max});
            entry_centerpoint = Point{x: R, y: 0., z: 0.};
            entry_sidepoint = Point{x: R-r2, y: 0., z: 0.};
            exit_centerpoint = Point{x: R*f64::cos(phi_max), y: R*f64::sin(phi_max), z: 0.};
            exit_sidepoint = Point{x: (R-r2)*f64::cos(phi_max), y: (R-r2)*f64::sin(phi_max), z: 0.};
        }
        let layer2 = (*center).entry_plane();
        let layer3 = (*center).exit_plane();

        let ent_norm = layer2.normal(entry_centerpoint);
        let exi_norm = layer3.normal(exit_centerpoint);

        //
        // Cone construction
        //

        //Parallel to the entrance Cone
        let cone_par_entry = Vector::new_from_points(entry_centerpoint, entry_sidepoint);
        //Parallel to the exit Cone
        let cone_par_exit = Vector::new_from_points(exit_centerpoint, exit_sidepoint);
        //Center of larger cone circle, entrance side
        let cone_far_ent = Point{x: entry_centerpoint.x+ent_norm.x*C,
                                    y: entry_centerpoint.y+ent_norm.y*C,
                                    z: entry_centerpoint.z+ent_norm.z*C};
        //Center of larger cone circle, exit side
        let cone_far_exi = Point{x: exit_centerpoint.x+exi_norm.x*C,
                                    y: exit_centerpoint.y+exi_norm.y*C,
                                    z: exit_centerpoint.z+exi_norm.z*C};
        //Side of larger cone circle, entrance side
        let cone_fside_ent = Point{x: cone_far_ent.x+cone_par_entry.x*r1,
                                    y: cone_far_ent.y+cone_par_entry.y*r1,
                                    z: cone_far_ent.z+cone_par_entry.z*r1};
        //Side of larger cone circle, exit side
        let cone_fside_exi = Point{x: cone_far_ent.x+cone_par_exit.x*r1,
                                    y: cone_far_ent.y+cone_par_exit.y*r1,
                                    z: cone_far_ent.z+cone_par_exit.z*r1};
        //Cone side vectors (the central vector is the plane's normal)
        let cone_svec_ent = Vector::new_from_points(entry_sidepoint, cone_fside_ent);
        let cone_svec_exi = Vector::new_from_points(exit_sidepoint, cone_fside_exi);
        //Angles of Cones (although, they should have the same angle...)
        let alpha_ent = f64::acos(ent_norm.dot_product(&cone_svec_ent)
                                /(ent_norm.length()*cone_svec_ent.length()));
        let alpha_exi = f64::acos(exi_norm.dot_product(&cone_svec_exi)
                                /(exi_norm.length()*cone_svec_exi.length()));
        //Apex points of Cones
        let apex_ent = ent_norm.intersect(&cone_svec_ent);
        let apex_exi = exi_norm.intersect(&cone_svec_exi);
        let mut entry_conevec = ent_norm;
        entry_conevec.p = apex_ent;
        let mut exit_conevec = exi_norm;
        exit_conevec.p = apex_exi;

        //Creating the Cones themselves
        let entry_connector = Cone{c: entry_conevec ,a: alpha_ent ,length: C};
        let exit_connector = Cone{c: exit_conevec ,a: alpha_exi ,length: C};

        //
        // Cylinder Construction
        //
        let mut entry_cvec = ent_norm;
        entry_cvec.p = cone_far_ent;
        let entry = Cylinder{c: entry_cvec,r: r1, length: A};
        let mut exit_cvec = exi_norm;
        exit_cvec.p = cone_far_exi;
        let exit = Cylinder{c: exit_cvec,r: r1, length: A};

        //restraint handling

        Arc::new(Wire{
            n1: n1,
            n2: n2,
            layer_entry: entry.exit_plane(), //This might not be intuitive:
            entry: entry,                        //Entry implies facing the inner area.
            layer1: entry.entry_plane(),
            entry_connector: entry_connector,
            layer2: layer2,
            layer3: layer3,
            exit_connector: exit_connector,
            layer4: exit.entry_plane(),
            exit: exit,
            layer_exit: exit.exit_plane(),
            center: center,})
    }
}
