//
// sim.rs
// Where the simulation itself takes place
//

extern crate rand;

use fiber::{Geometry3D,Torus,Cylinder,Cone};
use math::{Point,Vector,Plane};
use rand::Rng;
use std::sync::Arc;
use std::thread;
use std::f64;

struct Wire {
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
    center: Box<Geometry3D + Send + Sync>,
}

enum Progress {
    Nothing,
    Advance,
    Regress,
    Collide,
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

        let center: Box<Geometry3D + Send + Sync>;
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

        println!("{:?}", layer2);
        println!("{:?}", layer3);

        let ent_norm = layer2.normal(entry_centerpoint);
        let exi_norm = layer3.normal(exit_centerpoint);

        println!("normal to plane entry: {:?}", ent_norm);
        println!("normal to plane exit: {:?}", exi_norm);

        //
        // Cone construction
        //

        //Parallel to the entrance Cone
        let mut cone_par_entry = Vector::new_from_points(entry_centerpoint, entry_sidepoint);
        cone_par_entry.normalize();
        //Parallel to the exit Cone
        let mut cone_par_exit = Vector::new_from_points(exit_centerpoint, exit_sidepoint);
        cone_par_exit.normalize();
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
        let cone_fside_exi = Point{x: cone_far_exi.x+cone_par_exit.x*r1,
                                    y: cone_far_exi.y+cone_par_exit.y*r1,
                                    z: cone_far_exi.z+cone_par_exit.z*r1};
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

    /// An instance of a simulation
    /// Each thread runs this, and it is divided into 5 stages, defined by the planes.
    /// Starting from stage zero, a stage is advanced whenever the ray passes through a layer,
    /// in the order that they appear, but degrades if the ray passes through a layer it has
    /// already went through before.
    pub fn simulate(&self) -> Option<i32> {
        let mut stage = 0;
        let mut hits = 0; //Number of hits on inner tube
        let mut rng = rand::thread_rng();
        let w = &self;
        let start_point = w.entry.c.point_at_t_equals(w.entry.length);
        //At the current version, the positive starting direction is guaranteed to be y+.
        //There might be better ways to ensure that the random starting vector is correctly oriented.
        let mut ray = Vector{p: start_point,
                        x: rng.gen::<f64>(),
                        y: f64::abs(rng.gen::<f64>()),
                        z: rng.gen::<f64>()};
        ray.normalize();
        let mut prog: Progress;

        println!("The start point is :{:?}", start_point);
        println!("The entry's data is: {:?}", w.entry);

        while stage<5 {

            println!("Ray location: {:?}", ray);
            println!("At stage: {:?}", stage);

            match stage {
                -1 => {
                    //We passed through the entrance so we went back
                    return None;
                }
                0 => {
                    prog = Wire::next_location(w.layer_entry, &w.entry, w.layer1, &mut ray);
                },
                1 => {
                    prog = Wire::next_location(w.layer1, &w.entry_connector, w.layer2, &mut ray);
                },
                2 => {
                    prog = Wire::next_location(w.layer2, &(*w.center), w.layer3, &mut ray);
                },
                3 => {
                    prog = Wire::next_location(w.layer3, &w.exit_connector, w.layer4, &mut ray);
                },
                4 => {
                    prog = Wire::next_location(w.layer4, &w.exit, w.layer_exit, &mut ray);
                },
                _ => {
                    panic!("Invalid ray location. The stage reported is {:?}", stage);
                },
            }
            match prog {
                Progress::Nothing => {
                //The ray is in limbo
                    return None;
                },
                Progress::Advance => {
                    stage += 1;
                },
                Progress::Regress => {
                    stage -= 1;
                },
                Progress::Collide => {
                    if stage == 2 {
                        hits += 1;
                    }
                },
            }
        }
        if stage==5 {
            return Some(hits);
        }
        return None;
    }

    /// Moves the ray to its next location.
    /// Also changes its direction if there is a collision.
    fn next_location(layer_back: Plane, part: &Geometry3D, layer_front: Plane, ray: &mut Vector) -> Progress {
        //Negative values are illegal, so if it persists, there is an error.
        let mut t = -1.;
        let mut res = Progress::Nothing;
        let mut new_point = Point::new_blank();
        let backward_move = layer_back.intersect(*ray);
        let forward_move = layer_front.intersect(*ray);
        let hit_geometry = part.collision_point(*ray);
        match backward_move {
            Some(point) => {
                t=Vector::new_from_points(ray.p, point).length();
                new_point = point;
                res = Progress::Regress;
            }
            None => {},
        }
        match forward_move {
            Some(point) => {
                let tc=Vector::new_from_points(ray.p, point).length();
                if t==-1. {
                    t = tc;
                    new_point = point;
                    res = Progress::Advance;
                }
                else if tc<t {
                    t = tc;
                    new_point = point;
                    res = Progress::Advance;
                }
            }
            None => {},
        }
        match hit_geometry {
            Some(point) => {
                let tc=Vector::new_from_points(ray.p, point).length();
                if t==-1. {
                    //t = tc; //Apparently this value is never used, but I guess it is only for comparison...
                    new_point = point;
                    res = Progress::Collide;
                }
                else if tc<t {
                    //t = tc;
                    new_point = point;
                    res = Progress::Collide;
                }
            }
            None => {},
        }
        match res {
            Progress::Collide => {
                match part.normal(new_point) {
                    Some(normal) => {
                        let nray = ray.reflection(&normal, new_point);
                        ray.p = nray.p;
                        ray.x = nray.x;
                        ray.y = nray.y;
                        ray.z = nray.z;
                    },
                    None => {
                        panic!("Unable to find Normal vector at ray collision point.");
                    }
                }
            },
            _ => {
                ray.p=new_point;
            },
        }
        return res;
    }
}

pub fn simulation(out_len: f64, in_len: f64, cone_len: f64, n1: f64, n2: f64, out_rad: f64,
                    in_rad: f64, angle: f64, thread_amount: i32, output_to_console: bool)
                    -> Vec<Option<i32>> {
    let w = Wire::generate(out_len,in_len,cone_len,n1,n2,out_rad,in_rad,angle);
    let mut results : Vec<Option<i32>> = vec![];
    let mut children = vec![];
    for _ in 0..thread_amount {
        let wc = w.clone();
        children.push(thread::spawn(move || {
            wc.simulate()
        }));
    }
    for child in children {
        match child.join() {
            Ok(result) => {
                results.push(result);
                if output_to_console {
                    match result {
                        Some(amount) => {
                            println!("Success with {:?} hits!", amount);
                        }
                        None => {
                            println!("Failure occurred.");
                        }
                    }
                }
            },
            Err(_) => {
                results.push(None);
            },
        }
    }
    results
}
