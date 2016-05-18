//
// fiber.rs
// Defines all structs related to the structure of the fiber
//

extern crate roots;
use roots::Roots;
use roots::find_roots_quartic;

use math::Point;
use math::Vector;
use std::f64;

trait Geometry3D {
    fn point_is_on_surface(&self,p: Point) -> bool;
    fn check_collision(&self,v: Vector) -> bool;
    fn collision_point(&self,v: Vector) -> Option<Point>;
    fn normal(&self,p: Point) -> Option<Vector>;
}



//
// Torus type
// Describes a part of a Torus, anything between a 2D circle and half a torus.
//
    // Friendly reminder, these are point coordinates for a torus:
    // x = (R+r*f64::cos(theta))*f64::cos(phi);
    // y = (R+r*f64::cos(theta))*f64::sin(phi);
    // z = r*f64::sin(theta);
    //
#[derive(PartialEq, Copy, Clone, Debug)]
#[allow(non_snake_case)]
pub struct Torus {
    pub R: f64, // "Outer Radius", distance from torus center to a mid-point in the tube
    pub r: f64, // Radius of tube
    //pub theta: f64, //Angle, revolves around tube
    pub phi_max: f64, //Angle, revolves around whole toroid alongside the tube - maximum value
}

impl Torus{
    //Returns theta and phi angle values if the given point is on the torus surface
    //If the point is not on the torus surface, returns "None"
    fn theta_phi(&self,p: Point) -> Option<(f64,f64)> {
        let theta1 = f64::asin(p.z); //Possible solution for theta
        let theta2 = f64::consts::PI - theta1; //Another possible solution for theta
        let phi1 = f64::acos(p.x/(self.R+self.r*f64::cos(theta1))); //Phi for Theta 1
        let phi2 = f64::acos(p.x/(self.R+self.r*f64::cos(theta2))); //Phi for Theta 2
        if phi1<=self.phi_max && p.y == (self.R+self.r*f64::cos(theta1))*f64::sin(phi1) { //Check if y upholds for theta&phi set #1
            return Some((theta1,phi1));
        }
        else if phi2<=self.phi_max && p.y == (self.R+self.r*f64::cos(theta2))*f64::sin(phi2) { //Check if y upholds for theta&phi set #2
            return Some((theta2,phi2));
        }
        None
    }
}

impl Geometry3D for Torus{
    fn point_is_on_surface(&self,p: Point) -> bool {
        if self.theta_phi(p) == None {
            return false;
        }
        true
    }
    fn check_collision(&self,v: Vector) -> bool {
        unimplemented!()
    }
    #[allow(non_snake_case)]
    fn collision_point(&self,v: Vector) -> Option<Point> {
        let R2 = self.R*self.R;
        let a = Vector{p: Point{x: 0.,y: 0.,z: 0.},x: v.p.x, y: v.p.y, z: v.p.z};
        let K = a.dot_product(&a)-self.r*self.r-R2;
        let A = 4.*a.dot_product(&v);
        let B = 2.*(2.*a.dot_product(&v)*a.dot_product(&v)+K+2.*R2*v.z*v.z);
        let C = 4.*(K*a.dot_product(&v)+2.*R2*a.z*v.z);
        let D = K*K+4.*R2*(a.z*a.z-self.r*self.r);
        let t: f64;
        match find_roots_quartic(1., A, B, C, D) {
            Roots::Four(roots) => {
                let mut temp = f64::INFINITY;
                if roots[0]>0. && roots[0]<temp {
                    temp=roots[0];
                }
                if roots[1]>0. && roots[1]<temp {
                    temp=roots[1];
                }
                if roots[2]>0. && roots[2]<temp {
                    temp=roots[2];
                }
                if roots[3]>0. && roots[3]<temp {
                    temp=roots[3];
                }
                if temp==f64::INFINITY {
                    return None;
                }
                else {
                    t = temp;
                }
            },
            Roots::Three(roots) => {
                let mut temp = f64::INFINITY;
                if roots[0]>0. && roots[0]<temp {
                    temp=roots[0];
                }
                if roots[1]>0. && roots[1]<temp {
                    temp=roots[1];
                }
                if roots[2]>0. && roots[2]<temp {
                    temp=roots[2];
                }
                if temp==f64::INFINITY {
                    return None;
                }
                else {
                    t = temp;
                }
            },
            Roots::Two(roots) => {
                let mut temp = f64::INFINITY;
                if roots[0]>0. && roots[0]<temp {
                    temp=roots[0];
                }
                if roots[1]>0. && roots[1]<temp {
                    temp=roots[1];
                }
                if temp==f64::INFINITY {
                    return None;
                }
                else {
                    t = temp;
                }
            },
            Roots::One(roots) => {
                if roots[0]>0. {
                    t = roots[0];
                }
                else {
                    return None;
                }
            },
            _ => {
                return None;
            },
        }
        return Some(Point {x: v.p.x+v.x*t,y: v.p.y+v.y*t,z: v.p.z+v.z*t});
    }
    fn normal(&self,p: Point) -> Option<Vector> {
        let theta: f64;
        let phi: f64;
        match self.theta_phi(p) {
            Some(pair) => {
                theta = pair.0;
                phi = pair.1;
            },
            None => return None,
        }
        // tangent position with respect to R
        let tx = -f64::sin(phi);
        let ty = f64::cos(phi);
        let tz = 0.;
        // tangent position with respect to r
        let sx = f64::cos(phi)*(-f64::sin(theta));
        let sy = f64::sin(phi)*(-f64::sin(theta));
        let sz = f64::cos(theta);
        // cross_product
        let nx = ty*sz - tz*sy;
        let ny = tz*sx - tx*sz;
        let nz = tx*sy - ty*sx;
        // convert to vector, normalize, and return
        let mut v = Vector {p: p, x: nx, y: ny, z: nz};
        v.normalize();
        Some(v)
    }
}


//
// Cylinder typic
//

#[derive(PartialEq, Copy, Clone, Debug)]
struct Cylinder {
    c: Vector,  //"central" vector, alternatively, the "spine" along which the cylinder stretches.
    r: f64,     //the radius
}

impl Geometry3D for Cylinder{
    fn point_is_on_surface(&self,p: Point) -> bool {
        //let x = ( p.x - self.c.p.x - (self.c.x*p.x-self.c.p.x)*self.c.x).powi(2)-self.r.powi(2);
        unimplemented!()
    }
    fn check_collision(&self,v: Vector) -> bool {
        unimplemented!()
    }
    fn collision_point(&self,v: Vector) -> Option<Point> {
        unimplemented!()
    }
    fn normal(&self,p: Point) -> Option<Vector> {
        unimplemented!()
    }
}

//
// Cone typic
//

#[derive(PartialEq, Copy, Clone, Debug)]
struct Cone {
    c: Vector,  //"central" vector, alternatively, the "spine" along which the cylinder stretches.
    a: f64,     //Angle between central vector and cone surface
}

impl Geometry3D for Cone{
    fn point_is_on_surface(&self,p: Point) -> bool {
        unimplemented!()
    }
    fn check_collision(&self,v: Vector) -> bool {
        unimplemented!()
    }
    fn collision_point(&self,v: Vector) -> Option<Point> {
        unimplemented!()
    }
    fn normal(&self,p: Point) -> Option<Vector> {
        unimplemented!()
    }
}
