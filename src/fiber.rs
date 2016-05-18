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
        //This function is disorderly and can probably be simplified.
        //Somehow.
        let t: f64;
        let a = v.p.x;
        let b = self.c.p.x;
        let c = v.x;
        let g = self.c.x;
        let h = v.p.y;
        let k = self.c.p.y;
        let m = v.y;
        let n = self.c.y;
        let p = v.p.z;
        let r = self.c.p.z;
        let s = v.z;
        let w = self.c.z;
        let x = f64::cos(self.a);
        let to_root = f64::sqrt((2.*a*c*g.powi(2)-2.*b*c*g.powi(2)-2.*a*c*x*g.powi(2)
                                +2.*b*c*x*g.powi(2)-2.*h*m*x*g.powi(2)+2.*k*m*x*g.powi(2)
                                -2.*p*s*x*g.powi(2)+2.*r*s*x*g.powi(2)+2.*c*h*n*g-2.*c*k*n*g
                                +2.*a*m*n*g-2.*b*m*n*g+2.*c*p*w*g-2.*c*r*w*g+2.*a*s*w*g
                                -2.*b*s*w*g+2.*h*m*n.powi(2)-2.*k*m*n.powi(2)+2.*p*s*w.powi(2)
                                -2.*r*s*w.powi(2)+2.*m*n*p*w-2.*m*n*r*w+2.*h*n*s*w-2.*k*n*s*w
                                -2.*a*c*n.powi(2)*x+2.*b*c*n.powi(2)*x-2.*h*m*n.powi(2)*x
                                +2.*k*m*n.powi(2)*x-2.*a*c*w.powi(2)*x+2.*b*c*w.powi(2)*x
                                -2.*h*m*w.powi(2)*x+2.*k*m*w.powi(2)*x-2.*p*s*w.powi(2)*x
                                +2.*r*s*w.powi(2)*x-2.*n.powi(2)*p*s*x+2.*n.powi(2)*r*s*x).powi(2)
                                -4.*(a.powi(2)*g.powi(2)+b.powi(2)*g.powi(2)-2.*a*b*g.powi(2)
                                -a.powi(2)*x*g.powi(2)-b.powi(2)*x*g.powi(2)-h.powi(2)*x*g.powi(2)
                                -k.powi(2)*x*g.powi(2)-p.powi(2)*x*g.powi(2)-r.powi(2)*x*g.powi(2)
                                +2.*a*b*x*g.powi(2)+2.*h*k*x*g.powi(2)+2.*p*r*x*g.powi(2)
                                +2.*a*h*n*g-2.*b*h*n*g-2.*a*k*n*g+2.*b*k*n*g+2.*a*p*w*g-2.*b*p*w*g
                                -2.*a*r*w*g+2.*b*r*w*g+h.powi(2)*n.powi(2)+k.powi(2)*n.powi(2)
                                -2.*h*k*n.powi(2)+p.powi(2)*w.powi(2)+r.powi(2)*w.powi(2)
                                -2.*p*r*w.powi(2)+2.*h*n*p*w-2.*k*n*p*w-2.*h*n*r*w+2.*k*n*r*w
                                -a.powi(2)*n.powi(2)*x-b.powi(2)*n.powi(2)*x-h.powi(2)*n.powi(2)*x
                                -k.powi(2)*n.powi(2)*x+2.*a*b*n.powi(2)*x+2.*h*k*n.powi(2)*x
                                -n.powi(2)*p.powi(2)*x-n.powi(2)*r.powi(2)*x-a.powi(2)*w.powi(2)*x
                                -b.powi(2)*w.powi(2)*x-h.powi(2)*w.powi(2)*x-k.powi(2)*w.powi(2)*x
                                -p.powi(2)*w.powi(2)*x-r.powi(2)*w.powi(2)*x+2.*a*b*w.powi(2)*x
                                +2.*h*k*w.powi(2)*x+2.*p*r*w.powi(2)*x
                                +2.*n.powi(2)*p*r*x)*(c.powi(2)*g.powi(2)-c.powi(2)*x*g.powi(2)
                                -m.powi(2)*x*g.powi(2)-s.powi(2)*x*g.powi(2)+2.*c*m*n*g+2.*c*s*w*g
                                +m.powi(2)*n.powi(2)+s.powi(2)*w.powi(2)
                                +2.*m*n*s*w-c.powi(2)*n.powi(2)*x-m.powi(2)*n.powi(2)*x
                                -n.powi(2)*s.powi(2)*x-c.powi(2)*w.powi(2)*x-m.powi(2)*w.powi(2)*x
                                -s.powi(2)*w.powi(2)*x));
        let top = -2.*a*c*g.powi(2)+2.*b*c*g.powi(2)+2.*a*c*x*g.powi(2)-2.*b*c*x*g.powi(2)
                    +2.*h*m*x*g.powi(2)-2.*k*m*x*g.powi(2)+2.*p*s*x*g.powi(2)-2.*r*s*x*g.powi(2)
                    -2.*c*h*n*g+2.*c*k*n*g-2.*a*m*n*g+2.*b*m*n*g-2.*c*p*w*g+2.*c*r*w*g-2.*a*s*w*g
                    +2.*b*s*w*g-2.*h*m*n.powi(2)+2.*k*m*n.powi(2)-2.*p*s*w.powi(2)+2.*r*s*w.powi(2)
                    -2.*m*n*p*w+2.*m*n*r*w-2.*h*n*s*w+2.*k*n*s*w+2.*a*c*n.powi(2)*x
                    -2.*b*c*n.powi(2)*x+2.*h*m*n.powi(2)*x-2.*k*m*n.powi(2)*x+2.*a*c*w.powi(2)*x
                    -2.*b*c*w.powi(2)*x+2.*h*m*w.powi(2)*x-2.*k*m*w.powi(2)*x+2.*p*s*w.powi(2)*x
                    -2.*r*s*w.powi(2)*x+2.*n.powi(2)*p*s*x-2.*n.powi(2)*r*s*x;
        let bottom = 2.*(c.powi(2)*g.powi(2)-c.powi(2)*x*g.powi(2)-m.powi(2)*x*g.powi(2)
                    -s.powi(2)*x*g.powi(2)+2.*c*m*n*g+2.*c*s*w*g+m.powi(2)*n.powi(2)
                    +s.powi(2)*w.powi(2)+2.*m*n*s*w-c.powi(2)*n.powi(2)*x-m.powi(2)*n.powi(2)*x
                    -n.powi(2)*s.powi(2)*x-c.powi(2)*w.powi(2)*x-m.powi(2)*w.powi(2)*x
                    -s.powi(2)*w.powi(2)*x);

        let t1 = (top+to_root)/(bottom);
        let t2 = (top-to_root)/(bottom);

        if t1>0. {
            if t2>0. {
                if t1 < t2 {
                    t = t1;
                }
                else {
                    t = t2;
                }
            }
            else {
                t = t1;
            }
        }
        else if t2>0. {
            t = t2;
        }
        else {
            return None;
        }
        return Some(Point {x: v.p.x+v.x*t,y: v.p.y+v.y*t,z: v.p.z+v.z*t});
    }
    fn normal(&self,p: Point) -> Option<Vector> {
        unimplemented!()
    }
}
