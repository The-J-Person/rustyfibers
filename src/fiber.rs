//
// fiber.rs
// Defines all structs related to the structure of the fiber
//

extern crate roots;
use roots::Roots;
use roots::find_roots_quartic;

use math::{Point,Vector,Plane};
use std::f64;

pub trait Approximable {
    fn approx(&self, other: f64) -> bool;
}

impl Approximable for f64 {
    fn approx(&self, other: f64) -> bool {
        if self+f64::EPSILON>other && self-f64::EPSILON<other {
            return true;
        }
        false
    }
}

pub trait Geometry3D {
    fn point_is_on_surface(&self,p: Point) -> bool;
    fn check_collision(&self,v: Vector) -> bool;
    fn collision_point(&self,v: Vector) -> Option<Point>;
    fn normal(&self,p: Point) -> Option<Vector>;
    ///Returns the plane by which the ray can enter the structure for the first time.
    ///Used for detecting which segment of the fiber the ray is on.
    fn entry_plane(&self) -> Plane;
    ///Returns the plane by which the ray can exit the structure for the first time.
    fn exit_plane(&self) -> Plane;
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
    ///Returns theta and phi angle values if the given point is on the torus surface
    ///If the point is not on the torus surface, returns "None"
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
    #[allow(unused_variables)]
    fn check_collision(&self,v: Vector) -> bool {
        if self.collision_point(v) != None {
            return true;
        }
        false
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
        let mut roots_t = Vec::new();
        match find_roots_quartic(1., A, B, C, D) {
            Roots::Four(roots) => {
                roots_t.push(roots[0]);
                roots_t.push(roots[1]);
                roots_t.push(roots[2]);
                roots_t.push(roots[3]);
            },
            Roots::Three(roots) => {
                roots_t.push(roots[0]);
                roots_t.push(roots[1]);
                roots_t.push(roots[2]);
            },
            Roots::Two(roots) => {
                roots_t.push(roots[0]);
                roots_t.push(roots[1]);
            },
            Roots::One(roots) => {
                roots_t.push(roots[0]);
            },
            _ => {
                return None;
            },
        }
        let mut plausible: Vec<f64> = Vec::new();
        for r in roots_t {
            if r<=0.+f64::EPSILON || !self.point_is_on_surface(Point {x: v.p.x+v.x*r,y: v.p.y+v.y*r,z: v.p.z+v.z*r}) {
                continue;
            }
            plausible.push(r);
        }

        let mut t = f64::INFINITY;

        for r in plausible {
            if t>r {
                t = r;
            }
        }
        if t==f64::INFINITY {
            return None;
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
    ///Returns the plane by which the ray can enter the Torus for the first time.
    ///Used for detecting which segment of the fiber the ray is on.
    fn entry_plane(&self) -> Plane {
        let xb = (self.R+self.r*f64::cos(0.))*f64::cos(0.);
        let yb = (self.R+self.r*f64::cos(0.))*f64::sin(0.);
        let zb = self.r*f64::sin(0.);
        let xc = (self.R+self.r*f64::cos(f64::consts::PI/2.))*f64::cos(0.);
        let yc = (self.R+self.r*f64::cos(f64::consts::PI/2.))*f64::sin(0.);
        let zc = self.r*f64::sin(f64::consts::PI/2.);
        let xa = (self.R+self.r*f64::cos(f64::consts::PI))*f64::cos(0.);
        let ya = (self.R+self.r*f64::cos(f64::consts::PI))*f64::sin(0.);
        let za = self.r*f64::sin(f64::consts::PI);
        Plane{a: Point{x: xa, y: ya, z: za},
            b: Point{x: xb, y: yb, z: zb},
            c: Point{x: xc, y: yc, z: zc}}
    }
    ///Returns the plane by which the ray can exit the Torus for the first time.
    ///Used for detecting which segment of the fiber the ray is on.
    fn exit_plane(&self) -> Plane {
        let xc = (self.R+self.r*f64::cos(0.))*f64::cos(self.phi_max);
        let yc = (self.R+self.r*f64::cos(0.))*f64::sin(self.phi_max);
        let zc = self.r*f64::sin(0.);
        let xb = (self.R+self.r*f64::cos(f64::consts::PI/2.))*f64::cos(self.phi_max);
        let yb = (self.R+self.r*f64::cos(f64::consts::PI/2.))*f64::sin(self.phi_max);
        let zb = self.r*f64::sin(f64::consts::PI/2.);
        let xa = (self.R+self.r*f64::cos(f64::consts::PI))*f64::cos(self.phi_max);
        let ya = (self.R+self.r*f64::cos(f64::consts::PI))*f64::sin(self.phi_max);
        let za = self.r*f64::sin(f64::consts::PI);
        Plane{a: Point{x: xa, y: ya, z: za},
            b: Point{x: xb, y: yb, z: zb},
            c: Point{x: xc, y: yc, z: zc}}
    }
}


//
// Cylinder typic
//

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Cylinder {
    pub c: Vector,  //"central" vector, alternatively, the "spine" along which the cylinder stretches.
    pub r: f64,     //the radius
    pub length: f64,//endpoint
}

impl Geometry3D for Cylinder{
    #[allow(unused_variables)]
    fn point_is_on_surface(&self,p: Point) -> bool {
        //let x = ( p.x - self.c.p.x - (self.c.x*p.x-self.c.p.x)*self.c.x).powi(2)-self.r.powi(2);
        let cp  = self.c.cross_product(&Vector::new_from_points(p, self.c.p));
        let result = cp.length()/self.c.length();
        if result.approx(self.r) {
            return true;
        }
        false
    }
    #[allow(unused_variables)]
    fn check_collision(&self,v: Vector) -> bool {
        unimplemented!()
    }
    fn collision_point(&self,v: Vector) -> Option<Point> {
        //This function is disorderly and can probably be simplified.
        //Somehow.
        let t: f64;
        let x = self.c.x;
        let y = self.c.y;
        let z = self.c.z;
        let a = self.c.p.x;
        let b = v.p.x;
        let c = v.x;
        let h = self.c.p.y;
        let q = v.p.y;
        let j = v.y;
        let k = self.c.p.z;
        let m = v.p.z;
        let n = v.z;
        let r = self.r;
        let to_root = f64::sqrt((-2.*a*c*y.powi(2)-2.*a*c*z.powi(2)+2.*a*j*x*y+2.*a*n*x*z
                                +2.*b*c*y.powi(2)+2.*b*c*z.powi(2)-2.*b*j*x*y-2.*b*n*x*z
                                +2.*c*h*x*y+2.*c*k*x*z-2.*c*m*x*z-2.*c*q*x*y-2.*h*j*x.powi(2)
                                -2.*h*j*z.powi(2)+2.*h*n*y*z+2.*j*k*y*z-2.*j*m*y*z+2.*j*q*x.powi(2)
                                +2.*j*q*z.powi(2)-2.*k*n*x.powi(2)-2.*k*n*y.powi(2)
                                +2.*m*n*x.powi(2)+2.*m*n*y.powi(2)-2.*n*q*y*z).powi(2)
                                -4.*(c.powi(2)*y.powi(2)+c.powi(2)*z.powi(2)-2.*c*j*x*y
                                -2.*c*n*x*z+j.powi(2)*x.powi(2)+j.powi(2)*z.powi(2)
                                -2.*j*n*y*z+n.powi(2)*x.powi(2)
                                +n.powi(2)*y.powi(2))*(a.powi(2)*y.powi(2)+a.powi(2)*z.powi(2)
                                -2.*a*b*y.powi(2)-2.*a*b*z.powi(2)-2.*a*h*x*y-2.*a*k*x*z
                                +2.*a*m*x*z+2.*a*q*x*y+b.powi(2)*y.powi(2)+b.powi(2)*z.powi(2)
                                +2.*b*h*x*y+2.*b*k*x*z-2.*b*m*x*z-2.*b*q*x*y+h.powi(2)*x.powi(2)
                                +h.powi(2)*z.powi(2)-2.*h*k*y*z+2.*h*m*y*z-2.*h*q*x.powi(2)
                                -2.*h*q*z.powi(2)+k.powi(2)*x.powi(2)+k.powi(2)*y.powi(2)
                                -2.*k*m*x.powi(2)-2.*k*m*y.powi(2)+2.*k*q*y*z+m.powi(2)*x.powi(2)
                                +m.powi(2)*y.powi(2)-2.*m*q*y*z+q.powi(2)*x.powi(2)
                                +q.powi(2)*z.powi(2)-r.powi(2)));
        let t1 = (to_root+2.*a*c*y.powi(2)+2.*a*c*z.powi(2)-2.*a*j*x*y-2.*a*n*x*z-2.*b*c*y.powi(2)
                    -2.*b*c*z.powi(2)+2.*b*j*x*y+2.*b*n*x*z-2.*c*h*x*y-2.*c*k*x*z+2.*c*m*x*z
                    +2.*c*q*x*y+2.*h*j*x.powi(2)+2.*h*j*z.powi(2)-2.*h*n*y*z-2.*j*k*y*z+2.*j*m*y*z
                    -2.*j*q*x.powi(2)-2.*j*q*z.powi(2)+2.*k*n*x.powi(2)+2.*k*n*y.powi(2)
                    -2.*m*n*x.powi(2)-2.*m*n*y.powi(2)+2.*n*q*y*z)/(2.*(c.powi(2)*y.powi(2)
                    +c.powi(2)*z.powi(2)-2.*c*j*x*y-2.*c*n*x*z+j.powi(2)*x.powi(2)
                    +j.powi(2)*z.powi(2)-2.*j*n*y*z+n.powi(2)*x.powi(2)+n.powi(2)*y.powi(2)));
        let t2 = (-to_root+2.*a*c*y.powi(2)+2.*a*c*z.powi(2)-2.*a*j*x*y-2.*a*n*x*z-2.*b*c*y.powi(2)
                    -2.*b*c*z.powi(2)+2.*b*j*x*y+2.*b*n*x*z-2.*c*h*x*y-2.*c*k*x*z+2.*c*m*x*z
                    +2.*c*q*x*y+2.*h*j*x.powi(2)+2.*h*j*z.powi(2)-2.*h*n*y*z-2.*j*k*y*z+2.*j*m*y*z
                    -2.*j*q*x.powi(2)-2.*j*q*z.powi(2)+2.*k*n*x.powi(2)+2.*k*n*y.powi(2)
                    -2.*m*n*x.powi(2)-2.*m*n*y.powi(2)+2.*n*q*y*z)/(2.*(c.powi(2)*y.powi(2)
                    +c.powi(2)*z.powi(2)-2.*c*j*x*y-2.*c*n*x*z+j.powi(2)*x.powi(2)
                    +j.powi(2)*z.powi(2)-2.*j*n*y*z+n.powi(2)*x.powi(2)+n.powi(2)*y.powi(2)));
        if t1>0.+f64::EPSILON {
            if t2>0.+f64::EPSILON {
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
        else if t2>0.+f64::EPSILON {
            t = t2;
        }
        else {
            return None;
        }
        return Some(Point {x: v.p.x+v.x*t,y: v.p.y+v.y*t,z: v.p.z+v.z*t});
    }
    fn normal(&self,q: Point) -> Option<Vector> {
        //Assuming the point is on the surface...
        //We should actually test for that. TODO.
        let t = (q.x*self.c.x-self.c.p.x*self.c.x
                +q.y*self.c.y-self.c.p.y*self.c.y
                +q.z*self.c.z-self.c.p.z*self.c.z)
                /(self.c.x.powi(2)+self.c.y.powi(2)+self.c.z.powi(2));
        let mut v = Vector{p: q,
                            x: q.x-self.c.p.x-self.c.x*t,
                            y: q.y-self.c.p.y-self.c.y*t,
                            z: q.z-self.c.p.z-self.c.z*t};
        //sanity check, but commented out because I won't rely on it
        // if v.length()!=self.r {
        //     //If sanity test failed, return nothing
        //     return None;
        // }
        v.normalize();
        return Some(v);
    }
    fn entry_plane(&self) -> Plane {
        //We have the normal.
        //The plane can also be expressed as ax+by+cz=d,
        //where a,b,c are the normal.
        let d=self.c.x*self.c.p.x+self.c.y*self.c.p.y+self.c.z*self.c.p.z;
        let pa=self.c.p;
        let mut pb: Point;
        //With bad rounding, this is an error waiting to happen... TODO add epsilon?
        // For every flat plane, a point where two axes are zero is contained in it.
        if self.c.x!=0. {
            pb = Point{x: d/self.c.x, y: 0., z: 0.};
            //dirty hack
            if pb == pa {
                pb = Point{x: self.c.p.x, y: 1., z: 0.};
            }
        }
        else if self.c.y!=0. {
            pb = Point{x: 0., y: d/self.c.y, z: 0.};
            //dirty hack
            if pb == pa {
                pb = Point{x: 1., y: self.c.p.y, z: 0.};
            }
        }
        else {
            pb = Point{x: 0., y: 0., z: d/self.c.z};
            //dirty hack
            if pb == pa {
                pb = Point{x: 1., y: 0., z: self.c.p.z};
            }
        }



        let v = Vector{p: pa, x: pb.x-pa.x, y: pb.y-pa.y, z: pb.z-pa.z};
        //Cross product between a vector on the plane and the central vector,
        // gives us another point on the plane.
        let cp = v.cross_product(&self.c);
        let pc = Point{x: cp.p.x+cp.x, y: cp.p.y+cp.y, z: cp.p.z+cp.z};
        Plane{a: pa, b: pb, c: pc}
    }
    fn exit_plane(&self) -> Plane {
        let mut plane = self.entry_plane();
        plane.a.x += self.c.x*self.length;
        plane.a.y += self.c.y*self.length;
        plane.a.z += self.c.z*self.length;

        plane.b.x += self.c.x*self.length;
        plane.b.y += self.c.y*self.length;
        plane.b.z += self.c.z*self.length;

        plane.c.x += self.c.x*self.length;
        plane.c.y += self.c.y*self.length;
        plane.c.z += self.c.z*self.length;

        let temp = plane.b;
        plane.b = plane.c;
        plane.c = temp;

        return plane;
    }
}

//
// Cone typic
//

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Cone {
    pub c: Vector,  //"central" vector, its starting point being the cone's "pointy end"
    pub a: f64,     //Angle between central vector and cone surface
    pub length: f64,//endpoint
}

impl Geometry3D for Cone{
    #[allow(unused_variables)]
    fn point_is_on_surface(&self,p: Point) -> bool {
        unimplemented!()
    }
    #[allow(unused_variables)]
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
        let x = (f64::cos(self.a).powi(2) as f32) as f64; //Stupid workaround for rounding errors
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

        if t1>0.+f64::EPSILON {
            if t2>0.+f64::EPSILON {
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
        else if t2>0.+f64::EPSILON {
            t = t2;
        }
        else {
            return None;
        }
        return Some(Point {x: v.p.x+v.x*t,y: v.p.y+v.y*t,z: v.p.z+v.z*t});
    }
    fn normal(&self,p: Point) -> Option<Vector> {
        //Assuming point is on the surface
        let a = p.x;
        let b = self.c.p.x;
        let c = self.c.x;
        let f = p.y;
        let g = self.c.p.y;
        let h = self.c.y;
        let k = p.z;
        let m = self.c.p.z;
        let n = self.c.z;
        let x = (f64::cos(f64::consts::PI/2.-self.a).powi(2) as f32) as f64; //Weird rounding
        let t: f64;

        let to_root = f64::sqrt(x*(c.powi(2)+h.powi(2)+n.powi(2)-x)*(a.powi(2)*h.powi(2)+a.powi(2)*n.powi(2)-2.*a*b*h.powi(2)-2.*a*b*n.powi(2)-2.*a*c*f*h+2.*a*c*g*h-2.*a*c*k*n+2.*a*c*m*n+b.powi(2)*h.powi(2)+b.powi(2)*n.powi(2)+2.*b*c*f*h-2.*b*c*g*h+2.*b*c*k*n-2.*b*c*m*n+c.powi(2)*f.powi(2)-2.*c.powi(2)*f*g+c.powi(2)*g.powi(2)+c.powi(2)*k.powi(2)-2.*c.powi(2)*k*m+c.powi(2)*m.powi(2)+f.powi(2)*n.powi(2)-2.*f*g*n.powi(2)-2.*f*h*k*n+2.*f*h*m*n+g.powi(2)*n.powi(2)+2.*g*h*k*n-2.*g*h*m*n+h.powi(2)*k.powi(2)-2.*h.powi(2)*k*m+h.powi(2)*m.powi(2)));
        let top = a*c.powi(3)+a*c*h.powi(2)+a*c*n.powi(2)-a*c*x-b*c.powi(3)-b*c*h.powi(2)-b*c*n.powi(2)+b*c*x+c.powi(2)*f*h-c.powi(2)*g*h+c.powi(2)*k*n-c.powi(2)*m*n+f*h.powi(3)+f*h*n.powi(2)-f*h*x-g*h.powi(3)-g*h*n.powi(2)+g*h*x+h.powi(2)*k*n-h.powi(2)*m*n+k*n.powi(3)-k*n*x-m*n.powi(3)+m*n*x;
        let bottom = (c.powi(2)+h.powi(2)+n.powi(2))*(c.powi(2)+h.powi(2)+n.powi(2)-x);

        assert!(bottom!=0.,"A cone has a 90 degree angle!");

        let t1=(top+to_root)/bottom;
        let t2=(top-to_root)/bottom;

        // test area
        // println!("t1 = {:?}", t1);
        // println!("t2 = {:?}", t2);
        t=t1;
        // end test area

        let mut v = Vector{p: p,
                            x: p.x-self.c.p.x-self.c.x*t,
                            y: p.y-self.c.p.y-self.c.y*t,
                            z: p.z-self.c.p.z-self.c.z*t};
        v.normalize();
        return Some(v);
    }
    #[allow(unused_variables)]
    fn entry_plane(&self) -> Plane {
        unimplemented!()
    }
    #[allow(unused_variables)]
    fn exit_plane(&self) -> Plane {
        unimplemented!()
    }
}


//
// Tests
//

#[test]
fn torus_basic_collision() {
    let ray = Vector {p: Point{x: 0.,y: 8.,z: 0.},x: 0.,y: -1.,z: 0.};
    let torus = Torus {R: 2., r: 1., phi_max: f64::consts::PI};
    let point = torus.collision_point(ray);
    match point {
        Some(p) => {
            assert_eq!(p.x,0.);
            assert_eq!(p.y,3.);
            assert_eq!(p.z,0.);
        },
        _ => {
            panic!();
        },
    }
}

// #[test]
// fn torus_plane_generation() {
//     let torus = Torus {R: 2., r: 1., phi_max: f64::consts::PI/2.};
//     let gp = Point::new_blank();
//     let p1 = torus.entry_plane();
//     let p2 = torus.exit_plane();
//     let c1 = Plane{a: Point{x: 3.,y: 0.,z: 0.},b: Point{x: 2.,y: 0.,z: 1.},c: Point{x: 1.,y: 0.,z: 0.}};
//     let c2 = Plane{a: Point{x: -3.,y: 0.,z: 0.},b: Point{x: -2.,y: 0.,z: 1.},c: Point{x: -1.,y: 0.,z: 0.}};
//     // println!("{:?}", p1 );
//     // println!("{:?}", p2 );
//     println!("{:?}", p1.normal(gp) );
//     println!("{:?}", p2.normal(gp) );

//    panic!();
// }

#[test]
fn cylinder_basic_collision() {
    let ray = Vector {p: Point{x: 5.,y: 0.,z: 0.},x: -1.,y: 0.,z: 0.};
    let cylinder = Cylinder {c: Vector{p: Point{x: 0.,y: 0., z: 0.},x: 0., y:0., z: 1.}, r: 1., length: 4.};
    let point = cylinder.collision_point(ray);
    match point {
        Some(p) => {
            assert_eq!(p.x,1.);
            assert_eq!(p.y,0.);
            assert_eq!(p.z,0.);
        },
        _ => {
            panic!();
        },
    }
}

#[test]
fn cylinder_basic_normal() {
    let cylinder = Cylinder {c: Vector{p: Point{x: 0.,y: 0., z: 0.},x: 0., y:0., z: 1.}, r: 1., length: 4.};
    let point = Point{x: 1.,y: 0.,z: 0.};
    let v = cylinder.normal(point).unwrap() as Vector;
    assert_eq!(v.x,1.);
    assert_eq!(v.y,0.);
    assert_eq!(v.z,0.);
}

#[test]
fn cone_basic_collision() {
    let ray = Vector {p: Point{x: 0.,y: 5.,z: 1.},x: 0.,y: -1.,z: 0.};
    let cone = Cone {c: Vector{p: Point{x: 0.,y: 0., z: 0.},x: 0., y:0., z: 1.}, a: f64::consts::PI/4., length: 4.};
    let point = cone.collision_point(ray);
    match point {
        Some(p) => {
            println!("{:?}", p );
            assert_eq!(p.x,0.);
            assert_eq!(p.y,1.);
            assert_eq!(p.z,1.);
        },
        _ => {
            panic!();
        },
    }
}

// #[test]
// fn cone_basic_normal() {
//     let cone = Cone {c: Vector{p: Point{x: 0.,y: 0., z: 0.},x: 0., y:0., z: 1.}, a: f64::consts::PI/4., length: 4.};
//     let point = Point{x: 0.,y: 1.,z: 1.};
//     let v = cone.normal(point).unwrap() as Vector;;
//     println!("{:?}",v );
//     assert_eq!(v.x,0.);
//     assert_eq!(v.y,f64::cos(f64::consts::PI/4.));
//     assert_eq!(v.z,f64::cos(f64::consts::PI/4.));
// }
